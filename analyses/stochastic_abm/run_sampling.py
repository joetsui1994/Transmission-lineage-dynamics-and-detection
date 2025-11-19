#!/usr/bin/env python

import argparse
import yaml
import sys
import os
import numpy as np
import pandas as pd
from tqdm import tqdm


# ===============================
# Utility Functions
# ===============================

def validate_config(config):
    required = ['input_file', 'output_dir', 'num_draws']
    for k in required:
        if k not in config:
            raise ValueError(f"Config Error: Section 'sampling' missing '{k}'")

def generate_sampling_fractions(conf):
    """Generates the fine-grained grid of sampling fractions."""
    grid_conf = conf.get('sampling_grid', {})
    split = grid_conf.get('split_threshold', 0.05)
    n_low = grid_conf.get('num_divisions_lower', 30)
    n_high = grid_conf.get('num_divisions_upper', 30)
    
    lower = np.linspace(0, split, n_low)
    step = split / n_low if n_low > 0 else 0.01
    upper = np.linspace(split + step, 1.0, n_high)
    return np.unique(np.concatenate((lower, upper)))


# ===============================
# Main execution
# ===============================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run Combined Sampling Analysis")
    parser.add_argument('--config', type=str, required=True, help="Path to YAML config")
    parser.add_argument('--input', type=str, help="Override input file")
    parser.add_argument('--output', type=str, help="Override output dir")
    args = parser.parse_args()

    # load Config
    with open(args.config, 'r') as f:
        config = yaml.safe_load(f)

    input_file = args.input if args.input else config['input_file']
    output_dir = args.output if args.output else config['output_dir']
    num_draws = config['num_draws']
    rng = np.random.default_rng(seed=config.get('random_seed', None))

    # set up sampling grid (for detection probabilities)
    sfs_grid = generate_sampling_fractions(config)
    # temporal targets (for importation intensity over time)
    temporal_conf = config.get('temporal_analysis', {'enabled': False})
    sfs_temporal = set(temporal_conf.get('target_fractions', [])) if temporal_conf['enabled'] else set()

    # merge them so we can iterate only once, but flag which ones need temporal calculations
    all_sfs = np.unique(np.concatenate((sfs_grid, list(sfs_temporal))))

    print(f"loading data from: {input_file}")
    try:
        df = pd.read_csv(input_file, compression='gzip' if input_file.endswith('.gz') else None)
    except FileNotFoundError:
        print(f"error: input file not found.")
        sys.exit(1)

    # pre-filter
    infected_data = df[df['state'] >= 1].copy()

    # determine time-range for temporal analysis
    if temporal_conf['enabled']:
        valid_times = df[df['cluster_imported_time'] > -1]['cluster_imported_time']
        if not valid_times.empty:
            t_min, t_max = int(valid_times.min()), int(valid_times.max())
            time_range = np.arange(t_min, t_max + 1)
        else:
            time_range = np.array([])
    else:
        time_range = np.array([])

    # prepare output storage
    results_scalar = {sf: [] for sf in all_sfs}
    results_temporal = {sf: {t: [] for t in time_range} for sf in sfs_temporal}

    # group by experiment
    experiments = infected_data.groupby('experiment')

    for exp_id, df_exp in tqdm(experiments, desc="processing experiments"):
        true_num_clusters = df_exp['cluster'].nunique()
        if true_num_clusters == 0: continue
        
        num_infected = len(df_exp)
        indices = df_exp.index.to_numpy()
        
        # pre-extract cluster info to numpy for speed
        exp_clusters = df_exp['cluster'].to_numpy()
        exp_import_times = df_exp['cluster_imported_time'].to_numpy()

        for sf in all_sfs:
            num_sample = int(sf * num_infected)
            
            # Calculate Scalar Detection (Fast)
            if num_sample == 0:
                results_scalar[sf].extend([0.0] * num_draws)
                if sf in sfs_temporal:
                    for t in time_range: results_temporal[sf][t].extend([0]*num_draws)
                continue
            
            if num_sample >= num_infected:
                results_scalar[sf].extend([1.0] * num_draws)
                # full sampling: just count the actuals once and replicate
                if sf in sfs_temporal:
                    # count unique clusters per time
                    counts = df_exp.drop_duplicates('cluster')['cluster_imported_time'].value_counts()
                    for t in time_range:
                        results_temporal[sf][t].extend([counts.get(t, 0)] * num_draws)
                continue

            # random draws
            for _ in range(num_draws):
                sampled_idxs_rel = rng.choice(len(indices), num_sample, replace=False)
                
                # get the clusters found in this sample
                found_cluster_ids = exp_clusters[sampled_idxs_rel]
                unique_found_clusters, unique_idx = np.unique(found_cluster_ids, return_index=True)
                
                # scalar metric
                detection_rate = len(unique_found_clusters) / true_num_clusters
                results_scalar[sf].append(detection_rate)

                # temporal metric (conditional)
                if sf in sfs_temporal:                    
                    found_import_times = exp_import_times[sampled_idxs_rel][unique_idx]
                    valid_times_mask = found_import_times >= 0
                    counts = np.bincount(found_import_times[valid_times_mask].astype(int), minlength=time_range[-1]+1)
                    
                    for t in time_range:
                        if t < len(counts):
                            results_temporal[sf][t].append(counts[t])
                        else:
                            results_temporal[sf][t].append(0)


    # ===============================
    # Export results
    # ===============================

    os.makedirs(output_dir, exist_ok=True)

    # export scalar results
    print("exporting scalar lineage detection probabilities...")
    stats_rows = [] = []
    for sf in sorted(all_sfs):
        data = results_scalar[sf]
        if data:
            stats_rows.append({
                'sf': sf, 
                'median': np.median(data),
                'lw_95CI': np.percentile(data, 2.5),
                'up_95CI': np.percentile(data, 97.5)
            })
    pd.DataFrame(stats_rows).to_csv(os.path.join(output_dir, 'detection_probs.median_95CI.csv'), index=False)

    # export temporal results
    if temporal_conf['enabled']:
        print("exporting Temporal Detection Numbers...")
        temporal_rows = []
        for sf in sorted(list(sfs_temporal)):
            for t in time_range:
                data = results_temporal[sf][t]
                if data:
                    temporal_rows.append({
                        'sf': sf,
                        'time': t,
                        'median': np.median(data),
                        'lw_95CI': np.percentile(data, 2.5),
                        'up_95CI': np.percentile(data, 97.5)
                    })
        pd.DataFrame(temporal_rows).to_csv(os.path.join(output_dir, 'detection_numbers_by_time.csv'), index=False)