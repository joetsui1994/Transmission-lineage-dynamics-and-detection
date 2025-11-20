#!/usr/bin/env python

import argparse
import yaml
import sys
import os
import math as m
import numpy as np
from tqdm import tqdm


# ===============================
# Utility / Factory Functions
# ===============================

def validate_params(name, params, required_keys):
    missing = [key for key in required_keys if key not in params]
    if missing:
        raise ValueError(f"configuration error: '{name}' is missing required parameters: {missing}")

def build_importation_function(param_config):
    """
    Returns a function f(t) -> int (nuber of new imports at time t)
    based on the configuration dictionary.
    """
    if 'type' not in param_config:
        raise ValueError("configuration error: importation rate is missing 'type'.")
    
    func_type = param_config['type']
    
    if func_type == 'constant':
        if 'value' not in param_config:
            raise ValueError("configuration error: constant importation requires a 'value'.")
        val = int(param_config['value'])
        return lambda t: val

    elif func_type == 'exponential':
        p = param_config.get('params', {})
        validate_params('importation_rate', p, ['M0', 'rate'])
        M0 = p['M0']
        rate = p['rate']
        return lambda t: round(M0 * m.exp(rate * t))
        
    else:
        # catch unknown function types
        raise ValueError(f"configuration error: unknown function type '{func_type}'."
                         f"supported types are: constant and exponential.")

def get_lineage_size(import_time, t, r):
    return round(m.exp(r * (t - import_time)))

def generate_sampling_props(config):
    defaults = {
        'split_threshold': 1.0,
        'num_divisions_lower': 101,
        'num_divisions_upper': 0
    }
    
    # extract grid config, using defaults if missing
    grid = config.get('grid', defaults)
    
    split = grid.get('split_threshold', 1.0)
    n_low = grid.get('num_divisions_lower', 101)
    n_high = grid.get('num_divisions_upper', 0)
    
    props = []
    
    # lower section (0 to split)
    if n_low > 0:
        props.append(np.linspace(0, split, n_low + 1)) # +1 to include split
        
    # upper section (split to 1.0)
    if n_high > 0:
        start = split + ((1.0 - split) / n_high)
        props.append(np.linspace(start, 1.0, n_high))
        
    if not props:
        return np.array([0.0])
        
    # combine and ensure unique sorted
    return np.unique(np.concatenate(props))

# ===============================
# Main execution
# ===============================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run deterministic lineage simulation")
    parser.add_argument('--config', type=str, required=True, help="Path to YAML config file")
    args = parser.parse_args()

    # load config file
    print('loading config file...')
    with open(args.config, 'r') as f:
        config = yaml.safe_load(f)

    # extract parameters from config
    sim_conf = config['simulation']
    seed = sim_conf.get('random_seed', 42)  # default to 42 if not specified
    observe_time = sim_conf.get('observe_time', 40)  # default to 40 if not specified
    output_dir = sim_conf.get('output_dir', './results')  # default to ./results if not specified

    # set random seed for reproducibility
    np.random.seed(seed)

    # epidemiological Parameters
    params_conf = config['epi_parameters']
    r = params_conf['local_growth_rate']
    try:
        importation_func = build_importation_function(params_conf['importation_rate'])
    except ValueError as e:
        print(e)
        sys.exit(1)

    # generate sampling proportions from config
    sampling_conf = config.get('sampling_parameters', {})
    num_sampling_simulations = sampling_conf.get('num_sampling_simulations', 50)  # default to 50 if not specified
    sampling_props = generate_sampling_props(sampling_conf)

    # run simulation (evolve lineages deterministically)
    print('simulating outbreaks...')
    lineages_tracker = []
    for time in tqdm(range(observe_time + 1)):
        new_lineages_count = importation_func(time)
        lineage_size = get_lineage_size(time, observe_time, r)
        lineages_tracker += [(time, lineage_size)] * new_lineages_count

    # simulate sampling
    # iterate over lineages and assign start/end indices
    curr_i = 0
    for i, (time, size) in enumerate(lineages_tracker):
        lineages_tracker[i] = (time, size, curr_i, curr_i + size - 1)
        curr_i += size

    # create array of infected individuals labeled by lineage index
    labelled_infected = np.repeat(
        np.arange(len(lineages_tracker)), 
        [item[1] for item in lineages_tracker]
    )

    # convert tracker to numpy array for faster indexing
    lineages_tracker = np.array(lineages_tracker)

    ## useful quantities
    num_lineages = len(lineages_tracker)
    infected_population_size = len(labelled_infected)

    # storage for results
    sp_detection_probs = []
    sp_detection_probs_95CI = []
    sp_detected_import_times = []
    sp_detected_import_times_95CI = []

    # iterate over sampling proportions
    print('simulating sampling...')
    for sp in tqdm(sampling_props):
        detection_probs = []
        detection_probs_time = [[] for _ in range(observe_time + 1)]
        sample_size = int(infected_population_size * sp)

        for _ in range(num_sampling_simulations):
            if sp > 0.5: # sample those NOT included for efficiency
                sampled_infected = np.random.choice(infected_population_size, infected_population_size - sample_size, replace=False)
                sampled_lineages = np.unique(labelled_infected[~np.isin(np.arange(infected_population_size), sampled_infected)])
            else:
                sampled_infected = np.random.choice(infected_population_size, sample_size, replace=False)
                sampled_lineages = np.unique(labelled_infected[sampled_infected])

            # calculate detection probability
            detection_probs.append(len(sampled_lineages) / num_lineages)

            # get importation time of detected lineages
            if len(sampled_lineages) > 0:
                sampled_lineages_importation_times = lineages_tracker[sampled_lineages][:, 0]
                # count occurrences for each importation time
                counts = np.bincount(sampled_lineages_importation_times.astype(int), minlength=observe_time + 1)
                # only take up to observe_time (in case lineages started exactly at observe_time are tracked)
                counts = counts[:observe_time + 1]

                for t_idx, count in enumerate(counts):
                    detection_probs_time[t_idx].append(count)

            else:
                for t_idx in range(observe_time + 1):
                    detection_probs_time[t_idx].append(0)

        # aggregate results
        sp_detection_probs.append(np.median(detection_probs))
        sp_detection_probs_95CI.append((np.percentile(detection_probs, 2.5), np.percentile(detection_probs, 97.5)))

        detection_probs_time = np.array(detection_probs_time)
        sp_detected_import_times.append(np.median(detection_probs_time, axis=1))
        sp_detected_import_times_95CI.append((
            np.percentile(detection_probs_time, 2.5, axis=1),
            np.percentile(detection_probs_time, 97.5, axis=1)
        ))


    # ===============================
    # Export results
    # ===============================

    os.makedirs(output_dir, exist_ok=True)
    print("exporting results...")

    # export detection probs
    with open(os.path.join(output_dir, 'detection_probs.csv'), 'w') as outfile:
        outfile.write('sampling_prop,median,lw_0p025,up_0p975\n')
        for i, t in enumerate(sampling_props):
            outfile.write(f"{t},{sp_detection_probs[i]},{sp_detection_probs_95CI[i][0]},{sp_detection_probs_95CI[i][1]}\n")

    header_indices = [str(i) for i in range(observe_time + 1)]

    # export detection probs by importation time (median)
    with open(os.path.join(output_dir, 'detection_probs_time.csv'), 'w') as outfile:
        outfile.write('sampling_prop,' + ','.join(header_indices) + '\n')
        for i, t in enumerate(sampling_props):
            row_vals = ','.join([f"{x}" for x in sp_detected_import_times[i]])
            outfile.write(f"{t},{row_vals}\n")

    # export detection probs by importation time (lower CI)
    with open(os.path.join(output_dir, 'detection_probs_time.95CI_lw.csv'), 'w') as outfile:
        outfile.write('sampling_prop,' + ','.join(header_indices) + '\n')
        for i, t in enumerate(sampling_props):
            row_vals = ','.join([f"{x}" for x in sp_detected_import_times_95CI[i][0]])
            outfile.write(f"{t},{row_vals}\n")

    # export detection probs by importation time (upper CI)
    with open(os.path.join(output_dir, 'detection_probs_time.95CI_up.csv'), 'w') as outfile:
        outfile.write('sampling_prop,' + ','.join(header_indices) + '\n')
        for i, t in enumerate(sampling_props):
            row_vals = ','.join([f"{x}" for x in sp_detected_import_times_95CI[i][1]])
            outfile.write(f"{t},{row_vals}\n")