#!/usr/bin/env python

from tqdm import tqdm
import pandas as pd
import numpy as np
import json
import os


# specify path to simulated data
infile = './results/T40.csv.gz'
# read in simulated data
simulated_data = pd.read_csv(infile, compression='gzip')

# specify number of draws
num_draws = 200

# sampling proportions to calculate lineage detection rate
split_sf = 0.05 # more fine-grained sampling for sfs below split_sf
num_div_lw = 30 # number of divisions below split_sf
num_div_up = 30 # number of divisions above split_sf
sfs = np.concatenate((np.linspace(0, split_sf, num_div_lw),
                      np.linspace(split_sf + split_sf/num_div_lw, 1, num_div_up)))

# precompute variables for efficiency
num_experiments = len(simulated_data['experiment'].unique())
true_num_clusters = len(simulated_data['cluster'].unique())
infected_data = simulated_data[simulated_data['state'] >= 1]  # filter once for infected individuals

# prepare detection probs dictionary
sf_exp_detection_probs = {}
for sf in sfs:
    print(sf)
    exp_detection_probs = []

    for exp_i in tqdm(range(num_experiments)):
        # filter data for the current experiment
        df_exp = infected_data[infected_data['experiment'] == exp_i]
        num_infected = len(df_exp)
        num_sample = int(sf * num_infected)

        if num_sample == 0:
            continue  # skip if no individuals can be sampled

        # perform sampling and calculate detection probs
        for _ in range(num_draws):
            sampled_indices = np.random.choice(df_exp.index, num_sample, replace=False)
            sampled_clusters = df_exp.loc[sampled_indices, 'cluster'].nunique()
            detection_rate = sampled_clusters / true_num_clusters
            exp_detection_probs.append(detection_rate)

    # store results for the current sampling fraction
    sf_exp_detection_probs[sf] = exp_detection_probs

# for each sampling proportion, calculate median and 95% CI
medians = [np.median(sf_exp_detection_probs[sf]) if sf_exp_detection_probs[sf] else np.nan for sf in sfs]
lws_95CI = [np.percentile(sf_exp_detection_probs[sf], 2.5) if sf_exp_detection_probs[sf] else np.nan for sf in sfs]
ups_95CI = [np.percentile(sf_exp_detection_probs[sf], 97.5) if sf_exp_detection_probs[sf] else np.nan for sf in sfs]

# specify output folder
out_dir = './results'

# export detection probs
with open(os.path.join(out_dir, 'detection_probs.csv'), 'w+') as outfile:
    outfile.write(json.dumps(sf_exp_detection_probs))
# export detection probs (median, interquartile ranges)
with open(os.path.join(out_dir, 'detection_probs.median_95CI.csv'), 'w+') as outfile:
    outfile.write('sf,median,lw_95CI,up_95CI\n')
    outfile.write('\n'.join(['%f,%f,%f,%f' % (sf, medians[sf_i], lws_95CI[sf_i], ups_95CI[sf_i]) for sf_i, sf in enumerate(sfs)]))