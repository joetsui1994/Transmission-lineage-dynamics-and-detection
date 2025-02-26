#!/usr/bin/env python

from tqdm import tqdm
import numpy as np
import math as m


# ===============================
# Utility functions
# ===============================

# exponentially increasing/decreasing importation rate
def exponential_importation(t, M0, rate):
    return round(M0 * m.exp(rate * t))

# update lineage size assuming local exponential growth
def get_lineage_size(import_time, t, r):
    return round(m.exp(r * (t - import_time)))


# ===============================
# Simulation parameters
# ===============================

observe_time = 40
r = 0.25
importation_M0 = 100
importation_m = 0


# ===============================
# Simulate outbreak
# ===============================

# evolve lineages deterministically
lineages_tracker = []
for time in tqdm(range(observe_time + 1)):
    lineages_tracker += [(time, get_lineage_size(time, observe_time, r))] * exponential_importation(time, importation_M0, importation_m)


# ===============================
# Simulate sampling
# ===============================

# number of (sampling) simulations
num_sampling_simulations = 50

# iterate over lineages and assign each lineage start_i and end_i according to its index and size
curr_i = 0
for i, (time, size) in enumerate(lineages_tracker):
    lineages_tracker[i] = (time, size, curr_i, curr_i + size - 1)
    curr_i += size

# create an array of infected individuals labeled by their lineage index (for sampling)
labelled_infected = np.repeat(np.arange(len(lineages_tracker)), [size for _, size, _, _ in lineages_tracker])

# convert to numpy array
lineages_tracker = np.array(lineages_tracker)

# set up random generator
rng = np.random.default_rng(seed=0)

## useful quantities
num_lineages = len(lineages_tracker)
infected_population_size = len(labelled_infected)

# sampling proportions
sampling_props = [0.001 * i for i in range(101)]
# store detection probs (i.e. proportion of lineages detected)
sp_detection_probs = []
# store detection probs 95% CI
sp_detection_probs_95CI = []
# store detection probs by time of importation
sp_detected_import_times = []
# store detection probs by time of importation 95% CI
sp_detected_import_times_95CI = []

# iterate over sampling proportions
for sampling_prop in tqdm(sampling_props):
    detection_probs = []
    detection_probs_time = [[] for _ in range(observe_time)]
    sample_size = int(infected_population_size * sampling_prop)

    for _ in range(num_sampling_simulations):
        # to speed up the simulation, we sample those agents that are not included when sampling_prop > 0.5
        if sampling_prop > 0.5:
            sampled_infected = rng.choice(infected_population_size, infected_population_size - sample_size, replace=False)
            sampled_lineages = np.unique(labelled_infected[~np.isin(np.arange(infected_population_size), sampled_infected)])
        else:
            sampled_infected = rng.choice(infected_population_size, sample_size, replace=False)
            sampled_lineages = np.unique(labelled_infected[sampled_infected])

        # calculate detection prob
        detection_probs.append(len(sampled_lineages) / num_lineages)

        # get importation time of sampled lineages
        sampled_lineages_importation_times = lineages_tracker[sampled_lineages][:,0]
        # get count for each time
        for time in range(observe_time):
            detection_probs_time[time].append(np.sum(sampled_lineages_importation_times == time))

    sp_detection_probs.append(np.median(detection_probs))
    sp_detection_probs_95CI.append((np.percentile(detection_probs, 2.5), np.percentile(detection_probs, 97.5)))
    detection_probs_time = np.array(detection_probs_time)
    sp_detected_import_times.append(np.median(detection_probs_time, axis=1))
    sp_detected_import_times_95CI.append((np.percentile(detection_probs_time, 2.5, axis=1),
                                          np.percentile(detection_probs_time, 97.5, axis=1)))


# ===============================
# Export results
# ===============================

# specify output folder
output_folder = './results'

# export detection probs (median and 95% CI)
with open('./%s/detection_probs.csv' % output_folder, 'w') as outfile:
    outfile.write('sampling_prop,median,lw_0p025,up_0p975\n')
    outfile.write('\n'.join(['%f,%f,%f,%f' % (
        t, sp_detection_probs[i],
        sp_detection_probs_95CI[i][0], sp_detection_probs_95CI[i][1])
        for i, t in enumerate(sampling_props)]))

# export detection probs by importation time (median)
with open('./%s/detection_probs_time.csv' % output_folder, 'w') as outfile:
    outfile.write('sampling_prop,' + ','.join([str(i) for i in range(observe_time)]) + '\n')
    outfile.write('\n'.join([','.join(['%f' % t] + ['%f' % x for x in sp_detected_import_times[i]])
                             for i, t in enumerate(sampling_props)]))

# export detection probs by importation time (95% CI)
# lower bound
with open('./%s/detection_probs_time.95CI_lw.csv' % output_folder, 'w') as outfile:
    outfile.write('sampling_prop,' + ','.join([str(i) for i in range(observe_time)]) + '\n')
    outfile.write('\n'.join([','.join(['%f' % t] + ['%f' % x for x in sp_detected_import_times_95CI[i][0]])
                             for i, t in enumerate(sampling_props)]))
# upper bound
with open('./%s/detection_probs_time.95CI_up.csv' % output_folder, 'w') as outfile:
    outfile.write('sampling_prop,' + ','.join([str(i) for i in range(observe_time)]) + '\n')
    outfile.write('\n'.join([','.join(['%f' % t] + ['%f' % x for x in sp_detected_import_times_95CI[i][1]])
                             for i, t in enumerate(sampling_props)]))