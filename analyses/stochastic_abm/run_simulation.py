#!/usr/bin/env python

from tqdm import tqdm
import pandas as pd
import numpy as np
import argparse
import yaml
import sys
import os


# ===============================
# Parameters validation
# ===============================
def validate_params(func_name, params, required_keys):
    missing = [key for key in required_keys if key not in params]
    if missing:
        raise ValueError(
            f"configuration error: the '{func_name}' function is missing "
            f"required parameters: {missing}, please check your YAML file."
        )


# ===============================
# Function factory
# ===============================

def build_time_function(name, param_config):
    """
    Returns a function f(t) based on the configuration dictionary.
    """
    # check if 'type' is specified
    if 'type' not in param_config:
        raise ValueError(f"configuration error: '{name}' is missing the 'type' field.")

    func_type = param_config['type']
    
    # check if 'params' exists (unless type is constant and uses 'value')
    if func_type != 'constant' and 'params' not in param_config:
        raise ValueError(f"configuration error: '{name}' (type: {func_type}) requires a 'params' section.")

    if func_type == 'constant':
        if 'value' not in param_config:
             raise ValueError(f"configuration error: '{name}' (type: constant) requires a 'value'.")
        val = param_config['value']
        return lambda t: val
        
    elif func_type == 'sigmoid':
        p = param_config['params']
        # validate specific parameters for sigmoid
        validate_params(name, p, ['a', 'b', 'c', 'd'])
        
        # value constraints (e.g., growth rate must be positive)
        if p['b'] < 0:
             raise ValueError(f"configuration error: In '{name}', parameter 'b' (growth rate) must be positive.")
             
        return lambda t: sigmoid_func(t, p['a'], p['b'], p['c'], p['d'])
        
    elif func_type == 'exponential':
        p = param_config['params']
        validate_params(name, p, ['initial', 'growth'])
        return lambda t: int(p['initial'] * np.exp(p['growth'] * t))
        
    else:
        # catch unknown function types
        raise ValueError(f"configuration error: unknown function type '{func_type}' for '{name}'. "
                         f"supported types are: constant, sigmoid, exponential.")
        

# ===============================
# Utility functions
# ===============================

def init_df(num_agents, num_patient_zero):
    state = np.zeros(num_agents)
    infected = np.random.choice(num_agents, num_patient_zero, replace=False)
    state[infected] = 1

    # initialize the 'infected_by' column with NaN or a placeholder value
    infected_by = np.full(num_agents, -1)  # -1 indicates not infected or patient zero
    imported_time = np.full(num_agents, -1)  # -1 indicates not imported
    infected_time = np.full(num_agents, -1) # -1 indicates not infected or patient zero
    data = {
        'state': state,
        'infected_by': infected_by,
        'imported_time': imported_time,
        'infected_time': infected_time,
    }
    df = pd.DataFrame(data)
    return df

def import_cases(df, current_import_rate, time):
    # add entries to the dataframe for new cases
    state = np.ones(current_import_rate)
    infected_by = np.full(current_import_rate, -1)
    imported_time = np.full(current_import_rate, time)
    infected_time = np.full(current_import_rate, -1)
    data = { 'state': state, 'infected_by': infected_by, 'imported_time': imported_time, 'infected_time': infected_time }
    new_cases = pd.DataFrame(data)
    df = pd.concat([df, new_cases], ignore_index=True)
    return df

def infect(df, time, current_beta, current_kappa):
    num_agents = len(df)
    total_contacts = int(current_kappa * num_agents // 2)  # divide by 2 because each contact involves two agents

    # randomly select pairs of agents
    contact_pairs = np.random.choice(df.index, (total_contacts, 2), replace=True)

    # exclude self-contacts
    valid_contacts = contact_pairs[:, 0] != contact_pairs[:, 1]
    contact_pairs = contact_pairs[valid_contacts]

    # get states of both agents in each pair
    states_1 = df.loc[contact_pairs[:, 0], "state"].values
    states_2 = df.loc[contact_pairs[:, 1], "state"].values

    # identify where one is infected and the other is susceptible
    transmission_mask = (
        ((states_1 == 1) & (states_2 == 0)) |
        ((states_1 == 0) & (states_2 == 1))
    )

    # determine if infection occurs
    transmission_events = transmission_mask & (np.random.uniform(size=len(contact_pairs)) < current_beta)

    # identify newly infected agents
    infected_indices = []
    infectors = []
    # agents where the first is infected and second is susceptible
    mask_1_infected = transmission_events & (states_1 == 1) & (states_2 == 0)
    infected_indices.extend(contact_pairs[mask_1_infected, 1])
    infectors.extend(contact_pairs[mask_1_infected, 0])
    # agents where the second is infected and first is susceptible
    mask_2_infected = transmission_events & (states_1 == 0) & (states_2 == 1)
    infected_indices.extend(contact_pairs[mask_2_infected, 0])
    infectors.extend(contact_pairs[mask_2_infected, 1])

    # remove duplicates
    unique_infected_indices, indices = np.unique(infected_indices, return_index=True)
    unique_infectors = np.array(infectors)[indices]

    # update states and record infectors
    df.loc[unique_infected_indices, "state"] = 1
    df.loc[unique_infected_indices, "infected_by"] = unique_infectors
    df.loc[unique_infected_indices, "infected_time"] = time

def recover(df, current_gamma):
    ps = np.random.uniform(0, 1, len(df[df["state"] == 1]))
    change = np.array(ps <= current_gamma).astype(int)
    df.loc[df["state"] == 1, "state"] = 1 + change

def step(df, time, beta_t, kappa_t, gamma_t, import_t):
    # add new cases to df according to the import rate
    current_import_rate = import_t(time)
    df = import_cases(df, current_import_rate, time)

    # get contact rate at current time
    current_beta = beta_t(time)
    current_kappa = kappa_t(time)
    current_gamma = gamma_t(time)

    # infect contacts probabilisitically
    infect(df, time, current_beta, current_kappa)
    # recover infected agents probabilistically
    recover(df, current_gamma)

    return df

def sigmoid_func(x, a, b, c, d):
    """
    Returns a sigmoid function with parameters a, b, c and d
    a: time-shift (the greater a is, the later the growth starts)
    b: growth rate (the greater b is, the sharper the change)
    c: the initial value
    d: the final value
    """
    # check that b is positive
    if b < 0:
        raise ValueError("b must be positive")
    # check that both c and d are positive
    if c < 0 or d < 0:
        raise ValueError("both c and d must be positive")

    return c - (c - d) / (1 + np.exp(-b * (x - a)))

def find_index_cases(df):
    # find the index cases
    index_cases_indices = df.loc[(df.infected_by == -1) & (df.state != 0)].index.tolist()
    return index_cases_indices

def get_infector_infectees(df):
    infectee_infector = list(enumerate(df.infected_by.values))
    infector_infectees = {}
    for infectee, infector in infectee_infector:
        if infector in infector_infectees:
            infector_infectees[infector].append(infectee)
        else:
            infector_infectees[infector] = [infectee]
    return infector_infectees

def recursive_trace(infector, infector_infectees):
    new_infectors = infector_infectees.get(infector, [])
    infectees = [infector]
    for new_infector in new_infectors:
        infectees += recursive_trace(new_infector, infector_infectees)
    return infectees

def extract_clusters(df):
    # get index cases indices 
    index_cases = find_index_cases(df)
    # get infector_infectees
    infector_infectees = get_infector_infectees(df)
    # start with index cases
    clusters = []
    for index_case in index_cases:
        cluster = recursive_trace(index_case, infector_infectees)
        clusters.append(cluster)
    return clusters

def label_agents_by_cluster(df, clusters):
    # label agents by cluster
    agent_cluster = np.full(len(df), -1)
    for cluster_i, cluster in enumerate(clusters):
        agent_cluster[cluster] = cluster_i
    df["cluster"] = agent_cluster
    # add cluster_imported_time column, which is the importation time of the imported case in the cluster
    cluster_imported_time = np.full(len(df), -1)
    for cluster_i, cluster in enumerate(clusters):
        cluster_imported_time[cluster] = df.loc[cluster].query('imported_time > -1').imported_time.min()
    df["cluster_imported_time"] = cluster_imported_time
    return df


# ===============================
# Main execution
# ===============================

if __name__ == "__main__":
    # parse command line arguments
    parser = argparse.ArgumentParser(description="Run stochastic ABM simulation with viral importation.")
    parser.add_argument('--config', type=str, required=True, help="Path to YAML config file")
    args = parser.parse_args()

    # load config file
    print('loading config file...')
    with open(args.config, 'r') as f:
        config = yaml.safe_load(f)

    # extract parameters from config
    sim_conf = config['simulation']
    seed = sim_conf.get('seed', 42)  # default to 42 if not specified
    num_agents = sim_conf.get('num_agents', 10000)  # default to 10000 if not specified
    num_steps = sim_conf.get('num_steps', 100)  # default to 100 if not specified
    num_experiments = sim_conf.get('num_experiments', 10)  # default to 10 if not specified
    num_patient_zero = sim_conf.get('num_patient_zero', 0)  # default to 0 if not specified
    output_dir = sim_conf.get('output_dir', './results')  # default to ./results if not specified

    # set random seed for reproducibility
    np.random.seed(seed)

    # build time functions from config
    try:
        beta_t = build_time_function('transmission_prob', config['parameters']['beta'])
        kappa_t = build_time_function('contact_rate', config['parameters']['kappa'])
        gamma_t = build_time_function('recovery_rate', config['parameters']['gamma'])
        import_t = build_time_function('importation_rate', config['parameters']['importation'])        
    except ValueError as e:
        print(e)
        sys.exit(1)

    # run simulation
    print('simulating outbreaks...')
    exp_dfs = []
    exp_histories = []
    for exp_i in tqdm(range(num_experiments)):
        # initialize the experiment
        df = init_df(num_agents, num_patient_zero)
        # run
        history = []
        for i in tqdm(range(num_steps + 1), leave=False):
            df = step(df, i, beta_t, kappa_t, gamma_t, import_t)
            # get state counts
            history.append(dict(df["state"].value_counts()))

        # store the history of the experiment
        exp_histories.append(history)
        # store the final state of the experiment
        exp_dfs.append(df.copy())


# ===============================
# Export results
# ===============================

# extract clusters from each experiment
print('extracting clusters...')
exp_clusters = []
for exp_i in tqdm(range(num_experiments)):
    clusters = extract_clusters(exp_dfs[exp_i])
    exp_clusters.append(clusters)

# label agents by cluster (agents labelled -1 are not infected, so ignore them)
print('labelling agents by cluster...')
exp_dfs_clustered = []
for exp_i in tqdm(range(num_experiments)):
    df = label_agents_by_cluster(exp_dfs[exp_i], exp_clusters[exp_i])
    df = df[df.cluster != -1]
    exp_dfs_clustered.append(df)

# convert exp_dfs_clustered to a single dataframe, with columns for experiment number
exp_dfs_clustered_combined = pd.concat(exp_dfs_clustered, keys=range(num_experiments), names=["experiment"]).reset_index(level=0)
# export to csv (compressed), creating output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)
print('exporting results...')
exp_dfs_clustered_combined.to_csv("./%s/T%d.csv.gz" % (output_dir, num_steps), index=False, compression='gzip')