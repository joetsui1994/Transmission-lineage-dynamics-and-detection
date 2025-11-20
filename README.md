# Transmission lineage dynamics and the detection of viral importation in emerging epidemics

[Author list has been removed for double-blind review]

---

This repository contains data and scripts used to generate results
presented in [citation]. Please note that some scripts may need small adjustments to run depending on local setup.

## Abstract

_The accurate inference of pathogen movements between locations during an epidemic is crucial for measuring infectious disease spread and for informing effective control strategies. Phylogeographic methods can reconstruct historical patterns of disease dissemination by combining the evolutionary history of sampled pathogen genomes with geographic information. Despite a substantial expansion of pathogen genomics during and since the COVID-19 pandemic, only a small fraction of infections are typically sequenced, leading to an underestimation of the true intensity of viral importation. Here, we seek to understand the sampling processes underlying this underestimation. We show that the coupling of viral importation and local transmission dynamics can result in local transmission lineages with different size distributions, influencing the probability that individual viral importation events will be detected. Using analytical and simulation approaches, we demonstrate that both the proportion of importation events detected and the temporal patterns of inferred importation are highly sensitive to importation dynamics and local transmission parameters. Our findings highlight the importance of interpreting phylogeographic estimates in the context of outbreak conditions, particularly when comparing viral movements across time and among different epidemic settings. These insights are critical for improving the reliability of genomic epidemiology approaches in the design of public health responses._

## Repository usage and structure

The structure of this repository is shown below:

```
├── analyses
│   ├── COVID19_studies_plot.ipynb
│   ├── lineage_detection_prob_heatmap.ipynb
│   ├── deterministic_model
│   │   ├── cumulative_lineage_size_distribution_plot.ipynb
│   │   ├── inferred_importation_rates_plot.ipynb
│   │   ├── lineage_detection_probs_plot.ipynb
│   │   └── run_simulation.py
│   └── stochastic_abm
│       ├── inferred_importation_rates_plot.ipynb
│       ├── lineage_detection_probs_plot.ipynb
│       ├── run_sampling.py
│       ├── run_simulation.py
│       ├── sampling_config.yaml
│       └── simulation_config.yaml
├── data
│   ├── COVID-19_studies_extracted_data.tsv
│   └── world-administrative-boundaries
│       ├── world-administrative-boundaries.cpg
│       ├── world-administrative-boundaries.dbf
│       ├── world-administrative-boundaries.prj
│       ├── world-administrative-boundaries.shp
│       └── world-administrative-boundaries.shx
├── LICENSE.txt
├── README.md
└── gitignore.txt
```

## Stochastic agent-based model

### Running the model

Model parameters are specified through a YAML configuration file. To run a simulation, simply execute the `run_simulation.py` script located in the `analyses/stochastic_abm/` directory, providing the path to your configuration file as an argument (`--config`). For example:

```bash
python run_simulation.py --config config.yaml
```

This will run the simulation and save the results as a compressed CSV file to the output directory specified in your configuration file (default is `./results`).

#### Configuration file

The configuration file should be in YAML format and include the following sections:

```yaml
simulation:
  seed: <random_seed>
  num_agents: <number_of_agents>
  num_steps: <number_of_time_steps>
  num_experiments: <number_of_experiments>
  output_dir: <output_directory>
parameters:
  num_patient_zero: <number_of_initial_infected_agents>
  beta: <transmission_probability_parameters>
  kappa: <contact_rate_parameters>
  gamma: <recovery_rate_parameters>
  importation: <importation_rate_parameters>
```

##### General settings
- `<random_seed>`: Seed for random number generation to ensure reproducibility (integer).
- `<number_of_agents>`: Total population size (or number of agents) (integer).
- `<number_of_time_steps>`: Duration of the simulation in discrete time steps (integer).
- `<number_of_experiments>`: Number of independent simulation runs (integer).
- `<number_of_initial_infected_agents>`: Number of initially infected agents at the start of the simulation (integer).
- `<output_directory>`: Directory where simulation results will be saved.

##### Epidemiological parameters
- `<transmission_probability_parameters>`: Parameters defining the probability of infection per contact ($$\beta$$).
- `<contact_rate_parameters>`: Parameters defining the average number of contacts per agent per time step ($$\kappa$$).
- `<recovery_rate_parameters>`: Parameters defining the recovery rate per time step ($$\gamma$$).
- `<importation_rate_parameters>`: Parameters defining the rate at which new infections are imported into the population per day.

Each epidemiological parameter supports the following types, allowing for dynamic changes over time:
| Type | Description | Required Keys in `params` |
| :--- | :--- | :--- |
| **`constant`** | Value remains fixed throughout the simulation. | `value`: The constant value (float/int). |
| **`sigmoid`** | Transitions from an initial to a final value following a sigmoid trajectory (e.g., following the implementation/lifting of travel restrictions). | `a`: Time-shift (the greater it is, the later the growth/decay starts).<br>`b`: Growth rate (the greater it is, the sharper the change).<br>`c`: Initial value.<br>`d`: Final value. |
| **`exponential`** | Grows or decays exponentially over time. | `initial`: Value at $t=0$.<br>`growth`: Rate coefficient ($>0$ for growth, $<0$ for decay). |

### Simulating a phylogeographic importation analysis (with subsampling)

Once the stochastic ABM simulation is complete and you have located the output CSV file (e.g., `results/T10.csv.gz`), you can run the sampling analysis to simulate the process of phylogeographic importation analysis with subsampling.

To run the sampling analysis, execute the `run_sampling.py` script located in the `analyses/stochastic_abm/` directory, providing the path to your sampling configuration file:

```bash
python run_sampling.py --config sampling_config.yaml
```

This will run the sampling analysis and save the results to the output directory specified in your configuration file (default is `./results`). Specifically, it generates two key output files:

1. `detection_probs.median_95CI.csv`: Proportion (median, 95% CI) of total importation events detected at different sampling proportions beteween 0 and 1.
2. `detection_numbers_by_time.csv`: Inferred importation intensity over time at specified sampling proportions, assuming that the timing of importation of each detected lineage can be accurately recovered.

### Configuration file

The configuration file should be in YAML format and include the following sections:

```yaml
input_file: <path_to_simulation_output>
output_dir: <output_directory>
num_draws: <number_of_random_draws>
random_seed: <random_seed>
sampling_grid:
  split_threshold: <threshold_to_switch_between_fine_and_coarse_grid>
  num_divisions_lower: <number_of_divisions_below_threshold>
  num_divisions_upper: <number_of_divisions_above_threshold>
temporal_analysis:
  enabled: <true_or_false>
  target_proportions: [<list_of_sampling_proportions>]
```

- `<path_to_simulation_output>`: Path to the input simulation CSV file (e.g., `results/T10.csv.gz`).
- `<output_directory>`: Directory where sampling results will be saved.
- `<number_of_random_draws>`: Number of random draws per sampling fraction (integer).
- `<random_seed>`: Seed for random number generation to ensure reproducibility (integer).
- `<threshold_to_switch_between_fine_and_coarse_grid>`: Threshold to switch between fine and coarse sampling grid (float between 0 and 1).
- `<number_of_divisions_below_threshold>`: Number of divisions in the sampling grid below the threshold (integer).
- `<number_of_divisions_above_threshold>`: Number of divisions in the sampling grid above the threshold (integer).
- `<true_or_false>`: Whether to enable temporal analysis, i.e. to calculate inferred importation intensities over time (boolean).
- `<list_of_sampling_proportions>`: List of sampling proportions (between 0 and 1) for which to calculate inferred importation intensities over time.

### Visualisation

Once the importation analysis is complete, you can either analyse the output CSV files directly, or use the provided Jupyter notebooks in the `analyses/stochastic_abm/` directory to recreate the figures from the manuscript, specifically:

- `lineage_detection_probs_plot.ipynb`: To visualise the lineage detection probabilities (i.e. proportion of extant lineages sampled) at different sampling fractions.
- `inferred_importation_rates_plot.ipynb`: To visualise the inferred importation rates over time at different sampling fractions.

## Deterministic model with local exponential growth

### Running the model

Model parameters are specified through a YAML configuration file. To run a simulation, simply execute the `run_simulation.py` script located in the `analyses/deterministic_model/` directory, providing the path to your configuration file as an argument (`--config`). For example:

```bash
python run_simulation.py --config config.yaml
```

This will run both the outbreak simulation and importation analysis (with subsampling), and save the results as a compressed CSV file to the output directory specified in your configuration file (default is `./results`).

#### Configuration file

```yaml
simulation:
  random_seed: <random_seed>
  observe_time: <duration_in_days>
  output_dir: <output_directory>
sampling_parameters:
  num_sampling_simulations: <number_of_sampling_runs>
  grid:
    split_threshold: <sampling_fraction_threshold>
    num_divisions_lower: <number_of_divisions_below_threshold>
    num_divisions_upper: <number_of_divisions_above_threshold>
epi_parameters:
  local_growth_rate: <local_growth_rate>
  importation_rate:
    type: <constant_or_exponential>
    value/params: <function_specific_values>
```
- `<random_seed>`: Seed for random number generation to ensure reproducibility (integer).
- `<observe_time>`: Time at which samples are taken (and time horizon for which importation intensities are evaluated) (integer).
- `<output_directory>`: Directory where simulation results will be saved.
- `num_sampling_simulations`: Number of independent sampling runs to perform (integer).
- `<sampling_fraction_threshold>`: Threshold to switch between fine and coarse sampling grid (float between 0 and 1).
- `<number_of_divisions_below_threshold>`: Number of divisions in the sampling grid below the threshold (integer).
- `<number_of_divisions_above_threshold>`: Number of divisions in the sampling grid above the threshold (integer).
- `<local_growth_rate>`: Exponential growth rate of local outbreak (float).
- `importation_rate` supports `constant` (requires `value`) or `exponential` (requires `params.initial` and `params.growth`).

The model will output three CSV files to the specified output directory, specifically:

- `detection_probs.csv`: Median lineage detection probability and 95% CI bounds for each sampling proportion.
- `detection_probs_time.csv`: Median number of importation events inferred at each time point (up to `observe_time`) for each sampling proportion.
- `detection_probs_time.95CI_lw.csv` and `detection_probs_time.95CI_up.csv`: Lower and upper bounds of the 95% confidence interval for the number of importation events inferred at each time point (up to `observe_time`) for each sampling proportion.

### Visualisation

Once the importation analysis is complete, you can either analyse the output CSV files directly, or use the provided Jupyter notebooks in the `analyses/deterministic_model/` directory to recreate the figures from the manuscript, specifically:

- `lineage_detection_probs_plot.ipynb`: To visualise the lineage detection probabilities (i.e. proportion of extant lineages sampled) at different sampling fractions.
- `inferred_importation_rates_plot.ipynb`: To visualise the inferred importation rates (median and 95% CI) over time at different sampling fractions.

---