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
│       ├── run_sampling_simulation.py
│       ├── run_simulation.py
│       └── simulate_temporal_sampling.ipynb
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

## Stochastic agent-based simulation

Model parameters are specified through a YAML configuration file. To run a simulation, simply execute the `run_simulation.py` script located in the `analyses/stochastic_abm/` directory, providing the path to your configuration file as an argument (`--config`). For example:

```bash
python run_simulation.py --config config.yaml
```

This will run the simulation and save the results as a compressed CSV file to the output directory specified in your configuration file (default is `./results`).

### Configuration file format

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

#### General settings
- `<random_seed>`: Seed for random number generation to ensure reproducibility (integer).
- `<number_of_agents>`: Total population size (or number of agents) (integer).
- `<number_of_time_steps>`: Duration of the simulation in discrete time steps (integer).
- `<number_of_experiments>`: Number of independent simulation runs (integer).
- `<number_of_initial_infected_agents>`: Number of initially infected agents at the start of the simulation (integer).
- `<output_directory>`: Directory where simulation results will be saved.

#### Epidemiological parameters
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

---