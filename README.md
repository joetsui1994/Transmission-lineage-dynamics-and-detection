# Transmission lineage dynamics and the detection of viral importation in emerging epidemics

**Joseph L.-H. Tsui**, Prathyush Sambaturu, Rosario Evans Pena, Linus Too, Bernardo Gutierrez, Rhys Inward, **Moritz U. G. Kraemer**, **Louis du Plessis**, **Oliver G. Pybus**

---

This repository contains data and scripts used to generate results
presented in [citation]. Please note that some scripts may need small adjustments to run depending on local setup.

## Abstract

_The accurate inference of pathogen movements between locations during an epidemic is crucial for measuring infectious disease spread and for informing effective control strategies. Phylogeographic methods can reconstruct historical patterns of disease dissemination by combining the evolutionary history of sampled pathogen genomes with geographic information. Despite a substantial expansion of pathogen genomics during and since the COVID-19 pandemic, only a small fraction of infections are typically sequenced, leading to an underestimation of the true intensity of viral importation. Here, we seek to understand the sampling processes underlying this underestimation. We show that the coupling of viral importation and local transmission dynamics can result in local transmission lineages with different size distributions, influencing the probability that genomic surveillance will detect individual viral importation events. Using analytical and stochastic simulation approaches, we demonstrate that both the proportion of importation events detected and the temporal patterns of inferred importation are highly sensitive to importation dynamics and local transmission parameters. Our findings highlight the importance of interpreting phylogeographic estimates in the context of outbreak conditions, particularly when comparing viral movements across time and among different epidemic settings. These insights are critical for improving the reliability of genomic epidemiology approaches in the design of public health responses._

## Repository usage and structure

The structure of this repository is shown below:

```
├── LICENSE.txt
├── README.md
├── analyses
│   ├── COVID19_studies_plot.ipynb
│   ├── deterministic_model
│   │   ├── cumulative_lineage_size_distribution_plot.ipynb
│   │   ├── inferred_importation_rates_plot.ipynb
│   │   ├── lineage_detection_probs_plot.ipynb
│   │   └── run_simulation.py
│   ├── lineage_detection_prob_heatmap.ipynb
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
└── gitignore.txt
```

---