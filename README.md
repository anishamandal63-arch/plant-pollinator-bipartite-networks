# plant-pollinator-bipartite-networks
This repository contains data and R scripts used to analyze plant–pollinator
interaction networks in primary forest (PF) and invasive-dominated habitats (IDH)
from the Eastern Himalaya.

## Contents
### Data
- `pp_data.csv`  
  Raw interaction data used for all analyses.
  Columns include:
  - Forest type FT (PF / IDH)
  - Transect
  - Plant_species
  - Pollinator_species
  - Season (Pre-monsoon / Post-monsoon)
  - PN (Interaction_frequency)
  - Activity (pollinator's activity- Nectaring/Sitting/Basking)
### Code
- `analysis.R`  
  R script used to calculate network indices like nestedness and modularity, diversity indices, NMDS, PERMANOVA and to generate figures.
## Software
- R version 4.5.2
- Packages: vegan, bipartite, igraph, tidyverse
  
This repository is provided for peer review and reproducibility.
