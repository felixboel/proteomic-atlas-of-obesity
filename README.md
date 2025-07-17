# Multi-Tissue Proteomic Atlas of Obesity Regression in Mice

**Boel et al.**


## Overview

This repository contains the code for the interactive Shiny webtool accompanying the study “Multi-tissue proteomic atlas of obesity regression in mice” by Boel et al. The tool enables users to explore system-wide proteomic changes across 15 mouse tissues during obesity and its regression, integrating Reactome pathway data with differential expression metrics. 

You can access the hosted application here: https://felixboel.shinyapps.io/proteomic-atlas-of-obesity/ 


## Repository Contents

- **`app_launcher.R`**
  - **Function:** Launches the Shiny application locally.
  - **Description:** Installs missing R package dependencies, sources the main application logic, and starts the app.

- **`sub_folder/proteomics_explorer.R`**
  - **Function:** Function: Defines the UI and server logic for the Shiny app.
  - **Description:** Implements data loading, tree-based pathway navigation, interactive bubble plots, and pathway selection. Allows users to visualize differential protein regulation across tissues and timepoints.

- **`sub_folder/*.txt (data files will be uploaded upon publication)`**
  - **UniProt2Reactome_All_Levels.txt:** Maps proteins to Reactome pathways.
  - **ReactomePathways.txt:** Pathway metadata.
  - **ReactomePathwaysRelation.txt:** Defines the hierarchical structure of Reactome pathways.
  - **combined_all.txt:** Differential expression data across tissues and timepoints.
  
## Webtool Features
  - **Pathway-based exploration:** Navigate Reactome’s pathway hierarchy via an interactive tree.
  - **Tissue & timepoint context:** View protein-level changes across 4 timepoints (Ob, STR, MTR, LTR) in 15 tissues.
  - **Bubble plots:** Visualize the count and direction (Up/Down) of regulated proteins for any selected pathway..
  - **Custom filtering:** Set cutoffs for adjusted p-values and log fold-changes to refine the analysis..


## Installation & Running Locally

Ensure the following R packages are available:
```r
install.packages(c("shiny", "shinymanager", "shinyTree", "ggplot2", "plotly", "dplyr", "tidyr", "readr"))
```
To launch the app:
```r
source("app.R")
```


## Citation

If you use this webtool or underlying datasets in your research, please cite our publication:
**Boel, et al. Multi-tissue proteomic atlas of obesity regression in mice**
