# Lindera_SDM
Predicting habitat suitability of Korean Lindera as Tertiary relict plants under climate change scenarios.

# Lindera SDM under Climate Change

This repository contains R scripts used for species distribution modeling (SDM) and niche overlap analysis of four Lindera species under climate change scenarios.

## Data
The data used in this study are available at:
https://doi.org/xxxx

## Requirements
R version >= 4.0

Required packages:
biomod2  
raster  
terra  
dplyr  
ggplot2  

## Note
File paths in the scripts should be adjusted to the user’s local environment before execution.

## Scripts

01_data_preprocessing.R  
Preprocessing occurrence data (splitting CSV and preparing species occurrence records).

02_sdm_modeling.R  
Ensemble species distribution modeling using biomod2.

03_analysis_sensitivity_area.R  
Calculation of sensitivity indices, area change, and generation of related figures.

04_niche_overlap_analysis.R  
Niche overlap analysis (Schoener’s D) and visualization.
