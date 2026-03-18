# Lindera SDMs under Climate Change

This repository contains R scripts used for species distribution models (SDMs) and niche overlap analysis of four Lindera species under climate change scenarios.

The data used in this study are available at:
(https://doi.org/10.6084/m9.figshare.31796251)

## Requirements
R version >= 4.3.0

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
