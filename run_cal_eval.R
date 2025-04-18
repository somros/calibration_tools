#' Run Atlantis Model Evaluation
#' 
#' This script demonstrates how to use the atlantis model evaluation functions
#' to analyze and visualize model outputs.

# Load necessary libraries and source the functions
library(here)
library(tidyverse)
library(tidync)
library(ncdf4)
library(rbgm)

select <- dplyr::select

# Source the functions
source("calibration_eval_functions.R")

# Set the run number to evaluate
run_no <- 1961  # Change this to the run number you want to evaluate

# This will create plots and a summary report
# Only analyze biomass trends and weight-at-age with 15-year trend period
evaluate_atlantis_run(
  run_no = 1961,
  trend_yrs = 10,
  plot_biomass = TRUE,       # Skip standard biomass plots
  plot_naa = TRUE,           # Skip numbers-at-age plots
  plot_waa = TRUE,            # Include weight-at-age plots
  plot_biomass_trends = TRUE, # Include biomass trend analysis
  generate_report = TRUE      # Generate summary report
)
