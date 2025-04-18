# Atlantis Model Evaluation Tools

## Overview

This code provides a collection of R functions to evaluate and visualize outputs from Atlantis GOA ecosystem models. It helps assess model stability, identify potential issues, and summarize model behavior to support calibration and validation efforts.

The tools analyze key aspects of model performance including:
- Relative biomass over time
- Numbers at age
- Weight at age
- Biomass trends
- Model persistence (extinction checks)


## Usage

The master function has several parameters to customize the analysis:

```r
evaluate_atlantis_run(
  run_no = 1961,              # Run number
  trend_yrs = 10,             # Years for trend analysis
  plot_biomass = TRUE,        # Generate biomass plots
  plot_naa = TRUE,            # Generate numbers-at-age plots
  plot_waa = TRUE,            # Generate weight-at-age plots
  plot_biomass_trends = FALSE, # Generate biomass trend plots
  generate_report = TRUE,     # Generate summary report
  output_base_path = "../Parametrization/output_files/data/",
  data_path = "data/",
  bgm_file = "data/GOA_WGS84_V4_final.bgm",
  plot_output_dir = "plots/"
)
```

## Output Files

The code produces the following outputs in the run-specific directory:

### Plots
- `relbiom_verts_XXXX.png`: Relative biomass for vertebrate species
- `relbiom_pools_XXXX.png`: Relative biomass for biomass pools
- `relNAA_XXXX.png`: Relative numbers at age
- `relWAA_XXXX.png`: Relative weight at age
- `slope_naa_XXXX.png`: Trends in numbers at age
- `slope_waa_XXXX.png`: Trends in weight at age
- `biomass_trends_vertebrates_XXXX.png`: Biomass trends for vertebrates
- `biomass_trends_biomass_pools_XXXX.png`: Biomass trends for biomass pools
- `biomass_trends_all_XXXX.png`: Biomass trends for all groups

### Report
- `model_summary_XXXX.md`: Summary report with key findings

## Function Documentation

### Master Function

- `evaluate_atlantis_run()`: Main interface to run all analyses

### Analysis Functions

- `plot_relbiom()`: Plot relative biomass over time
- `analyze_biomass_trends()`: Analyze trends in biomass
- `extract_plot_naa()`: Extract and plot numbers at age
- `extract_plot_waa()`: Extract and plot weight at age
- `generate_summary_report()`: Create a summary report

### Helper Functions

- `setNA()`: Set values in boundary boxes to NA
- `collapse_array()`: Sum over depth layers in array

## Interpreting Results

### Relative Biomass Plots
These plots show the ratio of biomass to initial biomass over time:
- Blue line (1.0): Unchanged from initial conditions
- Green lines (0.5, 2.0): Moderate changes (±50%)
- Orange lines (0.25, 4.0): Large changes (±75%)
- Red lines (0.125, 8.0): Extreme changes (±87.5%)

### Trend Analysis
The trend plots show the annual percentage change over the specified period:
- Blue points: Decreasing trends
- Red points: Increasing trends
- Point size: Magnitude of change

### Summary Report
The report highlights potential issues:
- Extinction: Groups with biomass <1% of initial
- Extreme changes: Groups with biomass >8x or <0.125x initial
- Strong trends: Groups with >5% annual change
