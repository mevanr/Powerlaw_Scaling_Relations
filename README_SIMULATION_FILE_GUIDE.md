# Simulation Analysis File Guide

by Mevan Rajakaruna
1) PraxiModus Institute for Advanced Studies, Canada
2) University of Toronto, Mississauga, Canada  

This guide describes the folder and file structure used for the simulation analyses accompanying the manuscript on power-law scaling estimation under measurement error and habitat hierarchy.

The project contains two main simulation regimes:

1. **Single-habitat simulations**: homogeneous scaling data with measurement error.
2. **Multiple-habitat simulations**: hierarchical habitat-structured scaling data with habitat-specific intercepts and/or slopes.

The guide is intended to help readers, reviewers, and collaborators understand which files contain the raw simulation outputs, cleaned results, summary tables, figures, and analysis scripts.

---

## Top-level project structure

```text
Github_23June2026/
├── MetaAnalysis_multiple_habitats/
├── MetaAnalysis_single_habitat/
├── all_results_multiple_habitats.csv
├── all_results_single_habitats.csv
├── Github_Powerlaw_biomass_scaling_manuscript_23June2026.zip
├── Multiple_habitat_analytics.py
├── multiple_habitat_simu_and_estimation_code.R
├── multiple_habitat_simulation_and_estimations.zip
├── proposed_bayesian_model_estimation_on_multiple_habitats.R
├── Single_habitat_analytics.py
├── single_habitat_simu_and_estimation_code.R
└── single_habitat_simulation_and_estimations.zip
```

### Main files

| File or folder | Description |
|---|---|
| `single_habitat_simu_and_estimation_code.R` | R script used to generate single-habitat simulated datasets and apply the estimator set used in the manuscript. |
| `multiple_habitat_simu_and_estimation_code.R` | R script used to generate multiple-habitat/hierarchical simulated datasets and apply the estimator set used in the manuscript. |
| `proposed_bayesian_model_estimation_on_multiple_habitats.R` | R script for the Bayesian mixed errors-in-variables model evaluated for the multiple-habitat case. |
| `Single_habitat_analytics.py` | Python post-processing script for single-habitat results, including summary tables and manuscript-style figures. |
| `Multiple_habitat_analytics.py` | Python post-processing script for multiple-habitat results, including summary tables, heatmaps, ranking outputs, and figures. |
| `all_results_single_habitats.csv` | Combined single-habitat simulation results across scenarios and methods. |
| `all_results_multiple_habitats.csv` | Combined multiple-habitat simulation results across scenarios and methods. |
| `single_habitat_simulation_and_estimations.zip` | Compressed archive of the single-habitat simulation output folders. |
| `multiple_habitat_simulation_and_estimations.zip` | Compressed archive of the multiple-habitat simulation output folders. |
| `Github_Powerlaw_biomass_scaling_manuscript_23June2026.zip` | Full project archive prepared for GitHub submission, including scripts, results, and analysis folders. |

---

## Single-habitat simulation output folder

```text
single_habitat_simulation_and_estimations/
├── habitat_visualization_single_output_n36_mex20_mey40/
├── habitat_visualization_single_output_n36_mex30_mey30/
├── habitat_visualization_single_output_n36_mex40_mey20/
├── habitat_visualization_single_output_n72_mex20_mey40/
├── habitat_visualization_single_output_n72_mex30_mey30/
└── habitat_visualization_single_output_n72_mex40_mey20/
```

Each folder corresponds to one single-habitat simulation scenario.

### Folder-name convention

```text
habitat_visualization_single_output_n{sample_size}_mex{measurement_error_x}_mey{measurement_error_y}
```

For example:

```text
habitat_visualization_single_output_n36_mex20_mey40
```

means:

- `n36`: sample size = 36
- `mex20`: measurement error in the predictor/log-body-size variable = 20%
- `mey40`: measurement error in the response/log-biomass variable = 40%

### Files inside each single-habitat output folder

```text
all_results.csv
Figure1_boxplots.png
Figure2_bias.png
Figure3_rmse_heatmap.png
Figure4_ranking.png
intermediate_results.csv
method_ranking.csv
performance_metrics.csv
scenario_summary.csv
```

| File | Description |
|---|---|
| `all_results.csv` | Full estimator-level results for the scenario. This is the main output file for downstream meta-analysis. |
| `intermediate_results.csv` | Intermediate estimates saved during the simulation and model-fitting workflow. Useful for debugging or method-level inspection. |
| `performance_metrics.csv` | Summary performance metrics for each estimator, including bias and RMSE. |
| `method_ranking.csv` | Ranking of estimators according to performance criteria used in the analysis. |
| `scenario_summary.csv` | Scenario-level summary describing the simulation conditions and aggregate performance. |
| `Figure1_boxplots.png` | Boxplots comparing estimated scaling exponents across methods. |
| `Figure2_bias.png` | Bias comparison across methods. |
| `Figure3_rmse_heatmap.png` | RMSE heatmap across estimators and scenario settings. |
| `Figure4_ranking.png` | Visual ranking of estimator performance for the scenario. |

---

## Multiple-habitat simulation output folder

```text
multiple_habitat_simulation_and_estimations/
├── habitat_visualization_output_R_20_40_n36_6/
├── habitat_visualization_output_R_20_40_n36_12/
├── habitat_visualization_output_R_20_40_n72_6/
├── habitat_visualization_output_R_20_40_n72_12/
├── habitat_visualization_output_R_30_30_n36_6/
├── habitat_visualization_output_R_30_30_n36_12/
├── habitat_visualization_output_R_30_30_n72_6/
├── habitat_visualization_output_R_30_30_n72_12/
├── habitat_visualization_output_R_40_20_n36_6/
├── habitat_visualization_output_R_40_20_n36_12/
├── habitat_visualization_output_R_40_20_n72_6/
└── habitat_visualization_output_R_40_20_n72_12/
```

Each folder corresponds to one multiple-habitat simulation scenario.

### Folder-name convention

```text
habitat_visualization_output_R_{measurement_error_x}_{measurement_error_y}_n{sample_size}_{number_of_habitats}
```

For example:

```text
habitat_visualization_output_R_20_40_n36_6
```

means:

- `R_20_40`: predictor and response measurement-error setting = 20% and 40%
- `n36`: total or scenario-specific sample size = 36
- `6`: number of habitats = 6

Similarly:

```text
habitat_visualization_output_R_40_20_n72_12
```

means measurement-error setting 40%/20%, sample size 72, and 12 habitats.

### Files inside each multiple-habitat output folder

```text
complete_results.rds
Figure1_boxplots.png
Figure2_bias.png
Figure3_rmse_heatmap.png
Figure4_ranking.png
method_ranking.csv
results_clean.csv
results_raw.csv
session_info.txt
summary_statistics.csv
working_methods.csv
```

| File | Description |
|---|---|
| `complete_results.rds` | Complete R object containing scenario outputs, fitted estimates, and related objects saved from the R workflow. |
| `results_raw.csv` | Raw estimator outputs before cleaning or filtering. |
| `results_clean.csv` | Cleaned estimator-level results used for summary analysis. |
| `summary_statistics.csv` | Scenario-level summary statistics, including estimator performance metrics. |
| `method_ranking.csv` | Ranking of estimators under the scenario. |
| `working_methods.csv` | List of methods that successfully ran for the scenario. Useful for identifying convergence or model-fitting failures. |
| `session_info.txt` | R session information, package versions, and computational environment details for reproducibility. |
| `Figure1_boxplots.png` | Boxplots comparing estimated scaling exponents across methods. |
| `Figure2_bias.png` | Bias comparison across estimators. |
| `Figure3_rmse_heatmap.png` | RMSE heatmap across methods and settings. |
| `Figure4_ranking.png` | Visual ranking of estimator performance for the scenario. |

---

## Meta-analysis folders

The meta-analysis folders contain post-processed results aggregated across all simulation scenarios.

---

## `MetaAnalysis_single_habitat/`

This folder contains summary outputs for the single-habitat simulation regime. It is generated from the scenario-level single-habitat result folders and/or from `all_results_single_habitats.csv`.

Typical contents include estimator summaries, ranking tables, heatmaps, and figure files used to support the manuscript and supplementary material.

Use this folder when you need the final aggregated single-habitat results rather than the results from one individual scenario.

---

## `MetaAnalysis_multiple_habitats/`

```text
MetaAnalysis_multiple_habitats/
├── best_method_by_scenario.csv
├── bias_heatmap_exact_error_ratios_table.csv
├── bias_heatmap_exact_slopes_table.csv
├── decision_tree_best_method.png
├── decision_tree_rules.txt
├── heatmap_bias_rmse_exact_error_ratios.png
├── heatmap_bias_rmse_exact_slopes.png
├── hierarchical_heatmap_bias.png
├── hierarchical_heatmap_bias_table.csv
├── hierarchical_heatmap_rmse.png
├── hierarchical_heatmap_rmse_table.csv
├── long_metrics_rowwise.csv
├── overall_bias_rmse_barplots.png
├── overall_method_summary.csv
├── rmse_heatmap_exact_error_ratios_table.csv
├── rmse_heatmap_exact_slopes_table.csv
├── rmse_vs_habitats_and_sample_size.png
├── scenario_summary_metrics.csv
├── summary_by_error_ratio.csv
├── summary_by_slope.csv
└── summary_tables.csv
```

### Key multiple-habitat meta-analysis files

| File | Description |
|---|---|
| `overall_method_summary.csv` | Overall method-level performance summary across all multiple-habitat scenarios. |
| `scenario_summary_metrics.csv` | Scenario-level summary metrics used to compare methods across sample size, habitat number, and error structure. |
| `long_metrics_rowwise.csv` | Long-format metrics table suitable for plotting, filtering, or additional statistical analysis. |
| `summary_tables.csv` | Consolidated summary tables used for manuscript or supplementary reporting. |
| `best_method_by_scenario.csv` | Identifies the best-performing method within each scenario according to the selected performance criterion. |
| `summary_by_error_ratio.csv` | Method performance summarized by measurement-error ratio. |
| `summary_by_slope.csv` | Method performance summarized by true scaling exponent or slope scenario. |
| `hierarchical_heatmap_bias_table.csv` | Table underlying the hierarchical bias heatmap. |
| `hierarchical_heatmap_rmse_table.csv` | Table underlying the hierarchical RMSE heatmap. |
| `bias_heatmap_exact_error_ratios_table.csv` | Bias heatmap source table organized by exact error-ratio settings. |
| `rmse_heatmap_exact_error_ratios_table.csv` | RMSE heatmap source table organized by exact error-ratio settings. |
| `bias_heatmap_exact_slopes_table.csv` | Bias heatmap source table organized by exact slope settings. |
| `rmse_heatmap_exact_slopes_table.csv` | RMSE heatmap source table organized by exact slope settings. |
| `decision_tree_rules.txt` | Plain-text decision rules summarizing which methods perform best under different simulation conditions. |
| `decision_tree_best_method.png` | Decision-tree visualization of best-method selection across scenarios. |
| `heatmap_bias_rmse_exact_error_ratios.png` | Combined bias/RMSE heatmap by measurement-error ratio. |
| `heatmap_bias_rmse_exact_slopes.png` | Combined bias/RMSE heatmap by slope setting. |
| `hierarchical_heatmap_bias.png` | Bias heatmap for hierarchical/multiple-habitat settings. |
| `hierarchical_heatmap_rmse.png` | RMSE heatmap for hierarchical/multiple-habitat settings. |
| `overall_bias_rmse_barplots.png` | Overall bias and RMSE barplots across methods. |
| `rmse_vs_habitats_and_sample_size.png` | Figure showing how RMSE changes with number of habitats and sample size. |

---

## Recommended order for reproducing the analysis

### 1. Run single-habitat simulations

Run:

```text
single_habitat_simu_and_estimation_code.R
```

This generates the scenario-specific folders under:

```text
single_habitat_simulation_and_estimations/
```

The main scenario-level output is:

```text
all_results.csv
```

inside each scenario folder.

### 2. Run single-habitat post-processing

Run:

```text
Single_habitat_analytics.py
```

This combines scenario-level results, generates summary tables, and creates the single-habitat meta-analysis outputs.

### 3. Run multiple-habitat simulations

Run:

```text
multiple_habitat_simu_and_estimation_code.R
```

This generates the scenario-specific folders under:

```text
multiple_habitat_simulation_and_estimations/
```

The main scenario-level outputs are:

```text
results_raw.csv
results_clean.csv
summary_statistics.csv
complete_results.rds
```

### 4. Run multiple-habitat post-processing

Run:

```text
Multiple_habitat_analytics.py
```

This combines the multiple-habitat results, generates heatmaps, ranking outputs, decision-tree summaries, and manuscript-ready summary tables.

### 5. Run Bayesian mixed EIV analysis, if needed

Run:

```text
proposed_bayesian_model_estimation_on_multiple_habitats.R
```

This script evaluates the Bayesian mixed errors-in-variables model for the multiple-habitat/hierarchical case. Because Bayesian models may require longer computation time and diagnostic checking, these results should be interpreted together with convergence diagnostics such as R-hat, effective sample size, and divergent transitions.

---

## Main result files for manuscript reporting

For manuscript-level reporting, the most useful files are:

| Purpose | File |
|---|---|
| Combined single-habitat results | `all_results_single_habitats.csv` |
| Combined multiple-habitat results | `all_results_multiple_habitats.csv` |
| Single-habitat final summaries | `MetaAnalysis_single_habitat/` |
| Multiple-habitat final summaries | `MetaAnalysis_multiple_habitats/overall_method_summary.csv` |
| Multiple-habitat scenario summaries | `MetaAnalysis_multiple_habitats/scenario_summary_metrics.csv` |
| Multiple-habitat best methods by scenario | `MetaAnalysis_multiple_habitats/best_method_by_scenario.csv` |
| Multiple-habitat heatmap source tables | `MetaAnalysis_multiple_habitats/*heatmap*_table.csv` |
| Multiple-habitat manuscript figures | `MetaAnalysis_multiple_habitats/*.png` |

---

## Notes on interpretation

- The **single-habitat regime** evaluates estimator behavior under homogeneous power-law scaling with measurement error.
- The **multiple-habitat regime** evaluates estimator behavior when the data contain hierarchical habitat structure, including habitat-level variation in intercepts and/or slopes.
- Classical single-level estimators and global errors-in-variables corrections may perform differently across these two regimes.
- The scenario folders retain the full analysis trail from raw simulation output to cleaned results, rankings, summaries, and figures.
- The meta-analysis folders are intended for final manuscript reporting and should be treated as the primary location for aggregated results.

---

## Suggested GitHub citation note

If using this repository, please cite the associated manuscript:

> Rajakaruna M., Talmy D., and Rajakaruna H., When global errors-in-variables corrections fail: hierarchical structure induces a regime shift in power-law scaling estimation. Manuscript submitted to *Environmental and Ecological Statistics*.


---

## Reproducibility notes

- R session information is saved in each multiple-habitat scenario folder as `session_info.txt`.
- Scenario-level results are saved separately to avoid overwriting outputs from different sample-size, measurement-error, and habitat-number conditions.
- Compressed `.zip` archives are included for convenient transfer and GitHub upload.
- The `.rds` files preserve complete R objects and are useful for reloading results without rerunning the full simulation.
