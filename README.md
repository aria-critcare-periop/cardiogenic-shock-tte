# Cardiogenic shock phenotypes TTE

This repository contains the analysis code for the statistical analysis plan *Cardiogenic Shock Phenotypes Validation and Target Trial Emulation of Temporary Mechanical Circulatory Support Strategies*, version 1.2 (doi:10.5281/zenodo.22842184).

## Running the analysis

Run the full analysis with:

```bash
Rscript run_all.R
```

The code was developed and tested with R 4.3.1 using the following packages: `readxl`, `mgcv`, `mice`, `splines`, `cmprsk`, `ConsensusClusterPlus`, `cluster`, and `Rtsne`.

The path to the source workbook is defined in `R/00_settings.R`. Patient-level data are not included in this repository.

Tables, figures, and audit files are written to `results/`. Files containing patient-level information are written to `private/`.

Bootstrap replicates are saved individually in `private/bootstrap/`, allowing an interrupted analysis to resume from the last completed replicate. This folder should be deleted before rerunning the analysis after any change to the source data or the code.

## Repository structure

| File | Content |
|---|---|
| `R/00_settings.R` | Settings, adjustment set, labels, colors |
| `R/01_data.R` | Source workbook, strategies, time zero, outcomes |
| `R/02_phenotypes.R` | Consensus clustering and phenotype assignment |
| `R/03_imputation.R` | Multiple imputation |
| `R/04_models.R` | Outcome and hazard models, g-computation, propensity model, inference |
| `R/05_bootstrap.R` | Point estimates across imputations, center-stratified bootstrap |
| `R/06_phenotype_description.R` | Phenotype profiles, SCAI stage, prognostic models |
| `R/07_tables.R` | Tables |
| `R/08_figures.R` | Figures |
| `R/09_checks.R` | Numerical checks and audit files |

## Use of artificial intelligence tools

The analysis code was initially written without the use of artificial intelligence tools. ChatGPT Sol (OpenAI) and Claude Opus (Anthropic) were used to assist with code review, identify and correct errors, and improve the efficiency and clarity of the implementation. All suggested changes were reviewed and validated by the authors.
