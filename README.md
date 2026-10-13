# Vitamin D and Untreated Caries

This repository contains a consolidated R script for analyzing vitamin D intake and untreated caries among adults aged 60 years and older using NHANES 2015–2016 and 2017–March 2020 pre-pandemic data.

## File

`code.R` includes:

- NHANES data retrieval using the nhanesA package.
- Data preparation, participant selection, and variable definitions.
- Survey-weighted descriptive statistics and multivariable logistic regression.
- Restricted cubic splines and threshold effect analysis.
- Correlation analysis and LASSO feature selection.
- Training and evaluation of 15 machine learning models.
- Calibration curves and subgroup analyses.
- SHAP analysis for Logistic Regression.
- Figure and table exports.

## Requirements

R and the packages listed in `load_required_packages()` in `code.R` are required. The script reports any missing packages before starting.

Internet access is required to download NHANES data.

## Usage

From the directory containing `code.R`, run:

```bash
Rscript code.R --run
```

To specify output and download-cache directories:

```bash
Rscript code.R --run --output=results --cache=nhanes_cache
```

Without `--run`, the script does not execute the analysis.

## Outputs

Outputs include processed data, statistical results, model evaluations, SHAP values, figures, tables, and session information.

By default, results are saved in a timestamped folder under `code_results`.
