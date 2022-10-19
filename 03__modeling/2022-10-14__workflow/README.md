# Analysis workflow

## Data import

- `01__database/00__create_QFeatures_object.R`: import evidence.txt file into
R and created and `QFeature` object.

## Pre-processing

- `processing/00__imputation.R`: performs imputation using knn algoritm
implemented in `impute` R package.
- `processing/01__normalization.R`: performs normalization across samples
using the function `normalizeCyclicLoess` from `limma` R package.

## Exploratory analysis

- `2022-10-18__missing_values_inspection/missing_values_plots.R`: missing
values visualizations.
- `2022-10-19__visualizing_pre_processed_data`: normalization and comparisons
between log intensities abundance from original paper.
