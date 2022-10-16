# Analysis workflow

> This folder contains the scripts to reproduce the analysis

The steps are as follows:

1. `processing/01__normalization.R`: performs log intensity normalization.
2. `models/01__proteins_selection.R`: model the relationships between mean and variance of proteins intensity and then decompose the variance of protein abundance into technical and biological.
3. `models/02__proteins_correlation.R`: estimate and visualize using biplot the correlation between the highly variable proteins.
