#' 2022-10-04: Compute the correlation between proteins.
library(dplyr)

path_local <- "./03__modeling/2022-10-04__correlation_analysis"
path_data <- file.path(path_local, "/processing/processed_data/")
path_res <- file.path(path_local, "models/results")


data_wide_protein <- readRDS(file.path(path_data, "data_wide_protein.rds"))
n <- nrow(data_wide_protein)

# Correlation considering all groups --------------------------------------
cor_all <- cor(as.matrix(data_wide_protein[, -(1:2)]), method = "pearson")
dim(cor_all)
