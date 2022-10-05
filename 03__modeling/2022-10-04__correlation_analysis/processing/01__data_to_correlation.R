#' Organize the data to compute the correlation between proteins.
library(dplyr)

path_data <- "./01__database/"
path_export <- "./03__modeling/2022-10-04__correlation_analysis/processing/processed_data/"


# Import data -------------------------------------------------------------

data_all_protein <- read.csv(
  file.path(path_data, "processed_data/all_protein__pivotted.csv"))

head(data_all_protein)


# Pivoting the data -------------------------------------------------------

data_wide_protein <- data_all_protein |>
  select(protein__id, replicate, group, value) |>
  tidyr::pivot_wider(names_from = protein__id, values_from = value)
# Checking dimensions
all.equal(dim(data_wide_protein), c(15, 947))
data_wide_protein

saveRDS(object = data_wide_protein,
        file = file.path(path_export, "data_wide_protein.rds"))
