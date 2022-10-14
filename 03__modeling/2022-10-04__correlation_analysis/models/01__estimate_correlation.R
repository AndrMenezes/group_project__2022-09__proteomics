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
p <- dim(cor_all)[1L]
L <- lower.tri(cor_all)
rc_index <- which(L, arr.ind = TRUE)
data_cor <- tibble(
  protein__A = dimnames(cor_all)[[1L]][rc_index[, 1L]],
  protein__B = dimnames(cor_all)[[2L]][rc_index[, 2L]],
  correlation = cor_all[L]
)
choose(p, 2) == dim(data_cor)[1L]

# Computing the unbiased estimate and then filtering
data_cor <- data_cor |>
  mutate(
    correlation_unbiased = correlation * (1 + (1 - correlation^2) / 2 / n)) |>
  filter(abs(correlation_unbiased) >= 0.80) |>
  arrange(-abs(correlation_unbiased)) |>
  mutate(group = "All") |>
  select(group, everything())

# Correlation for each group ----------------------------------------------

splited_df <- split(data_wide_protein[, -(1:2)], f = data_wide_protein$group)
list_cor <- lapply(splited_df, function(x) {
  cor_all <- cor(as.matrix(x), method = "pearson")
  L <- lower.tri(cor_all)
  rc_index <- which(L, arr.ind = TRUE)
  data_cor <- tibble(
    protein__A = dimnames(cor_all)[[1L]][rc_index[, 1L]],
    protein__B = dimnames(cor_all)[[2L]][rc_index[, 2L]],
    correlation = cor_all[L])
  data_cor
})
names(list_cor) <- names(splited_df)
data_cor_grouped <- bind_rows(list_cor, .id = "group") |>
  mutate(
    correlation_unbiased = correlation * (1 + (1 - correlation^2) / 2 / n)) |>
  filter(abs(correlation_unbiased) >= 0.80) |>
  arrange(-abs(correlation_unbiased))


# Exporting ---------------------------------------------------------------
data_cor_all <- bind_rows(data_cor, data_cor_grouped)
write.csv(data_cor_all,
          file = file.path(path_res, "estimated_correlation.csv"),
          row.names = FALSE)
