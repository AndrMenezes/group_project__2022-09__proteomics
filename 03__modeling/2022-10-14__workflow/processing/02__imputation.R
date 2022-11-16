#' 2022-10-18: Imputation of missing values using knn from `impute::impute.knn`
#' function.
suppressMessages(library(QFeatures))

path_res <- "03__modeling/2022-10-14__workflow/processing/processed_data"

# Import data -------------------------------------------------------------
fts <- readRDS(file.path(path_res, "fts_processsed.rds"))


# Distribution of the number of peptide use to aggregate ------------------
plot(table(rowData(fts[["proteins_median"]])$.n))
aggcounts(fts[["proteins_median"]])
to_remove <- rowData(fts[["proteins_median"]])$.n <= 3
table(to_remove)
se_low_peptides <- fts[["proteins_median"]][to_remove, ]
range(rowMeans(aggcounts(se_low_peptides)))

# Removing proteins that have less than 3 peptides matches ----------------
fts[["proteins_median"]] <- fts[["proteins_median"]][!to_remove, ]
fts[["proteins_mean"]] <- fts[["proteins_mean"]][!to_remove, ]
fts[["proteins_sum"]] <- fts[["proteins_sum"]][!to_remove, ]

# Filtering proteins with more than 50% of missing values -----------------
x <- assay(fts[["proteins_mean"]])
th <- 0.50
to_remove <- rowMeans(is.na(x)) >= 0.5
table(to_remove)
prop.table(table(to_remove))

fts[["proteins_median"]] <- fts[["proteins_median"]][!to_remove, ]
fts[["proteins_mean"]] <- fts[["proteins_mean"]][!to_remove, ]
fts[["proteins_sum"]] <- fts[["proteins_sum"]][!to_remove, ]

head(assay(fts[["proteins_sum"]]))

# Log2 transformed --------------------------------------------------------
assay(fts[["proteins_median"]], "log2_intensity") <- log2(
  assay(fts[["proteins_median"]]))
assay(fts[["proteins_mean"]], "log2_intensity") <- log2(
  assay(fts[["proteins_mean"]]))
assay(fts[["proteins_sum"]], "log2_intensity") <- log2(
  assay(fts[["proteins_sum"]]))

# Imputation using k-nearest neighbor averaging ---------------------------
assay(fts[["proteins_median"]], "log2_imputed") <- impute::impute.knn(
  data = assay(fts[["proteins_median"]], "log2_intensity"))$data
assay(fts[["proteins_mean"]], "log2_imputed") <- impute::impute.knn(
  data = assay(fts[["proteins_mean"]], "log2_intensity"))$data
assay(fts[["proteins_sum"]], "log2_imputed") <- impute::impute.knn(
  data = assay(fts[["proteins_sum"]], "log2_intensity"))$data

# Saving the data ---------------------------------------------------------
saveRDS(object = fts, file = file.path(path_res, "fts_processsed.rds"))
