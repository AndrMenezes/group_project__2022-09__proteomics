#' 2022-10-18: Imputation of missing values using knn from `impute::impute.knn`
#' function.
suppressMessages(library(QFeatures))

path_data <- "./01__database/processed_data"
path_res <- "03__modeling/2022-10-14__workflow/processing/processed_data"

# Import data -------------------------------------------------------------
fts <- readRDS(file.path(path_data, "QFeature_obj.rds"))

# Aggregate at protein level ----------------------------------------------
fts <- aggregateFeatures(object = fts, i = "peptides", fcol = "protein",
                         name = "proteins", fun = colMedians, na.rm = TRUE)
# Distribution of the number of peptide use to aggregate
plot(table(rowData(fts[["proteins"]])$.n))
aggcounts(fts[["proteins"]])
to_remove <- rowData(fts[["proteins"]])$.n <= 3
table(to_remove)
se_low_peptides <- fts[["proteins"]][to_remove, ]
range(rowMeans(aggcounts(se_low_peptides)))

# Removing proteins that have less than 3 peptides matches
fts[["proteins"]] <- fts[["proteins"]][!to_remove, ]

# Filtering proteins with more than 50% of missing values -----------------
x <- assay(fts[["proteins"]])
th <- 0.50
to_remove <- rowMeans(is.na(x)) >= 0.5
table(to_remove)
prop.table(table(to_remove))
fts[["proteins"]] <- fts[["proteins"]][!to_remove, ]

# Imputation using k-nearest neighbor averaging ---------------------------
assay(fts[["proteins"]], "log2_intensity") <- log2(assay(fts[["proteins"]]))
logx_imputed <- impute::impute.knn(
  data = assay(fts[["proteins"]], "log2_intensity"))
assay(fts[["proteins"]], "log2_imputed") <- logx_imputed$data

# Saving the data ---------------------------------------------------------
saveRDS(object = fts, file = file.path(path_res, "fts_processsed.rds"))
