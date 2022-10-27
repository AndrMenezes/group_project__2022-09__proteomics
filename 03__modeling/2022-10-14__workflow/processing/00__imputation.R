#' 2022-10-18: Imputation of missing values using knn from `impute::impute.knn`
#' function.
suppressMessages(library(QFeatures))

path_data <- "./01__database/"
path_res <- "03__modeling/2022-10-14__workflow/processing/processed_data"

# Import data -------------------------------------------------------------
fts <- readRDS(file.path(path_data, "QFeature_obj.rds"))

# Aggregate at protein level ----------------------------------------------
fts <- aggregateFeatures(object = fts, i = "peptides", fcol = "protein",
                         name = "proteins", fun = colMedians, na.rm = TRUE)
# Distribution of the number of peptide use to aggregate
plot(table(rowData(fts[["proteins"]])$.n))

# Filtering proteins with more than 50% of missing values -----------------
x <- assay(fts[["proteins"]])
th <- 0.50
chosen <- rowMeans(is.na(x)) >= 0.5
table(chosen)
prop.table(table(chosen))
fts[["proteins"]] <- fts[["proteins"]][!chosen, ]

# Imputation using k-nearest neighbor averaging ---------------------------
assay(fts[["proteins"]], "log_intensity") <- log2(assay(fts[["proteins"]]))
logx_imputed <- impute::impute.knn(
  data = assay(fts[["proteins"]], "log_intensity"))
assay(fts[["proteins"]], "log_intensity_imputed") <- logx_imputed$data

# Saving the data ---------------------------------------------------------
saveRDS(object = fts, file = file.path(path_res, "fts_processsed.rds"))
