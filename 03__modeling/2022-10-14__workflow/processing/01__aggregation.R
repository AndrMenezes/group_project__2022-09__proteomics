#' 2022-11-15: Aggregation from peptides to protein.
suppressMessages(library(QFeatures))

path_data <- "./01__database/processed_data"
path_res <- "03__modeling/2022-10-14__workflow/processing/processed_data"

# Import data -------------------------------------------------------------
fts <- readRDS(file.path(path_data, "QFeature_obj.rds"))

# Aggregate at protein level ----------------------------------------------
fts <- aggregateFeatures(object = fts, i = "peptides", fcol = "protein",
                         name = "proteins_median", fun = colMedians,
                         na.rm = TRUE)
fts <- aggregateFeatures(object = fts, i = "peptides", fcol = "protein",
                         name = "proteins_sum", fun = colSums,
                         na.rm = TRUE)
fts <- aggregateFeatures(object = fts, i = "peptides", fcol = "protein",
                         name = "proteins_mean", fun = colMeans,
                         na.rm = TRUE)

fts[["proteins_sum"]] <- zeroIsNA(fts[["proteins_sum"]])

# Saving the data ---------------------------------------------------------
saveRDS(object = fts, file = file.path(path_res, "fts_processsed.rds"))
