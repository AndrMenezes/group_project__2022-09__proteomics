#' 2022-10-19: Perform normalization on the log-intensities at protein level.
#' Following the studied of Välikangas et. al (2016) we use the fast
#' normalization in `normalizeCyclicLoess` function from {limma} package.
rm(list = ls())
suppressMessages(library(limma))
suppressMessages(library(QFeatures))

path_data <- "./03__modeling/2022-10-14__workflow/processing/processed_data"

# Import data -------------------------------------------------------------
fts <- readRDS(file.path(path_data, "fts_processsed.rds"))
assayNames(fts[["proteins_median"]])[1L] <- "intensity"
assayNames(fts[["proteins_mean"]])[1L] <- "intensity"
assayNames(fts[["proteins_sum"]])[1L] <- "intensity"
head(assay(fts[["proteins_sum"]], 1))

# Performing the normalization --------------------------------------------
assay(fts[["proteins_median"]], "log2_normalized") <- normalizeCyclicLoess(
  x = assay(fts[["proteins_median"]], "log2_imputed"), method = "fast")
assay(fts[["proteins_mean"]], "log2_normalized") <- normalizeCyclicLoess(
  x = assay(fts[["proteins_mean"]], "log2_imputed"), method = "fast")
assay(fts[["proteins_sum"]], "log2_normalized") <- normalizeCyclicLoess(
  x = assay(fts[["proteins_sum"]], "log2_imputed"), method = "fast")

# Save SE object with normalized intensities ------------------------------
saveRDS(object = fts, file = file.path(path_data, "fts_processsed.rds"))
