#' 2022-10-19: Perform normalization on the log-intensities at protein level.
#' Following the studied of Välikangas et. al (2016) we use the fast
#' normalization in `normalizeCyclicLoess` function from {limma} package.
rm(list = ls())
suppressMessages(library(limma))
suppressMessages(library(QFeatures))

path_data <- "./03__modeling/2022-10-14__workflow/processing/processed_data"

# Import data -------------------------------------------------------------
fts <- readRDS(file.path(path_data, "fts_processsed.rds"))
se <- fts[["proteins"]]
assayNames(se)[1L] <- "intensity"
head(assay(se, 1))

# Performing the normalization --------------------------------------------
assay(se, "log2_normalized") <- normalizeCyclicLoess(
  x = assay(se, "log2_imputed"), method = "fast")
se
fts[["proteins"]] <- se

# Save SE object with normalized intensities ------------------------------
saveRDS(object = fts, file = file.path(path_data, "fts_processsed.rds"))
