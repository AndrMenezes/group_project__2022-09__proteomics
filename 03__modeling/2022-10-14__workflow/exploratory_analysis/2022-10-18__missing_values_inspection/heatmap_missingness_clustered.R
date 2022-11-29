#' 2022-10-18: Visualizations of missing values.
suppressMessages(library(QFeatures))
suppressMessages(library(gplots))
library(ggplot2)
library(cowplot)
theme_set(
  theme_cowplot() +
    background_grid() +
    theme(legend.position = "top")
)

path_data <- "./03__modeling/2022-10-14__workflow/processing/processed_data"
path_res <- "./03__modeling/2022-10-14__workflow/exploratory_analysis/2022-10-18__missing_values_inspection"

# Import data -------------------------------------------------------------
fts <- readRDS("./01__database/processed_data/QFeature_obj.rds")

fts <- aggregateFeatures(object = fts, i = "peptides", fcol = "protein",
                         name = "proteins", fun = colMeans,
                         na.rm = TRUE)

# Pct the missing values
mean(is.na(assay(fts[["peptides"]])))
mean(is.na(assay(fts[["proteins"]])))



# Heatmap with missing values at peptide and protein level ----------------
x_peptides <- assay(fts[["peptides"]])
colnames(x_peptides) <- colData(fts)$sample_names
x_peptides[!is.na(x_peptides)] <- 0
x_peptides[is.na(x_peptides)] <- 1
head(x_peptides)

x_proteins <- assay(fts[["proteins"]])
colnames(x_proteins) <- colData(fts)$sample_names
x_proteins[!is.na(x_proteins)] <- 0
x_proteins[is.na(x_proteins)] <- 1
head(x_proteins)


png(file.path(path_res, "results",  "heatmap_missing_values__peptides.png"),
    width = 1000)
heatmap.2(x_peptides[1:100, ], col = c("lightgray", "black"),
          scale = "none", dendrogram = "none",
          trace = "none", keysize = 0.5, key = FALSE)
title("Missing values at peptide level")
graphics.off()

png(file.path(path_res, "results",  "heatmap_missing_values__proteins.png"),
    width = 1000)
heatmap.2(x_proteins, col = c("lightgray", "black"),
          scale = "none", dendrogram = "none",
          trace = "none", keysize = 0.5, key = FALSE)
title("Missing values at proteins level")
graphics.off()
