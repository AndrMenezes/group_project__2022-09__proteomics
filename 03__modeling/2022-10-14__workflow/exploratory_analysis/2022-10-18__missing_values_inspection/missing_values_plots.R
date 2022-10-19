#' 2022-10-18: Visualizations of missing values.
suppressMessages(library(QFeatures))
library(ggplot2)
library(cowplot)
theme_set(
  theme_cowplot() +
    background_grid() +
    theme(legend.position = "top")
)

path_data <- "./01__database/"
path_res <- "./03__modeling/2022-10-14__workflow/exploratory_analysis/2022-10-18__missing_values_inspection"

# Import data -------------------------------------------------------------
fts <- readRDS(file.path(path_data, "QFeature_obj.rds"))

# Pct the missing values
mean(is.na(assay(fts[["peptides"]])))


# Inspect the missing values at peptides level ----------------------------
data_ <- as.data.frame(assay(fts[["peptides"]]))
colnames(data_) <- fts$sample_names
p <- naniar::vis_miss(data_) +
  ggtitle("Missing values at peptide level")
save_plot(
  filename = file.path(path_res, "results", "missing__peptides.png"),
  bg = "white", plot = p, base_height = 8)

# Aggregating at protein level and inspect the missing values -------------
fts <- aggregateFeatures(object = fts, i = "peptides", fcol = "protein",
                         name = "proteins", fun = colMedians, na.rm = TRUE)

data_ <- as.data.frame(assay(fts[["proteins"]]))
colnames(data_) <- fts$sample_names
p <- naniar::vis_miss(data_) +
  ggtitle("Missing values at protein level")
save_plot(
  filename = file.path(path_res, "results", "missing__protein.png"),
  bg = "white", plot = p, base_height = 8)
