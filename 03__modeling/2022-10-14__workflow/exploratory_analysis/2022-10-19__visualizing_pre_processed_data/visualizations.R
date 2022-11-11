#' 2022-10-18: Visualizations of normalization tecnique applied.
library(ggplot2)
library(cowplot)
theme_set(
  theme_cowplot() +
    background_grid() +
    theme(legend.position = "top")
)

path_data <- "./03__modeling/2022-10-14__workflow/processing/processed_data"
path_res <- "./03__modeling/2022-10-14__workflow/exploratory_analysis/2022-10-19__visualizing_pre_processed_data"

# Import data -------------------------------------------------------------
fts <- readRDS(file.path(path_data, "fts_processsed.rds"))
se_fiona <- readRDS("./01__database/processed_data/se_processed.rds")
rownames(se_fiona) <- rowData(se_fiona)$protein__id
se <- fts[["proteins"]]
colData(se) <- colData(fts)

# Organizing the data to compare pre and post normalization ---------------
pivotting_data <- function(se, i = "log2_imputed") {
  tb <- assay(se, i) |>
    t() |>
    dplyr::as_tibble() |>
    dplyr::mutate(replicate = se$replicate,
                  group = se$group,
                  variable = i) |>
    tidyr::pivot_longer(cols = -c(replicate, group, variable),
                        names_to = "proteins",
                        values_to = "value") |>
    dplyr::select(proteins, group, replicate, variable, value)
  tb
}

data_pivotted <- dplyr::bind_rows(
  pivotting_data(se = se, i = "intensity"),
  pivotting_data(se = se, i = "log2_imputed"),
  pivotting_data(se = se, i = "log2_normalized"),
  pivotting_data(se = se_fiona, i = "log_intensity")) |> 
  dplyr::mutate(variable = dplyr::case_when(
    variable == "log_intensity" ~ "Margalit et. al (2022)",
    variable == "intensity" ~ "Unprocessed intensity",
    variable == "log_intensity_imputed" ~ "Log2 intensity imputed",
    variable == "log_intensity_normalized" ~ "Log2 intensity normalized"),
    variable = forcats::fct_relevel(
      factor(variable), "Margalit et. al (2022)", "Unprocessed intensity",
      "Log2 intensity imputed"))

p_densities <- ggplot(data_pivotted, aes(x = value, fill = group,
                                         col = group)) +
  facet_wrap(~variable, scales = "free") +
  geom_density(alpha = 0.3) +
  geom_rug(show.legend = FALSE) +
  ggtitle("Densities of proteins abundance by variable and group") +
  labs(x = "Abundance", y = "Density", fill = "", col = "")
save_plot(filename = file.path(path_res, "results", "densities_comparison.png"),
          plot = p_densities, base_height = 7, bg = "white")

p_densities_appr <- data_pivotted |> 
  dplyr::filter(variable %in% c("Margalit et. al (2022)",
                         "Log2 intensity normalized"),
                proteins %in% rownames(se_fiona)) |> 
  ggplot(aes(x = value, fill = variable, col = variable)) +
  facet_wrap(~group) +
  geom_density(alpha = 0.4) +
  geom_rug(show.legend = FALSE) +
  ggtitle("Comparison of densities abundance between processing approaches") +
  labs(x = "Abundance", y = "Density", fill = "", col = "") +
  scale_x_continuous(breaks = scales::pretty_breaks(6)) +
  scale_y_continuous(breaks = scales::pretty_breaks(6))
save_plot(filename = file.path(path_res, "results",
                               "densities_comparison_processing.png"),
          plot = p_densities_appr, base_height = 7, bg = "white")



# PCA plot ----------------------------------------------------------------
pca_non_normalized <- scater::calculatePCA(
  x = se, exprs_values = "log2_imputed")
pca_normalized <- scater::calculatePCA(
  x = se, exprs_values = "log2_normalized")
pca_fiona <- scater::calculatePCA(
  x = se_fiona, exprs_values = "log_intensity")

class(pca_non_normalized)
str(pca_non_normalized)
attr(pca_non_normalized, "percentVar")
attr(pca_normalized, "percentVar")
attr(pca_fiona, "percentVar")
tb <- dplyr::as_tibble(pca_non_normalized) |>
  dplyr::mutate(sample = rownames(pca_non_normalized),
                group = colData(se)$group,
                measure = "Log2 intensity imputed") |>
  dplyr::bind_rows(
    dplyr::as_tibble(pca_normalized) |>
      dplyr::mutate(group = colData(se)$group,
                    measure = "Log2 intensity normalized")) |>
  dplyr::bind_rows(
    dplyr::as_tibble(pca_fiona) |>
      dplyr::mutate(group = colData(se_fiona)$group,
                    measure = "Margalit et. al (2022)"))

p_pca <- ggplot(tb, aes(x = PC1, y = PC2, col = group)) +
  facet_wrap(~measure, ncol = 1) +
  geom_point(size = 4) +
  labs(x = "PCA 1", y = "PCA 2", col = "") +
  scale_x_continuous(breaks = scales::pretty_breaks(8)) +
  scale_y_continuous(breaks = scales::pretty_breaks(8)) +
  ggtitle("PCA performed on imputed log2 intensity before and after normalization")
save_plot(filename = file.path(path_res, "results", "pca_comparison.png"),
          plot = p_pca, base_height = 7, bg = "white")

# UMAP plot ---------------------------------------------------------------

set.seed(666)
umap_non_normalized <- scater::calculateUMAP(
  x = se, exprs_values = "log2_imputed")
umap_normalized <- scater::calculateUMAP(
  x = se, exprs_values = "log2_normalized")
umap_fiona <- scater::calculateUMAP(
  x = se_fiona, exprs_values = "log_intensity")
colnames(umap_non_normalized) <- colnames(umap_normalized) <- colnames(
  umap_fiona) <- c("UMAP_1", "UMAP_2")


tb_umap <- dplyr::as_tibble(umap_non_normalized) |>
  dplyr::mutate(group = colData(se)$group,
                measure = "Log2 intensity imputed") |>
  dplyr::bind_rows(
    dplyr::as_tibble(umap_normalized) |>
      dplyr::mutate(group = colData(se)$group,
                    measure = "Log2 intensity normalized")) |>
  dplyr::bind_rows(
    dplyr::as_tibble(umap_fiona) |>
      dplyr::mutate(group = colData(se_fiona)$group,
                    measure = "Margalit et. al (2022)"))

p_umap <- ggplot(tb_umap, aes(x = UMAP_1, y = UMAP_2, col = group)) +
  facet_wrap(~measure, ncol = 1) +
  geom_point(size = 4) +
  labs(x = "UMAP 1", y = "UMAP 2", col = "") +
  scale_x_continuous(breaks = scales::pretty_breaks(8)) +
  scale_y_continuous(breaks = scales::pretty_breaks(8)) +
  ggtitle("UMAP performed on imputed log2 intensity before and after normalization")
save_plot(filename = file.path(path_res, "results", "umap_comparison.png"),
          plot = p_umap, base_height = 7, bg = "white")
