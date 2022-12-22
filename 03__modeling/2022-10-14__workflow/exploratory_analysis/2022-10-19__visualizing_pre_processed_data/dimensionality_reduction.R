#' 2022-11-29: Dimensionality reduction comparison.
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
se_fiona <- readRDS("./01__database/processed_data/se_margalit.rds")
rownames(se_fiona) <- rowData(se_fiona)$protein__id
se <- fts[["proteins_median"]]
colData(se) <- colData(fts)


# PCA plot ----------------------------------------------------------------
pca_normalized <- scater::calculatePCA(
  x = se, exprs_values = "log2_normalized")
pca_fiona <- scater::calculatePCA(
  x = se_fiona, exprs_values = "log2_intensity")

sum(attr(pca_normalized, "percentVar")[1:2])
sum(attr(pca_fiona, "percentVar")[1:2])
tb <- dplyr::as_tibble(pca_normalized) |>
  dplyr::mutate(group = colData(se)$group,
                replicate = colData(se)$replicate,
                measure = "Alternative") |>
  dplyr::bind_rows(
    dplyr::as_tibble(pca_fiona) |>
      dplyr::mutate(group = colData(se_fiona)$group,
                    replicate = colData(se_fiona)$replicate,
                    measure = "Margalit et. al (2022)")) |> 
  dplyr::mutate(measure = forcats::fct_relevel(
    measure, "Margalit et. al (2022)"))
ylim <- c(min(tb$PC2) - 0.4, max(tb$PC2) + 1.5)

p_pca <- ggplot(tb, aes(x = PC1, y = PC2, col = group)) +
  facet_wrap(~measure, ncol = 1) +
  geom_point(size = 4) +
  geom_text(aes(label = replicate), vjust = -0.9) +
  labs(x = "PCA 1", y = "PCA 2", col = "") +
  scale_x_continuous(breaks = scales::pretty_breaks(8)) +
  scale_y_continuous(breaks = scales::pretty_breaks(8), limits = ylim)
save_plot(filename = file.path(path_res, "results", "pca_comparison.png"),
          plot = p_pca, base_height = 7, bg = "white")


# UMAP plot ---------------------------------------------------------------

set.seed(666)
umap_normalized <- scater::calculateUMAP(
  x = se, exprs_values = "log2_normalized")
umap_fiona <- scater::calculateUMAP(
  x = se_fiona, exprs_values = "log2_intensity")
colnames(umap_normalized) <- colnames(umap_fiona) <- c("UMAP_1", "UMAP_2")


tb_umap <- dplyr::as_tibble(umap_normalized) |>
  dplyr::mutate(group = colData(se)$group,
                measure = "Log2 intensity normalized") |>
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
