#' 2022-10-18: Visualizations of normalization log2 intensities.
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
  pivotting_data(se = se_fiona, i = "log2_intensity")) |> 
  dplyr::mutate(variable = dplyr::case_when(
    variable == "log2_intensity" ~ "Margalit et. al (2022)",
    variable == "intensity" ~ "Unprocessed intensity",
    variable == "log2_imputed" ~ "Log2 intensity imputed",
    variable == "log2_normalized" ~ "Log2 intensity normalized"),
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

p_violin <- ggplot(data_pivotted, aes(x = group, y = value, fill = group,
                                         col = group)) +
  facet_wrap(~variable, scales = "free") +
  geom_violin(alpha = 0.6) +
  ggtitle("Violin plots of proteins abundance by variable and group") +
  labs(x = "Abundance", y = "Density", fill = "", col = "")
save_plot(filename = file.path(path_res, "results", "violins_comparison.png"),
          plot = p_violin, base_height = 7, bg = "white")

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
  scale_y_continuous(breaks = scales::pretty_breaks(6)) +
  colorspace::scale_fill_discrete_qualitative()
save_plot(filename = file.path(path_res, "results",
                               "densities_comparison_processing.png"),
          plot = p_densities_appr, base_height = 7, bg = "white")

