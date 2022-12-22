#' 2022-12-15: Visualizations of imputation and normalization steps.
library(ggplot2)
library(cowplot)
library(ggridges)
theme_set(
  theme_cowplot() +
    background_grid() +
    theme(legend.position = "top")
)

path_data <- "./03__modeling/2022-10-14__workflow/processing/processed_data"
path_res <- "./03__modeling/2022-10-14__workflow/exploratory_analysis/2022-11-12__exploring_missing_values"

# Import data -------------------------------------------------------------
fts <- readRDS(file.path(path_data, "fts_processsed.rds"))
se_fiona <- readRDS("./01__database/processed_data/se_margalit.rds")
rownames(se_fiona) <- rowData(se_fiona)$protein__id
se_median <- fts[["proteins_median"]]
se_mean <- fts[["proteins_mean"]]
se_sum <- fts[["proteins_sum"]]
colData(se_median) <- colData(se_mean) <- colData(se_sum) <- colData(fts)

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

pivotted_fiona_miss <- pivotting_data(se = se_fiona, i = "log2_intensity")
pivotted_fiona_impu <- pivotting_data(se = se_fiona, i = "log2_imputed")

pivotted_ours_miss__med <- pivotting_data(se = se_median, i = "log2_intensity")
pivotted_ours_impu__med <- pivotting_data(se = se_median, i = "log2_normalized")

pivotted_ours_miss__mea <- pivotting_data(se = se_mean, i = "log2_intensity")
pivotted_ours_impu__mea <- pivotting_data(se = se_mean, i = "log2_normalized")

pivotted_ours_miss__sum <- pivotting_data(se = se_sum, i = "log2_intensity")
pivotted_ours_impu__sum <- pivotting_data(se = se_sum, i = "log2_normalized")


join_pivotted <- function(data_miss, data_imputed) {
  data_miss |> 
    dplyr::rename(value_imp = value) |>
    dplyr::select(-variable) |> 
    dplyr::left_join(dplyr::rename(data_imputed, value_n_imp = value) |> 
                dplyr::select(-variable)) |> 
    dplyr::mutate(is_imputed = is.na(value_imp))
}

pivotted_fiona <- join_pivotted(
  data_miss = pivotted_fiona_miss, data_imputed = pivotted_fiona_impu)
pivotted_median <- join_pivotted(
  data_miss = pivotted_ours_miss__med, data_imputed = pivotted_ours_impu__med)
pivotted_mean <- join_pivotted(
  data_miss = pivotted_ours_miss__mea, data_imputed = pivotted_ours_impu__mea)
pivotted_sum <- join_pivotted(
  data_miss = pivotted_ours_miss__sum, data_imputed = pivotted_ours_impu__sum)

pivotted_all <- dplyr::bind_rows(
  dplyr::mutate(pivotted_fiona, method = "Margalit et. al (2022)"),
  dplyr::mutate(pivotted_median, method = "Median"),
  dplyr::mutate(pivotted_mean, method = "Mean"),
  dplyr::mutate(pivotted_sum, method = "Sum")) |>
  dplyr::mutate(method = forcats::fct_relevel(
    method, "Margalit et. al (2022)",  "Sum"
  ))


p_hist_all <- ggplot(pivotted_all, aes(x = value_n_imp, fill = is_imputed)) +
  facet_wrap(~ method) +
  geom_histogram(alpha = 0.5, bins = 50) +
  labs(x = "Log2 intensity", y = "Frequency", fill = "Is imputed?")
x11(); p_hist_all
save_plot(filename = file.path(path_res, "results", "hist_all.png"),
          base_height = 6, bg = "white", plot = p_hist_all)

p_hist_chosen <- pivotted_all |> 
  dplyr::filter(method %in% c("Median", "Margalit et. al (2022)")) |>
  dplyr::mutate(method = ifelse(method == "Median", "Alternative",
                                "Margalit et. al (2022)")) |>
  dplyr::mutate(method = forcats::fct_relevel(method, "Margalit et. al (2022)")) |> 
  ggplot(aes(x = value_n_imp, fill = is_imputed)) +
  facet_wrap(~ method, scales = "free_x") +
  geom_histogram(alpha = 0.5, bins = 50, col = "grey45") +
  labs(x = "Log2 intensity", y = "Frequency", fill = "Is imputed?") +
  colorspace::scale_fill_discrete_qualitative()
# x11(); p_hist_chosen
save_plot(filename = file.path(path_res, "results", "hist_chosen.png"),
          base_height = 6, bg = "white", plot = p_hist_chosen)

# Filtering out the 
x11()
pivotted_all |>
  dplyr::filter(proteins %in% rownames(se_fiona)) |> 
  ggplot(aes(x = value_n_imp, fill = is_imputed)) +
  facet_grid(~ method) +
  geom_histogram(alpha = 0.5, bins = 50)

# Comparing the sum intensity approaches
p_densities_all <- pivotted_all |> 
  ggplot(aes(x = value_n_imp, fill = method, col = method)) + 
  geom_density(alpha = 0.5) +
  scale_x_continuous(breaks = scales::pretty_breaks(8)) +
  scale_y_continuous(breaks = scales::pretty_breaks(6)) +
  labs(x = "Log2 intensity", y = "Density", fill = "", col = "")
save_plot(filename = file.path(path_res, "results", "densities_all.png"),
          base_height = 6, bg = "white", plot = p_densities_all)

p_densities_all_ridges <- pivotted_all |> 
  ggplot(aes(x = value_n_imp, y = method)) + 
  geom_density_ridges() +
  scale_x_continuous(breaks = scales::pretty_breaks(8)) +
  labs(x = "Log2 intensity", y = "")
save_plot(filename = file.path(path_res, "results", "densities_all_ridges.png"),
          base_height = 6, bg = "white", plot = p_densities_all_ridges)

p_densities_sum <- pivotted_all |> 
  dplyr::filter(method %in% c("Margalit et. al (2022)", "Sum"),
                proteins %in% rownames(se_fiona)) |> 
  ggplot(aes(x = value_n_imp, fill = method, col = method)) + 
  geom_density(alpha = 0.5) +
  geom_rug(show.legend = FALSE) +
  scale_x_continuous(breaks = scales::pretty_breaks(8)) +
  scale_y_continuous(breaks = scales::pretty_breaks(6)) +
  labs(x = "Log2 intensity", y = "Density", fill = "", col = "")
save_plot(filename = file.path(path_res, "results", "densities_sum.png"),
          base_height = 6, bg = "white", plot = p_densities_sum)

