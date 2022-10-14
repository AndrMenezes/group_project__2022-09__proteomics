#' 2022-10-14: Perform normalization on the log-intensities.
#' Following the studied of Välikangas et. al (2016) we use the fast
#' normalization in `normalizeCyclicLoess` function from {limma} package.
suppressMessages(library(SummarizedExperiment))
suppressMessages(library(HDF5Array))
suppressMessages(library(limma))
suppressMessages(library(scran))
library(ggridges)
library(magrittr)
library(ggplot2)
library(cowplot)
theme_set(
  theme_cowplot() +
    background_grid() +
    theme(legend.position = "top")
)



path_data <- "./01__database/"
path_res <- "./03__modeling/2022-10-14__workflow/processing"

se <- loadHDF5SummarizedExperiment(
  dir = file.path(path_data, "processed_data/se_obj"))
assay(se, "log_intensity_normalized") <- normalizeCyclicLoess(x = assay(se))
se


# Save SE object with normalized intensities in .h5 format ----------------
saveHDF5SummarizedExperiment(
  x = se, dir = file.path(path_res, "processed_data/se_obj"))


# Organizing the data to compare pre and post normalization ---------------

pivotting_data <- function(se, i = "log_intensity") {
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
  pivotting_data(se = se, i = "log_intensity"),
  pivotting_data(se = se, i = "log_intensity_normalized")) |> 
  mutate(is_normalized = variable == "log_intensity_normalized")

p_densities <- ggplot(data_pivotted, aes(x = value, fill = is_normalized,
                                         col = is_normalized)) +
  facet_wrap(~group) +
  geom_density(alpha = 0.4) +
  geom_rug(show.legend = FALSE) +
  ggtitle("Density of log-intensity abundance of all proteins with and without normalization") +
  labs(x = "Log2 intensity", y = "Density", fill = "Is normalized?", 
       col = "Is normalized?")

p_ecdfs <- ggplot(data_pivotted, aes(x = value, col = is_normalized)) +
  facet_wrap(~group) +
  stat_ecdf() +
  geom_rug(show.legend = FALSE) +
  ggtitle("Empirical cumulative distribution of log-intensity abundance of all proteins with and without normalization") +
  labs(x = "Log2 intensity", y = "Cumulative probability",
       col = "Is normalized?")

ggplot(data = filter(data_pivotted)) +
  aes(x = value, y = group) +
  geom_density_ridges(scale = 0.9) +
  facet_wrap(~is_normalized)


