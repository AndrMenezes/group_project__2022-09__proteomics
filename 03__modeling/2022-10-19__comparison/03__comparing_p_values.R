library(dplyr)
library(tidyr)
library(readxl)
library(ggplot2)
library(cowplot)

path_data <- "./03__modeling/2022-10-19__comparison"

project_data <- read_xlsx(file.path(path_data, "results/all_results.xlsx"))

# Setting up data for scatter plot -------------------------------------------
replicate_data <- project_data |>
  filter(Method == "Two Sample T Test") |>
  select(c("Treatment", "Method", "P_value"))

limma_data <- project_data |>
  filter(Method == "Moderated T Test") |>
  select(c("Treatment", "Method", "P_value"))

plotting <- data.frame(Treatment = replicate_data$Treatment,
                       replicate_p_value = replicate_data$P_value,
                       limma_p_value = limma_data$P_value)


# Scatter plot ---------------------------------------------------------------
scatter_plot <- ggplot(data = plotting,
                       aes(x = replicate_p_value, y = limma_p_value)) +
  geom_point(size = .25) +
  facet_wrap(~ Treatment, nrow = 2, ncol = 2) +
  geom_abline(slope = 1, intercept = 0) +
  ggtitle("Two Sample vs. Moderated T-Test P Values") +
  xlab("Two Sample P Values") +
  ylab("Moderated P Values") +
  theme_bw()

# Density plot ---------------------------------------------------------------
density_plot <- ggplot(data = project_data, 
                       aes(x = P_value, fill = Method)) +
  geom_density(alpha = .4) +
  facet_wrap(~ Treatment, nrow = 2, ncol = 2) +
  ggtitle("Two Sample vs. Moderated T-Test P Values") +
  xlab("Two Sample P Values") +
  ylab("Moderated P Values") +
  theme_bw()

# Exporting ------------------------------------------------------------------
save_plot(file.path(path_data, "results/Scatter_plot.png"), scatter_plot)
save_plot(file.path(path_data, "results/Density_plot.png"), density_plot)
