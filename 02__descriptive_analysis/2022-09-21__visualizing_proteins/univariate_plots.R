#' 2022-09-21: .
library(dplyr)
library(ggplot2)
library(cowplot)
theme_set(
  theme_cowplot(font_size = 16) +
    background_grid() +
    theme(text = element_text(size = 16))
)

path_data <- "./01__database/"
path_res <- "./02__descriptive_analysis/2022-09-21__visualizing_proteins/results"

data_all_protein <- read.csv(
  file.path(path_data, "processed_data/all_protein__pivotted.csv"))
head(data_all_protein)
glimpse(data_all_protein)

length(unique(data_all_protein$protein__id))
data_curr <- data_all_protein |> 
  filter(protein__id == "P04949")
data_curr$group <- relevel(factor(data_curr$group), ref = "Control")
mod_1 <- lm(value ~ group, data = data_curr)
summary(mod_1)


ggplot(data_all_protein, aes(x = value)) +
  facet_wrap(~group) +
  geom_density() +
  geom_rug()

data_summarized <- data_all_protein |> 
  group_by(protein__id, group) |> 
  summarise(avg = mean(value), sd = sd(value), .groups = "drop")

ggplot(data_summarized, aes(x = avg, y = sd)) +
  facet_wrap(~group) +
  geom_point(size = 1.2, alpha = 0.6) +
  geom_smooth(col = "blue")

antibiotics <- c("Ampicillin", "Cefotaxime", "Ciprofloxacin", "Impipenem")
k <- 1L
p_list <- list()
y_lim <- range(data_summarized$sd)
x_lim <- range(data_summarized$avg)
for (chosen in antibiotics) {
  data_curr <- data_summarized |> 
    filter(group %in% c(chosen, "Control")) |> 
    mutate(group = forcats::fct_relevel(group, "Control"))
  
  p_list[[k]] <- ggplot(data_curr, aes(x = avg, y = sd, col = group)) +
    geom_point(size = 1.8, alpha = 0.4) +
    geom_smooth(se = FALSE) +
    ggtitle(paste0("Relationship Mean-Std of each protein for Control and ",
                   chosen, " groups")) +
    colorspace::scale_color_discrete_qualitative() +
    labs(x = "Average of log intensity",
         y = "Standard deviation of log intensity",
         col = "") +
    scale_x_continuous(limits = x_lim, breaks = scales::pretty_breaks(6)) +
    scale_y_continuous(limits = y_lim, breaks = scales::pretty_breaks(6))
  k <- k + 1L
}
p_grid <- plot_grid(plotlist = p_list)
save_plot(
  filename = file.path(path_res, "scatter_avg_std.png"), plot = p_grid,
  base_height = 8, base_asp = 2.5, bg = "white")


