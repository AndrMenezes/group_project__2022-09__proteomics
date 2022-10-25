# 2022-10-25: Reducing proteins based off an effect size cut off
# Not fully completed yet - reasoning behind an effect size cut off needs to 
# be fully determined
library(dplyr)
library(readxl)
library(writexl)

path_data_1 <- "./03__modeling/2022-10-25__feature_selection/processing/"
path_data_2 <- "./03__modeling/2022-10-25__feature_selection/models/"

project_data <- read_xlsx(file.path(path_data_1, "processed_data/processed_data.xlsx"))
effect_sizes <- read_xlsx(file.path(path_data_2, "results/effect_sizes.xlsx"))

# Data manipulation -----------------------------------------------------------
# Step 1: Feature Selection
# Selecting proteins with large effect sizes: >= 0.8 or <= -.08)
# This is a TEMPORARY idea
reduced_proteins <- effect_sizes |>
  filter(abs(Effect) >= 0.8) |>
  arrange(Protein_id)

# Step 2: Match these proteins to the protein_data data set
normalized_data_reduced <- project_data |>
  filter(variable == "log_intensity_normalized") |>
  filter(proteins %in% reduced_proteins$Protein_id) |>
  arrange(proteins)

# Exporting -------------------------------------------------------------------
write_xlsx(normalized_data_reduced, 
           file.path(path_data_2, "results/normalised_data_reduced.xlsx"))
