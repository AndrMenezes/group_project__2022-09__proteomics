library(dplyr)
library(tidyr)
library(readxl)
library(openxlsx)

path_data_1 <- "./03__modeling/2022-09-30__Multiple_Comparisons/models/"
path_data_2 <- "./03__modeling/2022-10-14__multiple_comparisons/models/"

replicated_data <- read_xlsx(
  file.path(path_data_1, "results/replicated_results.xlsx"),
  sheet = "All proteins")

limma_data <- read_xlsx(
  file.path(path_data_2, "results/limma_results.xlsx"),
  sheet = "All proteins")

# Data manipulation------------------------------------------------------------
replicated_data <- replicated_data |>
  rename(Protein = "Protein_id", Gene = "Gene_name", logFC = "LFQ") |>
  mutate(Method = "Two Sample T Test") |>
  select(Control, Treatment, Protein, Gene, Method, logFC, Statistic, 
         CI_L, CI_R, P_value, FDR)

limma_data <- limma_data |>
  rename(Protein = "protein__name", Gene = "ID", CI_L = "CI.L", CI_R = "CI.R", 
         Statistic = "t", P_value = "P.Value", FDR = "adj.P.Val") |>
  mutate(Method = "Bayesian T Test") |>
  select(Control, Treatment, Protein, Gene, Method, logFC, Statistic, 
         CI_L, CI_R, P_value, FDR)

# Concatenating data ----------------------------------------------------------
all_results <- rbind(replicated_data, limma_data)
all_results <- all_results |>
  arrange(Protein)

# Exporting -------------------------------------------------------------------
export_path <- "./03__modeling/2022-10-19__comparison/"

write.xlsx(all_results, file.path(export_path, "results/all_results.xlsx"),
           rowNames = FALSE)