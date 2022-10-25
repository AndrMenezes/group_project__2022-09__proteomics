# 2022/10/25 - Calculating Hedges G effect sizes
library(dplyr)
library(tidyr)
library(effsize) # For calculating Hedges G
library(readxl)
library(writexl)

path_data <- "./03__modeling/2022-10-25__feature_selection/processing/"

# Note: processed_data.xlsx is on Jira (PROT-83)
project_data <- read_xlsx(file.path(path_data, "processed_data/processed_data.xlsx"))

# Data manipulation -----------------------------------------------------------
# Focusing only on normalized results
# Pivoting: rows = proteins // columns = samples
normalized_data <- project_data |>
  filter(variable == "log_intensity_normalized") |>
  arrange(proteins, group) |>
  select(-variable) |>
  pivot_wider(names_from = c(replicate, group),
              values_from = value)

# Creating data frames for storing output -------------------------------------
ampicillin_values <- data.frame(Protein_id = character(0),
                                Statistic = numeric(0),
                                logFC = numeric(0),
                                P_value = numeric(0),
                                Effect = numeric(0))

cefotaxime_values <- data.frame(Protein_id = character(0),
                                Statistic = numeric(0),
                                logFC = numeric(0),
                                P_value = numeric(0),
                                Effect = numeric(0))

ciprofloxacin_values <- data.frame(Protein_id = character(0),
                                   Statistic = numeric(0),
                                   logFC = numeric(0),
                                   P_value = numeric(0),
                                   Effect = numeric(0))

impipenem_values <- data.frame(Protein_id = character(0),
                               Statistic = numeric(0),
                               logFC = numeric(0),
                               P_value = numeric(0),
                               Effect = numeric(0))

# Performing T-Test and effect size calculation -------------------------------
for(i in 1:nrow(normalized_data)) {
  
  # Ampicillin
  ampicillin_values[i, 1] <- normalized_data$proteins[i]
  ampicillin_values[i, 2] <- t.test(normalized_data[i, 2:4], normalized_data[i, 11:13], var.equal = TRUE)$statistic
  ampicillin_values[i, 3] <- mean(as.numeric(normalized_data[i, 2:4]) - as.numeric(normalized_data[i, 11:13]))
  ampicillin_values[i, 4] <- t.test(normalized_data[i, 2:4], normalized_data[i, 11:13], var.equal = TRUE)$p.value
  ampicillin_values[i, 5] <- cohen.d(as.numeric(normalized_data[i, 2:4]), as.numeric(normalized_data[i, 11:13]), hedges = TRUE)$estimate
  
  # Cefotaxime
  cefotaxime_values[i, 1] <- normalized_data$proteins[i]
  cefotaxime_values[i, 2] <- t.test(normalized_data[i, 5:7], normalized_data[i, 11:13], var.equal = TRUE)$statistic
  cefotaxime_values[i, 3] <- mean(as.numeric(normalized_data[i, 5:7]) - as.numeric(normalized_data[i, 11:13]))
  cefotaxime_values[i, 4] <- t.test(normalized_data[i, 5:7], normalized_data[i, 11:13], var.equal = TRUE)$p.value
  cefotaxime_values[i, 5] <- cohen.d(as.numeric(normalized_data[i, 5:7]), as.numeric(normalized_data[i, 11:13]), hedges = TRUE)$estimate
  # Ciprofloxacin
  ciprofloxacin_values[i, 1] <- normalized_data$proteins[i]
  ciprofloxacin_values[i, 2] <- t.test(normalized_data[i, 8:10], normalized_data[i, 11:13], var.equal = TRUE)$statistic
  ciprofloxacin_values[i, 3] <- mean(as.numeric(normalized_data[i, 8:10]) - as.numeric(normalized_data[i, 11:13]))
  ciprofloxacin_values[i, 4] <- t.test(normalized_data[i, 8:10], normalized_data[i, 11:13], var.equal = TRUE)$p.value
  ciprofloxacin_values[i, 5] <- cohen.d(as.numeric(normalized_data[i, 8:10]), as.numeric(normalized_data[i, 11:13]), hedges = TRUE)$estimate
  
  # Impipenem
  impipenem_values[i, 1] <- normalized_data$proteins[i]
  impipenem_values[i, 2] <- t.test(normalized_data[i, 14:16], normalized_data[i, 11:13], var.equal = TRUE)$statistic
  impipenem_values[i, 3] <- mean(as.numeric(normalized_data[i, 14:16]) - as.numeric(normalized_data[i, 11:13]))
  impipenem_values[i, 3] <- t.test(normalized_data[i, 14:16], normalized_data[i, 11:13], var.equal = TRUE)$p.value
  impipenem_values[i, 4] <- cohen.d(as.numeric(normalized_data[i, 14:16]), as.numeric(normalized_data[i, 11:13]), hedges = TRUE)$estimate
}

# Adding in treatment type and binding ----------------------------------------
ampicillin_values <- mutate(ampicillin_values, Treatment = "Ampicillin")
cefotaxime_values <- mutate(cefotaxime_values, Treatment = "Cefotaxime")
ciprofloxacin_values <- mutate(cefotaxime_values, Treatment = "Ciprofloxacin")
impipenem_values <- mutate(cefotaxime_values, Treatment = "Impipenem")

effect_sizes <- rbind(ampicillin_values, cefotaxime_values, ciprofloxacin_values, impipenem_values)

effect_sizes <- effect_sizes|>
  select(Protein_id, Treatment, everything()) |>
  arrange(Protein_id)

# Exporting -------------------------------------------------------------------
export_path <- "./03__modeling/2022-10-25__feature_selection/models/"

write_xlsx(effect_sizes, file.path(export_path, "results/effect_sizes.xlsx"))