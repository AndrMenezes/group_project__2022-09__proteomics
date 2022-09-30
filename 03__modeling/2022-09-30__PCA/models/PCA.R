# 2022-09-30 Replicating PCA Analysis

suppressPackageStartupMessages(library(tidyverse))
library(readxl)
library(ggplot2)
library(cowplot)

# Reading in data file
path_data <- "./01__database/"

project_data <- read_xlsx(
  file.path(path_data, "processed_data/data_filter.xlsx"))

# Removing protein names
project_data <- subset(project_data, select = -protein__name)


# PCA Analysis --------------------------------------------------------------
# prcomp() assumes the data rows are samples and columns are proteins.
# This is not the case with the current data set, so we need to transpose it. 
# This is done with the t() function.

pca <- prcomp(t(project_data))

# Calculating PC variance
# We get the square of the standard deviations and convert to a %
pca_var <- pca$sdev^2
pca_var_percent <- round(pca_var / sum(pca_var) * 100, 1)

# Creating an antimicrobial variable (for ggplot colours and shapes)
antimicrobial <- c(rep("Ampicillin", 3), rep("Cefotaxime", 3), rep("Impipenem", 3), 
                   rep("Ciprofloxacin", 3), rep("Control", 3))

# Creating a dataframe
pca_data <- data.frame(antimicrobial = antimicrobial, 
                       X = pca$x[, 1], 
                       Y = pca$x[, 2])

# Converting antimicrobial column to a factor
pca_data$antimicrobial <- as.factor(pca_data$antimicrobial)

# Plotting  -----------------------------------------------------------------
pca_plot <- ggplot(data = pca_data, aes(x = X, y = Y)) +
  geom_point(aes(colour = antimicrobial, shape = antimicrobial), size = 5) +
  scale_shape_manual(values = c(1, 5, 0, 2, 11)) +
  xlab(paste("Component 1 (", pca_var_percent[1], "%)", sep = "")) +
  ylab(paste("Component 2 (", pca_var_percent[2], "%)", sep = "")) +
  theme_bw() +
  ggtitle("Replicated PCA Plot")

# Exporting -----------------------------------------------------------------
export_path <- "./03__modeling/2022-09-30__PCA/models/"

save_plot(file.path(export_path, "results/PCA_Plot.png"), pca_plot, 
          base_height = 8, base_aspect_ratio = 1.4)
