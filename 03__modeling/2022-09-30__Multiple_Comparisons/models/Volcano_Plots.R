# 2022-09-30 Replicating T-Tests and Volcano Plots
# TODO - Use cowplot to export graphs
#      - Label proteins

suppressPackageStartupMessages(library(tidyverse))
library(readxl)
library(ggplot2)
library(cowplot)

# Reading in data set
path_data <- "./01__database/"

project_data <- read_xlsx(
  file.path(path_data, "processed_data/data_filter.xlsx"))

project_data <- subset(project_data, select = -protein__name)

# view(project_data)


# Creating an empty data frame with column names
columns <- c("Ampicillin_p_value", "Cefotaxime_p_value", 
             "Imipenem_p_value", "Ciprofloxacin_p_value", 
             "Ampicillin_lfq", "Cefotaxime_lfq", 
             "Imipenem_lfq", "Ciprofloxacin_lfq")

analysis_values <- data.frame(matrix(nrow = 0, ncol = length(columns)))

colnames(analysis_values) <- columns

# Assigning p values and mean LFQ differences to dataframe
for(i in 1:nrow(project_data)){
  
  # Ampicillin
  p_amp <- t.test(project_data[i, 1:3], project_data[i, 13:15], var.equal = TRUE)
  lfq_amp <- mean(as.numeric(project_data[i, 1:3]) - as.numeric(project_data[i, 13:15]))
  
  # Cefotaxime
  p_cef <- t.test(project_data[i, 4:6], project_data[i, 13:15], var.equal = TRUE)
  lfq_cef <-  mean(as.numeric(project_data[i, 4:6]) - as.numeric(project_data[i, 13:15]))
  
  # Imipenem
  p_imp <- t.test(project_data[i, 7:9], project_data[i, 13:15], var.equal = TRUE)
  lfq_imp <-  mean(as.numeric(project_data[i, 7:9]) - as.numeric(project_data[i, 13:15]))
  
  # Ciprofloxacin
  p_cip <- t.test(project_data[i, 10:12], project_data[i, 13:15], var.equal = TRUE)
  lfq_cip <-  mean(as.numeric(project_data[i, 10:12]) - as.numeric(project_data[i, 13:15]))
  
  # Assigning the p values
  analysis_values[i, 1] <- -log10(p_amp$p.value)
  analysis_values[i, 2] <- -log10(p_cef$p.value)
  analysis_values[i, 3] <- -log10(p_imp$p.value)
  analysis_values[i, 4] <- -log10(p_cip$p.value)
  
  # Assigning the lfq differences
  analysis_values[i, 5] <- lfq_amp
  analysis_values[i, 6] <- lfq_cef
  analysis_values[i, 7] <- lfq_imp
  analysis_values[i, 8] <- lfq_cip
}

# view(analysis_values)


# Generating volcano plots --------------------------------------------------
# Ampicillin plot
ggplot(data = analysis_values, aes(x = Ampicillin_lfq, y = Ampicillin_p_value)) +
  geom_point(size = 2/5) +
  ggtitle("Control & Ampicillin") +
  xlab(expression("log"[2]*" Difference")) + 
  ylab(expression("-log"[10]* "p")) +
  geom_vline(xintercept=c(-0.5, 0.5), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") +
  theme_bw()

# Cefotaxime plot
ggplot(data = analysis_values, aes(x = Cefotaxime_lfq, y = Cefotaxime_p_value)) +
  geom_point(size = 2/5) +
  ggtitle("Control & Cefotaxime") +
  xlab(expression("log"[2]*" Difference")) + 
  ylab(expression("-log"[10]* "p")) +
  geom_vline(xintercept=c(-0.5, 0.5), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") +
  theme_bw()

# Imipenem plot
ggplot(data = analysis_values, aes(x = Imipenem_lfq, y = Imipenem_p_value)) +
  geom_point(size = 2/5) +
  ggtitle("Control & Imipenem") +
  xlab(expression("log"[2]*" Difference")) + 
  ylab(expression("-log"[10]* "p")) +
  geom_vline(xintercept=c(-0.5, 0.5), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") +
  theme_bw()

# Ciprofloxacin plot
ggplot(data = analysis_values, aes(x = Ciprofloxacin_lfq, y = Ciprofloxacin_p_value)) +
  geom_point(size = 2/5) +
  ggtitle("Control & Ciprofloxacin") +
  xlab(expression("log"[2]*" Difference")) + 
  ylab(expression("-log"[10]* "p")) +
  geom_vline(xintercept=c(-0.5, 0.5), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") +
  theme_bw()
