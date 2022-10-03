# 2022-09-30 Replicating T-Tests and Volcano Plots
# TODO - Label proteins

suppressPackageStartupMessages(library(tidyverse))
library(readxl)
library(ggplot2)
library(ggrepel)
library(cowplot)

# Reading in data set --------------------------------------------------------
path_data <- "./01__database/"

project_data <- read_xlsx(
  file.path(path_data, "processed_data/data_filter.xlsx"))

# view(project_data)


# Creating data frames --------------------------------------------------------
ampicillin_values <- data.frame(Protein_name = character(0),
                                Gene_name = character(0),
                                P_value = numeric(0),
                                LFQ= numeric(0))

cefotaxime_values <- data.frame(Protein_name = character(0),
                                Gene_name = character(0),
                                P_value = numeric(0),
                                LFQ= numeric(0))

imipenem_values <- data.frame(Protein_name = character(0),
                              Gene_name = character(0),
                              P_value = numeric(0),
                              LFQ= numeric(0))

ciprofloxacin_values <- data.frame(Protein_name = character(0),
                                   Gene_name = character(0),
                                   P_value = numeric(0),
                                   LFQ= numeric(0))


# Assigning p values and mean LFQ differences to dataframe --------------------
for(i in 1:nrow(project_data)){
  
  # Ampicillin
  p_amp <- t.test(project_data[i, 3:5], project_data[i, 15:17], var.equal = TRUE)
  lfq_amp <- mean(as.numeric(project_data[i, 3:5]) - as.numeric(project_data[i, 15:17]))
  
  # Cefotaxime
  p_cef <- t.test(project_data[i, 6:8], project_data[i, 15:17], var.equal = TRUE)
  lfq_cef <-  mean(as.numeric(project_data[i, 6:8]) - as.numeric(project_data[i, 15:17]))
  
  # Imipenem
  p_imp <- t.test(project_data[i, 9:11], project_data[i, 15:17], var.equal = TRUE)
  lfq_imp <-  mean(as.numeric(project_data[i, 9:11]) - as.numeric(project_data[i, 15:17]))
  
  # Ciprofloxacin
  p_cip <- t.test(project_data[i, 12:14], project_data[i, 15:17], var.equal = TRUE)
  lfq_cip <-  mean(as.numeric(project_data[i, 12:14]) - as.numeric(project_data[i, 15:17]))
  
  # Assigning Ampicillin values
  ampicillin_values[i, 1] <- project_data$protein__name[i]
  ampicillin_values[i, 2] <- project_data$gene__name[i]
  ampicillin_values[i, 3] <- -log10(p_amp$p.value)
  ampicillin_values[i, 4] <- lfq_amp
  
  # Assigning Cefotaxime values
  cefotaxime_values[i, 1] <- project_data$protein__name[i]
  cefotaxime_values[i, 2] <- project_data$gene__name[i]
  cefotaxime_values[i, 3] <- -log10(p_cef$p.value)
  cefotaxime_values[i, 4] <- lfq_cef
  
  # Assigning Imipenem values
  imipenem_values[i, 1] <- project_data$protein__name[i]
  imipenem_values[i, 2] <- project_data$gene__name[i]
  imipenem_values[i, 3] <- -log10(p_imp$p.value)
  imipenem_values[i, 4] <- lfq_imp
  
  # Assigning Ciprofloxacin values
  ciprofloxacin_values[i, 1] <- project_data$protein__name[i]
  ciprofloxacin_values[i, 2] <- project_data$gene__name[i]
  ciprofloxacin_values[i, 3] <- -log10(p_cip$p.value)
  ciprofloxacin_values[i, 4] <- lfq_cip
}


# Marking up & down regulation ------------------------------------------------
ampicillin_values <- ampicillin_values %>%
  mutate(
    Expression = case_when(LFQ >= 0.05 & P_value >= -log10(.05) ~ "Up Regulated",
                           LFQ <= .05 & P_value >= -log10(.05) ~ "Down Regulated",
                           TRUE ~ "Unchanged")
  )
cefotaxime_values <- cefotaxime_values %>%
  mutate(
    Expression = case_when(LFQ >= 0.05 & P_value >= -log10(.05) ~ "Up Regulated",
                           LFQ <= .05 & P_value >= -log10(.05) ~ "Down Regulated",
                           TRUE ~ "Unchanged")
  )
imipenem_values <- imipenem_values %>%
  mutate(
    Expression = case_when(LFQ >= 0.05 & P_value >= -log10(.05) ~ "Up Regulated",
                           LFQ <= .05 & P_value >= -log10(.05) ~ "Down Regulated",
                           TRUE ~ "Unchanged")
  )
ciprofloxacin_values <- ciprofloxacin_values %>%
  mutate(
    Expression = case_when(LFQ >= 0.05 & P_value >= -log10(.05) ~ "Up Regulated",
                           LFQ <= .05 & P_value >= -log10(.05) ~ "Down Regulated",
                           TRUE ~ "Unchanged")
  )


# Labelling differentially abundant proteins ----------------------------------
top_ampicillin <- bind_rows(
  ampicillin_values %>%
    filter(Expression == "Up Regulated") %>%
    arrange(desc(LFQ)) %>%
    head(10),
  ampicillin_values %>%
    filter(Expression == "Down Regulated") %>%
    arrange(LFQ) %>%
    head(10)
)
top_cefotaxime <- bind_rows(
  cefotaxime_values %>%
    filter(Expression == "Up Regulated") %>%
    arrange(desc(LFQ)) %>%
    head(10),
  cefotaxime_values %>%
    filter(Expression == "Down Regulated") %>%
    arrange(LFQ) %>%
    head(10)
)
top_imipenem <- bind_rows(
  imipenem_values %>%
    filter(Expression == "Up Regulated") %>%
    arrange(desc(LFQ)) %>%
    head(10),
  imipenem_values %>%
    filter(Expression == "Down Regulated") %>%
    arrange(LFQ) %>%
    head(10)
)
top_ciprofloxacin <- bind_rows(
  ciprofloxacin_values %>%
    filter(Expression == "Up Regulated") %>%
    arrange(desc(LFQ)) %>%
    head(10),
  ciprofloxacin_values %>%
    filter(Expression == "Down Regulated") %>%
    arrange(LFQ) %>%
    head(10)
)


# Generating volcano plots ----------------------------------------------------
# Ampicillin plot
ampicillin_plot <- ggplot(data = ampicillin_values, 
                          aes(x = LFQ, 
                              y = P_value,
                              col = Expression)) +
  geom_point(size = 2/5) +
  ggtitle("Control & Ampicillin") +
  xlab(expression("log"[2]*" Difference")) + 
  ylab(expression("-log"[10]* "p")) +
  geom_vline(xintercept=c(-0.5, 0.5), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") +
  theme_bw()

# Cefotaxime plot
cefotaxime_plot <- ggplot(data = cefotaxime_values, 
                          aes(x = LFQ, 
                              y = P_value,
                              col = Expression)) +
  geom_point(size = 2/5) +
  ggtitle("Control & Cefotaxime") +
  xlab(expression("log"[2]*" Difference")) + 
  ylab(expression("-log"[10]* "p")) +
  geom_vline(xintercept=c(-0.5, 0.5), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") +
  theme_bw()

# Imipenem plot
imipenem_plot <- ggplot(data = imipenem_values, 
                        aes(x = LFQ, 
                            y = P_value,
                            col = Expression)) +
  geom_point(size = 2/5) +
  ggtitle("Control & Imipenem") +
  xlab(expression("log"[2]*" Difference")) + 
  ylab(expression("-log"[10]* "p")) +
  geom_vline(xintercept=c(-0.5, 0.5), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") +
  theme_bw()

# Ciprofloxacin plot
ciprofloxacin_plot <- ggplot(data = ciprofloxacin_values, 
                             aes(x = LFQ,
                                 y = P_value,
                                 col = Expression)) +
  geom_point(size = 2/5) +
  ggtitle("Control & Ciprofloxacin") +
  xlab(expression("log"[2]*" Difference")) + 
  ylab(expression("-log"[10]* "p")) +
  geom_vline(xintercept=c(-0.5, 0.5), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") +
  theme_bw()

# Adding labels to plots -------------------------------------------
ampicillin_plot <- ampicillin_plot +
  geom_label_repel(data = top_ampicillin,
                   mapping = aes(x = LFQ, y = P_value, label = Gene_name),
                   size = 5)
cefotaxime_plot <- cefotaxime_plot +
  geom_label_repel(data = top_cefotaxime,
                   mapping = aes(x = LFQ, y = P_value, label = Gene_name),
                   size = 5)
imipenem_plot <- imipenem_plot +
  geom_label_repel(data = top_imipenem,
                   mapping = aes(x = LFQ, y = P_value, label = Gene_name),
                   size = 5)
ciprofloxacin_plot <- ciprofloxacin_plot +
  geom_label_repel(data = top_ciprofloxacin,
                   mapping = aes(x = LFQ, y = P_value, label = Gene_name),
                   size = 5)
# Exporting ------------------------------------------------------------------
export_path <- "./03__modeling/2022-09-30__Multiple_Comparisons/models"

save_plot(file.path(export_path, "results/Control & Ampicillin.png"),
          ampicillin_plot, base_height = 8, base_aspect_ratio = 1.4)
save_plot(file.path(export_path, "results/Control & Cefotaxime.png"),
          cefotaxime_plot, base_height = 8, base_aspect_ratio = 1.4)
save_plot(file.path(export_path, "results/Control & Imipenem.png"),
          imipenem_plot, base_height = 8, base_aspect_ratio = 1.4)
save_plot(file.path(export_path, "results/Control & Ciprofloxacin.png"),
          ciprofloxacin_plot, base_height = 8, base_aspect_ratio = 1.4)


