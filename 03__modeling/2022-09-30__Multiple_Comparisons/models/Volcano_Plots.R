# 2022-09-30 Replicating T-Tests and Volcano Plots
library(tidyverse)
library(readxl)
library(openxlsx)
library(ggplot2)
library(ggrepel)
library(cowplot)

path_data <- "./01__database/"

project_data <- read_csv(
  file.path(path_data, "processed_data/data_filter.csv"))

# Creating data frames for storing output -------------------------------------
ampicillin_values <- data.frame(Protein_id = character(0),
                                Gene_name = character(0),
                                Statistic = numeric(0),
                                CI_L = numeric(0),
                                CI_R = numeric(0),
                                P_value = numeric(0),
                                LFQ = numeric(0))

cefotaxime_values <- data.frame(Protein_id = character(0),
                                Gene_name = character(0),
                                Statistic = numeric(0),
                                CI_L = numeric(0),
                                CI_R = numeric(0),
                                P_value = numeric(0),
                                LFQ = numeric(0))

imipenem_values <- data.frame(Protein_id = character(0),
                              Gene_name = character(0),
                              Statistic = numeric(0),
                              CI_L = numeric(0),
                              CI_R = numeric(0),
                              P_value = numeric(0),
                              LFQ = numeric(0))

ciprofloxacin_values <- data.frame(Protein_id = character(0),
                                   Gene_name = character(0),
                                   Statistic = numeric(0),
                                   CI_L = numeric(0),
                                   CI_R = numeric(0),
                                   P_value = numeric(0),
                                   LFQ = numeric(0))

# Conducting T-Tests ----------------------------------------------------------
for(i in 1:nrow(project_data)){
  
  # Ampicillin
  p_amp <- t.test(project_data[i, 3:5], project_data[i, 15:17], 
                  var.equal = TRUE, conf.level = 0.95)
  lfq_amp <- mean(as.numeric(project_data[i, 3:5]) - as.numeric(project_data[i, 15:17]))
  
  # Cefotaxime
  p_cef <- t.test(project_data[i, 6:8], project_data[i, 15:17], 
                  var.equal = TRUE, conf.level = 0.95)
  lfq_cef <-  mean(as.numeric(project_data[i, 6:8]) - as.numeric(project_data[i, 15:17]))
  
  # Impipenem
  p_imp <- t.test(project_data[i, 9:11], project_data[i, 15:17], 
                  var.equal = TRUE, conf.level = 0.95)
  lfq_imp <-  mean(as.numeric(project_data[i, 9:11]) - as.numeric(project_data[i, 15:17]))
  
  # Ciprofloxacin
  p_cip <- t.test(project_data[i, 12:14], project_data[i, 15:17], 
                  var.equal = TRUE, conf.level = 0.95)
  lfq_cip <-  mean(as.numeric(project_data[i, 12:14]) - as.numeric(project_data[i, 15:17]))
  
  # Assigning Ampicillin values
  ampicillin_values[i, 1] <- project_data$protein__id[i]
  ampicillin_values[i, 2] <- project_data$gene__name[i]
  ampicillin_values[i, 3] <- p_amp$statistic
  ampicillin_values[i, 4] <- p_amp$conf.int[1]
  ampicillin_values[i, 5] <- p_amp$conf.int[2]
  ampicillin_values[i, 6] <- p_amp$p.value
  ampicillin_values[i, 7] <- lfq_amp
  
  # Assigning Cefotaxime values
  cefotaxime_values[i, 1] <- project_data$protein__id[i]
  cefotaxime_values[i, 2] <- project_data$gene__name[i]
  cefotaxime_values[i, 3] <- p_cef$statistic
  cefotaxime_values[i, 4] <- p_cef$conf.int[1]
  cefotaxime_values[i, 5] <- p_cef$conf.int[2]
  cefotaxime_values[i, 6] <- p_cef$p.value
  cefotaxime_values[i, 7] <- lfq_cef
  
  # Assigning Impipenem values
  imipenem_values[i, 1] <- project_data$protein__id[i]
  imipenem_values[i, 2] <- project_data$gene__name[i]
  imipenem_values[i, 3] <- p_imp$statistic
  imipenem_values[i, 4] <- p_imp$conf.int[1]
  imipenem_values[i, 5] <- p_imp$conf.int[2]
  imipenem_values[i, 6] <- p_imp$p.value
  imipenem_values[i, 7] <- lfq_imp
  
  # Assigning Ciprofloxacin values
  ciprofloxacin_values[i, 1] <- project_data$protein__id[i]
  ciprofloxacin_values[i, 2] <- project_data$gene__name[i]
  ciprofloxacin_values[i, 3] <- p_cip$statistic
  ciprofloxacin_values[i, 4] <- p_cip$conf.int[1]
  ciprofloxacin_values[i, 5] <- p_cip$conf.int[2]
  ciprofloxacin_values[i, 6] <- p_cip$p.value
  ciprofloxacin_values[i, 7] <- lfq_cip
}

# Marking up & down regulation & adding in treatment --------------------------
ampicillin_values <- ampicillin_values %>%
  mutate(
    Expression = case_when(LFQ >= 0.5 & -log10(P_value) >= -log10(.05) ~ "Up Regulated",
                           LFQ <= -0.5 & -log10(P_value) >= -log10(.05) ~ "Down Regulated",
                           TRUE ~ "Unchanged"),
    Control = "Control", Treatment = "Ampicillin"
  )
cefotaxime_values <- cefotaxime_values %>%
  mutate(
    Expression = case_when(LFQ >= 0.5 & -log10(P_value) >= -log10(.05) ~ "Up Regulated",
                           LFQ <= -0.5 & -log10(P_value) >= -log10(.05) ~ "Down Regulated",
                           TRUE ~ "Unchanged"),
    Control = "Control", Treatment = "Cefotaxime"
  )
imipenem_values <- imipenem_values %>%
  mutate(
    Expression = case_when(LFQ >= 0.5 & -log10(P_value) >= -log10(.05) ~ "Up Regulated",
                           LFQ <= -0.5 & -log10(P_value) >= -log10(.05) ~ "Down Regulated",
                           TRUE ~ "Unchanged"),
    Control = "Control", Treatment = "Impipenem"
  )
ciprofloxacin_values <- ciprofloxacin_values %>%
  mutate(
    Expression = case_when(LFQ >= 0.5 & -log10(P_value) >= -log10(.05) ~ "Up Regulated",
                           LFQ <= -0.5 & -log10(P_value) >= -log10(.05) ~ "Down Regulated",
                           TRUE ~ "Unchanged"),
    Control = "Control", Treatment = "Ciprofloxacin"
  )

# Computing adjusted p value --------------------------------------------------
FDR_amp <- p.adjust(ampicillin_values$P_value, method = "BH")
FDR_cef <- p.adjust(cefotaxime_values$P_value, method = "BH")
FDR_imp <- p.adjust(imipenem_values$P_value, method = "BH")
FDR_cip <- p.adjust(ciprofloxacin_values$P_value, method = "BH")

# Adding to data
ampicillin_values$FDR <- FDR_amp
cefotaxime_values$FDR <- FDR_cef
imipenem_values$FDR <- FDR_imp
ciprofloxacin_values$FDR <- FDR_cip

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

# Plotting --------------------------------------------------------------------
ampicillin_plot <- ggplot(data = ampicillin_values, 
                          aes(x = LFQ, 
                              y = -log10(P_value),
                              col = Expression)) +
  geom_point(size = 2/5) +
  ggtitle("Control Vs. Ampicillin") +
  xlab(expression("log"[2]*" Difference")) + 
  ylab(expression("-log"[10]* "p")) +
  scale_colour_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
  geom_vline(xintercept=c(-0.5, 0.5), linetype = "dashed", col = "black") +
  geom_hline(yintercept=-log10(0.05), linetype = "dashed", col = "black") +
  theme_bw()

cefotaxime_plot <- ggplot(data = cefotaxime_values, 
                          aes(x = LFQ, 
                              y = -log10(P_value),
                              col = Expression)) +
  geom_point(size = 2/5) +
  ggtitle("Control Vs.Cefotaxime") +
  xlab(expression("log"[2]*" Difference")) + 
  ylab(expression("-log"[10]* "p")) +
  scale_colour_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
  geom_vline(xintercept=c(-0.5, 0.5), linetype = "dashed", col = "black") +
  geom_hline(yintercept=-log10(0.05), linetype = "dashed", col = "black") +
  theme_bw()

imipenem_plot <- ggplot(data = imipenem_values, 
                        aes(x = LFQ, 
                            y = -log10(P_value),
                            col = Expression)) +
  geom_point(size = 2/5) +
  ggtitle("Control Vs.Imipenem") +
  xlab(expression("log"[2]*" Difference")) + 
  ylab(expression("-log"[10]* "p")) +
  scale_colour_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
  geom_vline(xintercept=c(-0.5, 0.5), linetype = "dashed", col = "black") +
  geom_hline(yintercept=-log10(0.05), linetype = "dashed", col = "black") +
  theme_bw()

ciprofloxacin_plot <- ggplot(data = ciprofloxacin_values, 
                             aes(x = LFQ,
                                 y = -log10(P_value),
                                 col = Expression)) +
  geom_point(size = 2/5) +
  ggtitle("Control Vs. Ciprofloxacin") +
  xlab(expression("log"[2]*" Difference")) + 
  ylab(expression("-log"[10]* "p")) +
  scale_colour_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
  geom_vline(xintercept=c(-0.5, 0.5), linetype = "dashed", col = "black") +
  geom_hline(yintercept=-log10(0.05), linetype = "dashed", col = "black") +
  theme_bw()

# Adding labels ---------------------------------------------------------------
ampicillin_plot <- ampicillin_plot +
  geom_text_repel(data = top_ampicillin,
                  mapping = aes(x = LFQ, y = -log10(P_value), label = Gene_name),
                  size = 2)
cefotaxime_plot <- cefotaxime_plot +
  geom_text_repel(data = top_cefotaxime,
                  mapping = aes(x = LFQ, y = -log10(P_value), label = Gene_name),
                  size = 2)
imipenem_plot <- imipenem_plot +
  geom_text_repel(data = top_imipenem,
                  mapping = aes(x = LFQ, y = -log10(P_value), label = Gene_name),
                  size = 2)
ciprofloxacin_plot <- ciprofloxacin_plot +
  geom_text_repel(data = top_ciprofloxacin,
                  mapping = aes(x = LFQ, y = -log10(P_value), label = Gene_name),
                  size = 2)

# Creating a combined data set -----------------------------------------------
all_proteins <- rbind(ampicillin_values, cefotaxime_values, 
                      imipenem_values, ciprofloxacin_values)

all_proteins <- all_proteins |>
  select(Control, Treatment, everything())

# Exporting -------------------------------------------------------------------
export_path <- "./03__modeling/2022-09-30__Multiple_Comparisons/models"

# Data
results <- list("All proteins" = all_proteins,
                "Ampicillin" = ampicillin_values,
                "Cefotaxime" = cefotaxime_values,
                "Impipenem" = imipenem_values,
                "Ciprofloxacin" = ciprofloxacin_values)

write.xlsx(results, file.path(export_path, "results/replicated_results.xlsx"))

save_plot(file.path(export_path, "results/Control & Ampicillin.png"),
          ampicillin_plot, base_height = 8, base_aspect_ratio = 1.4)
save_plot(file.path(export_path, "results/Control & Cefotaxime.png"),
          cefotaxime_plot, base_height = 8, base_aspect_ratio = 1.4)
save_plot(file.path(export_path, "results/Control & Imipenem.png"),
          imipenem_plot, base_height = 8, base_aspect_ratio = 1.4)
save_plot(file.path(export_path, "results/Control & Ciprofloxacin.png"),
          ciprofloxacin_plot, base_height = 8, base_aspect_ratio = 1.4)