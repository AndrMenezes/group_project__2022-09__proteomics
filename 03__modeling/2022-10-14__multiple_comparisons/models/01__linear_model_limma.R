library(dplyr)
library(tidyr)
library(openxlsx)
library(limma)
library(ggplot2)
library(ggrepel)
library(cowplot)

path_data <- "./01__database/"
project_data <- read.csv(
  file.path(path_data, "processed_data/data_pivot.csv"))

# Data manipulation -----------------------------------------------------------
# Pivoting data: rows = proteins // columns = samples
project_data <- project_data |>
  arrange(group) |>
  select(-c("variable")) |>
  pivot_wider(names_from = c(group, replicate),
              values_from = value)

# Saving protein ID (use for row names later)
id <- project_data$protein__id

# Removing protein ID column
project_data <- select(project_data, -c("protein__id"))

# Reordering columns
project_data <- project_data[, c(10, 11, 12, 1, 2, 3, 4, 5, 6, 
                                 13, 14, 15, 7, 8, 9)]

# Converting to matrix
project_data <- as.matrix(project_data)
rownames(project_data) <- id
colnames(project_data) <- NULL

# Creating design matrix ------------------------------------------------------
x <- factor(rep(c("Control", "Ampicillin", "Cefotaxime", "Impipenem",
                  "Ciprofloxacin"), each = 3),
            levels = c("Control", "Ampicillin",  "Cefotaxime", "Impipenem",
                       "Ciprofloxacin"))

design_matrix <- model.matrix(~ x)
colnames(design_matrix) <- levels(x)

# Fitting linear model --------------------------------------------------------
fit <- lmFit(project_data, design_matrix)
fit <- eBayes(fit)

# Storing limma results -------------------------------------------------------
ampicillin_control <- topTable(fit, coef = "Ampicillin", number = Inf,
                               sort.by = "none", confint = TRUE)
cefotaxime_control <- topTable(fit, coef = "Cefotaxime", number = Inf,
                               sort.by = "none", confint = TRUE)
impipenem_control <- topTable(fit, coef = "Impipenem", number = Inf,
                              sort.by = "none", confint = TRUE)
ciprofloxacin_control <- topTable(fit, coef = "Ciprofloxacin", number = Inf,
                                  sort.by = "none", confint = TRUE)

# Plotting --------------------------------------------------------------------
ampicillin_plot <- ggplot(ampicillin_control,
                          aes(x = logFC, y = -log10(P.Value))) +
  geom_point(size = .005) +
  ggtitle("Control vs. Ampicillin") +
  xlab(expression("log"[2]*" Fold Change")) + 
  ylab(expression("-log"[10]* "p")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", col = "black") +
  geom_hline(yintercept = -log10(0.05),linetype = "dashed", col = "black") +
  theme_bw()

cefotaxime_plot <- ggplot(cefotaxime_control,
                          aes(x = logFC, y = -log10(P.Value))) +
  geom_point(size = .005) +
  ggtitle("Control vs. Cefotaxime") +
  xlab(expression("log"[2]*" Fold Change")) + 
  ylab(expression("-log"[10]* "p")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", col = "black") +
  geom_hline(yintercept = -log10(0.05),linetype = "dashed", col = "black") +
  theme_bw()

impipenem_plot <- ggplot(impipenem_control,
                         aes(x = logFC, y = -log10(P.Value))) +
  geom_point(size = .005) +
  ggtitle("Control vs. Impipenem") +
  xlab(expression("log"[2]*" Fold Change")) + 
  ylab(expression("-log"[10]* "p")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", col = "black") +
  geom_hline(yintercept = -log10(0.05),linetype = "dashed", col = "black") +
  theme_bw()

ciprofloxacin_plot <- ggplot(ciprofloxacin_control,
                             aes(x = logFC, y = -log10(P.Value))) +
  geom_point(size = .005) +
  ggtitle("Control vs. Ciprofloxacin") +
  xlab(expression("log"[2]*" Fold Change")) + 
  ylab(expression("-log"[10]* "p")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", col = "black") +
  geom_hline(yintercept = -log10(0.05),linetype = "dashed", col = "black") +
  theme_bw()

# Labelling -------------------------------------------------------------------
# Selcting top 20 proteins for each treatment
top_ampicillin <- topTable(fit, coef = "Ampicillin", number = 20)
top_cefotaxime <- topTable(fit, coef = "Cefotaxime", number = 20)
top_impipenem <- topTable(fit, coef = "Impipenem", number = 20)
top_ciprofloxacin <- topTable(fit, coef = "Ciprofloxacin", number = 20)

# Adding labels
ampicillin_plot <- ampicillin_plot +
  geom_text_repel(data = top_ampicillin,
                  mapping = aes(x = logFC, 
                                y = -log10(P.Value), 
                                label = rownames(top_ampicillin)),
                  size = 2,
                  force = 2)

cefotaxime_plot <- cefotaxime_plot +
  geom_text_repel(data = top_cefotaxime,
                  mapping = aes(x = logFC, 
                                y = -log10(P.Value), 
                                label = rownames(top_cefotaxime)),
                  size = 2,
                  force = 3)

impipenem_plot <- impipenem_plot +
  geom_text_repel(data = top_impipenem,
                  mapping = aes(x = logFC, 
                                y = -log10(P.Value), 
                                label = rownames(top_impipenem)),
                  size = 2,
                  force = 2)

ciprofloxacin_plot <- ciprofloxacin_plot +
  geom_text_repel(data = top_ciprofloxacin,
                  mapping = aes(x = logFC, 
                                y = -log10(P.Value), 
                                label = rownames(top_ciprofloxacin)),
                  size = 2,
                  force = 2)

# Exporting -------------------------------------------------------------------
export_path <- "./03__modeling/2022-10-14__multiple_comparisons/models"

# Exporting data
sheets <- list("Ampicillin" = ampicillin_control,
               "Cefotaxime" = cefotaxime_control,
               "Impipenem" = impipenem_control,
               "Ciprofloxacin" = ciprofloxacin_control)

top_proteins <- list("Ampicillin" = top_ampicillin,
                     "Cefotaxime" = top_cefotaxime,
                     "Impipenem" = top_impipenem,
                     "Ciprofloxacin" = top_ciprofloxacin)

write.xlsx(sheets, file.path(export_path, "results/limma_results.xlsx"),
           rowNames = TRUE, colNames = TRUE)
write.xlsx(top_proteins, file.path(export_path, "results/top_proteins.xlsx"),
           rowNames = TRUE, colNames = TRUE)

# Exporting plots
save_plot(file.path(export_path, "results/Ampicillin_Control.png"),
          ampicillin_plot)
save_plot(file.path(export_path, "results/Cefotaxime_Control.png"),
          cefotaxime_plot)
save_plot(file.path(export_path, "results/Impipenem_Control.png"),
          impipenem_plot)
save_plot(file.path(export_path, "results/Ciprofloxacin_Control.png"),
          ciprofloxacin_plot)

