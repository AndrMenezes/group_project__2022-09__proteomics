library(dplyr)
library(tidyr)
library(openxlsx)
library(limma)
library(ggplot2)
library(cowplot)

path_data <- "./01__database/"
project_data <- read.csv(
  file.path(path_data, "processed_data/data_pivot.csv"))

# Data manipulation -----------------------------------------------------------
# Creating control vs. treatments // removing unnecessary columns
# Pivoting data: Rows = proteins // Columns = samples
ampicillin_control <- project_data |>
  filter(group == "Ampicillin" | group == "Control") |>
  select(-c("variable")) |>
  pivot_wider(names_from = c(group, replicate),
              values_from = value)

cefotaxime_control <- project_data |>
  filter(group == "Cefotaxime" | group == "Control") |>
  select(-c("variable")) |>
  pivot_wider(names_from = c(group, replicate),
              values_from = value)

impipenem_control <- project_data |>
  filter(group == "Impipenem" | group == "Control") |>
  select(-c("variable")) |>
  pivot_wider(names_from = c(group, replicate),
              values_from = value)

ciprofloxacin_control <- project_data |>
  filter(group == "Ciprofloxacin" | group == "Control") |>
  select(-c("variable")) |>
  pivot_wider(names_from = c(group, replicate),
              values_from = value)

# Saving protein IDs (use for row names later)
id <- ampicillin_control$protein__id

# Removing protein ID column
ampicillin_control <- select(ampicillin_control, -c("protein__id"))
cefotaxime_control <- select(cefotaxime_control, -c("protein__id"))
impipenem_control <- select(impipenem_control, -c("protein__id"))
ciprofloxacin_control <- select(ciprofloxacin_control, -c("protein__id"))

# Reordering columns: Control // Treatment
ampicillin_control <- ampicillin_control[, c(4, 5, 6, 1, 2, 3)]
cefotaxime_control <- cefotaxime_control[, c(4, 5, 6, 1, 2, 3)]
impipenem_control <- impipenem_control[, c(4, 5, 6, 1, 2, 3)]
ciprofloxacin_control <- ciprofloxacin_control[, c(4, 5, 6, 1, 2, 3)]

# Turning into matrix with just row names
ampicillin_control <- as.matrix(ampicillin_control)
rownames(ampicillin_control) <- id
colnames(ampicillin_control) <- NULL

cefotaxime_control <- as.matrix(cefotaxime_control)
rownames(cefotaxime_control) <- id
colnames(cefotaxime_control) <- NULL

impipenem_control <- as.matrix(impipenem_control)
rownames(impipenem_control) <- id
colnames(impipenem_control) <- NULL

ciprofloxacin_control <- as.matrix(ciprofloxacin_control)
rownames(ciprofloxacin_control) <- id
colnames(ciprofloxacin_control) <- NULL

# Creating Design Matrix ------------------------------------------------------
control_treatment <- factor(x = c(rep("Control", 3),
                                  rep("Treatment", 3)),
                            levels = c("Control", "Treatment"))

design <- model.matrix(~ control_treatment)
colnames(design) <- levels(control_treatment)

# Fitting Linear Model -------------------------------------------------------
# Ampicillin
ampicillin_fit <- lmFit(ampicillin_control, design)
ampicillin_fit <- eBayes(ampicillin_fit)

# Cefotaxime
cefotaxime_fit <- lmFit(cefotaxime_control, design)
cefotaxime_fit <- eBayes(cefotaxime_fit)

# Impipenem
impipenem_fit <- lmFit(impipenem_control, design)
impipenem_fit <- eBayes(impipenem_fit)

# Ciprofloxacin
ciprofloxacin_fit <- lmFit(ciprofloxacin_control, design)
ciprofloxacin_fit <- eBayes(ciprofloxacin_fit)

# Outputting Limma Results ----------------------------------------------------
ampicillin_results <- topTable(ampicillin_fit, coef = 2, number = Inf, 
                               sort.by = "none", confint = TRUE)
cefotaxime_results <- topTable(cefotaxime_fit, coef = 2, number = Inf, 
                               sort.by = "none", confint = TRUE)
impipenem_results <- topTable(impipenem_fit, coef = 2, number = Inf, 
                              sort.by = "none", confint = TRUE)
ciprofloxacin_results <- topTable(ciprofloxacin_fit, coef = 2, number = Inf, 
                                  sort.by = "none", confint = TRUE)

# Plotting --------------------------------------------------------------------
# Ampicillin
ampicillin_plot <- ggplot(ampicillin_results,
                          aes(x = logFC, y = -log10(P.Value))) +
  geom_point(size = .05) +
  ggtitle("Control vs. Ampicillin") +
  xlab(expression("log"[2]*" Fold Change")) + 
  ylab(expression("-log"[10]* "p")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", col = "black") +
  geom_hline(yintercept = -log10(0.05),linetype = "dashed", col = "black") +
  theme_bw()

# Cefotaxime
cefotaxime_plot <- ggplot(cefotaxime_results,
                          aes(x = logFC, y = -log10(P.Value))) +
  geom_point(size = .05) +
  ggtitle("Control vs. Cefotaxime") +
  xlab(expression("log"[2]*" Fold Change")) + 
  ylab(expression("-log"[10]* "p")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", col = "black") +
  geom_hline(yintercept = -log10(0.05),linetype = "dashed", col = "black") +
  theme_bw()

# Impipenem
impipenem_plot <- ggplot(impipenem_results,
                         aes(x = logFC, y = -log10(P.Value))) +
  geom_point(size = .05) +
  ggtitle("Control vs. Impipenem") +
  xlab(expression("log"[2]*" Fold Change")) + 
  ylab(expression("-log"[10]* "p")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", col = "black") +
  geom_hline(yintercept = -log10(0.05),linetype = "dashed", col = "black") +
  theme_bw()

# Ciprofloxacin
ciprofloxacin_plot <- ggplot(ciprofloxacin_results,
                             aes(x = logFC, y = -log10(P.Value))) +
  geom_point(size = .05) +
  ggtitle("Control vs. Ciprofloxacin") +
  xlab(expression("log"[2]*" Fold Change")) + 
  ylab(expression("-log"[10]* "p")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", col = "black") +
  geom_hline(yintercept = -log10(0.05),linetype = "dashed", col = "black") +
  theme_bw()

# Exporting --------------------------------------------------------------------
export_path <- "./03__modeling/2022-10-14__multiple_comparisons/models"

# Data
sheets <- list("Ampicillin" = ampicillin_results,
               "Cefotaxime" = cefotaxime_results,
               "Impipenem" = impipenem_results,
               "Ciprofloxacin" = ciprofloxacin_results)

write.xlsx(sheets, file.path(export_path, "results/Limma_results.xlsx"),
           rowNames = TRUE, colnames = TRUE)

# Plots
save_plot(file.path(export_path, "results/Ampicillin_Control.png"),
          ampicillin_plot)
save_plot(file.path(export_path, "results/Cefotaxime_Control.png"),
          cefotaxime_plot)
save_plot(file.path(export_path, "results/Impipenem_Control.png"),
          impipenem_plot)
save_plot(file.path(export_path, "results/Ciprofloxacin_Control.png"),
          ciprofloxacin_plot)

