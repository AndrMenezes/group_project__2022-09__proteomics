# 2022-11-02 Using Limma on the QFeatures Log Intensity Normalized Data
suppressMessages(library(QFeatures))
suppressMessages(library(dplyr))
library(tidyr)
library(limma)
library(ggplot2)
library(ggrepel)
library(cowplot)

path_data <- "./03__modeling/2022-10-14__workflow/processing/processed_data"

# Importing data --------------------------------------------------------------
# First we import the QFeatures Data (all of it)
fts <- readRDS(file.path(path_data, "fts_processed.rds"))

# We can then focus in on the normalized protein values
normalised_proteins <- assay(x = fts[["proteins"]], 
                             i = "log_intensity_normalized")
colnames(normalised_proteins) <- fts@colData@listData[["group"]]

# Reordering columns so control values are first
normalised_proteins <- normalised_proteins[, c(10, 11, 12, 1, 2, 3, 4, 
                                               5, 6, 13, 14, 15, 7, 8, 9)]

# Creating design matrix ------------------------------------------------------
x <- factor(rep(c("Control", "Ampicillin", "Cefotaxime", "Impipenem",
                  "Ciprofloxacin"), each = 3),
            levels = c("Control", "Ampicillin",  "Cefotaxime", "Impipenem",
                       "Ciprofloxacin"))

design_matrix <- model.matrix(~ x)
colnames(design_matrix) <- levels(x)

# Fitting linear model --------------------------------------------------------
fit <- lmFit(normalised_proteins, design_matrix)
fit <- eBayes(fit)

# Storing limma results -------------------------------------------------------
ampicillin_control <- topTable(fit, coef = "Ampicillin", number = Inf,
                               sort.by = "none")
cefotaxime_control <- topTable(fit, coef = "Cefotaxime", number = Inf,
                               sort.by = "none")
impipenem_control <- topTable(fit, coef = "Impipenem", number = Inf,
                              sort.by = "none")
ciprofloxacin_control <- topTable(fit, coef = "Ciprofloxacin", number = Inf,
                                  sort.by = "none")

# Selecting necessary columns and renaming them
amp <- ampicillin_control |>
  select(logFC, P.Value, adj.P.Val) |>
  rename("log_FC__ampicillin" = logFC, "P_Value__ampicillin" = P.Value, 
         "FDR__ampicillin" = adj.P.Val)

cef <- cefotaxime_control |>
  select(logFC, P.Value, adj.P.Val) |>
  rename("log_FC__cefotaxime" = logFC, "P_Value__cefotaxime" = P.Value, 
         "FDR__cefotaxime" = adj.P.Val)

imp <- impipenem_control |>
  select(logFC, P.Value, adj.P.Val) |>
  rename("log_FC__impipenem" = logFC, "P_Value__impipenem" = P.Value, 
         "FDR__impipenem" = adj.P.Val)

cip <- ciprofloxacin_control |>
  select(logFC, P.Value, adj.P.Val) |>
  rename("log_FC__ciprofloxacin" = logFC, "P_Value__ciprofloxacin" = P.Value, 
         "FDR__ciprofloxacin" = adj.P.Val)

# Creating data frame (use for QFeatures later)
limma_data <- data.frame(amp, cef, imp, cip)

# Filtering based on Hedge's G ------------------------------------------------
# We can use the previously calculated Hedge's G effect size to select and plot
# only the most variable proteins

# Read in the effect size data
effect_size <- data.frame(rowData(fts[["proteins"]]))

# Removing unnecessary columns and pivoting to long format
effect_size_long <- effect_size |>
  select(-c(.n)) |>
  pivot_longer(cols = -protein, names_to = "hedges_g", 
               values_to = "value") |>
  mutate(hedges_g = gsub("hedges_g__", "", hedges_g))

# Removing absolute values less than 0.5 
effect_size_filter <- effect_size_long |>
  filter(abs(value) > 0.5)

# Getting the index position of rows where the protein names are the same
x <- which(row.names(limma_data) %in% unique(effect_size_filter$protein))

# Only focusing on these indexed positions (these proteins had greater than > 0.5)
limma_data_reduced <- limma_data[x, ]


# Plotting --------------------------------------------------------------------
ampicillin_plot <- ggplot(limma_data_reduced,
                          aes(x = log_FC__ampicillin, 
                              y = -log10(P_Value__ampicillin))) +
  geom_point(size = .05) +
  ggtitle("Control vs. Ampicillin") +
  xlab(expression("log"[2]*" Fold Change")) + 
  ylab(expression("-log"[10]* "p")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", col = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", col = "black") +
  theme_bw()

cefotaxime_plot <- ggplot(limma_data_reduced,
                          aes(x = log_FC__cefotaxime, 
                              y = -log10(P_Value__cefotaxime))) +
  geom_point(size = .05) +
  ggtitle("Control vs. Cefotaxime") +
  xlab(expression("log"[2]*" Fold Change")) + 
  ylab(expression("-log"[10]* "p")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", col = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", col = "black") +
  theme_bw()

impipenem_plot <- ggplot(limma_data_reduced,
                         aes(x = log_FC__impipenem, 
                             y = -log10(P_Value__impipenem))) +
  geom_point(size = .05) +
  ggtitle("Control vs. Impipenem") +
  xlab(expression("log"[2]*" Fold Change")) + 
  ylab(expression("-log"[10]* "p")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", col = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", col = "black") +
  theme_bw()

ciprofloxacin_plot <- ggplot(limma_data_reduced,
                             aes(x = log_FC__ciprofloxacin, 
                                 y = -log10(P_Value__ciprofloxacin))) +
  geom_point(size = .05) +
  ggtitle("Control vs. Ciprofloxacin") +
  xlab(expression("log"[2]*" Fold Change")) + 
  ylab(expression("-log"[10]* "p")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", col = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", col = "black") +
  theme_bw()

# Adding labels ---------------------------------------------------------------
# Locating top 20 proteins
top_ampicillin <- topTable(fit, coef = "Ampicillin", number = 20, confint = TRUE)
top_cefotaxime <- topTable(fit, coef = "Cefotaxime", number = 20, confint = TRUE)
top_impipenem <- topTable(fit, coef = "Impipenem", number = 20, confint = TRUE)
top_ciprofloxacin <- topTable(fit, coef = "Ciprofloxacin", number = 20, confint = TRUE)

# Adding labels to plots
ampicillin_plot <- ampicillin_plot +
  geom_text_repel(data = top_ampicillin,
                  mapping = aes(x = logFC, y = -log10(P.Value), 
                                label = rownames(top_ampicillin)),
                  max.overlaps = Inf,
                  size = 2, force = 5)

cefotaxime_plot <- cefotaxime_plot +
  geom_text_repel(data = top_cefotaxime,
                  mapping = aes(x = logFC, y = -log10(P.Value), 
                                label = rownames(top_cefotaxime)),
                  max.overlaps = Inf,
                  size = 2, force = 5)

impipenem_plot <- impipenem_plot +
  geom_text_repel(data = top_impipenem,
                  mapping = aes(x = logFC, y = -log10(P.Value), 
                                label = rownames(top_impipenem)),
                  max.overlaps = Inf,
                  size = 2, force = 5)

ciprofloxacin_plot <- ciprofloxacin_plot +
  geom_text_repel(data = top_ciprofloxacin,
                  mapping = aes(x = logFC, y = -log10(P.Value), 
                                label = rownames(top_ciprofloxacin)), 
                  max.overlaps = Inf,
                  size = 2, force = 5)

# Integrating with QFeatures --------------------------------------------------
if (all.equal(rownames(fts[["proteins"]]), rownames(limma_data)))
  rowData(fts[["proteins"]]) <- cbind(rowData(fts[["proteins"]]), limma_data)

# Exporting -------------------------------------------------------------------
export_path <- "./03__modeling/2022-10-14__multiple_comparisons/models/QFeatures_data"

# Exporting data
saveRDS(object = fts, file = file.path(path_data, "fts_processed.rds"))

# Exporting plots
save_plot(file.path(export_path, "results/Ampicillin_Control.png"),
          ampicillin_plot)
save_plot(file.path(export_path, "results/Cefotaxime_Control.png"),
          cefotaxime_plot)
save_plot(file.path(export_path, "results/Impipenem_Control.png"),
          impipenem_plot)
save_plot(file.path(export_path, "results/Ciprofloxacin_Control.png"),
          ciprofloxacin_plot)
