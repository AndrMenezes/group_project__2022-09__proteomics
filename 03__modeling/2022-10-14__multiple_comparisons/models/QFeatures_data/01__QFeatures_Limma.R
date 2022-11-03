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

# Integrating with QFeatures --------------------------------------------------
if (all.equal(rownames(fts[["proteins"]]), rownames(limma_data)))
  rowData(fts[["proteins"]]) <- cbind(rowData(fts[["proteins"]]), limma_data)

# Exporting data
saveRDS(object = fts, file = file.path(path_data, "fts_processed.rds"))
