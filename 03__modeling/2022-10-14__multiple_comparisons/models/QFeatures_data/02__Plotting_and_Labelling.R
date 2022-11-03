# 03/11/2022 - Plotting the Limma Results for most variable proteins
suppressMessages(library(QFeatures))
suppressMessages(library(dplyr))
library(ggplot2)
library(ggrepel)
library(cowplot)

path_data <- "./03__modeling/2022-10-14__workflow/processing/processed_data"

fts <- readRDS(file.path(path_data, "fts_processed.rds"))

path_data <- "./03__modeling/2022-10-14__workflow/processing/processed_data"

# Selecting and plotting only the variable proteins ---------------------------
# Focusing on the RowData
all_results <- data.frame(rowData(fts[["proteins"]]))

# Filtering based on mean-varience and effect size
filtered_results <- all_results |>
  filter(bio > 0) |>
  filter(abs(hedges_g__ampicillin) > 0.5 | abs(hedges_g__cefotaxime) > 0.5 |
           abs(hedges_g__ciprofloxacin) > 0.5 | abs(hedges_g__impipenem) > 0.5)

# Selecting the top 20 most differentially abundant proteins per treatment
top_ampicillin <- bind_rows(
  filtered_results %>%
    filter(log_FC__ampicillin >= 0.5 & -log10(P_Value__ampicillin) >= -log10(0.05)) %>%
    arrange(desc(log_FC__ampicillin)) %>%
    head(10),
  filtered_results %>%
    filter(log_FC__ampicillin <= -0.5 & -log10(P_Value__ampicillin) >= -log10(0.05)) %>%
    arrange(log_FC__ampicillin) %>%
    head(10)  
)

top_cefotaxime <- bind_rows(
  filtered_results %>%
    filter(log_FC__cefotaxime >= 0.5 & -log10(P_Value__cefotaxime) >= -log10(0.05)) %>%
    arrange(desc(log_FC__cefotaxime)) %>%
    head(10),
  filtered_results %>%
    filter(log_FC__cefotaxime <= -0.5 & -log10(P_Value__cefotaxime) >= -log10(0.05)) %>%
    arrange(log_FC__cefotaxime) %>%
    head(10)  
)

top_impipenem <- bind_rows(
  filtered_results %>%
    filter(log_FC__impipenem >= 0.5 & -log10(P_Value__impipenem) >= -log10(0.05)) %>%
    arrange(desc(log_FC__impipenem)) %>%
    head(10),
  filtered_results %>%
    filter(log_FC__impipenem<= -0.5 & -log10(P_Value__impipenem) >= -log10(0.05)) %>%
    arrange(log_FC__impipenem) %>%
    head(10)  
)

top_ciprofloxacin <- bind_rows(
  filtered_results %>%
    filter(log_FC__ciprofloxacin >= 0.5 & -log10(P_Value__ciprofloxacin) >= -log10(0.05)) %>%
    arrange(desc(log_FC__ciprofloxacin)) %>%
    head(10),
  filtered_results %>%
    filter(log_FC__ciprofloxacin <= -0.5 & -log10(P_Value__ciprofloxacin) >= -log10(0.05)) %>%
    arrange(log_FC__ciprofloxacin) %>%
    head(10)  
)

# Plotting --------------------------------------------------------------------
ampicillin_plot <- ggplot(filtered_results,
                          aes(x = log_FC__ampicillin, 
                              y = -log10(P_Value__ampicillin))) +
  geom_point(size = .05) +
  ggtitle("Control vs. Ampicillin") +
  xlab(expression("log"[2]*" Fold Change")) + 
  ylab(expression("-log"[10]* "p")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", col = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", col = "black") +
  theme_bw()

cefotaxime_plot <- ggplot(filtered_results,
                          aes(x = log_FC__cefotaxime, 
                              y = -log10(P_Value__cefotaxime))) +
  geom_point(size = .05) +
  ggtitle("Control vs. Cefotaxime") +
  xlab(expression("log"[2]*" Fold Change")) + 
  ylab(expression("-log"[10]* "p")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", col = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", col = "black") +
  theme_bw()

impipenem_plot <- ggplot(filtered_results,
                         aes(x = log_FC__impipenem, 
                             y = -log10(P_Value__impipenem))) +
  geom_point(size = .05) +
  ggtitle("Control vs. Impipenem") +
  xlab(expression("log"[2]*" Fold Change")) + 
  ylab(expression("-log"[10]* "p")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", col = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", col = "black") +
  theme_bw()

ciprofloxacin_plot <- ggplot(filtered_results,
                             aes(x = log_FC__ciprofloxacin, 
                                 y = -log10(P_Value__ciprofloxacin))) +
  geom_point(size = .05) +
  ggtitle("Control vs. Ciprofloxacin") +
  xlab(expression("log"[2]*" Fold Change")) + 
  ylab(expression("-log"[10]* "p")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", col = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", col = "black") +
  theme_bw()

# Adding labels
ampicillin_plot <- ampicillin_plot +
  geom_text_repel(data = top_ampicillin,
                  mapping = aes(x = log_FC__ampicillin, y = -log10(P_Value__ampicillin),
                                label = protein),
                  max.overlaps = Inf, size = 2, force = 4)

cefotaxime_plot <- cefotaxime_plot +
  geom_text_repel(data = top_cefotaxime,
                  mapping = aes(x = log_FC__cefotaxime, y = -log10(P_Value__cefotaxime),
                                label = protein),
                  max.overlaps = Inf, size = 2, force = 4)

impipenem_plot <- impipenem_plot +
  geom_text_repel(data = top_impipenem,
                  mapping = aes(x = log_FC__impipenem, y = -log10(P_Value__impipenem),
                                label = protein),
                  max.overlaps = Inf, size = 2, force = 4)

ciprofloxacin_plot <- ciprofloxacin_plot +
  geom_text_repel(data = top_ciprofloxacin,
                  mapping = aes(x = log_FC__ciprofloxacin, y = -log10(P_Value__ciprofloxacin),
                                label = protein), 
                  max.overlaps = Inf, size = 2, force = 4)

# Exporting -------------------------------------------------------------------
export_path <- "./03__modeling/2022-10-14__multiple_comparisons/models/QFeatures_data"

# Exporting plots
save_plot(file.path(export_path, "results/Ampicillin_Control.png"),
          ampicillin_plot)
save_plot(file.path(export_path, "results/Cefotaxime_Control.png"),
          cefotaxime_plot)
save_plot(file.path(export_path, "results/Impipenem_Control.png"),
          impipenem_plot)
save_plot(file.path(export_path, "results/Ciprofloxacin_Control.png"),
          ciprofloxacin_plot)