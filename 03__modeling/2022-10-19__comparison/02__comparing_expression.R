library(dplyr)
library(readxl)
library(ggplot2)
library(ggrepel)
library(cowplot)

path_data <- "./03__modeling/2022-10-19__comparison"

combined_data <- read_xlsx(file.path(path_data, "results/all_results.xlsx"))

# Separating based on treatment and noting up/down regulation -----------------
ampicillin_data <- combined_data |>
  filter(Treatment == "Ampicillin") |>
  mutate(Expression = case_when(logFC >= 0.5 & -log10(P_value) >= -log10(.05) ~ "Up Regulated",
                                logFC <= -0.5 & -log10(P_value) >= -log10(.05) ~ "Down Regulated",
                                TRUE ~ "Unchanged"))
cefotaxime_data <- combined_data |>
  filter(Treatment == "Cefotaxime") |>
  mutate(Expression = case_when(logFC >= 0.5 & -log10(P_value) >= -log10(.05) ~ "Up Regulated",
                                logFC <= -0.5 & -log10(P_value) >= -log10(.05) ~ "Down Regulated",
                                TRUE ~ "Unchanged"))
impipenem_data <- combined_data |>
  filter(Treatment == "Impipenem") |>
  mutate(Expression = case_when(logFC >= 0.5 & -log10(P_value) >= -log10(.05) ~ "Up Regulated",
                                logFC <= -0.5 & -log10(P_value) >= -log10(.05) ~ "Down Regulated",
                                TRUE ~ "Unchanged"))
ciprofloxacin_data <- combined_data |>
  filter(Treatment == "Ciprofloxacin") |>
  mutate(Expression = case_when(logFC >= 0.5 & -log10(P_value) >= -log10(.05) ~ "Up Regulated",
                                logFC <= -0.5 & -log10(P_value) >= -log10(.05) ~ "Down Regulated",
                                TRUE ~ "Unchanged"))

# Labelling most abundant proteins for each treatment -------------------------
ampicillin_labels <- bind_rows(
  ampicillin_data %>%
    filter(Method == "Two Sample T Test" & Expression == "Up Regulated") %>%
    arrange(desc(logFC)) %>%
    head(10),
  ampicillin_data %>%
    filter(Method == "Two Sample T Test", Expression == "Down Regulated") %>%
    arrange(logFC) %>%
    head(10),
  ampicillin_data %>%
    filter(Method == "Bayesian T Test" & Expression == "Up Regulated") %>%
    arrange(desc(logFC)) %>%
    head(10),
  ampicillin_data %>%
    filter(Method == "Bayesian T Test", Expression == "Down Regulated") %>%
    arrange(logFC) %>%
    head(10)
)

cefotaxime_labels <- bind_rows(
  cefotaxime_data %>%
    filter(Method == "Two Sample T Test" & Expression == "Up Regulated") %>%
    arrange(desc(logFC)) %>%
    head(10),
  cefotaxime_data %>%
    filter(Method == "Two Sample T Test", Expression == "Down Regulated") %>%
    arrange(logFC) %>%
    head(10),
  cefotaxime_data %>%
    filter(Method == "Bayesian T Test" & Expression == "Up Regulated") %>%
    arrange(desc(logFC)) %>%
    head(10),
  cefotaxime_data %>%
    filter(Method == "Bayesian T Test", Expression == "Down Regulated") %>%
    arrange(logFC) %>%
    head(10)
)

impipenem_labels <- bind_rows(
  impipenem_data %>%
    filter(Method == "Two Sample T Test" & Expression == "Up Regulated") %>%
    arrange(desc(logFC)) %>%
    head(10),
  impipenem_data %>%
    filter(Method == "Two Sample T Test", Expression == "Down Regulated") %>%
    arrange(logFC) %>%
    head(10),
  impipenem_data %>%
    filter(Method == "Bayesian T Test" & Expression == "Up Regulated") %>%
    arrange(desc(logFC)) %>%
    head(10),
  impipenem_data %>%
    filter(Method == "Bayesian T Test", Expression == "Down Regulated") %>%
    arrange(logFC) %>%
    head(10)
)

ciprofloxacin_labels <- bind_rows(
  ciprofloxacin_data %>%
    filter(Method == "Two Sample T Test" & Expression == "Up Regulated") %>%
    arrange(desc(logFC)) %>%
    head(10),
  ciprofloxacin_data %>%
    filter(Method == "Two Sample T Test", Expression == "Down Regulated") %>%
    arrange(logFC) %>%
    head(10),
  ciprofloxacin_data %>%
    filter(Method == "Bayesian T Test" & Expression == "Up Regulated") %>%
    arrange(desc(logFC)) %>%
    head(10),
  ciprofloxacin_data %>%
    filter(Method == "Bayesian T Test", Expression == "Down Regulated") %>%
    arrange(logFC) %>%
    head(10)
)

# Plotting --------------------------------------------------------------------
ampicillin_plot <- ggplot(data = ampicillin_data,
                          aes(x = logFC, y = -log10(P_value), colour = Expression)) +
  facet_wrap(~ Method, ncol = 1) +
  ggtitle("Control Vs. Ampicillin") +
  geom_point(size = .0005) +
  xlab(expression("log"[2]*" Fold Change")) + 
  ylab(expression("-log"[10]* "p")) +
  scale_colour_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", col = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", col = "black") +
  theme_bw()

cefotaxime_plot <- ggplot(data = cefotaxime_data,
                          aes(x = logFC, y = -log10(P_value), colour = Expression)) +
  facet_wrap(~ Method, ncol = 1) +
  ggtitle("Control Vs. Cefotaxime") +
  geom_point(size = .0005) +
  xlab(expression("log"[2]*" Fold Change")) + 
  ylab(expression("-log"[10]* "p")) +
  scale_colour_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", col = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", col = "black") +
  theme_bw()

impipenem_plot <- ggplot(data = impipenem_data,
                         aes(x = logFC, y = -log10(P_value), colour = Expression)) +
  facet_wrap(~ Method, ncol = 1) +
  ggtitle("Control Vs. Impipenem") +
  geom_point(size = .0005) +
  xlab(expression("log"[2]*" Fold Change")) + 
  ylab(expression("-log"[10]* "p")) +
  scale_colour_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", col = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", col = "black") +
  theme_bw()

ciprofloxacin_plot <- ggplot(data = ciprofloxacin_data,
                             aes(x = logFC, y = -log10(P_value), colour = Expression)) +
  facet_wrap(~ Method, ncol = 1) +
  ggtitle("Control Vs. Ciprofloxacin") +
  geom_point(size = .0005) +
  xlab(expression("log"[2]*" Fold Change")) + 
  ylab(expression("-log"[10]* "p")) +
  scale_colour_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", col = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", col = "black") +
  theme_bw()

# Labelling
ampicillin_plot <- ampicillin_plot +
  geom_text_repel(data = ampicillin_labels,
                  mapping = aes(x = logFC, y = -log10(P_value), label = Gene),
                  max.overlaps = Inf, size = 2)

cefotaxime_plot <- cefotaxime_plot +
  geom_text_repel(data = cefotaxime_labels,
                  mapping = aes(x = logFC, y = -log10(P_value), label = Gene),
                  max.overlaps = Inf, size = 2)

impipenem_plot <- impipenem_plot +
  geom_text_repel(data = impipenem_labels,
                  mapping = aes(x = logFC, y = -log10(P_value), label = Gene),
                  max.overlaps = Inf, size = 2)

ciprofloxacin_plot <- ciprofloxacin_plot +
  geom_text_repel(data = ciprofloxacin_labels,
                  mapping = aes(x = logFC, y = -log10(P_value), label = Gene),
                  max.overlaps = Inf, size = 2)

# Exporting -------------------------------------------------------------------
export_path <- "./03__modeling/2022-10-19__comparison"

save_plot(file.path(export_path, "results/Ampicillin.png"), ampicillin_plot,
          base_height = 5)
save_plot(file.path(export_path, "results/Cefotaxime.png"), cefotaxime_plot,
          base_height = 5)
save_plot(file.path(export_path, "results/Impipenem.png"), impipenem_plot,
          base_height = 5)
save_plot(file.path(export_path, "results/Ciprofloxacin.png"), ciprofloxacin_plot,
          base_height = 5)

