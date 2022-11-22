#' 2022-11-04: Perform DE analysis using the processed data and limma.
library(ggplot2)
library(cowplot)
theme_set(
  theme_cowplot() +
    background_grid() +
    theme(legend.position = "top")
)

# Path -------------------------------------------------------------------

local_path <- "./03__modeling/2022-10-14__workflow"
path_data <- file.path(local_path, "processing", "processed_data")
path_res <- file.path(local_path, "models", "results")

# Importing data --------------------------------------------------------------
fts <- readRDS(file.path(path_data, "fts_processsed.rds"))
se_fiona <- readRDS("./01__database/processed_data/se_processed.rds")


# Limma analysis on processed data ----------------------------------------

# Creating the design matrix
g <- factor(colData(fts)$group)
g <- relevel(g, ref = "Control")
design_matrix <- model.matrix(~ g)
colnames(design_matrix) <- gsub("g", "", colnames(design_matrix))

fit <- limma::lmFit(
  object = assay(fts[["proteins_median"]], "log2_normalized"),
  design = design_matrix)
fit <- limma::eBayes(fit)
ampicillin <- limma::topTable(fit, coef = "Ampicillin", number = Inf,
                              sort.by = "none", confint = TRUE)
cefotaxime <- limma::topTable(fit, coef = "Cefotaxime", number = Inf,
                              sort.by = "none", confint = TRUE)
impipenem <- limma::topTable(fit, coef = "Impipenem", number = Inf,
                             sort.by = "none", confint = TRUE)
ciprofloxacin <- limma::topTable(fit, coef = "Ciprofloxacin", number = Inf,
                                 sort.by = "none", confint = TRUE)

# Appending
data_de_limma <- dplyr::bind_rows(
    dplyr::mutate(ampicillin, group = "Ampicillin",
                  protein = rownames(ampicillin)),
    dplyr::mutate(cefotaxime, group = "Cefotaxime",
                  protein = rownames(cefotaxime)),
    dplyr::mutate(impipenem, group = "Impipenem",
                  protein = rownames(impipenem)),
    dplyr::mutate(ciprofloxacin, group = "Ciprofloxacin",
                  protein = rownames(ciprofloxacin))) |> 
  dplyr::as_tibble() |> 
  dplyr::mutate(method = "limma") |> 
  dplyr::select(method, protein, group, logFC, P.Value, adj.P.Val) |> 
  dplyr::rename(log_fc = logFC, p_value = P.Value, fdr = adj.P.Val)

data_de_limma |> 
  dplyr::filter(abs(log_fc) > 0.5, p_value < 0.005)

# Fiona's Student t-test --------------------------------------------------
expr_mat <- assay(se_fiona)
list_student_t <- lapply(seq_len(nrow(expr_mat)), function(i) {
  x_contrl <- expr_mat[i, colData(se_fiona)$group == "Control"]
  x_ampici <- expr_mat[i, colData(se_fiona)$group == "Ampicillin"]
  x_cefota <- expr_mat[i, colData(se_fiona)$group == "Cefotaxime"]
  x_impipe <- expr_mat[i, colData(se_fiona)$group == "Impipenem"]
  x_ciprof <- expr_mat[i, colData(se_fiona)$group == "Ciprofloxacin"]
  tibble::tibble(
    protein = rowData(se_fiona)$protein__id[i],
    group = c("Ampicillin", "Cefotaxime", "Impipenem", "Ciprofloxacin"),
    log_fc = c(mean(x_ampici - x_contrl), mean(x_cefota - x_contrl),
               mean(x_impipe - x_contrl), mean(x_ciprof - x_contrl)),
    p_value = c(t.test(x_contrl, x_ampici, var.equal = TRUE)$p.value,
                t.test(x_contrl, x_cefota, var.equal = TRUE)$p.value,
                t.test(x_contrl, x_impipe, var.equal = TRUE)$p.value,
                t.test(x_contrl, x_ciprof, var.equal = TRUE)$p.value)
  )
})

data_de_fiona <- do.call(rbind, list_student_t) |> 
  dplyr::mutate(method = "student_t") |> 
  dplyr::group_by(group) |> 
  dplyr::mutate(fdr = p.adjust(p_value, method = "fdr")) |> 
  dplyr::ungroup() |> 
  dplyr::select(method, protein, group, log_fc, p_value, fdr)


# Appending results -------------------------------------------------------
data_de <- dplyr::bind_rows(data_de_limma, data_de_fiona) |> 
  dplyr::mutate(method_label = ifelse(method == "student_t",
                                      "Margalit et al. (2022)",
                                      "Processed data + Limma"))


# Visualizing the difference ----------------------------------------------
data_de_pivotted <- data_de |> 
  dplyr::group_by(group) |> 
  dplyr::mutate(rank_pvalue = rank(p_value)) |> 
  dplyr::ungroup() |> 
  dplyr::select(-c(method_label, fdr)) |> 
  tidyr::pivot_wider(names_from = method,
                     values_from = c(log_fc, p_value, rank_pvalue)) |> 
  dplyr::filter(!(is.na(log_fc_limma) | is.na(log_fc_student_t))) |> 
  dplyr::mutate(diff_log_fc = abs(log_fc_student_t) - abs(log_fc_limma),
                diff_p_value = p_value_student_t - p_value_limma,
                diff_rank = rank_pvalue_student_t - rank_pvalue_limma)

data_de_pivotted |> 
  dplyr::group_by(group) |> 
  dplyr::summarise(diff_log_fc_lt_0 = mean(diff_log_fc < 0),
                   diff_p_value_lt_0 = mean(diff_p_value < 0))


p_diff_log_fc <- ggplot(data_de_pivotted, aes(x = diff_log_fc)) +
  facet_wrap(~group) +
  geom_density() +
  geom_rug() +
  geom_vline(xintercept = 0, col = "red") +
  labs(x = "Difference in the estimated log fold change", y = "Density") +
  scale_x_continuous(breaks = scales::pretty_breaks(6)) +
  scale_y_continuous(breaks = scales::pretty_breaks(4))
save_plot(filename = file.path(path_res, "diff_log_fc.png"),
          plot = p_diff_log_fc, base_width = 8, bg = "white")

p_diff_p_value <- ggplot(data_de_pivotted, aes(x = diff_p_value)) +
  facet_wrap(~group) +
  geom_density() +
  geom_rug() +
  geom_vline(xintercept = 0, col = "red") +
  labs(x = "Difference in the p-value", y = "Density") +
  scale_x_continuous(breaks = scales::pretty_breaks(6)) +
  scale_y_continuous(breaks = scales::pretty_breaks(4))
save_plot(filename = file.path(path_res, "diff_p_value.png"),
          plot = p_diff_p_value, base_width = 8, bg = "white")


# Investigating the rank
ggplot(data_de_pivotted, aes(x = diff_rank)) +
  facet_wrap(~group) +
  geom_density() +
  geom_rug() +
  geom_vline(xintercept = 0, col = "red") +
  labs(x = "Difference in the rank of the p-value", y = "Density")


# Volcano plots -----------------------------------------------------------

volcano_plot <- function(data, cutoff_logfc = 0.5,
                         cutoff_pvalue = 0.005, top_de = 10) {
  data_labelled <- data |> 
    dplyr::filter(abs(log_fc) >= cutoff_logfc, p_value <= cutoff_pvalue) |> 
    dplyr::arrange(-abs(log_fc)) |> 
    head(top_de)
  p <- ggplot(data, aes(x = log_fc, y = -log10(p_value))) +
    geom_point(size = 2, shape = 21, col = "black", fill = "grey69") +
    geom_point(data = data_labelled, shape = 21, size = 2, fill = "red",
               colour = "black") +
    ggrepel::geom_text_repel(
      data = data_labelled, aes(label = protein), size = 4, color = "black") +
    geom_vline(xintercept = c(-cutoff_logfc, cutoff_logfc), linetype = "dashed",
               col = "black") +
    geom_hline(yintercept = -log10(cutoff_pvalue), linetype = "dashed",
               col = "black") +
    labs(x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) +
    scale_x_continuous(breaks = scales::pretty_breaks(8)) +
    scale_y_continuous(breaks = scales::pretty_breaks(6))
  
  p
}
p_a <- volcano_plot(data = dplyr::filter(data_de_limma, group == "Ampicillin"))
save_plot(filename = file.path(path_res, "volcano_plot__ampicillin.png"),
          plot = p_a, base_width = 6, bg = "white")

volcano_plot(data = dplyr::filter(data_de_limma, group == "Cefotaxime"))
volcano_plot(data = dplyr::filter(data_de_limma, group == "Impipenem"))
volcano_plot(data = dplyr::filter(data_de_limma, group == "Ciprofloxacin"))

# Comparing the total of DE proteins according to method ------------------
p_total_de_proteins <- data_de |> 
  dplyr::filter(p_value <= 0.005) |> 
  dplyr::group_by(method_label, group) |> 
  dplyr::count() |> 
  ggplot(aes(x = group, y = n, fill = method_label)) +
  geom_col(position = position_dodge(width = 0.9), alpha = 0.6,
           col = "black") +
  geom_text(aes(label = n), vjust = -0.05,
            position = position_dodge(width = 0.9), size = 6) +
  labs(x = "", y = "", fill = "") +
  scale_y_continuous(limits = c(0, 55)) +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        axis.text.y = element_blank())
save_plot(filename = file.path(path_res, "total_de_proteins.png"),
          plot = p_total_de_proteins, base_width = 8, bg = "white")


# Comparing the 3 most expressed proteins ---------------------------------

groups <- unique(colData(fts)$group)
groups <- groups[which(groups != "Control")]
list_de <- list()
j <- 1L
for (g in groups) {
  cat(g, "\n")

  most_de_fiona <- data_de_fiona |> 
    dplyr::filter(group == g, abs(log_fc) > 0.5) |> 
    dplyr::arrange(p_value) |> 
    head(3)
  
  most_de_limma <- data_de_limma |> 
    dplyr::filter(group == g, abs(log_fc) > 0.5) |> 
    dplyr::arrange(p_value) |> 
    head(3)

  data_curr <- data_de |> 
    dplyr::filter(group == g, protein %in% most_de_fiona$protein) |> 
    dplyr::mutate(most_de = "fiona") |> 
    dplyr::bind_rows(
      data_de |> 
        dplyr::filter(group == g, protein %in% most_de_limma$protein) |> 
        dplyr::mutate(most_de = "ours"))

  data_curr_f <- data_curr |> 
    dplyr::filter(most_de == "fiona")
  data_curr_o <- data_curr |> 
    dplyr::filter(most_de == "ours")
  
  # Limits
  y_lim_f <- c(min(data_curr_f$log_fc) - 0.25, max(data_curr_f$log_fc) + 0.25)
  y_lim_o <- c(min(data_curr_o$log_fc) - 0.25, max(data_curr_o$log_fc) + 0.25)
  
  # Plotting
  p_most_de_fiona <- ggplot(data_curr_f, aes(x = protein, y = log_fc,
                                             fill = method_label)) +
    geom_col(position = position_dodge(width = 0.9),
             alpha = 0.6, col = "black") +
    geom_hline(yintercept = 0, col = "black", size = 1.5) +
    geom_text(aes(label = paste0(round(log_fc, 2), " (",
                                 formatC(p_value, format = "e", digits = 1),
                                 ")")),
              position = position_dodge(width = 0.9), vjust = -0.2,
              size = 4) +
    labs(x = "Protein", y = "Log Fold Change", fill = "") +
    scale_y_continuous(breaks = scales::pretty_breaks(6), limits = y_lim_f) +
    ggtitle("3 most DE proteins by Margalit et al. (2022)")
  
  p_most_de_ours <- ggplot(data_curr_o, aes(x = protein, y = log_fc,
                                            fill = method_label)) +
    geom_col(position = position_dodge(width = 0.9),
             alpha = 0.6, col = "black") +
    geom_hline(yintercept = 0, col = "black", size = 1.5) +
    geom_text(aes(label = paste0(round(log_fc, 2), " (",
                                 formatC(p_value, format = "e", digits = 1),
                                 ")")),
              position = position_dodge(width = 0.9), vjust = -0.2,
              size = 4) +
    labs(x = "Protein", y = "Log Fold Change", fill = "") +
    scale_y_continuous(breaks = scales::pretty_breaks(6), limits = y_lim_o) +
    ggtitle("3 most DE proteins by Processed data + Limma")
  p_grided <- plot_grid(p_most_de_fiona, p_most_de_ours)
  fname <- paste0("most_de__", tolower(g), ".png")
  save_plot(filename = file.path(path_res, fname),
            plot = p_grided, base_width = 17, bg = "white")
  
  # Appending most DE proteins
  list_de[[j]] <- data.frame(
    group = g, fiona = most_de_fiona$protein, ours = most_de_limma$protein)
  j <- j + 1L

}
de_proteins <- do.call(rbind, list_de)



# Inspecting the raw data for some DE proteins ----------------------------
col_data <- dplyr::as_tibble(colData(fts)) |> 
  dplyr::mutate(colname = rownames(colData(fts)))

# Function to pivot the raw intensity matrix into a long-format tibble
pivot_qfeatures <- function(fts, chosen_proteins) {
  list_tbs <- lapply(seq_len(length(fts)), function(j) {
    level <- names(fts)[j] 
    se_curr <- fts[[j]]
    chosen_rows <- rowData(se_curr)$protein %in% chosen_proteins
    
    i_assay <- if (level == "proteins") "intensity" else 1L
    row_data_curr <- dplyr::as_tibble(
      rowData(se_curr)[chosen_rows, "protein", drop = FALSE])
    
    pivotted <- assay(se_curr[chosen_rows, ], i_assay) |> 
      dplyr::as_tibble() |> 
      dplyr::bind_cols(row_data_curr) |> 
      tidyr::pivot_longer(cols = -c("protein"), names_to = "colname") |> 
      dplyr::mutate(level = level)
    pivotted
  })
  do.call(rbind, list_tbs) 
}

for (g in groups) {
  
  cat(g, "\n")
  # Selecting the corresponding columns
  chosen_cols <- colData(fts)$group %in% c("Control", g)
  fts_curr <- fts[, chosen_cols]

  # Selecting the proteins
  most_de_fiona <- de_proteins[de_proteins$group == g, "fiona"]
  most_de_limma <- de_proteins[de_proteins$group == g, "ours"]
  
  # Get the intensity values pivoted in a long format for all levels
  pivotted_fiona <- pivot_qfeatures(fts = fts_curr,
                                    chosen_proteins = most_de_fiona)
  pivotted_limma <- pivot_qfeatures(fts = fts_curr,
                                    chosen_proteins = most_de_limma)
  
  # Join to get the columns data information
  pivotted_fiona <- pivotted_fiona |>
    dplyr::left_join(col_data, by = "colname") |> 
    dplyr::select(level, protein, colname, group, replicate, sample_names,
                  value) |>
    dplyr::mutate(level = forcats::fct_relevel(
      factor(level), "psms", "peptides", "proteins"),
      group = forcats::fct_relevel(group, g)) |>
    dplyr::filter(level != "psms")
  pivotted_limma <- pivotted_limma |> 
    dplyr::left_join(col_data, by = "colname") |> 
    dplyr::select(level, protein, colname, group, replicate, sample_names,
                  value) |>
    dplyr::mutate(level = forcats::fct_relevel(
      factor(level), "psms", "peptides", "proteins"),
      group = forcats::fct_relevel(group, g)) |>
    dplyr::filter(level != "psms")
  
  # Plotting
  p_fiona <- ggplot(pivotted_fiona, aes(x = group, y = log2(value),
                                        colour = group)) +
    geom_point(size = 4) +
    facet_grid(protein ~ level, scales = "free") +
    labs(x = "Group", y = "Log2 intensity", col = "") +
    ggtitle("Distribution of the 3 most DE by Margalit et al. (2022)")
  
  p_limma <- ggplot(pivotted_limma, aes(x = group, y = log2(value),
                                        colour = group)) +
    geom_point(size = 4) +
    facet_wrap(protein ~ level, scales = "free", ncol = 2) +
    labs(x = "Group", y = "Log2 intensity", col = "") +
    ggtitle("Distribution of the 3 most DE by  Processed data + Limma")
  
  # Saving
  fname <- paste0("intensity_distr_fiona__", g, ".png")
  save_plot(filename = file.path(path_res, fname), plot = p_fiona,
            base_height = 8, bg = "white")
  fname <- paste0("intensity_distr_limma__", g, ".png")
  save_plot(filename = file.path(path_res, fname), plot = p_limma,
            base_height = 8, bg = "white")
}
