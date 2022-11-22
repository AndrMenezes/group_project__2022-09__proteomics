



most_de_limma <- data_de_limma |> 
  dplyr::filter(group == "Ampicillin", abs(log_fc) > 0.5) |> 
  dplyr::arrange(p_value) |> 
  head(3)


chosen_proteins <- most_de_limma$protein



pivotted_levels <- do.call(rbind, list_tbs) |> 
  dplyr::left_join(col_data, by = "colname") |> 
  dplyr::select(level, protein, colname, group, replicate, sample_names,
                value) |>
  dplyr::mutate(level = forcats::fct_relevel(factor(level), "psms", "peptides",
                                             "proteins"))
x11()
ggplot(pivotted_levels, aes(x = sample_names, y = value, colour = group)) +
  geom_point(size = 4) +
  facet_wrap(protein ~ level, scales = "free_y")

ggplot(pivotted_levels, aes(x = group, y = log2(value), colour = group)) +
  geom_point(size = 4) +
  facet_wrap(protein ~ level, scales = "free_y")


longFormat(fts[most_de_limma$protein[1:2], chosen_cols], ) |> 
  dplyr::as_tibble() |> 
  dplyr::left_join(col_data, by = "colname") |> 
  dplyr::select(assay, rowname, replicate, sample_names, group, value) |>
  dplyr::mutate(assay = forcats::fct_relevel(factor(assay), "psms", "peptides", "proteins")) |> 
  ggplot(aes(x = sample_names, y = value, colour = group)) +
  geom_point(size = 4) +
  facet_wrap(~ assay, scales = "free")

assay(se[most_de_limma$protein[1], chosen_cols], "intensity")
rowData(se[most_de_limma$protein[1], ])

aux <- longFormat(fts[most_de_limma$protein[1], chosen_cols], ) |> 
  dplyr::as_tibble() |> 
  dplyr::left_join(col_data, by = "colname") |> 
  dplyr::select(assay, rowname, replicate, sample_names, group, value) |> 
  dplyr::filter(!is.na(value))

fts
se


chosen_cols <- colData(fts)$group %in% c("Control", "Impipenem")

id <- which(rowData(fts[["peptides"]])$protein == "P0A9P6")
assay(fts[["peptides"]][id, ])

table(rownames(fts[["proteins"]]) == most_de$protein[1])

se_most_de <- fts[["proteins"]]["P0A9P6", chosen_cols]
rowData(se_most_de)
x0 <- assay(se_most_de, "intensity")
x1 <- assay(se_most_de, "log2_intensity")
x2 <- assay(se_most_de, "log2_normalized")
mean(x1[1:3] - x1[-c(1:3)])
mean(x2[1:3] - x2[-c(1:3)])
rowData(se_most_de)

rownames(se_fiona) <- rowData(se_fiona)$protein__id
chosen_cols <- colData(se_fiona)$group %in% c("Control", "Impipenem")
x3 <- assay(se_fiona["P0A9P6", chosen_cols])
mean(x3[c(1:3)] - x3[-c(1:3)])

x <- rexp(, rate = 0.6)
median(x)
mean(x)

data_de_limma |> 
  dplyr::filter(protein == "P0A9Q5")

data_de_fiona 
data_de_limma 

a <- data_de_fiona |> 
  dplyr::filter(fdr <= 0.05) |> 
  dplyr::pull(protein)
data_de_limma[data_de_limma$protein %in% a, ] |> 
  dplyr::filter(group == "Ciprofloxacin") |> 
  dplyr::arrange(protein)
data_de_fiona |> 
  dplyr::filter(fdr <= 0.05) |> 
  dplyr::arrange(protein)

data_de_limma

data_de_limma[data_de_limma$protein == "P62554", ]
data_de_fiona[data_de_fiona$protein == "ACQ41974.1;P62554", ]
chosen <- rowData(fts[["psms"]])$protein == "ACQ41977.1"
table(chosen)
fts[["psms"]][chosen, ]


se

data_curr <- data_results_fiona |> 
  dplyr::filter(Method == "Two Sample T Test")

y_proteins <- unique(data_curr$Protein)
x_proteins <- unique(data_limma$protein)

length(x_proteins)
length(y_proteins)
length(intersect(x_proteins, y_proteins))
df_1 <- data_curr |> 
  dplyr::filter(Treatment == "Ampicillin",
                P_value <= 0.005, abs(logFC) > 0.5) |> 
  dplyr::arrange(P_value) |> 
  head(2)

df_2 <- data_limma |> 
  dplyr::filter(group == "Ampicillin", p_value <= 0.005, abs(log_fc) > 0.5) |> 
  dplyr::arrange(p_value)

df <- data_curr |> 
  dplyr::filter(Treatment == "Ampicillin", Protein %in% c("P0A9Q5", "P24202", "P37005", "P0AB18")) |> 
  dplyr::rename(method = Method, protein = Protein, log_fc = logFC,
                p_value = P_value, group = Treatment) |> 
  dplyr::select(method, protein, group, log_fc, p_value) |> 
  dplyr::bind_rows(
    data_limma |> 
      dplyr::filter(group == "Ampicillin",
                    protein %in% c("P0A9Q5", "P24202", "P37005", "P0AB18")) |> 
      dplyr::select(method, protein, group, log_fc, p_value))


rowData(se[df_2$protein, ])

chosen  <- rowData(fts[["peptides"]])$protein == "P0A9J8"
fts[["peptides"]][chosen, ]
limma::plotDensities(assay(fts[["peptides"]][chosen, ]), group =colData(fts)$group)

