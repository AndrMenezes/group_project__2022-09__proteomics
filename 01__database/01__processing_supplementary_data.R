#' 2022-09-21: Creating data set in both long and filtered format and a dictionary.
library(dplyr)
library(readxl)

path_data <- "./01__database/"

data_raw <- read_xlsx(
  file.path(path_data, "raw_data/PRE_POSTimputated_pEK499.xlsx"),
  sheet = "Suppl_2_Post-imputated dataset",
  range = "A2:Y947")

data_raw <- data_raw |> 
  rename(
    protein__name = `Protein names`,
    protein__id = `Majority protein IDs`,
    gene__name = `Gene names`,
    Ampicillin__1 = `LFQ intensity Ampicillin Rep 1`,
    Ampicillin__2 = `LFQ intensity Ampicillin Rep 2`,
    Ampicillin__3 = `LFQ intensity Ampicillin Rep 3`,
    Cefotaxime__1 = `LFQ intensity Cefotaxime Rep 1`,
    Cefotaxime__2 = `LFQ intensity Cefotaxime Rep 2`,
    Cefotaxime__3 = `LFQ intensity Cefotaxime Rep 3`,
    Impipenem__1 = `LFQ intensity Impipenem Rep 1`,
    Impipenem__2 = `LFQ intensity Impipenem Rep 2`,
    Impipenem__3 = `LFQ intensity Impipenem Rep 3`,
    Ciprofloxacin__1 = `LFQ intensity Ciprofloxacin Rep 1`,
    Ciprofloxacin__2 = `LFQ intensity Ciprofloxacin Rep 2`,
    Ciprofloxacin__3 = `LFQ intensity Ciprofloxacin Rep 3`,
    Control__1 = `LFQ intensity Control Rep 1`,
    Control__2 = `LFQ intensity Control Rep 2`,
    Control__3 = `LFQ intensity Control Rep 3`,
  ) |> 
  select(protein__name, protein__id, gene__name,
         Ampicillin__1:Control__3)

# Importing the not imputed data set
data_raw_not_imputed <- read_xlsx(
  file.path(path_data, "raw_data/PRE_POSTimputated_pEK499.xlsx"),
  sheet = "Suppl_1_Pre-imputated",
  range = "A2:Z947", na = "NaN")

data_raw_not_imputed <- data_raw_not_imputed |> 
  rename(
    protein__name = `Protein names`,
    protein__id = `Protein IDs`,
    gene__name = `Gene names`,
    Ampicillin__1 = `LFQ intensity Amp_Rep1`,
    Ampicillin__2 = `LFQ intensity Amp_Rep2`,
    Ampicillin__3 = `LFQ intensity Amp_Rep3`,
    Cefotaxime__1 = `LFQ intensity Cef_Rep1`,
    Cefotaxime__2 = `LFQ intensity Cef_Rep2`,
    Cefotaxime__3 = `LFQ intensity Cef_Rep3`,
    Impipenem__1 = `LFQ intensity Imp_Rep2...24`,
    Impipenem__2 = `LFQ intensity Imp_Rep2...25`,
    Impipenem__3 = `LFQ intensity Imp_Rep3`,
    Ciprofloxacin__1 = `LFQ intensity cip_Rep1`,
    Ciprofloxacin__2 = `LFQ intensity Cip_Rep2`,
    Ciprofloxacin__3 = `LFQ intensity Cip_Rep3`,
    Control__1 = `LFQ intensity Cont_Rep1`,
    Control__2 = `LFQ intensity Cont_Rep2`,
    Control__3 = `LFQ intensity Cont_Rep3`,
  ) |> 
  select(protein__name, protein__id, gene__name,
         Ampicillin__1:Impipenem__3) |> 
  mutate(protein__id = gsub(";.*", "", protein__id))

# Ordering the matrices ---------------------------------------------------
data_raw_not_imputed <- data_raw_not_imputed[
  match(data_raw$protein__id, data_raw_not_imputed$protein__id),
  colnames(data_raw)]

# Columns data referring the design of experiment -------------------------
samples_names <- colnames(data_raw)[-c(1:3)]
group_names <- gsub("__", "", gsub("[0-9]+", "", samples_names))
col_data <- S4Vectors::DataFrame(
  group = group_names,
  replicate = rep(1L:3L, length(unique(group_names))),
  row.names = paste0("Sample_", seq_len(length(group_names))))

# Creating a SummarizedExperiment for PSMs level --------------------------
row_data <- data_raw[, c(1:3)]
mat_imputed <- as.matrix(data_raw[, -c(1:3)])
colnames(mat_imputed) <- NULL
rownames(mat_imputed) <- data_raw$protein__id
mat_not_imputed <- as.matrix(data_raw_not_imputed[, -c(1:3)])
colnames(mat_not_imputed) <- NULL
rownames(mat_not_imputed) <- data_raw$protein__id

se <- SummarizedExperiment::SummarizedExperiment(
  assays = list(
    log2_intensity = mat_not_imputed,
    log2_imputed = mat_imputed),
  rowData = row_data, colData = col_data)

saveRDS(object = se,
        file = file.path(path_data, "processed_data", "se_margalit.rds"))
