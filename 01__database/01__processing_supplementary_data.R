#' 2022-09-21: Creating data set in both long and filtered format and a dictionary.
library(dplyr)
library(readxl)

path_data <- "./01__database/"

data_raw <- read_xlsx(
  file.path(path_data, "raw_data/raw_data.xlsx"),
  sheet = "All protein IDs",
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


# Columns data referring the design of experiment -------------------------
samples_names <- colnames(data_raw)[-c(1:3)]
group_names <- gsub("__", "", gsub("[0-9]+", "", samples_names))
col_data <- DataFrame(group = group_names,
                      replicate = rep(1L:3L, length(unique(group_names))),
                      row.names = paste0("Sample_",
                                         seq_len(length(group_names))))


# Creating a SummarizedExperiment for PSMs level --------------------------
row_data <- data_raw[, c(1:3)]
mat <- as.matrix(data_raw[, -c(1:3)])
colnames(mat) <- NULL
rownames(mat) <- data_raw$protein__id
se <- SummarizedExperiment(assays = list(log_intensity = mat),
                           rowData = row_data, colData = col_data)

saveRDS(object = se,
        file = file.path(path_data, "processed_data", "se_margalit.rds"))
