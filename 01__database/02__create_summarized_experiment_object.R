#' 2022-10-14: Creating a SummarizedExperiment object.
#' This helps in the reproducibility as well as to keep a pattern.
suppressMessages(library(SummarizedExperiment))
suppressMessages(library(HDF5Array))
library(readxl)

path_data <- "./01__database/"

data_raw <- read_xlsx(
  file.path(path_data, "raw_data/Supplementary Data Stat Sigs.xlsx"),
  sheet = "All protein IDs",
  range = "A2:Y947")
head(data_raw)
dplyr::glimpse(data_raw)

log_intensity <- as.matrix(data_raw[, 11:25])
colnames(log_intensity) <- NULL
row_data <- DataFrame(protein_name = data_raw$`Protein names`,
                      gene_name = data_raw$`Gene names`,
                      sequence_length = data_raw$`Sequence length`,
                      row.names =  data_raw$`Majority protein IDs`)
col_data <- DataFrame(replicate = rep(1L:3L, times = 5L),
                      group = rep(c("Ampicillin", "Cefotaxime", "Impipenem",
                                    "Ciprofloxacin", "Control"), each = 3),
                      row.names = paste0("Sample_", seq_len(ncol(log_intensity))))
se <- SummarizedExperiment(
  assays = list(log_intensity = log_intensity),
  rowData = row_data, colData = col_data)


# Save data using .h5 format ----------------------------------------------
saveHDF5SummarizedExperiment(se, dir = file.path(path_data, "processed_data/se_obj"))

# saveRDS(object = se,
#         file = file.path(path_data, "processed_data/se_object.rds"))
