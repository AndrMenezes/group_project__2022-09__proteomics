#' 2022-10-15 Read the evidence file and organize the data.
#' I'm using the QFeatures object to aggregate at peptides and protein levels
#' using the `colMedians()` function from {MatrixGenerics}
suppressMessages(library(QFeatures))
suppressMessages(library(HDF5Array))

path_data <- "./01__database/"

# Import the raw data (evidence) ------------------------------------------

data_raw <- read.delim(file.path(path_data, "raw_data/evidence.txt"))

cols <- list(
  charge = "Charge",
  sequence = 'Sequence',
  modified_sequence = 'Modified.sequence',
  modifications = 'Modifications',
  protein_group = 'Proteins',
  protein = 'Leading.Razor.Protein',
  experiment = 'Experiment',
  reverse = 'Reverse',
  contaminant = 'Contaminant',
  intensity = 'Intensity'
)
data_raw <- data_raw[, as.character(cols)]
colnames(data_raw) <- names(cols)
head(data_raw)


# Filters out contaminants and reverse sequences --------------------------
data_raw <- data_raw[
  which(data_raw$contaminant != '+' & data_raw$reverse != '+'),]

# Creating a unique id, since there are multiple peptides corresponding to the
# same sequence
chosen_cols <- c("modified_sequence", "protein", "experiment", "intensity")
data_raw <- data_raw[, chosen_cols] |> 
  dplyr::group_by(experiment, modified_sequence, protein) |> 
  dplyr::mutate(number_peptides = dplyr::n()) |> 
  dplyr::ungroup() |>
  dplyr::arrange(experiment, modified_sequence, protein)
data_raw$unique_id <- seq.int(1, nrow(data_raw))
head(data_raw)
dim(data_raw)


# Pivot to create the data matrix at PSMs level ---------------------------
pivotted_psms <- tidyr::pivot_wider(
  data = data_raw,
  id_cols = c(unique_id, number_peptides, modified_sequence, protein),
  names_from = experiment,
  values_from = intensity)

# Creating the QFeature object --------------------------------------------
se <- readQFeatures(table = pivotted_psms, ecol = 5:22, name = "psms")
assay(se[[1L]])
rowData(se[[1L]])

# Aggregate data at peptide level -----------------------------------------
# rows: peptides sequences and columns: samples
se <- aggregateFeatures(object = se, i = "psms", fcol = "modified_sequence",
                        name = "peptides", fun = colMedians, na.rm = TRUE)
se[["peptides"]]
head(assay(se[[2L]], 2))
rowData(se[[2L]])

# The same as:
# pivotted__peptides <- tidyr::pivot_wider(
#   data = data_raw,
#   id_cols = c(modified_sequence, protein),
#   names_from = experiment,
#   values_from = intensity,
#   values_fn = sum)

# Aggregate data at protein level -----------------------------------------
# rows: proteins and columns: samples
se <- aggregateFeatures(object = se, i = "peptides", fcol = "protein",
                        name = "proteins", fun = colMedians, na.rm = TRUE)
se
se[["proteins"]]
head(assay(se[[3L]], 1))
rowData(se[[3L]])


m <- assay(se[[3L]], 1)
head(m)
colMeans(is.na(m))
sum(rowMeans(is.na(m)) == 0) / nrow(m)


# Saving object as .h5 file -----------------------------------------------

saveHDF5SummarizedExperiment(se[[1L]], dir = file.path(
  path_data, "processed_data/psms_obj"))
saveHDF5SummarizedExperiment(se[[2L]], dir = file.path(
  path_data, "processed_data/peptides_obj"))
saveHDF5SummarizedExperiment(se[[3L]], dir = file.path(
  path_data, "processed_data/proteins_obj"))

