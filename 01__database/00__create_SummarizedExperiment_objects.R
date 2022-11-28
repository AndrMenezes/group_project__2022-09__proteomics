#' 2022-10-15: Organize the data into `SummarizedExperiment` objects.
#' To aggregate at peptides and protein levels we use the `colMedians()`
#' function from {MatrixGenerics}
suppressMessages(library(SummarizedExperiment))
suppressMessages(library(HDF5Array))

path_data <- "./01__database/"

# Import the raw data (evidence) ------------------------------------------

data_raw <- read.delim(file.path(path_data, "raw_data/evidence.txt"))

cols <- list(
  charge = "Charge",
  sequence = "Sequence",
  modified_sequence = "Modified.sequence",
  modifications = "Modifications",
  protein_group = "Proteins",
  protein = "Leading.Razor.Protein",
  experiment = "Experiment",
  reverse = "Reverse",
  contaminant = "Contaminant",
  intensity = "Intensity"
)
data_raw <- data_raw[, as.character(cols)]
colnames(data_raw) <- names(cols)
head(data_raw)


# Filters out contaminants and reverse sequences --------------------------
data_raw <- data_raw[
  which(data_raw$contaminant != "+" & data_raw$reverse != "+"), ]

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


# Columns data referring the design of experiment -------------------------
samples_names <- colnames(pivotted_psms)[-c(1:4)]
map_names <- list(
  "amp" = "Ampicillin",
  "cef" = "Cefotaxime",
  "imp" = "Imipenem",
  "cip" = "Ciprofloxacin",
  "cont" = "Control",
  "tet" = "tet"
)
group_names <- as.character(map_names[
  tolower(gsub("[0-9]+", "", samples_names))])
col_data <- DataFrame(group = group_names,
                      replicate = rep(1L:3L, length(map_names)),
                      row.names = paste0("Sample_",
                                         seq_len(length(group_names))))


# Creating a SummarizedExperiment for PSMs level --------------------------
row_data <- pivotted_psms[, c(1:4)]
m_psms <- as.matrix(pivotted_psms[, -c(1:4)])
colnames(m_psms) <- NULL
se_psms <- SummarizedExperiment(assays = list(intensity = m_psms),
                                rowData = row_data, colData = col_data)



# Aggregating at peptides level -------------------------------------------

m_peptide <- MsCoreUtils::aggregate_by_vector(
  x = assay(se_psms),
  INDEX = rowData(se_psms)$modified_sequence,
  FUN = colMedians, na.rm = TRUE)
head(m_peptide)
dim(m_peptide)

# Creating the reduced row DataFrame
row_data_peptide <- rowData(se_psms)[, -1]
row_data_peptide <- row_data_peptide[
  !duplicated(row_data_peptide$modified_sequence), ]
row_data_peptide <- row_data_peptide[
  match(rownames(m_peptide), row_data_peptide$modified_sequence), ]

# Sanity check of DataFrame order
stopifnot(all.equal(row_data_peptide$modified_sequence, rownames(m_peptide)))

# Creating the SE object
se_peptide <- SummarizedExperiment(
  assays = list(intensity = m_peptide),
  rowData = row_data_peptide, colData = col_data)

# Aggregating at proteins level -------------------------------------------

m_proteins <- MsCoreUtils::aggregate_by_vector(
  x = as.matrix(assay(se_peptide)),
  INDEX = rowData(se_peptide)$protein,
  FUN = colMedians, na.rm = TRUE)
head(m_proteins)
dim(m_proteins)

# Creating the reduced row DataFrame
row_data_proteins <- rowData(se_peptide) |>
  dplyr::as_tibble() |>
  dplyr::group_by(protein) |>
  dplyr::summarise(number_peptides = sum(number_peptides)) |>
  DataFrame()

# Sanity check of DataFrame order
stopifnot(all.equal(row_data_proteins$protein, rownames(m_proteins)))

# Creating the SE object
se_protein <- SummarizedExperiment(
  assays = list(intensity = DelayedArray(m_proteins)),
  rowData = row_data_proteins, colData = col_data)
se_protein


# Saving objects ----------------------------------------------------------

# This format provides on-disk representation of large data sets without the
# need to load them into memory.
saveHDF5SummarizedExperiment(se_psms, dir = file.path(
  path_data, "processed_data/psms"))
saveHDF5SummarizedExperiment(se_peptide, dir = file.path(
  path_data, "processed_data/peptides"))
saveHDF5SummarizedExperiment(se_protein, dir = file.path(
  path_data, "processed_data/proteins"))
