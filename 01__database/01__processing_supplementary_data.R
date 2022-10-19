#' 2022-09-21: Creating data set in both long and filtered format and a dictionary.
library(dplyr)
library(readxl)

path_data <- "./01__database/"

data_raw <- read_xlsx(
  file.path(path_data, "raw_data/raw_data.xlsx"),
  sheet = "All protein IDs",
  range = "A2:Y947")

# head(data_raw)
# glimpse(data_raw)


# Renaming --------------------------------------------------------------------
data_raw <- data_raw |> 
  rename(
    protein__name = `Protein names`,
    protein__id = `Majority protein IDs`,
    gene__name = `Gene names`,
    peptides = Peptides,
    sequence_coverage = `Sequence coverage [%]`,
    molecular_weight = `Mol. weight [kDa]`,
    score = Score,
    intensity = Intensity,
    count = `MS/MS count`,
    sequence_length = `Sequence length`,
    lqf__Ampicillin__1 = `LFQ intensity Ampicillin Rep 1`,
    lqf__Ampicillin__2 = `LFQ intensity Ampicillin Rep 2`,
    lqf__Ampicillin__3 = `LFQ intensity Ampicillin Rep 3`,
    lqf__Cefotaxime__1 = `LFQ intensity Cefotaxime Rep 1`,
    lqf__Cefotaxime__2 = `LFQ intensity Cefotaxime Rep 2`,
    lqf__Cefotaxime__3 = `LFQ intensity Cefotaxime Rep 3`,
    lqf__Impipenem__1 = `LFQ intensity Impipenem Rep 1`,
    lqf__Impipenem__2 = `LFQ intensity Impipenem Rep 2`,
    lqf__Impipenem__3 = `LFQ intensity Impipenem Rep 3`,
    lqf__Ciprofloxacin__1 = `LFQ intensity Ciprofloxacin Rep 1`,
    lqf__Ciprofloxacin__2 = `LFQ intensity Ciprofloxacin Rep 2`,
    lqf__Ciprofloxacin__3 = `LFQ intensity Ciprofloxacin Rep 3`,
    lqf__Control__1 = `LFQ intensity Control Rep 1`,
    lqf__Control__2 = `LFQ intensity Control Rep 2`,
    lqf__Control__3 = `LFQ intensity Control Rep 3`,
  )


# Pivoting --------------------------------------------------------------------
data_pivotted <- data_raw |> 
  select(-c(peptides, protein__name, sequence_coverage,
            molecular_weight, score, intensity, count, sequence_length)) |> 
  tidyr::pivot_longer(-c(protein__id, gene__name)) |> 
  tidyr::separate(col = name, into = c("variable", "group", "replicate"),
                  sep = "__") |> 
  select(protein__id, gene__name, group, replicate, variable, value)

# Adding values to new columns
data_pivotted$group <- rep(c("Ampicillin", "Cefotaxime", "Impipenem", 
                          "Ciprofloxacin", "Control"), each = 3, times = 945)
data_pivotted$replicate <- rep(c(1, 2, 3), times = 4725)
data_pivotted$variable <- rep("lfq", times = 14175)

# Filtering -------------------------------------------------------------------
data_filter <- data_raw |>
  select(
    protein__id,
    gene__name,
    lqf__Ampicillin__1,
    lqf__Ampicillin__2,
    lqf__Ampicillin__3,
    lqf__Cefotaxime__1,
    lqf__Cefotaxime__2,
    lqf__Cefotaxime__3,
    lqf__Impipenem__1,
    lqf__Impipenem__2,
    lqf__Impipenem__3,
    lqf__Ciprofloxacin__1,
    lqf__Ciprofloxacin__2,
    lqf__Ciprofloxacin__3,
    lqf__Control__1,
    lqf__Control__2,
    lqf__Control__3
  )

# Dictionary ------------------------------------------------------------------
data_dict <- data_raw |> 
  select(protein__id, protein__name, gene__name,
         peptides, sequence_coverage, molecular_weight, score, intensity,
         count, sequence_length)


# Exporting -------------------------------------------------------------------
write.csv(
  x = data_dict, file = file.path(path_data, "processed_data/dictionary.csv"),
  row.names = FALSE)

write.csv(
  x = data_pivotted,
  file = file.path(path_data, "processed_data/all_protein__pivotted.csv"),
  row.names = FALSE)

write.csv(
  x = data_filter,
  file = file.path(path_data, "processed_data/data_filter.csv"),
  row.names = FALSE
)
