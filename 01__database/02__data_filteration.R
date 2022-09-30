# Renaming columns and creating filtered data set
library(dplyr)
library(readxl)
library(writexl)

path_data <- "./01__database/"

data_raw <- read_xlsx(
  file.path(path_data, "raw_data/raw_data.xlsx"),
  sheet = "All protein IDs",
  range = "A2:Y947")

# head(data_raw)
# glimpse(data_raw)

# Renaming ----------------------------------------------------------------

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
    lqf_Ampicillin_1 = `LFQ intensity Ampicillin Rep 1`,
    lqf_Ampicillin_2 = `LFQ intensity Ampicillin Rep 2`,
    lqf_Ampicillin_3 = `LFQ intensity Ampicillin Rep 3`,
    lqf_Cefotaxime_1 = `LFQ intensity Cefotaxime Rep 1`,
    lqf_Cefotaxime_2 = `LFQ intensity Cefotaxime Rep 2`,
    lqf_Cefotaxime_3 = `LFQ intensity Cefotaxime Rep 3`,
    lqf_Impipenem_1 = `LFQ intensity Impipenem Rep 1`,
    lqf_Impipenem_2 = `LFQ intensity Impipenem Rep 2`,
    lqf_Impipenem_3 = `LFQ intensity Impipenem Rep 3`,
    lqf_Ciprofloxacin_1 = `LFQ intensity Ciprofloxacin Rep 1`,
    lqf_Ciprofloxacin_2 = `LFQ intensity Ciprofloxacin Rep 2`,
    lqf_Ciprofloxacin_3 = `LFQ intensity Ciprofloxacin Rep 3`,
    lqf_Control_1 = `LFQ intensity Control Rep 1`,
    lqf_Control_2 = `LFQ intensity Control Rep 2`,
    lqf_Control_3 = `LFQ intensity Control Rep 3`,
  )

# Filtering ----------------------------------------------------------------
data_filter <- data_raw |>
  select(
    protein__name,
    lqf_Ampicillin_1,
    lqf_Ampicillin_2,
    lqf_Ampicillin_3,
    lqf_Cefotaxime_1,
    lqf_Cefotaxime_2,
    lqf_Cefotaxime_3,
    lqf_Impipenem_1,
    lqf_Impipenem_2,
    lqf_Impipenem_3,
    lqf_Ciprofloxacin_1,
    lqf_Ciprofloxacin_2,
    lqf_Ciprofloxacin_3,
    lqf_Control_1,
    lqf_Control_2,
    lqf_Control_3
  )

# View(data_filter)

# Exporting ----------------------------------------------------------------

write_xlsx(
  x = data_filter, 
  file.path(path_data, "processed_data/data_filter.xlsx"))


