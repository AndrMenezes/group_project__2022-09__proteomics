suppressMessages(library(QFeatures))
library(dplyr)
library(tidyr)
library(effsize) # For calculating Hedges G

path_data <- "./03__modeling/2022-10-14__workflow/processing/processed_data"

# Importing data --------------------------------------------------------------
# First we import the QFeatures Data (all of it)
fts <- readRDS(file.path(path_data, "fts_processed.rds"))

# We can then focus at the protein level (fts[["proteins"]])
se_protein <- fts[["proteins"]]
colData(se_protein) <- colData(fts)

# Finally, we can focus at the normalised protein data
normalised_proteins <- assay(x = se_protein, i = "log_intensity_normalized")
colnames(normalised_proteins) <- se_protein@colData@listData[["group"]]

# Calculation of effect sizes -------------------------------------------------
ampicillin <- apply(normalised_proteins, 1, function(x)
  cohen.d(x[1:3], x[10:12], hedges = TRUE)$estimate
)
cefotaxime <- apply(normalised_proteins, 1, function(x)
  cohen.d(x[4:6], x[10:12], hedges = TRUE)$estimate
)

ciprofloxacin <- apply(normalised_proteins, 1, function(x)
  cohen.d(x[7:9], x[10:12], hedges = TRUE)$estimate
)

impipenem <- apply(normalised_proteins, 1, function(x)
  cohen.d(x[13:15], x[10:12], hedges = TRUE)$estimate
)

# Creating matrix with all values
effect_sizes <- data.frame(ampicillin, cefotaxime, ciprofloxacin, impipenem)
effect_sizes <- as.matrix(effect_sizes)