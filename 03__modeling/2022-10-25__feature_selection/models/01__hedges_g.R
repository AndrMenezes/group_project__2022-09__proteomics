suppressMessages(library(QFeatures))
library(dplyr)
library(tidyr)
library(effectsize)

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

# Calculating Hedges G --------------------------------------------------------
ampicillin <- apply(normalised_proteins, 1, function(x) 
  hedges_g(x[1:3], x[10:12])$Hedges_g)

cefotaxime <- apply(normalised_proteins, 1, function(x)
  hedges_g(x[4:6], x[10:12])$Hedges_g)

ciprofloxacin <- apply(normalised_proteins, 1, function(x)
  hedges_g(x[7:9], x[10:12])$Hedges_g)

impipenem <- apply(normalised_proteins, 1, function(x)
  hedges_g(x[13:15], x[10:12])$Hedges_g)

# Creating matrix with all values
effect_sizes <- data.frame(ampicillin, cefotaxime, ciprofloxacin, impipenem)
effect_sizes <- as.matrix(effect_sizes)

######################################
#       Hedges G Formula
######################################

# Hedges G = Cohen's D * (1 - (3 / 4 * (n1 + n1) - 9))

# Cohen's D   = (M1 - M2) / SD_pooled
# n1          = Group 1 Sample Size
# n2          = Group 2 Sample Size
# -----------------------------------------------------------------------------

# Cohen's D   = (M1 - M2) / SD Pooled
# M1          = Group 1 Mean
# M2          = Group 2 Mean
# SD Pooled   = Pooled and Weighted Standard Deviation 

# Pooled Standard Deviation formula
# SD Pooled = sqrt((sd1)^2 + (sd2)^2) / 2)

# sd1   = Group 1 Standard Deviation
# sd2   = Group 2 Standard Deviation
