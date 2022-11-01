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

# Hedges G = Cohen's D * (1 - (3 / (4 * (n1 + n1) - 9)))

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

######################################
#           Sanity Check
######################################
# Using Formula on Protein 1 for Ampicillin vs. Control ------------------------
# Getting means
m_1 <- mean(normalised_proteins[1, 1:3])
m_2 <- mean(normalised_proteins[1, 10:12])

# Getting standard deviations & squaring
sd_1 <- sd(normalised_proteins[1, 1:3])
sd_2 <- sd(normalised_proteins[1, 10:12])

sd_1 <- sd_1^2
sd_2 <- sd_2^2

# Calculating pooled standard deviation
pooled_sd <- sqrt((sd_1 + sd_2)/2)

# Calculating Cohen's D
D <- (m_1 - m_2) / pooled_sd

# Applying Hedges Correction
G <- D * (1 - (3/(4*(3+3) - 9)))

G

# Using Hedges G function on Protein 1 for Ampicillin vs. Control -------------
hedges_g(normalised_proteins[1, 1:3], normalised_proteins[1, 10:12])$Hedges_g
