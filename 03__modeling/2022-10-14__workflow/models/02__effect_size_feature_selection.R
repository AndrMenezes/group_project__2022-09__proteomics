#' 2022-10-14: Perform proteins selection based on Hedges' g effect size.
#' We use the `hedges_g` function from {effectsize} package.
suppressMessages(library(QFeatures))

local_path <- "./03__modeling/2022-10-14__workflow"
path_data <- file.path(local_path, "processing", "processed_data")
path_res <- file.path(local_path, "models", "results")

# Importing data --------------------------------------------------------------
fts <- readRDS(file.path(path_data, "fts_processsed.rds"))

# Calculating Hedges G --------------------------------------------------------
ampicillin <- apply(
  assay(fts[["proteins_median"]], i = "log2_normalized"), 1,
  function(z)  effectsize::hedges_g(z[1:3], z[10:12])$Hedges_g)

cefotaxime <- apply(
  assay(fts[["proteins_median"]], i = "log2_normalized"), 1,
  function(x) effectsize::hedges_g(x[4:6], x[10:12])$Hedges_g)

ciprofloxacin <- apply(
  assay(fts[["proteins_median"]], i = "log2_normalized"), 1,
  function(x) effectsize::hedges_g(x[7:9], x[10:12])$Hedges_g)

impipenem <- apply(
  assay(fts[["proteins_median"]], i = "log2_normalized"), 1,
  function(x) effectsize::hedges_g(x[13:15], x[10:12])$Hedges_g)

# Organize the results in a data frame
effect_sizes <- data.frame(ampicillin, cefotaxime, ciprofloxacin, impipenem)
colnames(effect_sizes) <- paste0("hedges_g__", colnames(effect_sizes))
head(effect_sizes)

# Concatenating
if (all.equal(rownames(fts[["proteins_median"]]), rownames(effect_sizes)))
  rowData(fts[["proteins_median"]]) <- cbind(rowData(fts[["proteins_median"]]), effect_sizes)


# Saving ------------------------------------------------------------------
saveRDS(object = fts, file = file.path(path_data, "fts_processsed.rds"))

