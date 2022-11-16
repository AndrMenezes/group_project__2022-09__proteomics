#' 2022-10-14: Perform proteins selection based on var-mean relationship.
#' We use the `modelGeneVar` function from {scran} package.
library(ggplot2)
library(cowplot)
theme_set(
  theme_cowplot() +
    background_grid() +
    theme(legend.position = "top")
)

# Path -------------------------------------------------------------------

local_path <- "./03__modeling/2022-10-14__workflow"
path_data <- file.path(local_path, "processing", "processed_data")
path_res <- file.path(local_path, "models", "results")


# Load data ---------------------------------------------------------------
se_margalit <- readRDS("./01__database/processed_data/se_margalit.rds")
fts <- readRDS(file.path(path_data, "fts_processsed.rds"))

# Modeling the relationship mean vs variance of proteins ------------------
dec_ours <- scran::modelGeneVar(x = fts[["proteins_median"]],
                                assay.type = "log2_normalized")
dec_margalit <- scran::modelGeneVar(x = se_margalit,
                                    assay.type = "log_intensity")

dec_ours[order(dec_ours$bio, decreasing = TRUE), ]
dec_margalit[order(dec_margalit$bio, decreasing = TRUE), ]

# Visualizing the model ---------------------------------------------------
plot_mean_variance <- function(dec) {
  ggplot(dplyr::as_tibble(dec), aes(x = mean, y = total)) +
    geom_point(size = 2, alpha = 0.5) +
    geom_line(aes(y = tech), size = 1.5, col = "blue") +
    scale_x_continuous(breaks = scales::pretty_breaks(8)) +
    scale_y_continuous(breaks = scales::pretty_breaks(8)) +
    labs(x = "Mean of the log(intensity) normalized",
         y = "Variance of the log(intensity) normalized")
}
plot_mean_variance(dec = dec_ours)
plot_mean_variance(dec = dec_margalit)

chosen <- ((dec_ours$mean < 22) | (dec_ours$mean > 29))
to_remove <- rownames(dec_ours)[chosen]
print(to_remove)
table(chosen)
rowData(fts[["proteins"]][to_remove, ])
assay(fts[["proteins"]][to_remove, ], "intensity")
assay(fts[["proteins"]][to_remove, ], "log2_normalized")

# Removing the two proteins with lower (greater) mean (m < 22 and m > 29) 
chosen <- !(rownames(fts[["proteins"]]) %in% to_remove)
fts[["proteins"]] <- fts[["proteins"]][chosen, ]

dec_ours_2 <- scran::modelGeneVar(x = fts[["proteins_median"]],
                                  assay.type = "log2_normalized")
plot_mean_variance(dec = dec_ours_2)

# Concatenate the decomposition into SE object ----------------------------
if (all.equal(rownames(se_margalit), rownames(dec_margalit)))
  rowData(se_margalit) <- cbind(rowData(se_margalit), dec_margalit)
if (all.equal(rownames(fts[["proteins_median"]]), rownames(dec_ours_2)))
  rowData(fts[["proteins_median"]]) <- cbind(
    rowData(fts[["proteins_median"]]), dec_ours_2)


# Save the objects --------------------------------------------------------
saveRDS(object = fts, file = file.path(path_data, "fts_processsed.rds"))
saveRDS(object = se_margalit, file = file.path(path_data, "se_margalit.rds"))
