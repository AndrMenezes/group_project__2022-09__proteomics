#' 2022-10-14: Perform proteins selection based on var-mean relationship.
#' We use the `ModelGeneVar` function from {scran} package.
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
dec_ours <- scran::modelGeneVar(x = fts[["proteins"]],
                                assay.type = "log_intensity_normalized")
dec_margalit <- scran::modelGeneVar(x = se_margalit,
                                    assay.type = "log_intensity")

dec_ours[order(dec_ours$bio, decreasing = TRUE), ]
dec_margalit[order(dec_margalit$bio, decreasing = TRUE), ]


# Concatenate the decomposition into SE object ----------------------------
if (all.equal(rownames(se_margalit), rownames(dec_margalit)))
  rowData(se_margalit) <- cbind(rowData(se_margalit), dec_margalit)
if (all.equal(rownames(fts[["proteins"]]), rownames(dec_ours)))
  rowData(fts[["proteins"]]) <- cbind(rowData(fts[["proteins"]]), dec_ours)

# Visualizing the model ---------------------------------------------------
ggplot(dplyr::as_tibble(rowData(fts[["proteins"]])),
       aes(x = mean, y = total)) +
  geom_point(size = 2, alpha = 0.5) +
  geom_line(aes(y = tech), size = 1.5, col = "blue") +
  scale_x_continuous(breaks = scales::pretty_breaks(8)) +
  scale_y_continuous(breaks = scales::pretty_breaks(8)) +
  labs(x = "Mean of the log(intensity) normalized",
       y = "Variance of the log(intensity)\n normalized")

chosen <- ((rowData(fts[["proteins"]])$mean < 22) |
             (rowData(fts[["proteins"]])$mean > 29))
table(chosen)
rowData(fts[["proteins"]][chosen, ])
assay(fts[["proteins"]][chosen, ], "log_intensity_normalized")

dec_ours_2 <- scran::modelGeneVar(x = fts[["proteins"]][!chosen, ],
                                  assay.type = "log_intensity_normalized")

ggplot(dplyr::as_tibble(dec_ours_2),
       aes(x = mean, y = total)) +
  geom_point(size = 2, alpha = 0.5) +
  geom_line(aes(y = tech), size = 1.5, col = "blue") +
  scale_x_continuous(breaks = scales::pretty_breaks(8)) +
  scale_y_continuous(breaks = scales::pretty_breaks(8)) +
  labs(x = "Mean of the log(intensity) normalized",
       y = "Variance of the log(intensity)\n normalized")

# Save the objects --------------------------------------------------------
saveRDS(object = fts, file = file.path(path_data, "fts_processsed.rds"))
saveRDS(object = se_margalit, file = file.path(path_data, "se_margalit.rds"))
