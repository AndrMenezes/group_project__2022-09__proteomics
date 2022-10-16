#' 2022-10-14: Perform proteins selection based on var-mean relationship.
#' We use the `ModelGeneVar` function from {scran} package.
suppressMessages(library(SummarizedExperiment))
suppressMessages(library(HDF5Array))
suppressMessages(library(scran))
library(magrittr)
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

se <- loadHDF5SummarizedExperiment(dir = file.path(path_data, "se_obj"))
se

# Modeling the relationship mean vs variance of proteins ------------------

dec <- modelGeneVar(x = se, assay.type = "log_intensity_normalized")
dec_ordered <- dec[order(dec$bio, decreasing = TRUE), ] 
dec_ordered[dec_ordered$bio > 0.01, ]
dec_ordered[dec_ordered$FDR <= 0.10, ]

# Visualizing the model ---------------------------------------------------

p <- ggplot(dplyr::as_tibble(dec), aes(x = mean, y = total)) +
  geom_point(size = 2, alpha = 0.5) +
  geom_line(aes(y = tech), size = 1.5, col = "blue") +
  scale_x_continuous(breaks = scales::pretty_breaks(8)) +
  scale_y_continuous(breaks = scales::pretty_breaks(8)) +
  labs(x = "Mean of the log(intensity) normalized",
       y = "Variance of the log(intensity)\n normalized")
save_plot(filename = file.path(path_res, "var_mean.png"), plot = p,
          base_width = 8, bg = "white")



# Concatenate the decomposition into SE object ----------------------------
if (all.equal(rownames(se), rownames(dec)))
  rowData(se) <- cbind(rowData(se), dec)


# Updating the saved SE object with new information -----------------------
quickResaveHDF5SummarizedExperiment(se)
