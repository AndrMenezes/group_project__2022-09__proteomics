#' 2022-10-14: Estimate and visualize the correlation between proteins.
#' Before compute the correlation we chose the highly variable proteins using
#' the analysis perform in script `01__proteins_selection.R`
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




# Estimating and visualizing the protein correlation ----------------------


biplot_correlation <- function(se, chosen_proteins, label = TRUE) {
  
  # Selecting and computing the correlation
  m__abundance <- assay(se[chosen_proteins, ], i = 2L)
  m__cor <- cor(as.matrix(t(m__abundance)))

  # Re-scaling the correlation matrix
  mu_row <- rowMeans(m__cor)
  mu_col <- colMeans(m__cor)
  mu <- mean(m__cor)
  E <- sweep(m__cor, 1, mu_row)
  E <- sweep(E, 2, mu_col)
  E <- E + mu
  
  # Spectral decomposition 
  dec <- eigen(E, symmetric = TRUE)
  kp <- max(abs(dec$values)) / min(abs(dec$values))
  l1 <- dec$values[1L]
  l2 <- dec$values[2L]
  tot <- sum(abs(dec$values))
  pct_retain <- 100 * (l1 + l2) / tot
  
  # Plotting
  data_dec <- data.frame(
    protein = rownames(m__cor),
    e1 = dec$vectors[, 1L], e2 = dec$vectors[, 2L])
  p <- ggplot(data = data_dec) +
    aes(x = e1, y = e2) +
    geom_point(size = 1) +
    labs(x = paste0("Component 1 (",  round(100 * l1 / tot, 1), "%)"),
         y = paste0("Component 2 (",  round(100 * l2 / tot, 1), "%)")) +
    scale_x_continuous(breaks = scales::pretty_breaks(6)) +
    scale_y_continuous(breaks = scales::pretty_breaks(6))
  if (label) 
    p <- p + geom_label(aes(label = protein))
  
  
  list(pct_retain = pct_retain, plot = p)
}

chosen__bio <- rowData(se)$bio > 0.01
chosen__fdr <- rowData(se)$FDR < 0.01
sum(chosen__bio) / nrow(se)
sum(chosen__fdr) / nrow(se)

lt__bio <- biplot_correlation(se = se, chosen_proteins = chosen__bio,
                              label = FALSE)
lt__fdr <- biplot_correlation(se = se, chosen_proteins = chosen__fdr,
                              label = FALSE)

# Saving the plots
p__bio <- lt__bio$plot +
  ggtitle("Biplot filtering proteins with estimated biological variance greater than 0.01",
          subtitle = paste0(round(lt__bio$pct_retain, 2), "% variation explained"))
save_plot(filename = file.path(path_res, "biplot__bio_selection.png"),
          plot = p__bio, base_width = 8, bg = "white")
p__fdr <- lt__fdr$plot +
  ggtitle("Biplot filtering proteins with FDR lesser than 0.01",
          subtitle = paste0(round(lt__fdr$pct_retain, 2), "% variation explained"))
save_plot(filename = file.path(path_res, "biplot__fdr_selection.png"),
          plot = p__fdr, base_width = 8, bg = "white")


lt__fdr <- biplot_correlation(se = se, chosen_proteins = chosen__fdr,
                              label = TRUE)
p__fdr <- lt__fdr$plot +
  ggtitle("Biplot filtering proteins with FDR lesser than 0.01",
          subtitle = paste0(round(lt__fdr$pct_retain, 2), "% variation explained"))
save_plot(filename = file.path(path_res, "biplot__fdr_selection_labeled.png"),
          plot = p__fdr, base_width = 18, base_height = 8, bg = "white")



