#' 2022-10-04: Compute the correlation between proteins.
rm(list = ls())
library(dplyr)
library(ggplot2)
library(cowplot)
theme_set(
  theme_cowplot() +
    background_grid()
)

path_local <- "./03__modeling/2022-10-04__correlation_analysis"
path_data <- file.path(path_local, "/processing/processed_data/")
path_res <- file.path(path_local, "models/results")


data_wide_protein <- readRDS(file.path(path_data, "data_wide_protein.rds"))
n <- nrow(data_wide_protein)

# Correlation considering all groups --------------------------------------
cor_all <- cor(as.matrix(data_wide_protein[, -(1:2)]), method = "pearson")
dim(cor_all)
p <- dim(cor_all)[1L]
L <- lower.tri(cor_all)
rc_index <- which(L, arr.ind = TRUE)
data_cor <- tibble(
  protein__A = dimnames(cor_all)[[1L]][rc_index[, 1L]],
  protein__B = dimnames(cor_all)[[2L]][rc_index[, 2L]],
  correlation = cor_all[L] 
)
choose(p, 2) == dim(data_cor)[1L]

# Half-normal plot --------------------------------------------------------
n <- nrow(data_wide_protein)
k <- choose(p, 2)
data_cur <- data_cor |> 
  as_tibble() |> 
  mutate(z_obs = abs(atanh(correlation))) |> 
  arrange(z_obs) |> 
  mutate(z_teo = qnorm((1:k + k - 1/8) / (2 * k + 0.5)),
         z_teo_s = z_teo * 1 / sqrt(n - 3))
p <- ggplot(data_cur, aes(x = z_teo, y = z_obs)) +
  geom_point(size = 1.3, alpha = 0.6) +
  geom_abline(slope = 1 / sqrt(n - 3), intercept = 0) +
  scale_x_continuous(breaks = scales::pretty_breaks(6)) +
  scale_y_continuous(breaks = scales::pretty_breaks(6)) +
  ggtitle(paste0("Half-normal plot of all ", prettyNum(k, big.mark = ","),
                 " correlation coefficients")) +
  labs(x = "Theoretical values", y = "|atanh(r)|")
save_plot(filename = file.path(path_res, "half_normal_plot__all.png"),
          plot = p, base_width = 8, bg = "white")



# Bi-plot -----------------------------------------------------------------

cor_all <- cor(as.matrix(data_wide_protein[, (3:23)]), method = "pearson")
dim(cor_all)
mu_row <- rowMeans(cor_all)
mu_col <- colMeans(cor_all)
mu <- mean(cor_all)
E <- sweep(cor_all, 1, mu_row)
E <- sweep(E, 2, mu_col)
E <- E + mu

dec <- eigen(E, symmetric = TRUE)
max(abs(dec$values)) / min(abs(dec$values))

plot(dec$vectors[, 1], dec$vectors[, 2])

