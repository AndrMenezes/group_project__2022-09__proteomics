rm(list = ls())
library(dplyr)
library(ggplot2)
library(cowplot)
theme_set(
  theme_cowplot() +
    background_grid()
)

res_path <- "./03__modeling/2022-09-29__simulation_studies/models/results"

###############################################################################
# True parameter values
n_samples <- c(3, 5, 10, 15, 20, 50)
v_mu <- c(0, 1, 3)
v_sd <- c(1, 2, 4)
p <- length(v_sd)
D <- diag(v_sd)
list_correlations <- list(
  scenario_1 = c(0.25, 0.25, 0.25),
  scenario_2 = c(0.50, 0.50, 0.50),
  scenario_3 = c(0.85, 0.85, 0.85)
)
M <- 5e4
R <- diag(p)
n_max <- max(n_samples)
###############################################################################
# Function to calculate the pearson correlation coefficient
compute_corr <- function(v_y) {
  p <- nrow(v_y)
  comb2v2 <- combn(p, 2)
  lt <- list()
  for (i in seq_len(ncol(comb2v2))) {
    
    tmp <- cor.test(v_y[comb2v2[, i][1L], ],
                    v_y[comb2v2[, i][2L], ],
                    method = "pearson",
                    alternative = "two.sided")
    lt[[i]] <- data.frame(
      comb = paste0(comb2v2[, i], collapse = " vs "),
      rho = tmp$estimate,
      pvalue = tmp$p.value,
      ci_lower = if (is.null(tmp$conf.int[1L])) 0 else tmp$conf.int[1L],
      ci_upper = if (is.null(tmp$conf.int[2L])) 0 else tmp$conf.int[2L])
  }
  out <- do.call(rbind, lt)
  rownames(out) <- NULL
  out
}

###############################################################################
# Perform the simulation study
if (!file.exists(file.path(res_path, "simulated_results.rds"))) {
  set.seed(6669)
  counter <- 1L
  list_results <- list()
  for (rho in list_correlations) {
    
    cat(rho, "\n")
    
    # Define the new correlation matrix
    R[lower.tri(R)] <- rho
    R[upper.tri(R)] <- rho
    # Compute covariance matrix from correlation matrix
    V <- (D %*% R) %*% D
    # Get Cholesky decomposition to generate MVN distribution
    P <- t(chol(V))
    
    # Compute correlation and p-value for different sample size
    ini <- proc.time()
    out <- replicate(n = M, expr = {
      z <- matrix(rnorm(n_max * p), nrow = p)
      y <- P %*% z + v_mu
      do.call(rbind, lapply(n_samples, function(j) {
        cbind(compute_corr(v_y = y[, seq_len(j)]), n_sample = j)
      }))
    }, simplify = FALSE)
    end <- proc.time()
    
    cat("\n Finished in", (end[3L] - ini[3L]) / 60, "minutes\n")  
    list_results[[counter]] <- do.call(rbind, out)
    counter <- counter + 1L
  }
  names(list_results) <- names(list_correlations)
  data_results <- bind_rows(list_results, .id = "scenario") |> 
    as_tibble()
  saveRDS(data_results, file = file.path(res_path, "simulated_results.rds"))
}
data_results <- readRDS(file = file.path(res_path, "simulated_results.rds"))

# Creating data.frame with the true values
data_true_values <- bind_rows(list_correlations) |> 
  mutate(comb = c("1 vs 2", "1 vs 3", "2 vs 3")) |> 
  tidyr::pivot_longer(cols = -comb, names_to = "scenario",
                      values_to = "true_rho")

# Merging to get the true values for each scenario
data_results <- data_results |> 
  left_join(data_true_values, by = c("scenario", "comb")) |> 
  mutate(n_sample = forcats::fct_relevel(
    factor(n_sample), "3", "5", "10", "15"))

###############################################################################
# Plotting the distribution of estimates
counter <- 1L
list_plots <- list()
for (chosen in unique(data_results$scenario)) {
  data_curr <- data_results |> 
    filter(scenario == chosen)
  
  aux <- list("scenario_1" = "Low (0.25)",
              "scenario_2" = "Moderate (0.50)",
              "scenario_3" = "High (0.85)")[[chosen]]
  
  list_plots[[counter]] <- ggplot(data_curr, aes(x = n_sample, y = rho,
                                                 group = n_sample)) +
    facet_wrap(~ comb) +
    geom_boxplot() +
    geom_hline(data = filter(data_true_values, scenario == chosen),
               aes(yintercept = true_rho), col = "red") +
    labs(x = "Sample size", y = expression(rho)) +
    ggtitle(expression("MC distribution of"~rho),
            paste0("Scenario: ", aux, " correlation")) +
    scale_y_continuous(breaks = scales::pretty_breaks(6))
  
  counter <- counter + 1L
}
p_grid <- plot_grid(plotlist = list_plots)
save_plot(filename = file.path(res_path, "boxplot_estimates.png"),
          plot = p_grid, bg = "white", base_height = 8)

###############################################################################
# Summarizing the data
data_summarized <- data_results |> 
  group_by(scenario, n_sample, comb) |> 
  summarise(mean_rho = mean(rho),
            sd_rho = sd(rho),
            relative_bias_rho = mean((rho - true_rho) / true_rho),
            true_rho = mean(true_rho),
            test_size = mean(pvalue <= 0.05),
            coverage = mean(true_rho >= ci_lower & true_rho <= ci_upper),
            .groups = "drop") |> 
  mutate(n_sample = forcats::fct_relevel(
           factor(n_sample), "3", "5", "10", "15"),
         comb_rho = paste0(comb, "  ", expression(rho), "=", true_rho))

for (chosen in unique(data_summarized$scenario)) {
  data_curr <- data_summarized |> 
    filter(scenario == chosen)
  
  p_mean <- ggplot(data_curr, aes(x = n_sample, y = mean_rho)) +
    facet_wrap(. ~ comb_rho, scales = "free_y") +
    geom_point(size = 3) +
    geom_hline(data = filter(data_true_values, scenario == chosen),
               aes(yintercept = true_rho), col = "red") +
    labs(x = "Sample size", y = "Mean") +
    ggtitle("Mean of the estimates correlation coeficient across the Monte Carlo replicates") +
    scale_y_continuous(breaks = scales::pretty_breaks(6))
  
  p_bias <- ggplot(data_curr, aes(x = n_sample, y = relative_bias_rho)) +
    facet_wrap(. ~ comb, scales = "free_y") +
    geom_point(size = 4) +
    geom_hline(yintercept = 0, col = "red") +
    labs(x = "Sample size", y = "Relative bias") +
    ggtitle("Relative bias of the estimates correlation coeficient across the Monte Carlo replicates") +
    scale_y_continuous(breaks = scales::pretty_breaks(6))

  p_sd <- ggplot(data_curr, aes(x = n_sample, y = sd_rho)) +
    facet_wrap(. ~ comb, scales = "free_y") +
    geom_point(size = 4) +
    labs(x = "Sample size", y = "Standard Deviation") +
    ggtitle("Standard deviation of the estimates correlation coeficient across the Monte Carlo replicates") +
    scale_y_continuous(breaks = scales::pretty_breaks(6))

  save_plot(
    filename = file.path(res_path, paste0("mean_estimate__", chosen, ".png")),
    plot = p_mean,
    base_width = 12,
    bg = "white")
  save_plot(
    filename = file.path(res_path, paste0("bias_estimate__", chosen, ".png")),
    plot = p_bias,
    base_width = 12,
    bg = "white")
  save_plot(
    filename = file.path(res_path, paste0("std_estimate__", chosen, ".png")),
    plot = p_sd,
    base_width = 12,
    bg = "white")
  
}
