rm(list = ls())
library(dplyr)
library(ggplot2)
library(cowplot)
library(parallel)
theme_set(
  theme_cowplot() +
    background_grid() +
    theme(legend.position = "top")
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
  n <- ncol(v_y)
  comb2v2 <- combn(p, 2)
  lt <- list()
  for (i in seq_len(ncol(comb2v2))) {
    tmp <- cor.test(v_y[comb2v2[, i][1L], ],
                    v_y[comb2v2[, i][2L], ],
                    method = "pearson",
                    alternative = "two.sided")
    r <- tmp$estimate
    r_adj <- r * (1 + (1 - r^2) / 2 / n)
    lt[[i]] <- data.frame(
      comb = paste0(comb2v2[, i], collapse = " vs "),
      rho = r,
      rho_adj = r_adj,
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
  list_results <- mclapply(X = list_correlations, FUN = function(rho) {

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
    # Tracking the scenarios
    write.table(
      x = data.frame(paste0(rho, collapse = ", "), (end[3L] - ini[3L]) / 60),
      file = file.path(res_path, "tracking.txt"),
      sep = "\t",
      append = TRUE,
      col.names = FALSE,
      row.names = FALSE,
      quote = FALSE)
    return(do.call(rbind, out))
  }, mc.cores = 3L)
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
    filter(scenario == chosen) |>
    select(comb, n_sample, rho, rho_adj) |>
    tidyr::pivot_longer(cols = -c(comb, n_sample))

  aux <- list("scenario_1" = "Low (0.25)",
              "scenario_2" = "Moderate (0.50)",
              "scenario_3" = "High (0.85)")[[chosen]]

  list_plots[[counter]] <- ggplot(data_curr, aes(x = n_sample, y = value,
                                                 fill = name)) +
    facet_wrap(~ comb) +
    geom_boxplot() +
    geom_hline(data = filter(data_true_values, scenario == chosen),
               aes(yintercept = true_rho), col = "red", size = 1) +
    labs(x = "Sample size", y = "Estimate", fill = "") +
    ggtitle(expression("MC distribution of" ~ rho),
            paste0("Scenario: ", aux, " correlation")) +
    scale_y_continuous(breaks = scales::pretty_breaks(6)) +
    colorspace::scale_fill_discrete_qualitative(
      labels = c("rho" = expression(rho), "rho_adj" = expression(rho[adj])))
  fname <- paste0("boxplot_estimates__", counter, ".png")
  save_plot(filename = file.path(res_path, fname),
            plot = list_plots[[counter]], bg = "white", base_height = 6)
  counter <- counter + 1L
}
p_grid <- plot_grid(plotlist = list_plots)
save_plot(filename = file.path(res_path, "boxplot_estimates.png"),
          plot = p_grid, bg = "white", base_height = 8)

###############################################################################
# Summarizing the data
data_summarized <- data_results |>
  select(scenario, comb, n_sample, true_rho, rho, rho_adj) |>
  tidyr::pivot_longer(cols = -c(scenario, comb, n_sample, true_rho),
                      names_to = "estimator",
                      values_to = "value") |>
  group_by(scenario, n_sample, comb, estimator) |>
  summarise(mean = mean(value),
            sd = sd(value),
            true_rho = mean(true_rho),
            relative_bias = mean((value - true_rho) / true_rho),
            .groups = "drop") |>
  mutate(n_sample = forcats::fct_relevel(
           factor(n_sample), "3", "5", "10", "15"),
         comb_rho = paste0(comb, "  ", expression(rho), "=", true_rho))

for (chosen in unique(data_summarized$scenario)) {

  aux <- list("scenario_1" = "Low (0.25)",
              "scenario_2" = "Moderate (0.50)",
              "scenario_3" = "High (0.85)")[[chosen]]

  cat(chosen, "\n")
  data_curr <- data_summarized |>
    filter(scenario == chosen)

  p_mean <- ggplot(data_curr, aes(x = n_sample, y = mean, col = estimator,
                                  group = estimator)) +
    facet_wrap(. ~ comb) +
    geom_point(size = 3) +
    geom_line(linetype = "dotted") +
    geom_hline(data = filter(data_true_values, scenario == chosen),
               aes(yintercept = true_rho), col = "red") +
    labs(x = "Sample size", y = "Mean", col = "") +
    ggtitle("Mean of the estimates correlation coeficient across the Monte Carlo replicates",
            paste0("Scenario: ", aux, " correlation")) +
    scale_y_continuous(breaks = scales::pretty_breaks(6)) +
    colorspace::scale_color_discrete_qualitative(
      labels = c("rho" = expression(rho), "rho_adj" = expression(rho[adj])))

  p_bias <- ggplot(data_curr, aes(x = n_sample, y = relative_bias,
                                  col = estimator, group = estimator)) +
    facet_wrap(. ~ comb) +
    geom_point(size = 4) +
    geom_line(linetype = "dotted") +
    geom_hline(yintercept = 0, col = "red") +
    labs(x = "Sample size", y = "Relative bias") +
    ggtitle("Relative bias of the estimates correlation coeficient across the Monte Carlo replicates",
            paste0("Scenario: ", aux, " correlation")) +
    scale_y_continuous(breaks = scales::pretty_breaks(6)) +
    colorspace::scale_color_discrete_qualitative(
      labels = c("rho" = expression(rho), "rho_adj" = expression(rho[adj])))

  p_sd <- ggplot(data_curr, aes(x = n_sample, y = sd,
                                col = estimator, group = estimator)) +
    facet_wrap(. ~ comb) +
    geom_point(size = 4) +
    geom_line(linetype = "dotted") +
    labs(x = "Sample size", y = "Standard Deviation") +
    ggtitle("Standard deviation of the estimates correlation coeficient across the Monte Carlo replicates",
            paste0("Scenario: ", aux, " correlation")) +
    scale_y_continuous(breaks = scales::pretty_breaks(6)) +
    colorspace::scale_color_discrete_qualitative(
      labels = c("rho" = expression(rho), "rho_adj" = expression(rho[adj])))

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
