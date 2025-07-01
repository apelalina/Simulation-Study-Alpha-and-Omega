
# Calculate Omega using Lavaan -------------------------------------------------

lavaan_omega <- function(Y) {

  items <- colnames(Y)
  n_items <- length(items)

  # Create loadings string: first fixed to 1, others labeled l2, l3, ...
  loadings_labels <- c(paste0("l", 1:n_items, "*"))
  loadings_string <- paste0(loadings_labels, items, collapse = " + ")

  # Residual variances string: e1*item1 ~~ item1, e2*item2 ~~ item2, ...
  res_vars_string <- paste0(items, " ~~ e", 1:n_items, "*", items, collapse = "\n")

  # Omega formula components: sum of loadings terms (1 + l2 + l3 + ...)
  loading_terms <- paste0(paste0("l", 1:n_items), collapse = " + ")
  omega_formula <- paste0("Omega := (", loading_terms, ")^2 / ((", loading_terms, ")^2 + (", paste0("e", 1:n_items, collapse = " + "), "))")

  # Build full model string
  model <- paste0("
    EM1 =~ ", loadings_string, "
    ", res_vars_string, "
    EM1 ~~ EM1
    ", omega_formula
  )

  # Fit the CFA model
  fit <- lavaan::cfa(model, as.data.frame(Y), std.lv = TRUE)

  # Extract omega estimate from defined parameters
  par_est <- lavaan::parameterEstimates(fit)
  omega_val <- par_est$est[par_est$lhs == "Omega" & par_est$op == ":="]

  return(par_est$est[par_est$lhs == "Omega" & par_est$op == ":="])
}



# Helper functions for Bootstrapping -------------------------------------------

alpha_boot <- function(data, indices) {
  result <- tryCatch({
    psych::alpha(data[indices, ], discrete = FALSE, check.keys = FALSE)$total$raw_alpha
  }, error = function(e) NA)
  return(result)
}

omega_boot <- function(data, indices) {
  lavaan_omega(data[indices, ])
}


# Simulation Function ----------------------------------------------------------

simulate_reliability_trial <- function(loading_type, n, n_items) {

  # Loadings
  if (loading_type == "tau_equivalent") {
    loadings <- rep(0.7, n_items)
  } else if (loading_type == "low_var") {
    loadings <- seq(0.6, 0.8, length.out = n_items)
  } else if (loading_type == "medium_var") {
    loadings <- seq(0.4, 0.9, length.out = n_items)
  } else if (loading_type == "high_var") {
    loadings <- seq(0.2, 0.99, length.out = n_items)
  }

  # choose error variances such that the variance is 1
  error_vars <- 1 - loadings^2

  # Simulate Responses
  eta <- rnorm(n, mean = 0, sd = 1) # Latent Variable

  Y <- sapply(1:n_items, function(j) {
    eta * loadings[j] + rnorm(n, sd = sqrt(error_vars[j]))
  })

  colnames(Y) <- paste0("item", 1:ncol(Y))

  # Observed Cronbachs Alpha (with CI)
  alpha_obs <- psych::alpha(Y, discrete = FALSE, check.keys = FALSE)$total$raw_alpha
  bootstrapped_alpha <- boot(data = Y, statistic = alpha_boot, R = N_BOOT)
  alpha_ci <- boot.ci(bootstrapped_alpha, type = "perc")$percent[4:5]


  # Observed McDonalds Omega (with CI)
  omega_obs <- lavaan_omega(Y)
  bootstrapped_omega <- boot(data = Y, statistic = omega_boot, R = N_BOOT)
  omega_ci <- boot.ci(bootstrapped_omega, type = "perc")$percent[4:5]

  # True Omega
  true_omega <- (sum(loadings))^2 / ((sum(loadings))^2 + sum(error_vars))

  return(
    c(
      "alpha" = alpha_obs,
      "alpha_ci_lower" = alpha_ci[1],
      "alpha_ci_upper" = alpha_ci[2],
      "omega" = omega_obs,
      "omega_ci_lower" = omega_ci[1],
      "omega_ci_upper" = omega_ci[2],
      "true_omega" = true_omega,
      "coverage_alpha" = if (!anyNA(alpha_ci)) (true_omega >= alpha_ci[1] && true_omega <= alpha_ci[2]) else NA,
      "coverage_omega" = if (!anyNA(omega_ci)) (true_omega >= omega_ci[1] && true_omega <= omega_ci[2]) else NA,
      "bias_alpha" = true_omega - alpha_obs,
      "bias_omega" = true_omega - omega_obs
    )
  )
}

