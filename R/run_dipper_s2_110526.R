lib_path <- "/projappl/project_2017487/project_rpackages_r452"
.libPaths(c(lib_path, .libPaths()))

library(tidyverse)
library(cmdstanr)
library(DiPPER)

set_cmdstan_path("/appl/soft/math/r-env/452-stan/cmdstan-2.38.0")

c_args <- commandArgs(trailingOnly = TRUE)
data_file <- paste0("data_090625_counts_", c_args[1], ".rds")
save_file <- paste0("res_dipper_s2_110526_", c_args[1], ".rds")

# Load data
load(data_file)
datasets <- keep(
  data_10, ~ str_detect(.$meta$data_id[1], "original|null")
)

# Run DiPPER on datasets -------------------------------------------------------
res_list <- list()

for(j in 1:length(datasets)) {

  # Count data
  counts_raw <- datasets[[j]]$assay
  
  # Filter out taxa that are present in ALL samples (100% prevalence)
  taxa_to_keep <- colSums(counts_raw > 0) < nrow(counts_raw)
  counts <- counts_raw[, taxa_to_keep]
  
  # Meta data
  meta <- datasets[[j]]$meta
  
  # Run DiPPER
  n_chains <- 4
  n_sampling <- 1000
  n_warmup <- 1000
  
  st <- Sys.time()

  fit_dipper <- DiPPER(
    meta = meta,
    abund_matrix = t(counts),
    formula = ~ group,
    symmetric = F,
    assay_type = "counts",
    pa_threshold = 0,
    min_present = 0,
    min_absent = 0,
    chains = n_chains,
    cores = n_chains,
    iter_warmup = n_warmup,
    iter_sampling = n_sampling,
    seed = 1,
    run_diagnostics = F,
    print_progress = 0,
    prior_tau_sd = 0.5,
    prior_nu_sd = 0.05
  )
  
  et <- Sys.time()

  # Extract the raw CmdStanMCMC object
  fit <- fit_dipper$stanfit
  
  # Convergence diagnostics and hyper parameters -------------------------------
  all_summaries <- fit$summary()
  
  diags <- data.frame(
    method = "DiPPER_s2",
    data_id = meta$data_id[1],
    total_divergences = sum(fit$diagnostic_summary()$num_divergent),
    max_rhat = max(all_summaries$rhat, na.rm = TRUE),
    min_ess_bulk = min(all_summaries$ess_bulk, na.rm = TRUE),
    runtime = as.numeric(difftime(et, st, units = "mins"))
  )
  
  diags_beta <- fit$summary(variables = "beta") |>
    select(variable, rhat, ess_bulk, ess_tail) |>
    mutate(taxon = colnames(counts)) |>
    select(-variable)
  
  res_hyper = fit$summary(
    variables = c("nu", "tau"),
    est = ~ quantile(.x, probs = 0.500),
    lwr95 = ~ quantile(.x, probs = 0.050),
    upr95 = ~ quantile(.x, probs = 0.950),
    rhat = posterior::rhat,
    ess_bulk = posterior::ess_bulk) |>
    select(variable, est = "50%", lwr90 = "5%", upr90 = "95%",
           rhat, ess_bulk) |>
    mutate(data_id = meta$data_id[1],
           method = "DiPPER_s2")
  
  # Posterior quantities for differential prevalence estimates -----------------
  beta_draws <- fit$draws(variables = "beta", format = "matrix")
  stopifnot(ncol(beta_draws) == ncol(counts))
  colnames(beta_draws) <- colnames(counts)
  
  res_beta <- beta_draws |>
    as.data.frame() |>
    pivot_longer(cols = everything(), names_to = "taxon") |>
    group_by(taxon) |>
    summarize(
      est = quantile(value, .500, na.rm = TRUE),
      lwr95 = quantile(value, .025, na.rm = TRUE),
      lwr90 = quantile(value, .050, na.rm = TRUE),
      lwr80 = quantile(value, .100, na.rm = TRUE),
      upr80 = quantile(value, .900, na.rm = TRUE),
      upr90 = quantile(value, .950, na.rm = TRUE),
      upr95 = quantile(value, .975, na.rm = TRUE),
      prob_low = mean(value < 0, na.rm = TRUE),
      prob_hi = mean(value > 0, na.rm = TRUE),
      .groups = "drop"
    ) |>
    full_join(data.frame(taxon = colnames(counts_raw)), by = "taxon") |>
    left_join(diags_beta, by = "taxon") |>
    mutate(method = "DiPPER_s2",
           data_id = meta$data_id[1])
  
  # Summary of results and model diagnostics -----------------------------------
  res_list[[j]] <- list( 
    diags = diags,
    res_hyper = res_hyper,
    res_beta = res_beta
  )
  
  cat(paste0("---- Data number ", j, " from set ", c_args[1],
             " (", meta$data_id[1], ") ----\n"))
  cat(paste0(
    "Divergences: ", diags$total_divergences,
    "\nMax R-hat: ", round(diags$max_rhat, 3),
    "\nMin ESS bulk: ", round(diags$min_ess_bulk, 0),
    "\nMin beta ESS tail: ", round(min(diags_beta$ess_tail, na.rm = TRUE), 0),
    "\nRuntime: ", floor(diags$runtime), "min ",
    round((diags$runtime %% 1) * 60), "s\n"
  ))
  
  res_beta |> summarize(
    n_sign_05 = sum(lwr95 > 0 | upr95 < 0, na.rm = TRUE),
    n_sign_10 = sum(lwr90 > 0 | upr90 < 0, na.rm = TRUE),
    n_sign_20 = sum(lwr80 > 0 | upr80 < 0, na.rm = TRUE)
  ) |> print()
  
  # Save results
  save(res_list, file = save_file)
}