library(tidyverse)
library(TreeSummarizedExperiment)

pure_maaslin2 <- function(counts, meta, formula){
  
  #(TSS) Normalize counts to relative abundances
  rel_abs_0 <- counts / rowSums(counts)
  
  #Replace zeros by 1/2 * the smallest non_zero rel abundance for each taxon
  replace_zeros <- function(y){replace(y, y == 0, min(y[y > 0]) / 2)}
  
  rel_abs <- apply(rel_abs_0, 2, replace_zeros)
  
  #Log_2 transform relative abundances
  log_rel_abs <- log2(rel_abs)
  
  
  #Create model matrix
  X <- model.matrix(formula, meta)
  
  #Function to run linear model for response y and predictors x
  run_lm <- function(y, x){
    RcppEigen::fastLmPure(x, y)
  }
  
  #Run a linear model for each taxon
  lm_results <- apply(log_rel_abs, 2, run_lm, x = X)
  
  
  #Extract regression coefficients, standard errors and degrees of freedom
  coefs <- lm_results |> purrr::map_dfr(~ .$coefficients) |> as.matrix()
  ses <- lm_results |> purrr::map_dfr(~ .$se) |> as.matrix() |> t()
  
  #Calculate p values and confidence intervals
  df <- nrow(counts) - ncol(X)
  ts <- coefs / ses
  
  p_values <- 2 * (1 - pt(abs(ts), df = df))
  
  ci_lwr <- coefs - qt(.975, df = df) * ses
  ci_upr <- coefs + qt(.975, df = df) * ses
  
  
  #Combine the results
  results <- data.frame(
    taxon = rep(colnames(counts), ncol(X) - 1),
    variable = rep(colnames(X)[-1], each = ncol(counts)),
    estimate = as.numeric(coefs[, -1]),
    p_value = as.numeric(p_values[, -1]),
    ci_lwr = as.numeric(ci_lwr[, -1]),
    ci_upr = as.numeric(ci_upr[, -1]))
  
  return(results)
}


# MaAsLin2 ---------------------------------------------------------------------
run_maaslin2 <- function(rel_abs, meta, fm = ~ group){

  obj <- pure_maaslin2(counts = rel_abs,
                       meta = meta,
                       formula = fm)
  
  res <- obj %>%
    filter(variable == 'groupcase') %>%
    mutate(method = paste0('MaAsLin2'),
           q = p.adjust(p_value, method = 'BH')) %>%
    select(taxon, est = estimate, p = p_value, q, method,
           lwr = ci_lwr, upr = ci_upr)
  
  return(res %>% mutate(data_id = meta$data_id[1]))
}


#Load data and run DAA----------------------------------------------------------
load('data_090625.rds')

meta_list <- tses |> map(~ colData(.) |> as.data.frame())
assay_list <- tses |> 
  map(~ SummarizedExperiment::assay(., 'relabundance') |> t())
rm(tses)

res_maaslin2 <- map2(assay_list, meta_list, ~ run_maaslin2(.x, .y),
                     .progress = T)


#Save results and SessionInfo---------------------------------------------------
date <- format(Sys.Date(), "%d%m%y")

save(res_maaslin2, file = paste0('res_maaslin2_', date, '.rds'))

si <- capture.output(sessionInfo())
writeLines(si, con = paste0("session_info_maaslin2_", date, ".txt"))