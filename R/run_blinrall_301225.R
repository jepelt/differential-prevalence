library(tidyverse)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

c_args <- commandArgs(trailingOnly = TRUE)
data_file <- paste0('data_090625_counts_', c_args[1], '.rds')
save_file <- paste0('res_blinrall_301225_', c_args[1], '.rds')

# Function to impute zeros with half the minimum non-zero value
impute_zeros <- function(x){
  min_nonzero <- min(x[x > 0])
  x[x == 0] <- min_nonzero / 2
  return(x)
}


#Data---------------------------------------------------------------------------
load(data_file)

datasets <- data_10
d_inds <- 1:length(datasets)

mod <- stan_model('stan_blinrall_301225.stan')


res_list <- list()
for(j in 1:length(d_inds)){

  counts <- datasets[[d_inds[j]]]$assay
  meta <- datasets[[d_inds[j]]]$meta
  
  # Transform to relative abundances
  rel_abs <- counts / rowSums(counts)
  stopifnot(all(round(rowSums(rel_abs), 8) == 1))
  
  # Impute zeros in the relative abundance data
  rel_abs_imputed <- apply(rel_abs, 2, impute_zeros)
  
  # Log transform the relative abundances
  log_rel_abs <- log(rel_abs_imputed)
  
  # Standardize the log-transformed data for each taxon
  log_rel_abs_std <- scale(log_rel_abs, center = TRUE, scale = TRUE)
  
  stopifnot(round(colMeans(log_rel_abs_std), 6) == 0)
  stopifnot(round(Rfast::colVars(log_rel_abs_std), 6) == 1)

  
  #Prepare data ----------------------------------------------------------------
  
  taxa <- colnames(log_rel_abs_std)
  n_taxa <- length(taxa)
  
  x <- ifelse(meta$group == 'case', 1, 0)

  #Data object for Stan
  d_stan <- list(N = nrow(counts),   #Number of samples
                 K = n_taxa,         #Number of taxa in P/A analysis
                 y = t(log_rel_abs_std), #Transformed relative abundance matrix
                 x = x               #Group variable
  ) 
  
  
  #Run Stan---------------------------------------------------------------------
  n_chains <- 4
  n_iter <- 3000
  n_warmup <- 1000
  
  st <- Sys.time()
  fit <- sampling(mod,
                  data = d_stan,
                  chains = n_chains,
                  iter = n_iter,
                  warmup = n_warmup,
                  cores = 4,
                  seed = 1)
  et <- Sys.time()
  
  mod_summary <- summary(fit)$summary
  div_trans <- sum(get_divergent_iterations(fit))
  
  
  #Calculate posterior quantities-----------------------------------------------
  res_pa <- extract(fit, pars = 'beta')$beta %>% 
    as.data.frame %>% 
    rename_all(~ taxa) %>% 
    mutate(iter = 1:(n_chains * (n_iter - n_warmup))) %>% 
    pivot_longer(cols = -iter) %>% 
    select(taxon = name, iter, pa = value)
  
  res <- res_pa %>% 
    group_by(taxon) %>% 
    summarize(est_pa = quantile(pa, .500, na.rm = T),
              lwr95_pa = quantile(pa, .025, na.rm = T),
              lwr90_pa = quantile(pa, .050, na.rm = T),
              lwr80_pa = quantile(pa, .100, na.rm = T),
              upr80_pa = quantile(pa, .900, na.rm = T),
              upr90_pa = quantile(pa, .950, na.rm = T),
              upr95_pa = quantile(pa, .975, na.rm = T),
              pl_pa = mean(pa < 0, na.rm = T),
              ph_pa = mean(pa > 0, na.rm = T)) %>% 
    
    mutate(s95 = lwr95_pa > 0 | upr95_pa < 0,
           s90 = lwr90_pa > 0 | upr90_pa < 0,
           s80 = lwr80_pa > 0 | upr80_pa < 0) %>% 
    
    left_join(tibble(taxon = colnames(log_rel_abs_std),
                     N = nrow(log_rel_abs_std),
                     N_zeros = colSums(counts < .5)),
              by = 'taxon') %>%
    mutate(method = 'blinrall',
           specs = paste(n_chains, n_iter, n_warmup),
           N_pos = N - N_zeros,
           prevl = N_pos / N,
           data_id = meta$data_id[1],
           st = st,
           et = et)
  
  res_list[[j]] <- list(res = res,
                        mod_summary = mod_summary,
                        div_trans = div_trans)
  
  print(paste('----- Data number', j, '-----'))
  print(et)
  print(meta$data_id[1])
  print(paste('Divergent transitions', div_trans))
  print(res %>% summarize(s95 = sum(s95),
                          s90 = sum(s90),
                          s80 = sum(s80)))
  
  save(res_list, file = save_file)
}
