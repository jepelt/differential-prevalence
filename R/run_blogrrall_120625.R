library(tidyverse)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

c_args <- commandArgs(trailingOnly = TRUE)
data_file <- paste0('data_090625_', c_args[1], '.rds')
save_file <- paste0('res_blogrrall_120625_', c_args[1], '.rds')


#Data---------------------------------------------------------------------------
load(data_file)

datasets <- data_10
d_inds <- 1:length(datasets)

mod <- stan_model('stan_blogrrall_120625.stan')


res_list <- list()
for(j in 1:length(d_inds)){
  
  pa <- datasets[[d_inds[j]]]$assay
  meta <- datasets[[d_inds[j]]]$meta
  
  
  #Prepare data ----------------------------------------------------------------
  
  taxa_pa <- colnames(pa)
  n_taxa_pa <- length(taxa_pa)
  
  x <- ifelse(meta$group == 'case', 1, 0)
  reads <- log10(meta$reads) %>% scale(scale = F) %>% as.numeric()
  
  #Data object for Stan
  d_stan <- list(N = nrow(pa),   #Number of samples
                 K = n_taxa_pa,  #Number of taxa in P/A analysis
                 y = t(pa),      #Presence/absence matrix
                 x = x,          #Group variable
                 reads = reads
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
    rename_all(~ taxa_pa) %>% 
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
    
    left_join(tibble(taxon = colnames(pa),
                     N = nrow(pa),
                     N_zeros = colSums(pa < .5)),
              by = 'taxon') %>%
    mutate(method = 'blogrrall',
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
