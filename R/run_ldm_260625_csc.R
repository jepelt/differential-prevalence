library(tidyverse)
library(TreeSummarizedExperiment)

.libPaths(c("/projappl/project_2014622/project_rpackages_r442", .libPaths()))  
library(LDM)


c_args <- commandArgs(trailingOnly = TRUE)
data_file <- paste0('data_090625_counts_', c_args[1], '.rds')
save_file <- paste0('res_ldm_260625_', c_args[1], '.rds')


# Function to run LDM-----------------------------------------------------------
run_ldm <- function(counts, meta, fm = ~ group, n_rarefy = 5){

  covs <- paste(all.vars(fm)[-1], collapse = ' + ')
  
  if(covs == ""){
    obj <- LDM::ldm(counts ~ group,
                    data = meta,
                    comp.anal = F,
                    binary = T,
                    n.rarefy = n_rarefy,
                    verbose = F,
                    n.cores = 1)
  }else{
    f <- paste0('counts | (', covs, ') ~ group') %>% as.formula
    obj <- LDM::ldm(f,
                    data = meta,
                    comp.anal = F,
                    binary = T,
                    n.rarefy = n_rarefy,
                    verbose = F,
                    n.cores = 1)
  }
  
  stopifnot(rownames(obj$beta) %in% c('groupcase', 'groupcontrol'))
  stopifnot(names(obj$q.otu.pa) == names(obj$beta))
  stopifnot(names(obj$q.otu.pa) == names(obj$p.otu.omni))
  
  res <- tibble(taxon = names(obj$q.otu.pa),
                name = rownames(obj$beta),
                dir = ifelse(name == 'groupcontrol', -1, 1),
                est_abs = obj$beta[1, ] %>% as.numeric,
                est = dir * est_abs,
                p = obj$p.otu.pa %>% as.numeric,
                q = obj$q.otu.pa %>% as.numeric,
                p_global = obj$p.global.pa[1] %>% as.numeric,
                method = paste0('LDM (', n_rarefy, ')')) %>%
    select(taxon, est, p, q, method)
  
  return(res %>% mutate(data_id = meta$data_id[1]))
}


#Load data and run DAA----------------------------------------------------------
load(data_file)

assay_list <- data_10 |> map(~ .$assay)
meta_list <- data_10 |> map(~ .$meta)

set.seed(1)
res_ldm <- list()
for(i in 1:length(assay_list)){
  counts <- assay_list[[i]]
  meta <- meta_list[[i]]
  res_ldm[[i]] <- run_ldm(counts, meta, n_rarefy = 'all')
  print(paste(c_args[1], i))
}

save(res_ldm, file = save_file)
