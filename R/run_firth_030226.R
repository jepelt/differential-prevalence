library(tidyverse)
library(TreeSummarizedExperiment)


# Firth logistic regression------------------------------------------------------------
run_firth <- function(assay, meta, fm = ~ group + log_reads){

  meta$log_reads <- log10(meta$reads) |> scale(scale = F) |> as.numeric()

  X <- model.matrix(fm, meta)
  
  res <- apply(assay, 2, dpa_firth, x = X) %>% 
      t() %>%
      as.data.frame() %>%
      dplyr::rename(est = 1, p0 = 2, lwr = 3, upr = 4) %>%
      mutate(prevl_1 = colMeans(assay) == 1,
             p = ifelse(prevl_1, NA, p0),
             q = p.adjust(p, method = 'BH'),
             method = 'Firth') %>% 
      rownames_to_column('taxon') %>%
      select(taxon, est, p, q, method, lwr, upr)
  
  return(res %>% mutate(data_id = meta$data_id[1]))
}


dpa_firth <- function(y, x){

  data <- cbind(y = y, x) %>% as.data.frame() %>% select(- '(Intercept)')
  m <- logistf::logistf(y ~ ., data = data,
                        alpha = .10,
                        control= logistf::logistf.control(maxit = 1000))
  
  c(est = m$coefficients['groupcase'],
    p = m$prob['groupcase'],
    lwr = m$ci.lower['groupcase'],
    upr = m$ci.upper['groupcase'])
}


# Load data and run the analysis -----------------------------------------------
load('data_090625.rds')

meta_list <- tses |> map(~ colData(.) |> as.data.frame())
assay_list <- tses |> 
  map(~ SummarizedExperiment::assay(., 'pa') |> t())

future::plan(future::multisession, workers = 7)


res_firth <- furrr::future_map2(assay_list,
                                meta_list,
                                ~ run_firth(.x, .y),
                                   .options = furrr::furrr_options(seed = 1),
                                   .progress = T)



#Save results and SessionInfo---------------------------------------------------
date <- format(Sys.Date(), "%d%m%y")

save(res_firth, file = paste0('res_firth_', date, '.rds'))

si <- capture.output(sessionInfo())
writeLines(si, con = paste0("session_info_firth_", date, ".txt"))
