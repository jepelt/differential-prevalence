library(tidyverse)
library(TreeSummarizedExperiment)


# Wald logistic regression------------------------------------------------------
run_wald <- function(assay, meta, fm = ~ group + log_reads){

  meta$log_reads <- log10(meta$reads) |> scale(scale = F) |> as.numeric()
  
  # Find edge cases
  pa_case <- assay[meta$group == 'case', ]
  pa_control <- assay[meta$group == 'control', ]
  is_edge <- colMeans(pa_case) == 0 | colMeans(pa_control) == 0
  
  X <- model.matrix(fm, meta)
  
  res <- apply(assay, 2, dpa_wald, x = X) %>% 
    bind_rows(.id = 'taxon') %>%
    mutate(prevl_1 = colMeans(assay) == 1,
           is_edge = is_edge,
           p = ifelse(prevl_1, NA, p_wald),
           q = p.adjust(p, method = 'BH'),
           method = 'Wald') %>% 
    select(taxon, est, se, p, q, method, w, is_edge)
  
  return(res %>% mutate(data_id = meta$data_id[1]))
}


dpa_wald <- function(y, x){
  
  warn_msg <- 'None'
  m1 <- withCallingHandlers(
    expr = fastglm::fastglm(x, y, family = binomial()),
    warning = function(w) {
      warn_msg <<- w$message
      invokeRestart("muffleWarning")
    }
  ) 
  
  cf <- coef(summary(m1))
  est <- as.numeric(cf['groupcase', 1])
  se <- as.numeric(cf['groupcase', 2])
  p_wald <- as.numeric(cf['groupcase', 4])
  
  list(est = est, se = se, p_wald = p_wald, w = warn_msg)
}


# Load data and run the analysis -----------------------------------------------
load('data_090625.rds')

meta_list <- tses |> map(~ colData(.) |> as.data.frame())
assay_list <- tses |> 
  map(~ SummarizedExperiment::assay(., 'pa') |> t())

future::plan(future::multisession, workers = 7)


res_wald <- furrr::future_map2(assay_list,
                                meta_list,
                                ~ run_wald(.x, .y),
                                   .options = furrr::furrr_options(seed = 1),
                                   .progress = T)



#Save results and SessionInfo---------------------------------------------------
date <- format(Sys.Date(), "%d%m%y")

save(res_wald, file = paste0('res_wald_', date, '.rds'))

si <- capture.output(sessionInfo())
writeLines(si, con = paste0("session_info_wald_", date, ".txt"))
