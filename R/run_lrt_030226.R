library(tidyverse)
library(TreeSummarizedExperiment)


# Likelihood ratio test logistic regression ------------------------------------
run_lrt <- function(assay, meta, fm = ~ group + log_reads){

  meta$log_reads <- log10(meta$reads) |> scale(scale = F) |> as.numeric()
  
  X <- model.matrix(fm, meta)
  
  res <- apply(assay, 2, dpa_lrt, x = X) %>% 
    bind_rows(.id = 'taxon') %>%
    mutate(prevl_1 = colMeans(assay) == 1,
           l = 2 * (l1 - l0),
           p0 = pchisq(l, df = 1, lower.tail = F),
           p = ifelse(prevl_1, NA, p0),
           q = p.adjust(p, method = 'BH'),
           method = 'LRT') %>% 
    select(taxon, est, p, q, method, w0, w1)
  
  return(res %>% mutate(data_id = meta$data_id[1]))
}


dpa_lrt <- function(y, x){
  
  x0 <- if(ncol(x) > 2){x[, -2]}else{matrix(rep(1, length(y)), ncol = 1)}
  
  warn_msg_0 <- warn_msg_1 <- 'None'
  
  m0 <- withCallingHandlers(
    expr = fastglm::fastglm(x0, y, family = binomial()),
    warning = function(w) {
      warn_msg_0 <<- w$message
      invokeRestart("muffleWarning")
    }
  )
  
  m1 <- withCallingHandlers(
    expr = fastglm::fastglm(x, y, family = binomial()),
    warning = function(w) {
      warn_msg_1 <<- w$message
      invokeRestart("muffleWarning")
    }
  ) 
  
  l0 <- as.numeric(logLik(m0))
  l1 <- as.numeric(logLik(m1))
  
  cf <- coef(summary(m1))
  est <- as.numeric(cf['groupcase', 1])

  list(l0 = l0, l1 = l1, est = est, w0 = warn_msg_0, w1 = warn_msg_1)
}


# Load data and run the analysis -----------------------------------------------
load('data_090625.rds')

meta_list <- tses |> map(~ colData(.) |> as.data.frame())
assay_list <- tses |> 
  map(~ SummarizedExperiment::assay(., 'pa') |> t())

future::plan(future::multisession, workers = 7)


res_lrt <- furrr::future_map2(assay_list,
                                meta_list,
                                ~ run_lrt(.x, .y),
                                   .options = furrr::furrr_options(seed = 1),
                                   .progress = T)



#Save results and SessionInfo---------------------------------------------------
date <- format(Sys.Date(), "%d%m%y")

save(res_lrt, file = paste0('res_lrt_', date, '.rds'))

si <- capture.output(sessionInfo())
writeLines(si, con = paste0("session_info_lrt_", date, ".txt"))
