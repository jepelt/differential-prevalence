library(tidyverse)
library(TreeSummarizedExperiment)

run_maaslin3 <- function(counts, meta, fm = ~ g + reads){

  meta$g <- meta$group
  fmc <- paste(paste(fm), collapse = ' ')
  
  obj <- maaslin3::maaslin3(input_data = counts,
                            input_metadata = meta,
                            formula = fmc,
                            output = 'output',
                            min_prevalence = 0,
                            plot_associations = F,
                            plot_summary_plot = F,
                            verbosity = "ERROR")
  
  res <- obj$fit_data_prevalence$results %>%
    filter(metadata == 'g' & value == 'case') %>%
    mutate(p = pval_individual,
           q = p.adjust(p, method = 'BH'),
           method = paste0('MaAsLin3-DP')) |> 
    select(taxon = feature, est = coef, p, q, method, error)
  
  return(res %>% mutate(data_id = meta$data_id[1]))
}


#Load data and run DAA----------------------------------------------------------
load('data_090625.rds')

meta_list <- tses |> map(~ colData(.) |> as.data.frame())
assay_list <- tses |> 
  map(~ SummarizedExperiment::assay(., 'counts') |> t())

future::plan(future::multisession, workers = 7)


res_maaslin3 <- furrr::future_map2(
  assay_list, meta_list, ~ run_maaslin3(.x, .y, fm = ~ g + reads),
  .options = furrr::furrr_options(seed = 1),
  .progress = T)


#Save results and SessionInfo---------------------------------------------------
date <- format(Sys.Date(), "%d%m%y")

save(res_maaslin3, file = paste0('res_maaslin3_', date, '.rds'))

si <- capture.output(sessionInfo())
writeLines(si, con = paste0("session_info_maaslin3_", date, ".txt"))
