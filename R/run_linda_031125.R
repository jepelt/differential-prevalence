library(tidyverse)
library(TreeSummarizedExperiment)


#LinDA--------------------------------------------------------------------------
run_linda <- function(rel_abs, meta){

  obj <- LinDA::linda(otu.tab = rel_abs |> t(),
                      meta = meta,
                      formula = '~ group',
                      type = 'proportion')
  
  res <- obj$output$groupcase %>%
    rownames_to_column('taxon') %>%
    mutate(method = paste0('LinDA')) %>%
    select(taxon, est = log2FoldChange, p = pvalue, q = padj, method,
           se = lfcSE, df)
  
  return(res %>% mutate(data_id = meta$data_id[1]))
}


#Load data and run DAA----------------------------------------------------------
load('data_090625.rds')

meta_list <- tses |> map(~ colData(.) |> as.data.frame())
assay_list <- tses |> 
  map(~ SummarizedExperiment::assay(., 'relabundance') |> t())
rm(tses)

res_linda <- map2(assay_list, meta_list, ~ run_linda(.x, .y),
                     .progress = T)


#Save results and SessionInfo---------------------------------------------------
date <- format(Sys.Date(), "%d%m%y")

save(res_linda, file = paste0('res_linda_', date, '.rds'))

si <- capture.output(sessionInfo())
writeLines(si, con = paste0("session_info_linda_", date, ".txt"))