library(tidyverse)

load('data_090625.rds')

data_list <- 1:length(tses) |> as.list()
for(i in 1:length(tses)){
  
  data_list[[i]]$meta <- SummarizedExperiment::colData(tses[[i]]) |> 
    as.data.frame()
  
  data_list[[i]]$assay <- SummarizedExperiment::assay(tses[[i]], 'pa') |> t()
  
}

base_inds <- seq(0, 1170, 10)
for(i in 1:10){
 
  inds <- base_inds + i
  inds <- inds[inds <= 1175]
  
  data_10 <- data_list[inds]
  
  file_name <- paste0('data_090625_', i, '.rds')
  save(data_10, file = file_name)
}
