library(tidyverse)
library(TreeSummarizedExperiment)

load('data_090625.rds')

mt_list <- list()

for(i in 1:length(tses)){
  
  counts <- assay(tses[[i]], 'counts')
  rel_abs <- assay(tses[[i]], 'relabundance')
  meta <- colData(tses[[i]]) |> as.data.frame()
  
  mt <- tibble(taxon = rownames(counts),
               data_id = meta$data_id[1],
               prevl = rowMeans(counts > 0),
               prevl_control = rowMeans(counts[, meta$group == 'control'] > 0),
               prevl_case = rowMeans(counts[, meta$group == 'case'] > 0),
               edge = prevl_control == 0 | prevl_case == 0,
               mean_ra = rowMeans(rel_abs),
               md_ra = Rfast::rowMedians(rel_abs),
               min_ra = rowMins(rel_abs),
               max_ra = rowMaxs(rel_abs))
  
  mt_list[[i]] <- mt
}

meta_taxa <- bind_rows(mt_list)

save(meta_taxa, file = 'data_meta_taxa_090625.rds')
