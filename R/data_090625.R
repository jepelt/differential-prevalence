library(tidyverse)
library(TreeSummarizedExperiment)

load('data_tses_16s_030625.rds')
load('data_tses_sg_030625.rds')

tses_both <- c(tses_16s, tses_sg)


# Tidy colData------------------------------------------------------------------
tses_raw <- list()

for(i in 1:length(tses_both)){
  
  tse <- tses_both[[i]]
  
  colData(tse) <- colData(tse) |> 
    as.data.frame() |> 
    mutate(reads = ifelse(seq_type == '16s',
                          colSums(assay(tse)),
                          number_reads)) |>
    select(seq_type, study, body_site, disease = dis, reads,
           group, bmi, age, sex) |> 
    DataFrame()
  
  if(tse$seq_type[1] == '16s'){
    tse <- tse[, tse$reads >= 200]
  }else{
    tse <- tse[, tse$reads >= 100000]
  }
    
  print(min(tse$reads))
  
  tses_raw[[i]] <- tse 
}

rm(tses_16s, tses_sg, tses_both)

# 85 Original datasets----------------------------------------------------------
tses_original <- list()
meta_original_list <- list()

for(i in 1:length(tses_raw)){
  
  tse <- tses_raw[[i]]
  
  tse$data_id <- paste('original', tse$seq_type, tse$study, sep = '_')
  
  # Taxa with at least 4 non-zero counts are included
  n <- ncol(tse)
  prevl <- 4 / n - 10 ^ -6
  
  stopifnot(tse$seq_type[1] %in% c('16s', 'sg'))
  
  if(tse$seq_type[1] == '16s'){
    tse <- mia::subsetByPrevalent(tse, rank = "genus", prevalence = prevl)
  }else{
    tse <- mia::subsetByPrevalent(tse, rank = "species", prevalence = prevl)
  }
  
  tse <- tse |> 
    mia::transformAssay(assay.type = "counts", method = "relabundance")
  
  tse <- tse |> 
    mia::transformAssay(assay.type = "counts", method = "pa")
  
  tses_original[[i]] <- tse
  
  meta <- tibble(
    data_id = tse$data_id[1],
    seq_type = tse$seq_type[1],
    study = tse$study[1],
    body_site = tse$body_site[1],
    disease = tse$disease[1],
    n = ncol(tse),
    n_case = sum(tse$group == 'case'),
    n_control = sum(tse$group == 'control'),
    n_taxa = nrow(tse),
    min_reads = min(tse$reads, na.rm = TRUE),
    median_reads = median(tse$reads, na.rm = TRUE),
    max_reads = max(tse$reads, na.rm = TRUE)
  )
  
  meta_original_list[[i]] <- meta
}


# Null datasets-----------------------------------------------------------------

# Distinct studies with at least 20 controls
distinct_studies <- 
  tses_raw |> 
  keep(~ sum(.$group == 'control') >= 20) |> 
  map_dfr(~ colData(.) |> as.data.frame() |> dplyr::slice(1)) |> 
  remove_rownames() |> 
  mutate(study0 = str_replace(tolower(study),
                              paste0('_', tolower(disease)),
                              '')) |> 
  distinct(study0, .keep_all = T) |> 
  filter(!study %in% c('ChngKR_2016_1', 'ChngKR_2016_2')) |> 
  pull(study)

tses_controls <- 
  tses_raw |>
  keep(~ .$study[1] %in% distinct_studies) |> 
  map(~ .[, .$group == 'control'])


tses_null <- list()
meta_null_list <- list()
k <- 0

set.seed(1)
for(i in 1:length(tses_controls)){
  
  tse <- tses_controls[[i]]
  
  for(j in 1:10){

    equal_sizes <- sample(c(F, T), 1)
    n <- ncol(tse)
    
    if(equal_sizes){
      n_controls <- ceiling(n / 2)
    }else{
      n_controls <- sample(10:(n - 10), size = 1)
    }
    
    n_cases <- n - n_controls
    
    tse_null <- tse
    
    tse_null$group <- 
      sample(c(rep('control', n_controls), rep('case', n_cases))) |> 
      factor(levels = c('control', 'case'))
    
    tse_null$data_id <- 
      paste('null', tse_null$seq_type[1], tse_null$study[1], j, sep = '_')
    
    # Taxa with at least 4 non-zero counts are included
    n <- ncol(tse_null)
    prevl <- 4 / n - 10 ^ -6
    
    if(tse_null$seq_type[1] == '16s'){
      tse_null <-
        mia::subsetByPrevalent(tse_null, rank = "genus", prevalence = prevl)
    }else{
      tse_null <-
        mia::subsetByPrevalent(tse_null, rank = "species", prevalence = prevl)
    }
    
    tse_null <- tse_null |> 
      mia::transformAssay(assay.type = "counts", method = "relabundance")
    
    tse_null <- tse_null |> 
      mia::transformAssay(assay.type = "counts", method = "pa")
    
    k <- k + 1
    tses_null[[k]] <- tse_null
    
    meta <- tibble(
      data_id = tse_null$data_id[1],
      seq_type = tse_null$seq_type[1],
      study = tse_null$study[1],
      body_site = tse_null$body_site[1],
      disease = tse_null$disease[1],
      n = ncol(tse_null),
      n_case = sum(tse_null$group == 'case'),
      n_control = sum(tse_null$group == 'control'),
      n_taxa = nrow(tse_null),
      min_reads = min(tse_null$reads, na.rm = TRUE),
      median_reads = median(tse_null$reads, na.rm = TRUE),
      max_reads = max(tse_null$reads, na.rm = TRUE)
    )
    
    meta_null_list[[k]] <- meta
  }
}

rm(tses_controls)


# Split datasets----------------------------------------------------------------

# Studies with at least 20 cases and controls
tses_splittable <- 
  tses_raw |> 
  keep(~ sum(.$group == 'control') >= 20 & 
         sum(.$group == 'case') >= 20)

tses_split <- list()
meta_split_list <- list()
k <- 0

set.seed(1)
for(i in 1:length(tses_splittable)){

  tse <- tses_splittable[[i]]
  
  for(j in 1:5){

    split_inds <- tibble(group = tse$group) |> 
      group_by(group) |> 
      mutate(split_ind = sample(rep(1:2, 5000)[1:n()])) |>
      ungroup()
    
    tse$split_ind <- split_inds$split_ind
    
    tse_splits <- mia::splitOn(tse, by = 'samples', group = 'split_ind')
    
    for(split_ind in 1:2){
      
      tse_split <- tse_splits@listData[[split_ind]]
      
      tse_split$data_id <- 
        paste('split', tse$seq_type[1], tse$study[1], j, split_ind, sep = '_')
      
      # Taxa with at least 4 non-zero counts are included
      n <- ncol(tse_split)
      prevl <- 4 / n - 10 ^ -6
      
      if(tse_split$seq_type[1] == '16s'){
        tse_split <-
          mia::subsetByPrevalent(tse_split, rank = "genus", prevalence = prevl)
      }else{
        tse_split <-
          mia::subsetByPrevalent(tse_split, rank = "species", prevalence = prevl)
      }
      
      tse_split <- tse_split |> 
        mia::transformAssay(assay.type = "counts", method = "relabundance")
      
      tse_split <- tse_split |> 
        mia::transformAssay(assay.type = "counts", method = "pa")
      
      k <- k + 1
      tses_split[[k]] <- tse_split
      
      meta <- tibble(
        data_id = tse_split$data_id[1],
        seq_type = tse_split$seq_type[1],
        study = tse_split$study[1],
        body_site = tse_split$body_site[1],
        disease = tse_split$disease[1],
        n = ncol(tse_split),
        n_case = sum(tse_split$group == 'case'),
        n_control = sum(tse_split$group == 'control'),
        n_taxa = nrow(tse_split),
        min_reads = min(tse_split$reads, na.rm = TRUE),
        median_reads = median(tse_split$reads, na.rm = TRUE),
        max_reads = max(tse_split$reads, na.rm = TRUE)
      )
      
      meta_split_list[[k]] <- meta
    }
  }
}

rm(tses_splittable)


# Combine and save results. Print SessionInfo-----------------------------------
date <- format(Sys.Date(), "%d%m%y")

tses <- c(tses_original, tses_null, tses_split)
metas <- bind_rows(meta_original_list, meta_null_list, meta_split_list)

save(tses, file = 'data_090625.rds')
save(metas, file = 'meta_090625.rds')

si <- capture.output(sessionInfo())
writeLines(si, con = paste0("session_info_data_", date, ".txt"))
