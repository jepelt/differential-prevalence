library(tidyverse)


# Import datasets---------------------------------------------------------------
studies <- tibble(folder = list.files(path = 'Data/Duvallet')) %>% 
  filter(str_detect(folder, 'results') & !str_detect(folder, 'tar\\.gz')) %>% 
  filter(!str_detect(folder, 'crc_zhu|ob_escobar')) %>% #No metadata
  mutate(folder = str_replace(folder, '_results', '')) %>% 
  pull(folder)


datasets <- list()
for(i in 1:length(studies)){

  # Import count data-----------------------------------------------------------
  counts_file <- paste0('Data/Duvallet/', studies[i], '_results/RDP/',
                        studies[i], '.otu_table.100.denovo.rdp_assigned')
  
  counts_raw <- read.table(file = counts_file, sep = '\t', header = T)
  
  counts <- counts_raw |> 
    column_to_rownames('X') |> 
    rename_all(~ str_replace_all(., '\\.', '-')) |> 
    rename_all(~ case_when(. == "Ademona1-2065" ~ "Adenoma1-2065",
                           . == "ADenoma12-2799" ~ "Adenoma12-2799",
                           . == "HG113" ~ "XHG113",
                           . == "HG237" ~ "XHG237",
                           . == "HG216" ~ "XHG216",
                           TRUE ~ .)) |> 
    rename_all(~ paste0(., '_', studies[i])) |> 
    as.data.frame()
  
  tax <- counts_raw |>
    mutate(id = X) |>
    select(X, id) |> 
    column_to_rownames('X') |>  
    separate(id, into = c('kingdom', 'phylum', 'class', 'order',
                          'family', 'genus', 'species', 'otu'),
             sep = ';')
  
  
  # Import metadata-------------------------------------------------------------
  meta_file <- paste0('Data/Duvallet/', studies[i], '_results/', studies[i],
                      '.metadata.txt')
  
  header <- ifelse(studies[i] %in% c("autism_kb",
                                     "cdi_youngster",
                                     "crc_xiang",
                                     "crc_zhao",
                                     "hiv_noguerajulian",
                                     "ibd_engstrand_maxee",
                                     "ibd_huttenhower",
                                     "ra_littman",
                                     "t1d_mejialeon"), F, T)
  
  meta_raw <- read.table(file = meta_file,
                         sep = '\t',
                         header = header,
                         fileEncoding = "latin1")
  
  meta <- meta_raw |> 
    dplyr::rename(id = 1) %>% 
    mutate(study = studies[i],
           id = str_replace_all(as.character(id), '\\.', '-'),
           id = ifelse(study %in% c("crc_baxter", "hiv_dinh", "ibd_alm",
                                    "ibd_huttenhower", "ra_littman"),
                       paste0('X', id), id),
           id = ifelse(study == "ob_zupancic" & str_sub(id, 1, 1) != 'p',
                       paste0('X', id), id),
           id = paste0(id, '_', studies[i])) |> 
    column_to_rownames('id')
  
  
  datasets[[i]] <- list(counts = counts,
                        tax = tax,
                        meta = meta)
}


#Curate metadata----------------------------------------------------------------
m1 <- datasets[[1]]$meta %>% 
  mutate(dis = 'ASD',
         study = study,
         bmi = NA,
         age = NA,
         sex = factor(host_sex_s, levels = c('female', 'male')),
         group = factor(DiseaseState, levels = c('H', 'ASD'),
                        labels = c('control', 'case'))) %>% 
  select(dis, study, bmi, age, sex, group)

m2 <- datasets[[2]]$meta %>% 
  mutate(dis = 'ASD',
         study = study,
         bmi = NA,
         age = NA,
         sex = NA,
         group = factor(V6, levels = c('H', 'ASD'),
                        labels = c('control', 'case'))) %>% 
  select(dis, study, bmi, age, sex, group)

m3_1 <- datasets[[3]]$meta %>% 
  mutate(dis = 'CDI',
         study = paste0(study, '_cdi'),
         bmi = NA,
         age = age,
         sex = factor(gender, levels = c('F', 'M'),
                      labels = c('female', 'male')),
         group = factor(DiseaseState, levels = c('H', 'CDI'),
                        labels = c('control', 'case'))) %>% 
  filter(!is.na(group)) %>% 
  select(dis, study, bmi, age, sex, group)

m3_2 <- datasets[[3]]$meta %>% 
  mutate(dis = 'diarrhea',
         study = paste0(study, '_diarrhea'),
         bmi = NA,
         age = age,
         sex = factor(gender, levels = c('F', 'M'),
                      labels = c('female', 'male')),
         group = factor(DiseaseState, levels = c('H', 'nonCDI'),
                        labels = c('control', 'case'))) %>% 
  filter(!is.na(group)) %>% 
  select(dis, study, bmi, age, sex, group)

m4 <- datasets[[4]]$meta %>% 
  mutate(dis = 'CDI',
         study = study,
         bmi = NA,
         age = NA,
         sex = NA,
         group = factor(DiseaseState, levels = c('H', 'CDI'),
                        labels = c('control', 'case'))) %>% 
  select(dis, study, bmi, age, sex, group)

m5 <- datasets[[5]]$meta %>% 
  filter(V25 == 0 & V23 %in% c('H', 'CDI')) %>% 
  mutate(dis = 'CDI',
         study = study,
         bmi = NA,
         age = NA,
         sex = NA,
         group = factor(V23, levels = c('H', 'CDI'),
                        labels = c('control', 'case'))) %>% 
  select(dis, study, bmi, age, sex, group)

m6_1 <- datasets[[6]]$meta %>% 
  mutate(dis = 'adenoma',
         study = paste0(study, '_adenoma'),
         bmi = as.numeric(BMI_s) %>% if_else(. < 10, NA, ., NA),
         age = Age_s,
         sex = as.factor(Gender_s),
         group = factor(DiseaseState, levels = c('H', 'nonCRC'),
                        labels = c('control', 'case'))) %>% 
  filter(!is.na(group)) %>% 
  select(dis, study, bmi, age, sex, group)

m6_2 <- datasets[[6]]$meta %>% 
  mutate(dis = 'CRC',
         study = paste0(study, '_crc'),
         bmi = as.numeric(BMI_s) %>% if_else(. < 10, NA, ., NA),
         age = Age_s,
         sex = as.factor(Gender_s),
         group = factor(DiseaseState, levels = c('H', 'CRC'),
                        labels = c('control', 'case'))) %>% 
  filter(!is.na(group)) %>% 
  select(dis, study, bmi, age, sex, group)

m7 <- datasets[[7]]$meta %>% 
  mutate(dis = 'CRC',
         study = study,
         bmi = NA,
         age = NA,
         sex = NA,
         group = factor(V3, levels = c('H', 'CRC'),
                        labels = c('control', 'case'))) %>% 
  select(dis, study, bmi, age, sex, group)

m8_1 <- datasets[[8]]$meta %>% 
  mutate(dis = 'adenoma',
         study = paste0(study, '_adenoma'),
         bmi = NA,
         age = as.numeric(str_sub(age, 1, 2)),
         sex = factor(str_sub(sex, 1, 1), levels = c('f', 'm'),
                      labels = c('female', 'male')),
         group = factor(DiseaseState, levels = c('H', 'nonCRC'),
                        labels = c('control', 'case'))) %>% 
  filter(!is.na(group)) %>% 
  select(dis, study, bmi, age, sex, group)

m8_2 <- datasets[[8]]$meta %>% 
  mutate(dis = 'CRC',
         study = paste0(study, '_CRC'),
         bmi = NA,
         age = as.numeric(str_sub(age, 1, 2)),
         sex = factor(str_sub(sex, 1, 1), levels = c('f', 'm'),
                      labels = c('female', 'male')),
         group = factor(DiseaseState, levels = c('H', 'CRC'),
                        labels = c('control', 'case'))) %>% 
  filter(!is.na(group)) %>% 
  select(dis, study, bmi, age, sex, group)

m9 <- datasets[[9]]$meta %>% 
  mutate(dis = 'CRC',
         study = study,
         bmi = BMI..kg.m..,
         age = Age..years.,
         sex = factor(Gender, levels = c('F', 'M'),
                      labels = c('female', 'male')),
         group = factor(DiseaseState, levels = c('H', 'CRC'),
                        labels = c('control', 'case'))) %>% 
  filter(!is.na(group)) %>% 
  select(dis, study, bmi, age, sex, group)

#Very low library sizes!
m10 <- datasets[[10]]$meta %>% 
  mutate(dis = 'CRC',
         study = study,
         bmi = NA,
         age = NA,
         sex = NA,
         group = factor(V2, levels = c('H', 'CRC'),
                        labels = c('control', 'case'))) %>% 
  filter(!is.na(group)) %>% 
  select(dis, study, bmi, age, sex, group)

m11 <- datasets[[11]]$meta %>% 
  mutate(dis = 'EDD',
         study = study,
         bmi = NA,
         age = NA,
         sex = NA,
         group = factor(DiseaseState, levels = c('H', 'EDD'),
                        labels = c('control', 'case'))) %>% 
  filter(!is.na(group)) %>% 
  select(dis, study, bmi, age, sex, group)

m12 <- datasets[[12]]$meta %>% 
  mutate(dis = 'HIV',
         study = study,
         bmi = as.numeric(BMI_s),
         age = as.numeric(age_s),
         sex = as.factor(tolower(sex_s)),
         group = factor(DiseaseState, levels = c('H', 'HIV'),
                        labels = c('control', 'case'))) %>% 
  filter(!is.na(group)) %>% 
  select(dis, study, bmi, age, sex, group)

#For samples with ls < 100 (after taxa filtering)
m13 <- datasets[[13]]$meta %>% 
  filter(time_point == 1) %>% 
  mutate(dis = 'HIV',
         study = study,
         bmi = body_mass_index,
         age = age,
         sex = NA,
         group = factor(DiseaseState, levels = c('H', 'HIV'),
                        labels = c('control', 'case'))) %>% 
  filter(!is.na(group)) %>% 
  select(dis, study, bmi, age, sex, group)

#A few low readcounts (60, 102, 150, 242...)
m14 <- datasets[[14]]$meta %>% 
  filter(V66 %in% c('BCN0', 'STK')) %>% 
  mutate(dis = 'HIV',
         study = study,
         bmi = NA,
         age = NA,
         sex = as.factor(V34),
         group = factor(V64, levels = c('H', 'HIV'),
                        labels = c('control', 'case'))) %>% 
  filter(!is.na(group)) %>% 
  select(dis, study, bmi, age, sex, group)

m15 <- datasets[[15]]$meta %>% 
  mutate(dis = 'IBD',
         study = study,
         bmi = NA,
         age = age,
         sex = as.factor(gender),
         group = case_when(DiseaseState == 'nonIBD' ~ 'control',
                           DiseaseState %in% c('UC', 'CD') ~ 'case',
                           T ~ NA),
         group = factor(group, levels = c('control', 'case'))) %>% 
  filter(!is.na(group)) %>% 
  select(dis, study, bmi, age, sex, group)

#Smallest ls: 68  90 130 143 145 191
m16 <- datasets[[16]]$meta %>% 
  filter(V4 == 'stool') %>% 
  mutate(dis = 'IBD',
         study = study,
         bmi = NA,
         age = V12,
         sex = factor(V10, levels = c('F', 'M'),
                      labels = c('female', 'male')),
         group = case_when(V6 == 'H' ~ 'control',
                           V6 %in% c('UC', 'CD') ~ 'case',
                           T ~ NA),
         group = factor(group, levels = c('control', 'case'))) %>% 
  filter(!is.na(group)) %>% 
  select(dis, study, bmi, age, sex, group)

m17 <- datasets[[17]]$meta %>% 
  mutate(dis = 'IBD',
         study = study,
         bmi = NA,
         age = NA,
         sex = NA,
         group = case_when(DiseaseState == 'nonIBD' ~ 'control',
                           DiseaseState %in% c('UC', 'CD') ~ 'case',
                           T ~ NA),
         group = factor(group, levels = c('control', 'case'))) %>% 
  filter(!is.na(group)) %>% 
  select(dis, study, bmi, age, sex, group)

m18 <- datasets[[18]]$meta %>% 
  filter(V3 == 'stool') %>% 
  mutate(dis = 'IBD',
         study = study,
         bmi = NA,
         age = NA,
         sex = factor(V15, levels = c('F', 'M'),
                      labels = c('female', 'male')),
         group = case_when(V5 == 'H' ~ 'control',
                           V5 %in% c('UC', 'CD') ~ 'case',
                           T ~ NA),
         group = factor(group, levels = c('control', 'case'))) %>% 
  filter(!is.na(group)) %>% 
  select(dis, study, bmi, age, sex, group)

m19_1 <- datasets[[19]]$meta %>% 
  mutate(dis = 'CIRR',
         study = paste0(study, '_CIRR'),
         bmi = NA,
         age = NA,
         sex = NA,
         group = factor(DiseaseState, levels = c('H', 'CIRR'),
                        labels = c('control', 'case'))) %>% 
  filter(!is.na(group)) %>% 
  select(dis, study, bmi, age, sex, group)

m19_2 <- datasets[[19]]$meta %>% 
  mutate(dis = 'MHE',
         study = paste0(study, '_MHE'),
         bmi = NA,
         age = NA,
         sex = NA,
         group = factor(DiseaseState, levels = c('H', 'MHE'),
                        labels = c('control', 'case'))) %>% 
  filter(!is.na(group)) %>% 
  select(dis, study, bmi, age, sex, group)

m20 <- datasets[[20]]$meta %>% 
  filter(Status_s == 'Baseline') %>% 
  mutate(dis = 'NASH',
         study = study,
         bmi = NA,
         age = NA,
         sex = NA,
         group = factor(DiseaseState, levels = c('H', 'NASH'),
                        labels = c('control', 'case'))) %>% 
  filter(!is.na(group)) %>% 
  select(dis, study, bmi, age, sex, group)

m21_1 <- datasets[[21]]$meta %>% 
  mutate(dis = 'NASH',
         study = paste0(study, '_NASH'),
         bmi = NA,
         age = age,
         sex = factor(sex, levels = c('F', 'M'),
                      labels = c('female', 'male')),
         group = factor(DiseaseState, levels = c('H', 'NASH'),
                        labels = c('control', 'case'))) %>% 
  filter(!is.na(group)) %>% 
  select(dis, study, bmi, age, sex, group)

m21_2 <- datasets[[21]]$meta %>% 
  mutate(dis = 'OB',
         study = paste0(study, '_OB'),
         bmi = NA,
         age = age,
         sex = factor(sex, levels = c('F', 'M'),
                      labels = c('female', 'male')),
         group = factor(DiseaseState, levels = c('H', 'nonNASH-OB'),
                        labels = c('control', 'case'))) %>% 
  filter(!is.na(group)) %>% 
  select(dis, study, bmi, age, sex, group)

m22_1 <- datasets[[22]]$meta %>% 
  filter(n_sample == 0) %>% 
  mutate(dis = 'OB',
         study = paste0(study, '_OB'),
         bmi = body_mass_index,
         age = age,
         sex = NA,
         group = factor(DiseaseState, levels = c('H', 'OB'),
                        labels = c('control', 'case'))) %>% 
  filter(!is.na(group)) %>% 
  select(dis, study, bmi, age, sex, group)

m22_2 <- datasets[[22]]$meta %>% 
  filter(n_sample == 0) %>% 
  mutate(dis = 'OW',
         study = paste0(study, '_OW'),
         bmi = body_mass_index,
         age = age,
         sex = NA,
         group = factor(DiseaseState, levels = c('H', 'OW'),
                        labels = c('control', 'case'))) %>% 
  filter(!is.na(group)) %>% 
  select(dis, study, bmi, age, sex, group)

m23_1 <- datasets[[23]]$meta %>% 
  mutate(dis = 'OB',
         study = paste0(study, '_OB'),
         bmi = NA,
         age = NA,
         sex = NA,
         group = factor(DiseaseState, levels = c('H', 'OB'),
                        labels = c('control', 'case'))) %>% 
  filter(!is.na(group)) %>% 
  select(dis, study, bmi, age, sex, group)

m23_2 <- datasets[[23]]$meta %>% 
  mutate(dis = 'OW',
         study = paste0(study, '_OW'),
         bmi = NA,
         age = NA,
         sex = NA,
         group = factor(DiseaseState, levels = c('H', 'OW'),
                        labels = c('control', 'case'))) %>% 
  filter(!is.na(group)) %>% 
  select(dis, study, bmi, age, sex, group)

m24_1 <- datasets[[24]]$meta %>% 
  mutate(dis = 'OB',
         study = paste0(study, '_OB'),
         bmi = BMI,
         age = age_at_visit,
         sex = factor(sex, levels = c('F', 'M'),
                      labels = c('female', 'male')),
         group = factor(DetailedDiseaseState, levels = c('H', 'OB'),
                        labels = c('control', 'case'))) %>% 
  filter(!is.na(group)) %>% 
  select(dis, study, bmi, age, sex, group)

m24_2 <- datasets[[24]]$meta %>% 
  mutate(dis = 'OW',
         study = paste0(study, '_OW'),
         bmi = BMI,
         age = age_at_visit,
         sex = factor(sex, levels = c('F', 'M'),
                      labels = c('female', 'male')),
         group = factor(DetailedDiseaseState, levels = c('H', 'OW'),
                        labels = c('control', 'case'))) %>% 
  filter(!is.na(group)) %>% 
  select(dis, study, bmi, age, sex, group)

m25 <- datasets[[25]]$meta %>% 
  rownames_to_column('id') %>%
  group_by(submitted_subject_id_s) %>% 
  arrange(visit_number) %>% 
  filter(row_number() == 1) %>% 
  ungroup() %>% 
  mutate(dis = 'OB',
         study = study,
         bmi = NA,
         age = NA,
         sex = as.factor(sex_s),
         group = factor(DiseaseState, levels = c('H', 'OB'),
                        labels = c('control', 'case'))) %>% 
  filter(!is.na(group)) %>% 
  column_to_rownames('id') %>%
  select(dis, study, bmi, age, sex, group)

m26 <- datasets[[26]]$meta %>% 
  mutate(dis = 'parkinson',
         study = study,
         bmi = NA,
         age = NA,
         sex = NA,
         group = factor(DiseaseState, levels = c('H', 'PAR'),
                        labels = c('control', 'case'))) %>% 
  filter(!is.na(group)) %>% 
  select(dis, study, bmi, age, sex, group)

m27_1 <- datasets[[27]]$meta %>% 
  mutate(dis = 'NORA',
         study = paste0(study, '_NORA'),
         bmi = NA,
         age = NA,
         sex = NA,
         group = factor(V6,
                        levels = c('Healthy', 'New onset rheumatoid arthritis'),
                        labels = c('control', 'case'))) %>% 
  filter(!is.na(group)) %>% 
  select(dis, study, bmi, age, sex, group)

m27_2 <- datasets[[27]]$meta %>% 
  mutate(dis = 'CRA',
         study = paste0(study, '_CRA'),
         bmi = NA,
         age = NA,
         sex = NA,
         group = factor(V6,
                        levels = c('Healthy',
                                   'Treated chronic rheumatoid arthritis'),
                        labels = c('control', 'case'))) %>% 
  filter(!is.na(group)) %>% 
  select(dis, study, bmi, age, sex, group)

m27_3 <- datasets[[27]]$meta %>% 
  mutate(dis = 'PSA',
         study = paste0(study, '_PSA'),
         bmi = NA,
         age = NA,
         sex = NA,
         group = factor(V6,
                        levels = c('Healthy',
                                   'Psoriatic arthritis'),
                        labels = c('control', 'case'))) %>% 
  filter(!is.na(group)) %>% 
  select(dis, study, bmi, age, sex, group)

m28 <- datasets[[28]]$meta %>% 
  mutate(dis = 'T1D',
         study = study,
         bmi = NA,
         age = Age,
         sex = factor(Gender, levels = c('F', 'M'),
                      labels = c('female', 'male')),
         group = factor(DiseaseState, levels = c('H', 'T1D'),
                        labels = c('control', 'case'))) %>% 
  filter(!is.na(group)) %>% 
  select(dis, study, bmi, age, sex, group)

m29 <- datasets[[29]]$meta %>% 
  mutate(dis = 'T1D',
         study = study,
         bmi = NA,
         age = NA,
         sex = factor(V19, levels = c('f', 'm'),
                      labels = c('female', 'male')),
         group = factor(V52, levels = c('H', 'T1D'),
                        labels = c('control', 'case'))) %>% 
  filter(!is.na(group)) %>% 
  select(dis, study, bmi, age, sex, group)


# Create TSE objects------------------------------------------------------------
metas <- list(list(m1),
              list(m2),
              list(m3_1, m3_2),
              list(m4),
              list(m5),
              list(m6_1, m6_2),
              list(m7),
              list(m8_1, m8_2),
              list(m9),
              list(m10),
              list(m11),
              list(m12),
              list(m13),
              list(m14),
              list(m15),
              list(m16),
              list(m17),
              list(m18),
              list(m19_1, m19_2),
              list(m20),
              list(m21_1, m21_2),
              list(m22_1, m22_2),
              list(m23_1, m23_2),
              list(m24_1, m24_2),
              list(m25),
              list(m26),
              list(m27_1, m27_2, m27_3),
              list(m28),
              list(m29))

tses_16s <- list()
k <- 0
for(i in 1:length(datasets)){
  
  counts_0 <- datasets[[i]]$counts
  tax <- datasets[[i]]$tax
  
  for(j in 1:length(metas[[i]])){

    meta_0 <- metas[[i]][[j]]
    sample_ids <- intersect(rownames(meta_0), colnames(counts_0))
    
    counts <- counts_0[, sample_ids]
    meta <- meta_0[sample_ids, ]
    meta$body_site <- 'stool'
    meta$seq_type <- '16s'
    
    print(c(ncol(counts), nrow(meta)))
    
    stopifnot(rownames(meta) == colnames(counts))
    stopifnot(rownames(counts) == rownames(tax))
    
    k <- k + 1
    tses_16s[[k]] <- TreeSummarizedExperiment::TreeSummarizedExperiment(
      assays = list(counts = counts),
      colData = S4Vectors::DataFrame(meta),
      rowData = S4Vectors::DataFrame(tax))
  }
}


#Save results and SessionInfo---------------------------------------------------
date <- format(Sys.Date(), "%d%m%y")

save(tses_16s, file = 'data_tses_16s_030625.rds')

si <- capture.output(sessionInfo())
writeLines(si, con = paste0("session_info_data_tses_16s_", date, ".txt"))
