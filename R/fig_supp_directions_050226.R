library(tidyverse)
library(patchwork)

load('res_firth_030226.rds')
load('data_expl_studies.rds')
load('meta_090625.rds')
load('data_meta_taxa_090625.rds')


results_firth <- res_firth |>  
  bind_rows() |> 
  filter(str_detect(data_id, 'original')) |>
  dplyr::left_join(metas, by = c('data_id')) |>
  dplyr::left_join(meta_taxa, by = c('data_id', 'taxon')) |>
  mutate(disease = ifelse(disease %in% c('EDD', 'CDI', 'diarrhea'),
                          'Diarrhea',
                          disease),
         dir = sign(est)) |> 
  filter(body_site == 'stool') |> 
  select(method, data_id, taxon, q, dir, seq_type,
         disease, study, median_reads, n)


sgn_level <- 0.10


res_dir <- results_firth |> 
  group_by(seq_type, disease, data_id) |> 
  filter(!is.na(q)) |>
  summarize(dir_pos = sum(dir > 0 & q < sgn_level),
            dir_neg = sum(dir < 0 & q < sgn_level),
            n = first(n),
            .groups = 'drop') |> 
  ungroup() |> 
  mutate(total = dir_pos + dir_neg,
         prop_pos = dir_pos / total) |>
  filter(total >= 5) |>
  mutate(disease_clean = ifelse(disease == toupper(disease), 
                                disease, 
                                stringr::str_to_sentence(disease)),
         seq_type = toupper(seq_type),
         label_text = paste0(disease_clean, " (", seq_type, ", N = ", n, ")")) |> 
  pivot_longer(cols = c(dir_pos, dir_neg),
               names_to = 'direction',
               values_to = 'count') |> 
  mutate(
    row_id = fct_reorder(data_id, prop_pos),
    direction = recode(
      direction,
      'dir_pos' = 'Positive associations\n(Features significanlty more\nprevalent in cases)',
      'dir_neg' = 'Negative associations\n(Features significantly more\nprevalent in controls)'),
    prop = count / total
  )


y_labels <- setNames(res_dir$label_text, res_dir$data_id)

(p <- ggplot(res_dir, aes(x = direction, y = row_id)) +
  geom_tile(aes(fill = prop), color = "gray80") +
  geom_text(aes(label = count, 
                color = prop > 0.5), 
            size = 3) +
  scale_fill_gradient(low = "white", high = "black", limits = c(0, 1)) +
  scale_color_manual(values = c("black", "white"), guide = "none") +
  scale_y_discrete(labels = y_labels) + 
  labs(x = 'Direction of significant differential prevalence results',
       y = "Dataset",
       fill = 'Proportion') +
  theme_minimal(base_size = 10) +
  theme(axis.title.y = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 10),
        panel.grid = element_blank())
)


date <- format(Sys.Date(), "%d%m%y")
ggsave(plot = p, filename = paste0('fig_supp_directions_', date, '.png'),
       width = 170, height = 180, dpi = 300, unit = 'mm',
       bg = 'white')

