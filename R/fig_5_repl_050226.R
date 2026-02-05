library(tidyverse)
library(ggrepel)
library(ggbeeswarm)
library(patchwork)
library(latex2exp)


load('res_blogrrall_120625.rds')
load('res_wald_030226.rds')
load('res_lrt_030226.rds')
load('res_firth_030226.rds')
load('res_maaslin3_030226.rds')
load('res_ldm_260625.rds')
load('res_maaslin2_250126.rds')
load('res_linda_031125.rds')

load('data_expl_studies.rds')
load('meta_090625.rds')
load('data_meta_taxa_090625.rds')


results_dipper <- res_blogrrall |>  
  map(~ .$res) |> 
  bind_rows() |> 
  mutate(q = 2 * pmin(pl_pa, ph_pa),
         dir = sign(est_pa),
         method = 'DiPPER') |> 
  select(method, data_id, taxon, q, dir)

results_wald <- res_wald |> 
  bind_rows() |> 
  mutate(p = ifelse(is_edge, NA, p)) |>
  group_by(data_id) |>
  mutate(q = p.adjust(p, method = 'fdr')) |>
  ungroup() |> 
  mutate(q = ifelse(is.na(q), 1, q),
         dir = sign(est)) |> 
  select(method, data_id, taxon, q, dir)

results_lrt <- res_lrt |>
  bind_rows() |>
  mutate(q = ifelse(is.na(q), 1, q),
         dir = sign(est)) |>
  select(method, data_id, taxon, q, dir)

results_firth <- res_firth |>  
  bind_rows() |> 
  mutate(q = ifelse(is.na(q), 1, q),
         dir = sign(est)) |> 
  select(method, data_id, taxon, q, dir)

results_maaslin3 <- res_maaslin3 |>
  bind_rows() |>
  mutate(q = ifelse(is.na(q), 1, q),
         dir = sign(est)) |> 
  select(method, data_id, taxon, q, dir)

results_ldm <- res_ldm_all |>
  bind_rows() |>
  left_join(results_firth |> select(data_id, taxon, dir),
            by = c('data_id', 'taxon')) |>
  mutate(method = 'LDM-DP') |>
  select(method, data_id, taxon, q, dir)

results_maaslin2 <- res_maaslin2 |> 
  bind_rows() |> 
  mutate(dir = sign(est),
         method = 'MaAsLin2 (DA)') |> 
  select(method, data_id, taxon, q, dir)

results_linda <- res_linda |> 
  bind_rows() |> 
  mutate(dir = sign(est),
         method = 'LinDA (DA)') |> 
  select(method, data_id, taxon, q, dir)


res <- bind_rows(results_dipper,
                 results_wald,
                 results_lrt,
                 results_firth,
                 results_maaslin3,
                 results_ldm,
                 results_maaslin2,
                 results_linda) |>
  dplyr::left_join(metas, by = c('data_id')) |>
  dplyr::left_join(meta_taxa, by = c('data_id', 'taxon')) |>
  mutate(disease = ifelse(disease == 'EDD', 'CDI', disease),
         method = as.factor(method)) |> 
  filter(body_site == 'stool') |> 
  filter(str_detect(data_id, 'original|null')) |>
  select(method, data_id, taxon, q, dir, seq_type,
         disease, study, median_reads, n, n_control, n_case, n_taxa)

rm(list = ls(pattern = 'res_'))
rm(metas, meta_taxa)

colors <- c('DiPPER' = '#CD534CFF',
            'Wald' = '#EFC000FF',
            'LRT' = '#8F7700FF',
            'Firth' = '#6A3D9AFF',
            'MaAsLin3-DP' = '#0073C2FF',
            'LDM-DP' = '#003C67FF',
            'MaAsLin2 (DA)' = 'black',
            'LinDA (DA)' = '#868686FF')



# Definition of replicated and conflicting results------------------------------

dd <- tibble(taxon = paste('Feature', c('A', 'B')),
             est_1 = c(1.0,  0.9),
             est_2 = c(1.5, -1.0)) |> 
  pivot_longer(cols = c(est_1, est_2), values_to = 'est') |> 
  mutate(Study = factor(name, levels = c('est_2', 'est_1'),
                        labels = c('Study 2', 'Study 1')),
         Feature = fct_rev(taxon),
         lwr = est - ifelse(Study == 'Study 1', 0.60, 0.60),
         upr = est + ifelse(Study == 'Study 1', 0.60, 0.60))

dt <- tibble(Feature = paste('Feature', c('A', 'B')),
             est = c(2.4, 2.4),
             lbl = c('Replicated', 'Conflicting'))


# Number of replicated at alpha = .10 ------------------------------------------
expl_studies_2 <- expl_studies |> 
  str_replace('whole', 'original') |> 
  str_replace('16S', '16s') |> 
  str_replace_all(' ', '_')

resl_sep_1 <- list()

for(i in 1:1){
  sl <- .10
  
  res_sep1 <- res |> 
    filter(str_detect(data_id, 'original')) |> 
    filter(data_id %in% expl_studies_2 & !is.na(q)) |> 
    mutate(n1 = n + median_reads / 10 ^ 10,
           sgn1 = q < sl) |> 
    select(method, seq_type, disease, taxon, study1 = study, n1, dir1 = dir, sgn1)
  
  res_sep2 <- res |> 
    filter(str_detect(data_id, 'original')) |> 
    filter(!is.na(q) & n_case >= 10 & n_control >= 10) |> 
    mutate(n2 = n + median_reads / 10 ^ 10,
           sgn2 = q < sl) |> 
    select(method, seq_type, disease, taxon, study2 = study, n2, dir2 = dir, sgn2)
  
  res_sepc <- dplyr::inner_join(res_sep1, res_sep2,
                                by = join_by(method, seq_type, disease, taxon, n1 < n2)) |> 
    mutate(match = sgn1 & sgn2 & dir1 == dir2,
           confl = sgn1 & sgn2 & dir1 == -dir2,
           either = sgn1 | sgn2)
  
  resl_sep_1[[i]] <- res_sepc |> 
    group_by(method, study1, study2) |> 
    summarize(m = sum(match),
              c = sum(confl),
              e = sum(either),
              sl = sl)
  
  print(i)
}


q75 <- function(x) quantile(x, .75, na.rm = T)

bind_rows(resl_sep_1) |> 
  ungroup() |> 
  group_by(study1, study2) |>
  summarize(max_m = max(m)) |> 
  ungroup() |> 
  summarize(m = mean(max_m > 0),
            s = sum(max_m > 0))

res_sep_1 <- bind_rows(resl_sep_1) |> 
  ungroup() |> 
  mutate(method = fct_reorder(method, m, .fun = q75))


# Replicated & conflicting "RO curves" -----------------------------------------
sls <- seq(.00, 1, .01)

resl_sep_2 <- list()

for(i in 1:length(sls)){
  sl <- sls[i]
  
  res_sep1 <- res |> 
    filter(str_detect(data_id, 'original')) |> 
    filter(data_id %in% expl_studies_2 & !is.na(q)) |> 
    mutate(n1 = n + median_reads / 10 ^ 10,
           sgn1 = q < sl) |> 
    select(method, seq_type, disease, taxon, study1 = study, n1, dir1 = dir, sgn1)
  
  res_sep2 <- res |> 
    filter(str_detect(data_id, 'original')) |> 
    filter(!is.na(q) & n_case >= 10 & n_control >= 10) |> 
    mutate(n2 = n + median_reads / 10 ^ 10,
           sgn2 = q < sl) |> 
    select(method, seq_type, disease, taxon, study2 = study, n2, dir2 = dir, sgn2)
  
  res_sepc <- dplyr::inner_join(res_sep1, res_sep2,
                                by = join_by(method, seq_type, disease, taxon, n1 < n2)) |> 
    mutate(match = sgn1 & sgn2 & dir1 == dir2,
           confl = sgn1 & sgn2 & dir1 == -dir2)
  
  resl_sep_2[[i]] <- res_sepc |> 
    group_by(method) |> 
    summarize(m = sum(match),
              c = sum(confl),
              sl = sl)
  
  print(i)
}

n_data_pairs <- res_sepc |> 
  distinct(study1, study2) |> 
  nrow()
stopifnot(n_data_pairs == 110)


res_sep_2 <- bind_rows(resl_sep_2) |> 
  mutate(analysis = 'sep',
         sl = round(sl, 2),
         sl_lab = case_when(sl == .05 ~ '0.05',
                            sl == .10 ~ '0.10',
                            sl == .20 ~ '0.20',
                            T ~ NA),
         method = factor(
           method,
           levels = c('Wald', 'MaAsLin3-DP', 'LinDA (DA)', 'LDM-DP',
                      'Firth', 'MaAsLin2 (DA)', 'LRT', 'DiPPER')))


# Subfigures--------------------------------------------------------------------

or_breaks <- c(0.2, 0.5, 1, 2, 5)

p1 <- ggplot(dd, aes(est, Feature, color = Study)) +
  geom_point(size = 2.0,
             position = position_dodge(width = .50)) +
  geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0,
                 position = position_dodge(width = 0.50),
                 linewidth = 0.8) +
  geom_text(data = dt, mapping = aes(label = lbl), color = 'black',
            size = 2.8, hjust = 0) +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  scale_color_manual(values = c('orange', 'grey40'),
                     guide = guide_legend(reverse = T, nrow = 2)) +
  scale_x_continuous(
    breaks = log(or_breaks),
    labels = or_breaks
  ) +
  labs(x = 'Differential prevalence (OR)',
       y = 'Feature', 
       color = 'Study') +
  coord_cartesian(xlim = c(-1.9, 3.9)) +
  
  theme_classic(base_size = 10) +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size = 8, hjust = 0.0),
        legend.title = element_blank(),
        legend.position = 'inside',
        legend.position.inside = c(0.80, 1.25),
        legend.spacing.y = unit(3, 'mm'),
        legend.key.width = unit(6.0, 'mm'),
        legend.key.height = unit(2.0, 'mm'),
        legend.box.background = element_rect(colour = "black"),
        legend.text = element_text(size = 8),
        plot.margin = margin(t = 30, r = 15, b = 5.5, l = 2.5))


log10p <- function(x) log10(x + 1)

(p2 <- ggplot(res_sep_1, aes(x = method, y = log10(m + 1), color = method)) +
    ggbeeswarm::geom_quasirandom(size = 1.5, alpha = .3, width = 0.4,
                                 bandwidth = 0.9) +
    geom_boxplot(alpha = .2) +
    scale_y_continuous(breaks = log10p(c(0, 1, 3, 10, 30)),
                       labels = c('0', '1', '3', '10', '30')) +
    # geom_point(data = resm_sep, aes(y = m), size = 5, color = 'black') +
    scale_color_manual(values = colors) +
    labs(y = "Number of replicated findings\non each pair of studies (α = 0.10)") +
    theme_classic(base_size = 10) +
    theme(legend.position = 'none',
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
          plot.margin = margin(t = 15, r = 15, b = 10, l = 5.5))
)


(p3 <- ggplot(res_sep_2, aes(x = c, y = m, color = method, shape = sl_lab)) +
    geom_point(size = 4, alpha = .5) +
    geom_line(aes(group = method), linewidth = 1, alpha = .5) +
    coord_cartesian(xlim = c(0, 120), ylim = c(0, 1050)) +
    scale_color_manual(values = colors,
                       guide = guide_legend(reverse = T)) +
    scale_shape_manual(values = c(16, 17, 15, 3),
                       na.translate = F,
                       guide = guide_legend(reverse = T)) +
    labs(y = 'Total number of replicated findings',
         x = 'Total number of conflicting findings',
         color = '',
         shape = 'Significance\nlevel (α)') +
    theme_classic(base_size = 10) +
    theme(legend.position = 'right',
          plot.margin = margin(t = 15, r = 5.5, b = 25, l = 10))
)


# Combine subfigures------------------------------------------------------------
p4 <- ggpubr::ggarrange(p1, p2, common.legend = F, heights = c(.26, .74),
                        ncol = 1, nrow = 2, labels = c('a', 'b'))

(p <- ggpubr::ggarrange(p4, p3, ncol = 2, nrow = 1, widths = c(.4, .6),
                        labels = c('', 'c')))

date <- format(Sys.Date(), "%d%m%y")
ggsave(plot = p, filename = paste0('fig_5_repl_', date, '.png'),
       width = 170, height = 140, dpi = 300, unit = 'mm',
       bg = 'white')


#-------------------------------------------------------------------------------
(res_fig_3 <- res_sep_1 |> 
  group_by(method, sl) |> 
  summarize(mean = mean(m),
            median = median(m),
            q75 = quantile(m, .75))
)

writexl::write_xlsx(res_fig_3, paste0('res_fig_5_', date, '.xlsx'))
