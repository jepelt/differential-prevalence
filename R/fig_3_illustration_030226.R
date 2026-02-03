library(tidyverse)
library(patchwork)

load('res_blogrrall_120625.rds')
load('res_wald_030226.rds')
load('meta_090625.rds')
load('data_meta_taxa_090625.rds')

orig_studies <- metas |> 
  filter(str_detect(data_id, 'original')) |> 
  distinct(data_id, .keep_all = TRUE)

null_studies <- metas |> 
  filter(str_detect(data_id, 'null')) |> 
  distinct(data_id, .keep_all = TRUE) |> 
  filter(n_case == 30 & n_control == 31)


#-------------------------------------------------------------------------------
sl <- .10

null_data <- 'null_sg_ZellerG_2014_CRC_3'
original_data <- 'original_sg_GuptaA_2019'

res_freq <- res_wald |> 
  bind_rows() |> 
  filter(data_id %in% c(null_data, original_data)) |> 
  mutate(method = 'Frequentist (Wald)\n(unadjusted)',
         lwr = est - qnorm(.95) * se,
         upr = est + qnorm(.95) * se) |> 
  select(method, data_id, taxon, est, se, p, lwr, upr)

res_bonf <- res_freq |> 
  group_by(data_id) |> 
  mutate(n_tests = sum(!is.na(p))) |> 
  ungroup() |>
  mutate(method = 'Frequentist (Wald)\n(multiplicity adjusted)',
         sl_bonf = sl / n_tests,
         qnorm_bonf = qnorm(1 - sl_bonf / 2),
         lwr = est - qnorm_bonf * se,
         upr = est + qnorm_bonf * se,
         p = p.adjust(p, method = 'bonferroni')) |> 
  select(method, data_id, taxon, est, se, p, lwr, upr)

res_bayes <- res_blogrrall |> 
  map(~ .$res) |> 
  bind_rows() |> 
  filter(data_id %in% c(null_data, original_data)) |> 
  mutate(method = 'BDPE') |> 
  select(method, data_id, taxon, est = est_pa, lwr = lwr90_pa, upr = upr90_pa)

res_freq |> group_by(data_id) |> summarise(n = n())
res_bonf |> group_by(data_id) |> summarise(n = n())
res_bayes |> group_by(data_id) |> summarise(n = n())

setdiff(res_bayes$taxon, res_bonf$taxon)


#-------------------------------------------------------------------------------
res_null <- bind_rows(res_freq, res_bonf, res_bayes) |> 
  left_join(meta_taxa, by = c('data_id', 'taxon')) |>
  left_join(metas, by = 'data_id') |>
  filter(data_id == null_data) |> 
  filter(taxon != 'g__' & !str_detect(taxon, 'g___')) |> 
  filter(prevl != 1) |> 
  mutate(method = factor(method, 
                         levels = c('BDPE',
                                    'Frequentist (Wald)\n(unadjusted)', 
                                    'Frequentist (Wald)\n(multiplicity adjusted)')),
         taxon = str_replace(taxon, 'g__', ''),
         Result = ifelse(lwr > 0 | upr < 0,
                         'Excludes OR = 1\n(significant)',
                         'Includes OR = 1\n(non-significant)'),
         Result = factor(Result,
                         levels = c('Includes OR = 1\n(non-significant)',
                                    'Excludes OR = 1\n(significant)'))) |> 
  arrange(taxon, method) |> 
  dplyr::slice(1:75) |> 
  filter(abs(est) < 10)


or_breaks <- c(0.01, 0.1, 1, 10, 100)

p1 <- ggplot(res_null, aes(x = est, y = taxon, color = Result)) +
  geom_vline(xintercept = 0, linetype = '22', color = 'grey0') +
  geom_point(size = 1) +
  geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.0) +
  facet_grid(. ~ method) +
  scale_x_continuous(
    breaks = log(or_breaks),
    labels = or_breaks 
  ) +
  labs(x = 'Differential prevalence (OR)') +
  scale_y_discrete(limits = rev) +
  coord_cartesian(xlim = c(-5.5, 5.5)) +
  scale_color_manual(
    values = c('Excludes OR = 1\n(significant)' = '#FA7413',
               'Includes OR = 1\n(non-significant)' = 'grey40')) +
  theme_classic(base_size = 10) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(size = 7),
        strip.text.x = element_text(size = 9),
        panel.grid.minor = element_blank(),
        legend.position = 'top',
        legend.title = element_blank(),
        legend.text = element_text(size = 9),
        panel.border = element_rect(colour = "grey0",
                                    fill = NA,
                                    linewidth = 0.5),
        axis.line = element_line(color = 'grey0'),
        axis.ticks = element_line(color = "grey0"),
        plot.margin = margin(t = 5, r = 5, b = 10, l = 5))


#-------------------------------------------------------------------------------
res_orig <- bind_rows(res_freq, res_bonf, res_bayes) |> 
  left_join(meta_taxa, by = c('data_id', 'taxon')) |>
  left_join(metas, by = 'data_id') |>
  filter(data_id == original_data) |> 
  filter(taxon != 'g__' & !str_detect(taxon, 'g___')) |> 
  filter(prevl != 1) |> 
  mutate(method = factor(method, 
                         levels = c('BDPE',
                                    'Frequentist (Wald)\n(unadjusted)', 
                                    'Frequentist (Wald)\n(multiplicity adjusted)')),
         taxon = str_replace(taxon, 'g__', ''),
         Result = ifelse(lwr > 0 | upr < 0,
                         'Excludes OR = 1\n(significant)',
                         'Includes OR = 1\n(non-significant)'),
         Result = factor(Result,
                         levels = c('Includes OR = 1\n(non-significant)',
                                    'Excludes OR = 1\n(significant)'))) |> 
  arrange(taxon, method) |> 
  dplyr::slice(1:75)

res_orig_1 <- res_orig |> filter(abs(est) < 10)

res_orig_2 <- res_orig |> 
  filter(abs(est) > 10) |> 
  mutate(est = -0.3)


(p2 <- ggplot(res_orig_1, aes(x = est, y = taxon, color = Result)) +
    geom_vline(xintercept = 0, linetype = '22', color = 'grey0') +
    geom_point(size = 1) +
    geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.0) +
    geom_text(data = res_orig_2, aes(x = est, y = taxon),
              label = 'N/A', color = 'grey40', size = 2.0,
              vjust = 0.5, hjust = 1.0, fontface = 'bold') +
    facet_grid(. ~ method) +
    scale_x_continuous(
      breaks = log(or_breaks),
      labels = or_breaks 
    ) +
    labs(x = 'Differential prevalence (OR)') +
    scale_y_discrete(limits = rev) +
    coord_cartesian(xlim = c(-5.5, 5.5)) +
    scale_color_manual(
      values = c('Excludes OR = 1\n(significant)' = '#FA7413',
                 'Includes OR = 1\n(non-significant)' = 'grey40')) +
    theme_classic(base_size = 10) +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_text(size = 7),
          strip.text.x = element_text(size = 9),
          panel.grid.minor = element_blank(),
          legend.position = 'top',
          legend.title = element_blank(),
          legend.text = element_text(size = 9),
          panel.border = element_rect(colour = "grey0",
                                      fill = NA,
                                      linewidth = 0.5))
)


#-------------------------------------------------------------------------------
(p <- (p1 / p2) +
   plot_layout(nrow = 2, ncol = 1, guides = 'collect') +
   plot_annotation(tag_levels = "a") &
   theme(legend.position = 'bottom',
         legend.title = element_blank(),
         legend.text = element_text(size = 9),
         plot.tag = element_text(size = 14, face = "bold"))
)

date <- format(Sys.Date(), "%d%m%y")
ggsave(filename = paste0('fig_3_illustration_', date,'.png'),
       plot = p, width = 170, height = 180, dpi = 300, unit = 'mm',
       bg = 'white')
