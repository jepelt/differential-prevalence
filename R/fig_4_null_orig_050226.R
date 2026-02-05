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

sl <- 0.10


# Null data results ------------------------------------------------------------
glm_ci <- function(x, n) {
  prop.test(x, n, conf.level = 0.90)$conf.int |> as.numeric()
}

res_null <- res |> 
  filter(str_detect(data_id, 'null')) |> 
  mutate(sgn = q < sl) |> 
  group_by(method, data_id) |> 
  summarize(s = sum(sgn),
            t = n(),
            .groups = 'drop') |> 
  group_by(method) |>
  summarize(fdr = mean(s > .5),
            ci_lower = glm_ci(sum(s > .5), n()) [1],
            ci_upper = glm_ci(sum(s > .5), n()) [2],
            ms = mean(s),
            n_datasets = n(),
            .groups = 'drop')


res_null_w <- res |> 
  filter(str_detect(data_id, 'null')) |> 
  filter(!is.na(q)) |>
  mutate(sgn = q < sl) |> 
  group_by(method, data_id) |> 
  summarize(s = sum(sgn),
            .groups = 'drop') |> 
  ungroup() |> 
  pivot_wider(names_from = method, values_from = s)

res_null_wm <- res |> 
  filter(str_detect(data_id, 'null')) |> 
  filter(!is.na(q)) |>
  mutate(sgn = q < sl) |> 
  group_by(method, data_id) |> 
  summarize(m = round(mean(sgn), 3),
            .groups = 'drop') |> 
  ungroup() |> 
  pivot_wider(names_from = method, values_from = m)
  


# Original datasets results ----------------------------------------------------
res_orig <- res |> 
  filter(str_detect(data_id, 'original')) |> 
  mutate(sgn = q < sl) |> 
  group_by(data_id, method) |> 
  summarize(s0 = sum(sgn), .groups = 'drop')

res_orig_summary <- res_orig |> 
  group_by(method) |> 
  summarize(md_s0 = median(s0),
            m_s0 = mean(s0),
            q75_s0 = quantile(s0, .75),
            .groups = 'drop')


# Create figure ----------------------------------------------------------------

# Panel a: Scatter plot
d_scatter <- left_join(res_null, res_orig_summary, by = "method")

(p_scatter <- ggplot(d_scatter, aes(x = fdr, y = md_s0, color = method)) +
    # geom_vline(xintercept = 0, linetype = 'solid', color = 'gray40',
    #            linewidth = .5) +
    # geom_hline(yintercept = 0, linetype = 'solid', color = 'gray40',
    #            linewidth = .5) +
    geom_vline(xintercept = 0.10, linetype = 'dashed', color = 'gray40') +
    geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper),
                   height = 0, 
                   linewidth = .7, 
                   alpha = 0.7) +
    geom_point(size = 4, color = "white", alpha = 1) +
    geom_point(size = 4, alpha = .7) +
    geom_text_repel(aes(label = method),
                    fontface = 'plain',
                    size = 3.5, 
                    box.padding = 0.6,   
                    point.padding = 0.5,
                    force = 10,          
                    min.segment.length = 0,
                    show.legend = FALSE,
                    color = "black") +
    scale_color_manual(values = colors) +
    scale_x_continuous(breaks = seq(0, 0.25, 0.05), limits = c(0, 0.22)) +
    scale_y_continuous(breaks = seq(0, 30, 2)) +
    labs(x = TeX("Null data error rate ($\\lambda$)"),
         y = 'Median number of significant findings') +
    theme_classic(base_size = 11) +
    theme(legend.position = 'none',
          axis.title = element_text(size = 11,
                                      margin = margin(t = 0, r = 1, b = 0, l = 0)),
          axis.text = element_text(size = 10))
)



# Panel b: Distribution of significant findings
md_q75 <- function(x) {
  median(x) + .001 * quantile(x, 0.75)
}

d_distr <- res_orig |> 
  ungroup() |> 
  mutate(method = fct_reorder(method, s0, .fun = md_q75))


log10p <- function(x) log10(x + 1)

(p_distr <- ggplot(d_distr, aes(method, log10p(s0), color = method)) + 
    ggbeeswarm::geom_quasirandom(size = 1.2, alpha = .3, width = 0.3,
                                 bandwidth = 0.5) +
    geom_boxplot(alpha = .0, outlier.shape = NA, width = 0.6) +
    
    scale_y_continuous(breaks = log10p(c(0, 1, 3, 10, 30, 100, 300)),
                       labels = c('0', '1', '3', '10', '30', '100', '300')) +
    scale_color_manual(values = colors) +
    labs(y = 'Number of significant findings') +
    theme_classic(base_size = 11) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
          axis.title.y = element_text(size = 11,
                                      margin = margin(t = 0, r = 1, b = 0, l = 0)),
          axis.text.y = element_text(size = 10),
          legend.position = 'none')
)


# Combine panels
(p_final <- free(p_scatter, type = "label") + free(p_distr, type = "label") +
    plot_layout(widths = c(1.0, 0.75)) +
    plot_annotation(tag_levels = 'a') &
    theme(plot.tag = element_text(face = "bold"),
          plot.margin = margin(t = 1, r = 1, b = -4, l = 1, unit = "mm"))
)

date <- format(Sys.Date(), "%d%m%y")
ggsave(plot = p_final,
       filename = paste0('fig_4_null_orig_', date, '.png'),
       width = 170, height = 120, dpi = 300, unit = 'mm', 
       bg = 'white')


#-------------------------------------------------------------------------------
(res_fig_2 <- res_orig_summary |> 
   left_join(res_null, by = 'method')
)

writexl::write_xlsx(res_fig_2, paste0('res_fig_null_orig_', date, '.xlsx'))
