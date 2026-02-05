library(tidyverse)
library(patchwork)
library(logistf)

# Function to run logistic regression for presence/absence data
run_dpa <- function(df) {
  
  df <- mutate(df, group = factor(group, levels = c("Control", "Case")))
  
  check_separation <- df |>
    group_by(group) |>
    summarise(is_constant = n_distinct(pa) == 1, .groups = "drop")
  is_edge <- any(check_separation$is_constant)
  
  fit <- logistf(pa ~ group, data = df, alpha = alpha)
  
  list(
    est = coef(fit)["groupCase"],
    lower = fit$ci.lower["groupCase"],
    upper = fit$ci.upper["groupCase"],
    significant = fit$prob["groupCase"] < alpha,
    edge_case = is_edge
  )
}

# Data generation and parameters -----------------------------------------------
fs_axis_text  <- 9
fs_axis_title <- 10
fs_strip_text <- 12
fs_lgnd_text  <- 9
fs_lgnd_title <- 10
key_size      <- 0.78
alpha         <- 0.20 

n_per_group <- 5
n_taxa <- 5

d <- tibble(
  id = rep(1:(n_per_group * 2), n_taxa),
  group = rep(c(rep("Control", n_per_group), rep("Case", n_per_group)),
              n_taxa),
  taxon = rep(paste("Feature", LETTERS[1:n_taxa]), each = n_per_group * 2),
  
  pa = c(1, 1, 1, 1, 0,   0, 1, 0, 0, 0, # Feature A
         1, 0, 0, 1, 1,   1, 1, 1, 1, 0, # Feature B
         0, 0, 0, 0, 0,   0, 1, 1, 0, 1, # Feature C
         0, 1, 1, 0, 0,   0, 0, 0, 1, 0, # Feature D
         1, 1, 1, 1, 1,   1, 1, 1, 0, 1) # Feature E
) |> 
  mutate(taxon = fct_rev(taxon))

res <- d |> 
  group_by(taxon) |> 
  reframe(
    as_tibble(run_dpa(pick(everything())))
  )


# A: P/A matrix figure ---------------------------------------------------------

dA <- d |> 
  mutate(sample = factor(id, labels = paste("Subject", 1:(2 * n_per_group))),
         group = factor(group, levels = c("Control", "Case")),
         pa = factor(pa),
         taxon = factor(taxon, levels = paste("Feature", LETTERS[1:n_taxa])),
         taxon = fct_rev(taxon))

(pA <- ggplot(dA, aes(x = sample, y = taxon, fill = pa)) +
    geom_tile(color = "white", size = 1.5) +
    facet_grid(~ group, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = c("0" = "grey80", "1" = "#2C7BB6"),
                      name = "Observation", 
                      labels = c("Absent", "Present")) +
    labs(y = "Taxonomic features", x = "Samples/subjects") +
    theme_classic() +
    theme(
      axis.title.x = element_text(size = fs_axis_title, color = 'grey30'),
      axis.title.y = element_text(size = fs_axis_title, color = 'grey30',
                                  margin = margin(r = 10)),
      axis.text.x = element_text(size = fs_axis_text, angle = 90, vjust = 0.5,
                                 hjust = 1, color = "grey30",
                                 margin = margin(t = 0, unit = "pt")),
      axis.text.y = element_text(size = fs_axis_text, color = "grey30"),
      strip.background = element_blank(),
      strip.text = element_text(size = fs_strip_text, face = "plain",
                                color = "black", vjust = 1),
      panel.spacing = unit(0.5, "lines"),
      legend.title = element_text(size = fs_lgnd_title),
      legend.text = element_text(size = fs_lgnd_text),
      legend.position = "bottom",
      legend.key.size = unit(key_size, "cm"), ,
      legend.box.spacing = unit(0, "pt"),
      legend.margin = margin(t = 0, b = 0, unit = "pt"),
      legend.box.margin = margin(t = -5),
      axis.line = element_blank(),
      axis.ticks = element_blank() 
    )
)


# B: DPA figure ----------------------------------------------------------------

dB <- res |> 
  mutate(taxon = factor(taxon, levels = paste("Feature", LETTERS[1:n_taxa])),
         taxon = fct_rev(taxon),
         sig_class = ifelse(significant, "Significant", "Non-significant"))

(pB <- ggplot(dB, aes(x = exp(est), y = taxon)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey60") +
    geom_errorbarh(aes(xmin = exp(lower), xmax = exp(upper), color = edge_case),
                   height = 0.0, size = 0.6) +
    geom_point(aes(color = edge_case), size = 1.5) +
    geom_text(data = subset(dB, edge_case), 
              aes(label = "?"), 
              vjust = -0.4,
              color = "grey50", 
              size = 5, 
              fontface = "bold") +
    
    scale_x_log10(breaks = c(0.1, 1, 10), 
                  labels = c("0.1", "1", "10")) +
    scale_color_manual(values = c("TRUE" = "grey50", "FALSE" = "black")) +
    labs(x = "Differential prevlence (OR)", y = NULL) +
    theme_classic() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),
      axis.text.x = element_text(size = fs_axis_text, color = "grey30"),
      axis.title.x = element_text(size = fs_axis_title, 
                                  margin = margin(t = 12, unit = "pt")),
      strip.background = element_blank(),
      strip.text = element_text(size = fs_strip_text, face = "plain", 
                                color = "black", vjust = 1),
      legend.position = "none"
    ))


# Combine panels A and B -------------------------------------------------------

(p_combined <- pA + free(pB, type = 'label') + 
   plot_layout(widths = c(1.4, 1)) +
   plot_annotation(tag_levels = 'a') & 
   theme(
     plot.tag = element_text(face = "bold"),
     plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "mm") 
   )
)

date <- format(Sys.Date(), "%d%m%y")
ggsave(plot = p_combined, file = paste0("fig_1_schematic_", date, ".png"),
       width = 170, height = 83, unit = 'mm', dpi = 300,
       bg = 'white')
