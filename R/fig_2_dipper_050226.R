library(tidyverse)
library(patchwork)
library(ggforce)
library(latex2exp)

# ───────────────────────────────────────────────────────────────
# 0   Global parameters
# ───────────────────────────────────────────────────────────────
BASE_PT    <- 9          # base font of everything (pt)
READS_REL  <- 0.80       # alpha or beta_reads text / normal node text

LW_DENS    <- 0.70       # density curves
LW_VLINE   <- 0.30       # helper v-lines
LW_ARROW   <- 0.30       # arrows
LW_CIRCLE  <- 0.30       # node borders

TIP_S      <- 0.15       # arrow tip size (cm)
NODE_R     <- 0.25       # circle radius (data units)
NODE_R2    <- 0.18       # circle radius for alpha and reads (data units)

PNG_W_CM   <- 170        # final width of the exported figure (mm)
PNG_H_CM   <- 95         # final width of the exported figure (mm)
DPI        <- 300        # resolution


# ───────────────────────────────────────────────────────────────
# 1   Helpers
# ───────────────────────────────────────────────────────────────
dald <- function(x, mu, tau, nu){
  h <- ifelse(x < mu, (1 - nu)*(mu - x), nu*(x - mu))
  2*nu*(1 - nu)/tau * exp(-(2/tau)*h)
}
dhn  <- function(x, mu, sd){
  ifelse(x < mu, 0,
         sqrt(2/pi)/sd * exp(-(x - mu)^2/(2*sd^2)))
}
pt2mm <- function(pt) pt / ggplot2::.pt

TXT_NODE  <- pt2mm(BASE_PT + 2)
TXT_READS <- TXT_NODE * READS_REL
TXT_DOTS  <- pt2mm(BASE_PT + 6)
TXT_ANN   <- pt2mm(BASE_PT)


# ───────────────────────────────────────────────────────────────
# 2   PRIOR PANELS
# ───────────────────────────────────────────────────────────────
len <- 1000

d_tau <- tibble(x = seq(0, 4, length.out = len), y = dhn(x, 0, 1))
max_tau <- max(d_tau$y) + .06
p_tau <- ggplot(d_tau, aes(x, y)) +
  geom_vline(xintercept = 0, linewidth = LW_VLINE) +
  geom_line(colour = "grey0", linewidth = LW_DENS, alpha = .4) +
  scale_x_continuous(limits = c(0, 4.2), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, max_tau)) +
  labs(x = expression(tau[0])) +
  theme_classic(base_size = BASE_PT) +
  theme(axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y  = element_blank(),
        plot.margin = margin(t = 0, r = 5, b = 0, l = 0))

d_nu <- tibble(x = seq(0, 1, length.out = len), y = dald(x, .5, .05, .5))
max_nu <- max(d_nu$y) + .7
p_nu <- ggplot(d_nu, aes(x, y)) +
  geom_vline(xintercept = 0.5, linetype = '33',
             colour = "grey0", linewidth = LW_VLINE) +
  geom_vline(xintercept = c(0, 1), linewidth = LW_VLINE) +
  geom_line(colour = "grey0", linewidth = LW_DENS, alpha = .4) +
  scale_x_continuous(breaks = c(0, .5, 1), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, max_nu)) +
  labs(x = expression(nu[0])) +
  theme_classic(base_size = BASE_PT) +
  theme(axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y  = element_blank(),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 5))


d_al <- expand_grid(x   = seq(-2.5, 2.5, length.out = 3 * len),
                    tau = c(.2, 1.5),
                    nu  = c(.5, .3)) |>
  mutate(y   = dald(x, 0, tau, nu),
         grp = interaction(tau, nu),
         tau = factor(tau),
         nu = fct_rev(factor(nu)))

max_al <- max(d_al$y) + .2
tau_cols <- c("0.2" = "#67ACBC",
              "1.5" = "#F8931F")

p_al <- ggplot(d_al, aes(x, y,
                         group  = grp,
                         colour = tau,
                         linetype = nu)) +
  geom_vline(xintercept = 0, linetype = "solid",
             colour = "grey0", linewidth = LW_VLINE) +
  geom_line(linewidth = LW_DENS, alpha = .8) +
  scale_colour_manual(values = tau_cols,
                      name    = expression(tau[0])) +
  scale_linetype_manual(values = c("solid", "21"),
                        name    = expression(nu[0])) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, max_al)) +
  labs(x = TeX("$\\beta_{j}$ (log(OR))")) +
  theme_classic(base_size = BASE_PT) +
  theme(axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y  = element_blank(),
        legend.position.inside = c(.9, .6),
        legend.position = "inside",
        legend.title = element_text(margin = margin(b = 0)),
        legend.key.height = unit(0.30, "cm"),
        legend.key.width = unit(0.35, "cm"),
        legend.spacing.y = unit(-0.10, "cm"))


# ───────────────────────────────────────────────────────────────
# 3   DIAGRAM
# ───────────────────────────────────────────────────────────────
h1 <- 1
h2 <- 2
h3 <- 3
h4 <- 5

shift <- .5
nodes <- tribble(
  ~name,             ~x,   ~y,
  "tau[0]", 1.5 + shift / 2, h4,
  "nu[0]",  3.5 + shift / 2, h4,
  "beta[1]",         1,   h3,
  "beta[2]",         2,   h3,
  "beta[3]",         3,   h3,
  "beta[K]",    4 + shift, h3,
  "alpha[1]",       0.38, h2,
  "beta[1]^cdots",     0.78, h2,
  "alpha[2]",       1.38, h2,
  "beta[2]^cdots",     1.78, h2,
  "alpha[3]",       2.38, h2,
  "beta[3]^cdots",     2.78, h2,
  "alpha[K]",   3.38 + shift, h2,
  "beta[K]^cdots", 3.78 + shift, h2,
  "bold(y)[1]",      1,   h1,
  "bold(y)[2]",      2,   h1,
  "bold(y)[3]",      3,   h1,
  "bold(y)[K]", 4 + shift, h1
) |> mutate(r = ifelse(str_detect(name, "alpha|cdots"), NODE_R2, NODE_R),
            f_small = str_detect(name, "alpha|cdots"))

nodes <- nodes |>
  mutate(
    name = case_when(
      str_detect(name, "^y\\[i") ~
        # turn   y[i1]  → bold(y)[plain(1)]
        str_replace(name,
                    "^y\\[i(.+)\\]$",
                    "bold(y)[plain(\\1)]"),
      TRUE ~ name
    )
  )

edges <- tribble(
  ~from,                 ~to,
  rep("tau[0]", 4),      c("beta[1]","beta[2]","beta[3]","beta[K]"),
  rep("nu[0]",  4),      c("beta[1]","beta[2]","beta[3]","beta[K]"),
  c("beta[1]","beta[2]","beta[3]","beta[K]"),
  c("bold(y)[1]","bold(y)[2]","bold(y)[3]","bold(y)[K]"),
  c("alpha[1]","alpha[2]","alpha[3]","alpha[K]"),
  c("bold(y)[1]","bold(y)[2]","bold(y)[3]","bold(y)[K]"),
  c("beta[1]^cdots","beta[2]^cdots",
    "beta[3]^cdots","beta[K]^cdots"),
  c("bold(y)[1]","bold(y)[2]","bold(y)[3]","bold(y)[K]")
) |> unnest(c(from, to))

edges <- edges |>
  left_join(nodes |> select(name, x, y, r),
            by = c("from" = "name")) |>
  rename(x_from = x, y_from = y, r_from = r) |>
  left_join(nodes |> select(name, x, y, r),
            by = c("to" = "name")) |>
  rename(x_to = x, y_to = y, r_to = r) |>
  mutate(dx = x_to - x_from,
         dy = y_to - y_from,
         L  = sqrt(dx^2 + dy^2),
         ux = dx / L,
         uy = dy / L,
         x_start = x_from + ux * r_from,
         y_start = y_from + uy * r_from,
         x_end   = x_to   - ux * r_to * 1.10,
         y_end   = y_to   - uy * r_to * 1.10)

diagram <- ggplot() +
  geom_segment(data = edges,
               aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
               arrow = arrow(length = unit(TIP_S, "cm"),
                             type = "open"),
               linewidth = LW_ARROW) +
  ggforce::geom_circle(data = nodes,
                       aes(x0 = x, y0 = y, r = r),
                       fill = "white", linewidth = LW_CIRCLE) +
  geom_text(data = nodes |> filter(!f_small),
            aes(x, y, label = name),
            size = TXT_NODE, parse = TRUE) +
  geom_text(data = nodes |> filter(f_small),
            aes(x, y, label = name),
            size = TXT_READS, parse = TRUE) +
  geom_text(aes(x = 3.5 + shift / 2, y = c(3.1, 1.1), label = "..."),
            size = TXT_DOTS) +
  geom_text(aes(x = -0.8, y = h4, label = "Hyperparameters"),
            hjust = 0, size = TXT_ANN) +
  geom_text(aes(x = -0.8, y = h3, label = "Differential\nprevalence\nparameters"),
            hjust = 0, size = TXT_ANN) +
  geom_text(aes(x = -0.8, y = h2, label = "Nuisance\nparameters"),
            hjust = 0, size = TXT_ANN * READS_REL) +
  geom_text(aes(x = -0.8, y = h1, label = "Presence/\nabsence\ndata"),
            hjust = 0, size = TXT_ANN) +
  # lims(x = c(0, 5), y = c(1.8, 4.2)) +
  coord_fixed(clip = "off") +
  theme_void(base_size = BASE_PT) +
  theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0))


# ───────────────────────────────────────────────────────────────
# 4   PANEL LAYOUT  (patchwork)
# ───────────────────────────────────────────────────────────────
top_right <- (p_tau | p_nu) + plot_layout(ncol = 2, widths = c(1, 1))

right_col <- (top_right / p_al) + plot_layout(heights = c(1, 1.4))

combined  <- (diagram | right_col) + 
  plot_layout(widths = c(3.5, 2)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 14, face = "bold"))


# ───────────────────────────────────────────────────────────────
# 5   EXPORT
# ───────────────────────────────────────────────────────────────
date <- format(Sys.Date(), "%d%m%y")
ggsave(paste0("fig_2_dipper_", date, ".png"),
       combined,
       width  = PNG_W_CM,
       height = PNG_H_CM,
       units  = "mm",
       dpi    = DPI,
       bg     = "white")

