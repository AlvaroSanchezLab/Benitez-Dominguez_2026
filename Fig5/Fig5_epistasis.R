# =============================================================================
# Title:   Fig 5. Pairwise interactions are context-dependent
# Author:  Belén Benítez Domínguez (IBFG-CSIC, Salamanca)
# GitHub:  https://github.com/AlvaroSanchezLab/Benitez-Dominguez_2026
# =============================================================================

# install.packages(c("ggplot2", "patchwork", "svglite"))
library(ggplot2)
library(patchwork)

# Paths -----------------------------------------------------------------------

# fig_dir <- "."
fig_dir <- "C:/Users/USUARIO/Downloads/github_figures/Fig5"

data_dir <- file.path(fig_dir, "data")
out_dir  <- file.path(fig_dir, "output")
dir.create(out_dir, showWarnings = FALSE)

colors_it <- c(it1 = "#E29CA8", it2 = "#B5485D", it3 = "#7B1E2B")

txt_axis_title <- 11
txt_axis       <- 9
txt_legend     <- 9
txt_tag        <- 12
txt_note       <- 9 / .pt


# 1. Read the data ------------------------------------------------------------

lactic_it1 <- read.csv(file.path(data_dir, "it1", "lactic_it1_data.csv"),
                       colClasses = c(Environment = "character"))
lactic_it3 <- read.csv(file.path(data_dir, "it3", "lactic_it3_data.csv"),
                       colClasses = c(Environment = "character"))


# 2. Functions ----------------------------------------------------------------

# diamond_nodes (Methods, Eq. M6)
# A community that did not finish the fermentation is not in the table: F = 0.
# A community in two plate positions: both positions are averaged.
diamond_nodes <- function(d, A, B, background) {
  node <- function(a, b) {
    rows <- d[[A]] == a & d[[B]] == b
    for (s in names(background)) rows <- rows & d[[s]] == background[[s]]
    if (sum(rows) == 0) return(c(F = 0, se = 0))
    se_each <- d$SD_lactic[rows] / sqrt(3)
    c(F = mean(d$Mean_lactic[rows]), se = sqrt(sum(se_each^2)) / sum(rows))
  }
  n <- rbind(ab = node(0, 0), Ab = node(1, 0), aB = node(0, 1), AB = node(1, 1))
  additive <- n["Ab", "F"] + n["aB", "F"] - n["ab", "F"]
  list(A = A, B = B, nodes = n, additive = additive, eps = n["AB", "F"] - additive)
}

# plot_diamond: one diamond, measured in colour and additive expectation in grey
plot_diamond <- function(dn, color, round_label, x_labels) {
  n <- dn$nodes
  nodes <- data.frame(x = c(0, 1, 1, 2), F = n[, "F"], se = n[, "se"])
  measured <- data.frame(x = c(0, 0, 1, 1), xend = c(1, 1, 2, 2),
                         y = n[c("ab", "ab", "Ab", "aB"), "F"],
                         yend = n[c("Ab", "aB", "AB", "AB"), "F"])
  additive <- data.frame(x = c(1, 1), xend = c(2, 2),
                         y = n[c("Ab", "aB"), "F"], yend = dn$additive)
  gap <- c(n["AB", "F"], dn$additive)

  p <- ggplot() +
    geom_segment(data = measured, aes(x = x, y = y, xend = xend, yend = yend, color = round_label),
                 linewidth = 0.6) +
    geom_segment(data = additive, aes(x = x, y = y, xend = xend, yend = yend,
                                      linetype = "Null (additive)"),
                 color = "grey60", linewidth = 0.6) +
    geom_errorbar(data = nodes, aes(x = x, ymin = F - se, ymax = F + se), width = 0.06, linewidth = 0.4, color = color) +
    geom_point(data = nodes, aes(x = x, y = F), shape = 21, fill = "white", color = color,
               size = 2, stroke = 0.8) +
    geom_point(aes(x = 2, y = dn$additive), shape = 21, fill = "white", color = "grey55",
               size = 2, stroke = 0.8) +
    annotate("text", x = 2.2, y = mean(gap), label = "\u03b5", color = color,
             fontface = "italic", size = txt_note, hjust = 0) +
    annotate("text", x = 2.45, y = mean(gap), label = sprintf("%+.1f", dn$eps),
             color = "grey20", fontface = "bold", size = txt_note, hjust = 0)

  if (abs(diff(gap)) > 0.08 * diff(range(c(n[, "F"], dn$additive)))) {
    p <- p + annotate("segment", x = 2.14, xend = 2.14, y = min(gap), yend = max(gap),
                      arrow = arrow(ends = "both", length = unit(0.12, "cm"), type = "closed"),
                      color = color, linewidth = 0.5)
  }

  p +
    scale_color_manual(values = setNames(color, round_label), name = NULL) +
    scale_linetype_manual(values = c("Null (additive)" = "dashed"), name = NULL) +
    scale_x_continuous(breaks = 0:2, labels = x_labels, limits = c(-0.3, 2.8)) +
    labs(x = NULL, y = "Function F") +
    theme_minimal(base_size = txt_axis) +
    theme(plot.title    = element_text(face = "bold", size = txt_axis_title),
          plot.subtitle = element_text(size = txt_axis, color = "grey30"),
          plot.caption  = element_text(size = txt_axis, color = "grey40", hjust = 0),
          plot.tag      = element_text(face = "bold", size = txt_tag),
          axis.title    = element_text(face = "bold", size = txt_axis_title),
          axis.text     = element_text(size = txt_axis, color = "black"),
          panel.grid.minor = element_blank(),
          legend.position = "bottom",
          legend.text = element_text(size = txt_legend))
}

x_default <- c("0\nboth\nlow", "1\none\nhigh", "2\nboth\nhigh")

# strain_label: strain name next to a node
strain_label <- function(x, y, label, hjust, vjust = 0.5) {
  annotate("text", x = x, y = y, label = label, hjust = hjust, vjust = vjust,
           fontface = "bold", color = "grey20", size = txt_note)
}


# 3. Panels -------------------------------------------------------------------

# Panel A: Lt24 x SC2, the largest pairwise term of round 1 (Fig 4D), in the
# reference background (all other strains at the low dose)
nA <- diamond_nodes(lactic_it1, "Lt24", "SC2",
                    background = c(Mp21 = 0, Td10 = 0, SC3 = 0, SC1 = 0))
pA <- plot_diamond(nA, colors_it[["it1"]], "it1",
                   x_labels = c("both low\n10\u00b3", "one high\n10\u2075", "both high\n10\u2075")) +
  strain_label(0.92, 0.07, "Lt24", hjust = 1) +
  strain_label(1.08, 0.07, "SC2", hjust = 0) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5)) +
  coord_cartesian(ylim = c(-0.05, 0.6)) +
  labs(title = "Lt24 \u00d7 SC2", subtitle = "Synergy (non-finishing background)", tag = "A",
       caption = "Background: Mp21, Td10, SC3, SC1 at 10\u00b3")

# Panel B: the same pair in round 3
nB <- diamond_nodes(lactic_it3, "Lt24", "SC2",
                    background = c(Mp21 = 1, Td10 = 0, SC3 = 0, SC1 = 1))
pB <- plot_diamond(nB, colors_it[["it3"]], "it3", x_labels = x_default) +
  strain_label(0.74, nB$nodes["Ab", "F"], "Lt24 10\u2077", hjust = 1) +
  strain_label(1.26, nB$nodes["aB", "F"], "SC2 10\u2074", hjust = 0) +
  scale_y_continuous(breaks = c(5, 10)) +
  coord_cartesian(ylim = c(2, 11.7)) +
  labs(title = "Lt24 \u00d7 SC2", subtitle = "Antagonism (finishing background)", tag = "B",
       caption = "Background: Mp21 10\u2074 \u00b7 SC3 10\u2076 \u00b7 SC1 10\u2074\nTd10 absent")

# Panel C: Lt24 x SC3 with SC2 and SC1 low. Changing only the dose of Mp21
# changes the sign of epsilon (left: Mp21 low, right: Mp21 high)
x_C <- c("both low\n10\u2076", "one high\n10\u2077", "both high\n10\u2077")
nC_low  <- diamond_nodes(lactic_it3, "Lt24", "SC3",
                         background = c(Mp21 = 0, Td10 = 0, SC2 = 0, SC1 = 0))
nC_high <- diamond_nodes(lactic_it3, "Lt24", "SC3",
                         background = c(Mp21 = 1, Td10 = 0, SC2 = 0, SC1 = 0))

pC_low <- plot_diamond(nC_low, colors_it[["it3"]], "it3", x_labels = x_C) +
  strain_label(0.82, nC_low$nodes["Ab", "F"], "Lt24", hjust = 1, vjust = 0.3) +
  strain_label(1.18, nC_low$nodes["aB", "F"], "SC3", hjust = 0, vjust = 0.3) +
  coord_cartesian(ylim = c(1.8, 8.9)) +
  labs(title = "Lt24 \u00d7 SC3", subtitle = "Background finishes \u00b7 Mp21 low (10\u00b3)", tag = "C",
       caption = "Background: Mp21 10\u00b3 \u00b7 SC2 10\u00b3 \u00b7 SC1 10\u00b3\nTd10 absent")

pC_high <- plot_diamond(nC_high, colors_it[["it3"]], "it3", x_labels = x_C) +
  strain_label(0.82, nC_high$nodes["Ab", "F"], "Lt24", hjust = 1, vjust = 0.3) +
  strain_label(1.18, nC_high$nodes["aB", "F"], "SC3", hjust = 0, vjust = 1.6) +
  coord_cartesian(ylim = c(1.8, 8.9)) +
  labs(title = "Lt24 \u00d7 SC3", subtitle = "Background finishes \u00b7 Mp21 high (10\u2074)",
       caption = "Background: Mp21 10\u2074 \u00b7 SC2 10\u00b3 \u00b7 SC1 10\u00b3\nTd10 absent")

# Panel D: the diamond whose top node is the winning consortium of round 3
nD <- diamond_nodes(lactic_it3, "Mp21", "Lt24",
                    background = c(Td10 = 0, SC3 = 1, SC2 = 1, SC1 = 1))
pD <- plot_diamond(nD, colors_it[["it3"]], "it3", x_labels = x_default) +
  strain_label(0.74, nD$nodes["Ab", "F"], "Mp21 10\u2074", hjust = 1) +
  strain_label(1.26, nD$nodes["aB", "F"], "Lt24 10\u2077", hjust = 0) +
  annotate("text", x = 2, y = nD$nodes["AB", "F"], label = "winner", vjust = -1.5,
           fontface = "italic", color = "#2166ac", size = txt_note) +
  coord_cartesian(ylim = c(1.6, 16)) +
  labs(subtitle = "Winning consortium", tag = "D",
       caption = "Background: SC3 10\u2077 \u00b7 SC2 10\u2074 \u00b7 SC1 10\u2074\nTd10 absent")


# 4. Assemble and save --------------------------------------------------------
fig5 <- (pA | pB) / (pC_low | pC_high) / (pD | plot_spacer()) +
  plot_layout(guides = "collect") +
  plot_annotation(caption = "Inoculation dose in CFU/mL",
                  theme = theme(legend.position = "bottom",
                                legend.text  = element_text(size = txt_legend),
                                plot.caption = element_text(size = txt_axis)))

ggsave(file.path(out_dir, "Fig5_epistasis.svg"), fig5, width = 6.3, height = 7.9, bg = "white")
ggsave(file.path(out_dir, "Fig5_epistasis.pdf"), fig5, width = 6.3, height = 7.9, device = cairo_pdf)
ggsave(file.path(out_dir, "Fig5_epistasis.png"), fig5, width = 6.3, height = 7.9, bg = "white", dpi = 300)
