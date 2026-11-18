# =============================================================================
# Title:   Fig 4. Lactic acid production is interactive
# Author:  Belén Benítez Domínguez (IBFG-CSIC, Salamanca)
# GitHub:  https://github.com/AlvaroSanchezLab/Benitez-Dominguez_2026
# =============================================================================

# install.packages(c("ggplot2", "patchwork", "scales", "svglite"))
library(ggplot2)
library(patchwork)
library(scales)

# Paths -----------------------------------------------------------------------

# fig_dir <- "."
fig_dir <- "C:/Users/USUARIO/Downloads/github_figures/Fig4"

data_dir <- file.path(fig_dir, "data")
out_dir  <- file.path(fig_dir, "output")
dir.create(out_dir, showWarnings = FALSE)

strains   <- c("Mp21", "Lt24", "Td10", "SC3", "SC2", "SC1")
rounds    <- c("it1", "it2", "it3")
colors_it <- c(it1 = "#E29CA8", it2 = "#B5485D", it3 = "#7B1E2B")

txt_axis_title <- 11
txt_axis       <- 9
txt_strip      <- 10
txt_legend     <- 9
txt_tag        <- 12
txt_note       <- 9 / .pt


# 1. Read the data ------------------------------------------------------------

hsfr   <- list()
lactic <- list()
for (it in rounds) {
  hsfr[[it]]   <- read.csv(file.path(data_dir, it, paste0("hsfr_", it, "_data.csv")),
                           colClasses = c(Environment = "character"))
  lactic[[it]] <- read.csv(file.path(data_dir, it, paste0("lactic_", it, "_data.csv")),
                           colClasses = c(Environment = "character"))
  hsfr[[it]]$func   <- hsfr[[it]]$Mean_sugar
  lactic[[it]]$func <- lactic[[it]]$Mean_lactic
}


# 2. Models (Methods, Eqs M1-M7) ----------------------------------------------
# In:  hsfr and lactic (both are needed: panel E shows the other function in grey)
# Out: r2_table   R2 of the null, additive and additive-plus-pairwise models (panel C)
#      term_table fraction of variance of each term (panel D)
#      pair_table normalised interaction I_ij of each pair (panel E)

all_data <- list("HS-FR" = hsfr, "Lactic acid" = lactic)

r2_table   <- data.frame()
term_table <- data.frame()
pair_table <- data.frame()

for (fun in names(all_data)) {
  for (it in rounds) {
    d <- all_data[[fun]][[it]]

    varying <- strains[sapply(strains, function(s) length(unique(d[[s]])) > 1)]

    fit_add  <- lm(as.formula(paste("func ~", paste(varying, collapse = " + "))), data = d)
    fit_pair <- lm(as.formula(paste("func ~ (", paste(varying, collapse = " + "), ")^2")), data = d)

    # Panel C
    r2_table <- rbind(r2_table, data.frame(
      fun = fun, iteration = it, order = 0:2,
      r2 = c(0, summary(fit_add)$r.squared, summary(fit_pair)$r.squared)))

    # Panel D
    a <- anova(fit_pair)
    frac <- a[["Sum Sq"]] / sum(a[["Sum Sq"]])
    keep <- rownames(a) != "Residuals"
    term_table <- rbind(term_table, data.frame(
      fun = fun, iteration = it,
      term  = rownames(a)[keep],
      order = ifelse(grepl(":", rownames(a)[keep]), 2, 1),
      frac  = frac[keep]))

    # Panel E
    b <- coef(fit_pair)
    for (pair in combn(varying, 2, simplify = FALSE)) {
      b_ij  <- b[[paste(pair, collapse = ":")]]
      F_add <- b[["(Intercept)"]] + b[[pair[1]]] + b[[pair[2]]]
      pair_table <- rbind(pair_table, data.frame(
        fun = fun, iteration = it,
        pair = paste(pair, collapse = " x "),
        I_ij = b_ij / F_add))
    }
  }
}

pair_table <- pair_table[!is.na(pair_table$I_ij), ]


# 3. Panels -------------------------------------------------------------------
fig_theme <- theme_minimal(base_size = txt_axis) +
  theme(axis.text        = element_text(size = txt_axis, color = "black"),
        axis.title       = element_text(size = txt_axis_title, face = "bold"),
        axis.title.x     = element_text(margin = margin(t = 5)),
        axis.title.y     = element_text(margin = margin(r = 5)),
        strip.text       = element_text(size = txt_strip, face = "bold"),
        strip.background = element_rect(fill = "grey95", color = NA),
        legend.text      = element_text(size = txt_legend),
        plot.tag         = element_text(size = txt_tag, face = "bold"),
        panel.grid.minor = element_blank())

lactic_cols <- c("produced_lactic_R1", "produced_lactic_R2", "produced_lactic_R3")

# Panel A: distribution of lactic acid produced in all replicate wells of each round
dens <- data.frame()
for (it in rounds) {
  dens <- rbind(dens, data.frame(iteration = it,
                                 value = unlist(lactic[[it]][lactic_cols])))
}
dens$iteration <- factor(dens$iteration, levels = c("it3", "it2", "it1"))

pA <- ggplot(dens, aes(x = value, fill = iteration, color = iteration)) +
  geom_density(alpha = 0.45, linewidth = 0.55, bw = 0.42) +
  facet_wrap(~ iteration, ncol = 1, strip.position = "right") +
  scale_fill_manual(values = colors_it) +
  scale_color_manual(values = colors_it) +
  scale_x_continuous(limits = c(0, 14), breaks = seq(0, 14, 2)) +
  scale_y_continuous(breaks = c(0, 0.5), labels = c("0", "0.5")) +
  labs(x = "Produced lactic acid (g/L)", y = "Density", tag = "A") +
  fig_theme +
  theme(legend.position = "none",
        panel.spacing = unit(0.6, "lines"))

# Panel B: replicate wells of the winning community of each round (highest
# mean minus SD). In round 3 Td10 is absent at both doses, so the winner sits
# in two plate positions and both are used.
win <- data.frame()
for (it in rounds) {
  d <- lactic[[it]]
  best <- d[which.max(d$Mean_lactic - d$SD_lactic), ]
  same <- if (it == "it3") setdiff(strains, "Td10") else strains
  d <- d[apply(d[same], 1, function(r) all(r == best[same])), ]
  win <- rbind(win, data.frame(iteration = it,
                               x = match(it, rounds),
                               value = unlist(d[lactic_cols])))
}
win_mean <- aggregate(value ~ iteration + x, data = win, FUN = mean)

pB <- ggplot() +
  geom_line(data = win_mean, aes(x = x, y = value),
            linetype = "dashed", color = "grey20", linewidth = 0.4) +
  geom_point(data = win_mean, aes(x = x, y = value, color = iteration), size = 3) +
  geom_point(data = win, aes(x = x, y = value, color = iteration),
             position = position_jitter(width = 0.05, height = 0, seed = 1),
             shape = 21, size = 1.25, stroke = 0.5) +
  scale_color_manual(values = colors_it, guide = "none") +
  scale_x_continuous(breaks = 1:3, labels = rounds, limits = c(0.5, 3.6)) +
  labs(x = "Iteration", y = "Produced\nlactic acid (g/L)", tag = "B") +
  fig_theme +
  theme(panel.grid.major.x = element_blank(),
        axis.line.x = element_line(color = "grey30", linewidth = 0.3))

# Panel C: R² of the null, additive and additive-plus-pairwise models
d_C <- r2_table[r2_table$fun == "Lactic acid", ]

pC <- ggplot(d_C, aes(x = order, y = r2, color = iteration)) +
  geom_hline(yintercept = 1, linetype = "dotted", color = "grey60") +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1.6) +
  geom_text(data = d_C[d_C$order == 2, ], aes(label = iteration),
            hjust = -0.35, size = txt_note, fontface = "bold") +
  scale_color_manual(values = colors_it, guide = "none") +
  scale_x_continuous(breaks = 0:2, labels = c("null", "additive", "additive\n+ pairwise"),
                     expand = expansion(mult = c(0.08, 0.26))) +
  scale_y_continuous(labels = percent, limits = c(0, 1.05)) +
  labs(x = "Maximum\ninteraction order", y = "Cumulative R\u00b2", tag = "C") +
  fig_theme

# Panel D: fraction of variance explained by the six largest terms of each round
d_D <- data.frame()
for (it in rounds) {
  x <- term_table[term_table$fun == "Lactic acid" & term_table$iteration == it, ]
  d_D <- rbind(d_D, head(x[order(-x$frac), ], 6))
}
d_D$bar <- reorder(paste(d_D$iteration, d_D$term, sep = "@"), d_D$frac)
d_D$order <- factor(d_D$order, levels = c(1, 2))

pD <- ggplot(d_D, aes(x = frac, y = bar, fill = order)) +
  geom_col(width = 0.75) +
  facet_wrap(~ iteration, ncol = 3, scales = "free_y") +
  scale_y_discrete(labels = function(z) sub(".*@", "", z)) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1"),
                     expand = expansion(mult = c(0.02, 0.04))) +
  scale_fill_manual(values = c("1" = "#B2182B", "2" = "#EF8A62"),
                    labels = c("1" = "Additive effect", "2" = "Pairwise interaction"),
                    name = NULL) +
  labs(x = "Fraction of variance explained", y = NULL, tag = "D") +
  fig_theme +
  theme(legend.position = "bottom",
        legend.key.size = unit(0.3, "cm"),
        panel.grid.major.y = element_blank())

# Panel E: interaction of each pair relative to its additive expectation
# (lactic acid in colour, HS-FR in grey)
pair_table$fun <- factor(pair_table$fun, levels = c("Lactic acid", "HS-FR"))

pE <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey55") +
  geom_point(data = pair_table[pair_table$fun == "HS-FR", ],
             aes(x = iteration, y = I_ij, shape = fun),
             position = position_jitter(width = 0.18, height = 0, seed = 1),
             color = "grey75", alpha = 0.7, size = 1.6) +
  geom_point(data = pair_table[pair_table$fun == "Lactic acid", ],
             aes(x = iteration, y = I_ij, color = iteration, shape = fun),
             position = position_jitter(width = 0.18, height = 0, seed = 1),
             alpha = 0.85, size = 2) +
  scale_color_manual(values = colors_it, guide = "none") +
  scale_shape_manual(values = c(16, 16), breaks = c("Lactic acid", "HS-FR"), name = NULL,
                     guide = guide_legend(override.aes = list(
                       color = c(colors_it[["it3"]], "grey65"), size = 2, alpha = 1))) +
  scale_y_continuous(transform = pseudo_log_trans(sigma = 0.1),
                     breaks = c(-10, -1, 0, 1, 10)) +
  labs(x = "Iteration", y = expression(bold("\u03b5 / F")[bold("add")]), tag = "E") +
  fig_theme +
  theme(legend.position = "bottom",
        legend.direction = "vertical",
        panel.grid.major.x = element_blank())


# 4. Assemble and save --------------------------------------------------------
top_row    <- pA + plot_spacer() + pB + plot_spacer() + pC +
  plot_layout(widths = c(3.5, 0.1, 3.6, 0.1, 4.7))
bottom_row <- pD + plot_spacer() + pE +
  plot_layout(widths = c(8.8, 0.3, 3.35))

fig4 <- wrap_elements(full = top_row) / wrap_elements(full = bottom_row) +
  plot_layout(heights = c(1.3, 1.1))

ggsave(file.path(out_dir, "Fig4_lactic_acid.svg"), fig4, width = 7.5, height = 6.2, bg = "white")
ggsave(file.path(out_dir, "Fig4_lactic_acid.pdf"), fig4, width = 7.5, height = 6.2, device = cairo_pdf)
ggsave(file.path(out_dir, "Fig4_lactic_acid.png"), fig4, width = 7.5, height = 6.2, bg = "white", dpi = 300)
