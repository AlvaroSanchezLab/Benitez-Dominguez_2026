# =============================================================================
# Title:   Fig 2. Two target functions at contrasting levels of functional
#          complexity, and the monoculture performance that motivates
#          selecting communities (bottom panels)
# Author:  Belén Benítez Domínguez (IBFG-CSIC, Salamanca)
# GitHub:  https://github.com/AlvaroSanchezLab/Benitez-Dominguez_2026
# =============================================================================

# install.packages(c("ggplot2", "ggh4x", "svglite"))
library(ggplot2)
library(ggh4x)   # one colour per strain in the facet strips

# Paths -----------------------------------------------------------------------

# fig_dir <- "."
fig_dir <- "C:/Users/USUARIO/Downloads/github_figures/Fig2"

data_dir <- file.path(fig_dir, "data")
out_dir  <- file.path(fig_dir, "output")
dir.create(out_dir, showWarnings = FALSE)

strains <- c("SC1", "SC2", "SC3", "Lt24", "Td10", "Mp21")
colors_strain <- c(SC1 = "#E1B12C", SC2 = "#E67E22", SC3 = "#C0392B",
                   Lt24 = "#3498DB", Td10 = "#8E44AD", Mp21 = "#2DD154")
colors_ref <- c("Theoretical maximum (total sugar)" = "grey25",
                "Best it1 consortium" = "#E29CA8",
                "Best it3 consortium" = "#7B1E2B",
                "Unfinished" = "#C0392B")

total_sugar <- c(HSFR = 300, SGM = 200)
max_residual <- 5   # g/L of residual sugar to count a fermentation as finished

txt_axis_title   <- 16
txt_axis         <- 16
txt_strip        <- 16
txt_legend       <- 11
txt_legend_title <- 12


# 1. Read the data ------------------------------------------------------------

# One row per well
read_monoc <- function(file) {
  x <- read.csv(file.path(data_dir, "monocultures", file))
  data.frame(strain = rep(x$strain, 3),
             dose   = rep(x$dose, 3),
             rep    = rep(1:3, each = nrow(x)),
             value  = c(x$R1, x$R2, x$R3))
}

hsfr     <- read_monoc("hsfr_sugar_consumed.csv")
lactic   <- read_monoc("sgm_lactic_acid.csv")
residual <- read_monoc("sgm_residual_sugar.csv")

well <- function(d) paste(d$strain, d$dose, d$rep)
lactic$finished <- residual$value[match(well(lactic), well(residual))] <= max_residual


# 2. Reference lines ----------------------------------------------------------
# Function of the winning consortium of rounds 1 and 3, chosen as in Figs 3B
# and 4B

winner_mean <- function(d, mean_col, sd_col, rep_cols, ignore = NULL) {
  best <- d[which.max(d[[mean_col]] - d[[sd_col]]), ]
  same <- setdiff(strains, ignore)
  d <- d[apply(d[same], 1, function(r) all(r == best[same])), ]
  mean(unlist(d[rep_cols]))
}

read_round <- function(fun, it) {
  read.csv(file.path(data_dir, it, paste0(fun, "_", it, "_data.csv")),
           colClasses = c(Environment = "character"))
}

sugar_cols  <- c("consumed_sugar_R1", "consumed_sugar_R2", "consumed_sugar_R3")
lactic_cols <- c("produced_lactic_R1", "produced_lactic_R2", "produced_lactic_R3")

ref_lines <- data.frame(
  must = c("HSFR", "HSFR", "HSFR", "SGM", "SGM"),
  ref  = c("Theoretical maximum (total sugar)", "Best it1 consortium", "Best it3 consortium",
           "Best it1 consortium", "Best it3 consortium"),
  y    = c(total_sugar[["HSFR"]],
           winner_mean(read_round("hsfr", "it1"), "Mean_sugar", "SD_sugar", sugar_cols),
           winner_mean(read_round("hsfr", "it3"), "Mean_sugar", "SD_sugar", sugar_cols),
           winner_mean(read_round("lactic", "it1"), "Mean_lactic", "SD_lactic", lactic_cols),
           winner_mean(read_round("lactic", "it3"), "Mean_lactic", "SD_lactic", lactic_cols,
                       ignore = "Td10")))
ref_lines$linetype <- ifelse(grepl("Theoretical", ref_lines$ref), "solid", "dashed")


# 3. Panels -------------------------------------------------------------------

fig_theme <- theme_bw() +
  theme(axis.text          = element_text(size = txt_axis),
        axis.title         = element_text(size = txt_axis_title),
        strip.text         = element_text(size = txt_strip, face = "bold", color = "white"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank(),
        panel.spacing      = unit(0.5, "lines"),
        legend.position    = "none")

dose_labels <- c(expression(10^3), "", expression(10^5), "", expression(10^7))

plot_monoc <- function(d, must, y_title, y_max) {
  d$strain <- factor(d$strain, levels = strains)
  d$dose   <- factor(log10(d$dose))
  lines    <- ref_lines[ref_lines$must == must, ]

  p <- ggplot(d, aes(x = dose, y = value)) +
    geom_boxplot(fill = "#0B6E6E", color = "grey30", width = 0.7,
                 linewidth = 0.3, outlier.size = 0.4) +
    geom_hline(data = lines, aes(yintercept = y, color = ref, linetype = linetype),
               linewidth = 1.2)

  if ("finished" %in% names(d)) {
    p <- p + geom_boxplot(data = d[!d$finished, ], aes(color = "Unfinished"),
                          fill = NA, width = 0.7, linewidth = 0.8, outlier.shape = NA)
  }

  p +
    facet_wrap2(~ strain, nrow = 2,
                strip = strip_themed(background_x = elem_list_rect(fill = colors_strain))) +
    scale_color_manual(values = colors_ref) +
    scale_linetype_identity() +
    scale_x_discrete(labels = dose_labels) +
    expand_limits(y = c(0, y_max)) +
    labs(x = "Inoculum dose (CFU/mL)", y = y_title) +
    fig_theme
}

p_hsfr <- plot_monoc(hsfr, "HSFR", "Sugar consumed (g/L)", 315)
p_sgm  <- plot_monoc(lactic, "SGM", "L-lactic acid (g/L)", 13)


# 4. Legends ------------------------------------------------------------------
# Each legend is saved on its own, to place it when the figure is assembled.

get_legend <- function(p) {
  g <- ggplotGrob(p)
  g$grobs[[grep("guide-box", g$layout$name)[1]]]
}

legend_theme <- theme_bw(base_size = 14) +
  theme(legend.title = element_text(size = txt_legend_title, face = "bold"),
        legend.text  = element_text(size = txt_legend),
        legend.key   = element_blank())

legend_lines <- function(keys, types, title) {
  d <- data.frame(key = keys)
  p <- ggplot(d) +
    geom_hline(aes(yintercept = 1, color = key, linetype = key), linewidth = 1.2) +
    scale_color_manual(values = colors_ref[keys], name = title) +
    scale_linetype_manual(values = types, name = title) +
    legend_theme
  get_legend(p)
}

leg_consortia   <- legend_lines(c("Best it1 consortium", "Best it3 consortium"),
                                c("dashed", "dashed"), "Consortium thresholds")
leg_theoretical <- legend_lines("Theoretical maximum (total sugar)", "solid",
                                "Theoretical limit")

p <- ggplot(data.frame(x = 1, y = 1), aes(x, y, color = "Unfinished")) +
  geom_boxplot(fill = NA, linewidth = 0.8) +
  scale_color_manual(values = colors_ref["Unfinished"], name = "Fermentation state") +
  legend_theme
leg_unfinished <- get_legend(p)


# 5. Save ---------------------------------------------------------------------

ggsave(file.path(out_dir, "Fig2_monoculture_HSFR.svg"), p_hsfr, width = 8, height = 6)
ggsave(file.path(out_dir, "Fig2_monoculture_SGM.svg"),  p_sgm,  width = 8, height = 6)
ggsave(file.path(out_dir, "Fig2_legend_consortia.svg"),   leg_consortia,   width = 3.5, height = 1.8)
ggsave(file.path(out_dir, "Fig2_legend_unfinished.svg"),  leg_unfinished,  width = 3.0, height = 1.5)
ggsave(file.path(out_dir, "Fig2_legend_theoretical.svg"), leg_theoretical, width = 4.0, height = 1.5)
