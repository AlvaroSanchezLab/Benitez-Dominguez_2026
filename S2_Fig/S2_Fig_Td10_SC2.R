# =============================================================================
# Title:   S2 Fig. Interaction between T. delbrueckii (Td10) and S. cerevisiae
#          SC2 near the high-sugar optimum (third round)
# Author:  Belén Benítez Domínguez (IBFG-CSIC, Salamanca)
# GitHub:  https://github.com/AlvaroSanchezLab/Benitez-Dominguez_2026
# =============================================================================

# install.packages(c("ggplot2", "svglite"))
library(ggplot2)

# Paths -----------------------------------------------------------------------

# fig_dir <- "."
fig_dir <- "C:/Users/USUARIO/Downloads/github_figures/S2_Fig"

data_dir <- file.path(fig_dir, "data")
out_dir  <- file.path(fig_dir, "output")
dir.create(out_dir, showWarnings = FALSE)

col_it3 <- "#7B1E2B"

txt_axis_title <- 11
txt_axis       <- 9
txt_legend     <- 9
txt_note       <- 9 / .pt


# 1. Read the data ------------------------------------------------------------
hsfr_it3 <- read.csv(file.path(data_dir, "it3", "hsfr_it3_data.csv"),
                     colClasses = c(Environment = "character"))


# 2. The four communities of the diamond (Methods, Eq. M6) --------------------
# Background of the round 3 winner: Mp21 and Lt24 low, SC3 and SC1 high
background <- hsfr_it3$Mp21 == 0 & hsfr_it3$Lt24 == 0 & hsfr_it3$SC3 == 1 & hsfr_it3$SC1 == 1
d <- hsfr_it3[background, ]

ab <- d[d$SC2 == 0 & d$Td10 == 0, ]   # both low
Ab <- d[d$SC2 == 1 & d$Td10 == 0, ]   # SC2 high (the round 3 winner)
aB <- d[d$SC2 == 0 & d$Td10 == 1, ]   # Td10 high
AB <- d[d$SC2 == 1 & d$Td10 == 1, ]   # both high

nodes <- data.frame(x  = c(0, 1, 1, 2),
                    F  = c(ab$Mean_sugar, Ab$Mean_sugar, aB$Mean_sugar, AB$Mean_sugar),
                    se = c(ab$SD_sugar, Ab$SD_sugar, aB$SD_sugar, AB$SD_sugar) / sqrt(3))

additive <- Ab$Mean_sugar + aB$Mean_sugar - ab$Mean_sugar
eps      <- AB$Mean_sugar - additive

# The expectation exceeds the 300 g/L available: drawn at 300, epsilon uncapped
additive_drawn <- min(additive, 300)


# 3. Plot and save ------------------------------------------------------------
measured <- data.frame(x = c(0, 0, 1, 1), xend = c(1, 1, 2, 2),
                       y    = nodes$F[c(1, 1, 2, 3)],
                       yend = nodes$F[c(2, 3, 4, 4)])
expected <- data.frame(x = c(1, 1), xend = c(2, 2),
                       y = nodes$F[c(2, 3)], yend = additive_drawn)
gap <- c(AB$Mean_sugar, additive_drawn)

s2 <- ggplot() +
  geom_segment(data = measured, aes(x = x, y = y, xend = xend, yend = yend, color = "it3"),
               linewidth = 0.6) +
  geom_segment(data = expected, aes(x = x, y = y, xend = xend, yend = yend,
                                    linetype = "Null (additive)"),
               color = "grey60", linewidth = 0.6) +
  geom_errorbar(data = nodes, aes(x = x, ymin = F - se, ymax = F + se), width = 0.06, linewidth = 0.4, color = col_it3) +
  geom_point(data = nodes, aes(x = x, y = F), shape = 21, fill = "white", color = col_it3,
             size = 2, stroke = 0.8) +
  geom_point(aes(x = 2, y = additive_drawn), shape = 21, fill = "white", color = "grey55",
             size = 2, stroke = 0.8) +
  annotate("segment", x = 2.14, xend = 2.14, y = min(gap), yend = max(gap),
           arrow = arrow(ends = "both", length = unit(0.12, "cm"), type = "closed"),
           color = col_it3, linewidth = 0.5) +
  annotate("text", x = 2.2, y = mean(gap), label = "\u03b5", color = col_it3,
           fontface = "italic", size = txt_note, hjust = 0) +
  annotate("text", x = 2.45, y = mean(gap), label = sprintf("%+.1f", eps),
           color = "grey20", fontface = "bold", size = txt_note, hjust = 0) +
  annotate("text", x = 0.74, y = Ab$Mean_sugar, label = "SC2 10\u2074", hjust = 1,
           fontface = "bold", color = "grey20", size = txt_note) +
  annotate("text", x = 1.26, y = aB$Mean_sugar, label = "Td10 10\u2077", hjust = 0,
           fontface = "bold", color = "grey20", size = txt_note) +
  annotate("text", x = 1, y = Ab$Mean_sugar, label = "winner", vjust = -1.1,
           fontface = "italic", color = "#2166ac", size = txt_note) +
  scale_color_manual(values = c(it3 = col_it3), name = NULL) +
  scale_linetype_manual(values = c("Null (additive)" = "dashed"), name = NULL) +
  scale_x_continuous(breaks = 0:2, labels = c("0\nboth\nlow", "1\none\nhigh", "2\nboth\nhigh"),
                     limits = c(-0.3, 2.8)) +
  scale_y_continuous(breaks = c(200, 250, 300)) +
  coord_cartesian(ylim = c(195, 300)) +
  labs(title = "HS-FR", subtitle = "Optimum keeps Td10 at 10\u2076",
       caption = "Background: Mp21 10\u2074 \u00b7 Lt24 10\u2074 \u00b7 SC3 10\u2077 \u00b7 SC1 10\u2074",
       x = NULL, y = "Consumed sugar (g/L)") +
  theme_minimal(base_size = txt_axis) +
  theme(plot.title    = element_text(face = "bold", size = txt_axis_title),
        plot.subtitle = element_text(size = txt_axis, color = "grey30"),
        plot.caption  = element_text(size = txt_axis, color = "grey40", hjust = 0),
        axis.title    = element_text(face = "bold", size = txt_axis_title),
        axis.text     = element_text(size = txt_axis, color = "black"),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = txt_legend))

ggsave(file.path(out_dir, "S2_Fig_Td10_SC2.svg"), s2, width = 4.4, height = 4.4, bg = "white")
ggsave(file.path(out_dir, "S2_Fig_Td10_SC2.pdf"), s2, width = 4.4, height = 4.4, device = cairo_pdf)
ggsave(file.path(out_dir, "S2_Fig_Td10_SC2.png"), s2, width = 4.4, height = 4.4, bg = "white", dpi = 300)
