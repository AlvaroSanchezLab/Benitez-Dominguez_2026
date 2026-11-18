# =============================================================================
# Title:   Fig A in S1 Methods. Contribution of each interaction order compared
#          with experimental noise
# Author:  Belén Benítez Domínguez (IBFG-CSIC, Salamanca)
# GitHub:  https://github.com/AlvaroSanchezLab/Benitez-Dominguez_2026
# =============================================================================

# install.packages(c("ggplot2", "scales", "svglite"))
library(ggplot2)
library(scales)

# Paths -----------------------------------------------------------------------

# fig_dir <- "."
fig_dir <- "C:/Users/USUARIO/Downloads/github_figures/S1_Methods"

data_dir <- file.path(fig_dir, "data")
out_dir  <- file.path(fig_dir, "output")
dir.create(out_dir, showWarnings = FALSE)

strains <- c("Mp21", "Lt24", "Td10", "SC3", "SC2", "SC1")
colors_it <- c(it1 = "#E29CA8", it2 = "#B5485D", it3 = "#7B1E2B")
n_sim   <- 1000
set.seed(123)

txt_axis_title <- 11
txt_axis       <- 9
txt_strip      <- 10
txt_legend     <- 9


# 1. Complete landscapes to analyse -------------------------------------------

files <- data.frame(
  panel  = c("HS-FR it1", "HS-FR it2", "HS-FR it3", "Lactic acid it3"),
  folder = c("it1", "it2", "it3", "it3"),
  file   = c("hsfr_it1_data.csv", "hsfr_it2_data.csv", "hsfr_it3_data.csv",
             "lactic_it3_data.csv"),
  value  = c("consumed_sugar", "consumed_sugar", "consumed_sugar", "produced_lactic"))


# 2. Function -----------------------------------------------------------------

# variance_by_order (S1 Methods, Eqs 1-3)
variance_by_order <- function(F, X) {
  N <- ncol(X)
  V <- numeric(N)
  for (m in 1:(2^N - 1)) {
    S <- which(bitwAnd(m, 2^(0:(N - 1))) > 0)       # strains in term S
    beta_S <- mean(F * apply(X[, S, drop = FALSE], 1, prod))
    V[length(S)] <- V[length(S)] + beta_S^2
  }
  V
}


# 3. Observed landscapes and noise-only null ----------------------------------

order_table <- data.frame()

for (i in 1:nrow(files)) {
  d <- read.csv(file.path(data_dir, files$folder[i], files$file[i]),
                colClasses = c(Environment = "character"))
  varying <- strains[sapply(strains, function(s) length(unique(d[[s]])) > 1)]

  wells <- list()
  for (r in 1:nrow(d)) {
    key <- paste(d[r, varying], collapse = "")
    wells[[key]] <- c(wells[[key]], unlist(d[r, paste0(files$value[i], "_R", 1:3)]))
  }
  X <- 2 * do.call(rbind, lapply(strsplit(names(wells), ""), as.numeric)) - 1
  F <- sapply(wells, mean)

  V_obs <- variance_by_order(F, X)
  V_tot <- sum(V_obs)

  V_null <- replicate(n_sim, {
    F_null <- sapply(wells, function(w) {
      res <- (w - mean(w)) * sqrt(length(w) / (length(w) - 1))
      mean(sample(res, length(w), replace = TRUE))
    })
    variance_by_order(F_null, X)
  })

  order_table <- rbind(order_table, data.frame(
    panel    = files$panel[i],
    order    = seq_along(V_obs),
    observed = V_obs / V_tot,
    null_95  = apply(V_null, 1, quantile, 0.95) / V_tot))
}

order_table$significant <- order_table$observed > order_table$null_95
print(order_table, row.names = FALSE, digits = 3)


# 4. Plot and save ------------------------------------------------------------
# Line: observed fraction of variance of each order. Shaded band: range of the
# same fraction in the noise-only landscapes, up to their 95th percentile.
order_table$iteration <- sub(".* ", "", order_table$panel)
order_table$panel <- factor(order_table$panel, levels = files$panel)

figA <- ggplot(order_table, aes(x = order)) +
  geom_ribbon(aes(ymin = 0, ymax = null_95, fill = "Noise (95th percentile)"), alpha = 0.6) +
  geom_line(aes(y = observed, color = iteration), linewidth = 0.6) +
  geom_point(aes(y = observed, color = iteration), size = 1.6) +
  facet_wrap(~ panel, nrow = 1, scales = "free_x") +
  scale_color_manual(values = colors_it, guide = "none") +
  scale_fill_manual(values = "#F9C9AE", name = NULL) +
  scale_x_continuous(breaks = 1:6) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  labs(x = expression(bold("Interaction order") ~ bolditalic(k)),
       y = expression(bold(V(k) / V[tot]))) +
  theme_minimal(base_size = txt_axis) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = txt_legend),
        strip.text = element_text(face = "bold", size = txt_strip),
        strip.background = element_rect(fill = "grey95", color = NA),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = txt_axis_title),
        axis.text = element_text(size = txt_axis, color = "black"))

ggsave(file.path(out_dir, "S1_Methods_FigA_WH.svg"), figA, width = 6.3, height = 2.7, bg = "white")
ggsave(file.path(out_dir, "S1_Methods_FigA_WH.pdf"), figA, width = 6.3, height = 2.7, device = cairo_pdf)
ggsave(file.path(out_dir, "S1_Methods_FigA_WH.png"), figA, width = 6.3, height = 2.7, bg = "white", dpi = 300)
