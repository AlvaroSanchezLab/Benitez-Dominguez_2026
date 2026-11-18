# =============================================================================
# Title:   S1 Fig. Contribution of pairwise interactions compared with
#          experimental noise
# Author:  Belén Benítez Domínguez (IBFG-CSIC, Salamanca)
# GitHub:  https://github.com/AlvaroSanchezLab/Benitez-Dominguez_2026
# =============================================================================

# install.packages(c("ggplot2", "scales", "svglite"))
library(ggplot2)
library(scales)

# Paths -----------------------------------------------------------------------

# fig_dir <- "."
fig_dir <- "C:/Users/USUARIO/Downloads/github_figures/S1_Fig"

data_dir <- file.path(fig_dir, "data")
out_dir  <- file.path(fig_dir, "output")
dir.create(out_dir, showWarnings = FALSE)

strains <- c("Mp21", "Lt24", "Td10", "SC3", "SC2", "SC1")
rounds  <- c("it1", "it2", "it3")
n_sim   <- 1000
set.seed(123)

txt_axis_title <- 11
txt_axis       <- 9
txt_strip      <- 10
txt_legend     <- 9


# 1. Files to analyse ---------------------------------------------------------
files <- data.frame(
  fun       = rep(c("Lactic acid", "HS-FR"), each = 3),
  iteration = rep(rounds, 2),
  file      = c("lactic_it1_data.csv", "lactic_it2_data.csv", "lactic_it3_data.csv",
                "hsfr_it1_data.csv", "hsfr_it2_data.csv", "hsfr_it3_data.csv"))


# 2. Pairwise terms against the additive-plus-noise null (Methods, Eq. M4) ----

null_table <- data.frame()

for (i in 1:nrow(files)) {
  d <- read.csv(file.path(data_dir, files$iteration[i], files$file[i]),
                colClasses = c(Environment = "character"))
  rep_cols <- grep("_R[123]$", names(d), value = TRUE)
  d$func <- rowMeans(d[rep_cols])

  se <- apply(d[rep_cols], 1, sd) / sqrt(3)

  varying <- strains[sapply(strains, function(s) length(unique(d[[s]])) > 1)]
  f_add  <- as.formula(paste("func ~", paste(varying, collapse = " + ")))
  f_pair <- as.formula(paste("func ~ (", paste(varying, collapse = " + "), ")^2"))

  fit_add <- lm(f_add, data = d)
  observed <- summary(lm(f_pair, data = d))$r.squared - summary(fit_add)$r.squared

  simulated <- numeric(n_sim)
  for (k in 1:n_sim) {
    d_sim <- d
    d_sim$func <- fitted(fit_add) + rnorm(nrow(d), 0, se)
    simulated[k] <- summary(lm(f_pair, data = d_sim))$r.squared -
                    summary(lm(f_add,  data = d_sim))$r.squared
  }

  null_table <- rbind(null_table, data.frame(
    fun       = files$fun[i],
    iteration = files$iteration[i],
    observed  = observed,
    null_95   = quantile(simulated, 0.95),
    p_value   = mean(simulated >= observed)))
}

null_table$significant <- null_table$p_value < 0.05
null_table$fun <- factor(null_table$fun, levels = c("HS-FR", "Lactic acid"))
print(null_table, row.names = FALSE)


# 3. Plot and save ------------------------------------------------------------
s1 <- ggplot(null_table, aes(x = iteration)) +
  geom_col(aes(y = observed, fill = significant), width = 0.6) +
  geom_segment(aes(x = as.numeric(factor(iteration)) - 0.35,
                   xend = as.numeric(factor(iteration)) + 0.35,
                   y = null_95, yend = null_95),
               color = "red", linewidth = 0.6) +
  facet_wrap(~ fun, ncol = 2) +
  scale_fill_manual(values = c("TRUE" = "#1a9641", "FALSE" = "grey70"),
                    labels = c("TRUE" = "p < 0.05", "FALSE" = "n.s."),
                    name = "Pairwise") +
  scale_y_continuous(labels = percent) +
  labs(x = "Iteration", y = "Marginal R\u00b2 (pairwise)") +
  theme_minimal(base_size = txt_axis) +
  theme(legend.position = "bottom",
        legend.title = element_text(face = "bold", size = txt_legend),
        legend.text = element_text(size = txt_legend),
        strip.text = element_text(face = "bold", size = txt_strip),
        strip.background = element_rect(fill = "grey95", color = NA),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = txt_axis_title, face = "bold"),
        axis.text = element_text(size = txt_axis, color = "black"))

ggsave(file.path(out_dir, "S1_Fig_null.svg"), s1, width = 5.2, height = 3, bg = "white")
ggsave(file.path(out_dir, "S1_Fig_null.pdf"), s1, width = 5.2, height = 3, device = cairo_pdf)
ggsave(file.path(out_dir, "S1_Fig_null.png"), s1, width = 5.2, height = 3, bg = "white", dpi = 300)
