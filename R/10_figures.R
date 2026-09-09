# 10_figures.R
# -----------------------------------------------------------------------------
# Manuscript figures.
#
#   Fig_Survival : proportion of seedlings alive at each census, by origin,
#                  with 95% confidence intervals.
#   Fig_BLUPs    : maternal-family conditional modes (BLUPs) from the
#                  cultivar-origin survival GLMM, showing among-family
#                  variation in survival.
#
# Inputs : data/derived/survival_by_census.csv
#          data/derived/blups_cultivar_survival.csv   (from 05_survival.R)
# Outputs: media/fig_survival.png / .pdf
#          media/fig_blups.png  / .pdf
# -----------------------------------------------------------------------------

source("R/00_setup.R")
suppressPackageStartupMessages({ library(readr); library(dplyr); library(ggplot2) })

theme_ms <- theme_classic(base_size = 11) +
  theme(legend.position = "top", legend.title = element_blank(),
        axis.text = element_text(colour = "black"))

origin_cols <- c(native = "#2c7fb8", cultivar = "#d95f0e")

# ---- Fig_Survival --------------------------------------------------------

sbc <- read_csv(file.path(paths$derived, "survival_by_census.csv"), show_col_types = FALSE) %>%
  mutate(
    date   = as.Date(date),
    origin = factor(origin, levels = c("native", "cultivar")),
    lo = pmax(0, prop_alive - 1.96 * se),
    hi = pmin(1, prop_alive + 1.96 * se)
  )

fig_survival <- ggplot(sbc, aes(date, prop_alive, colour = origin, fill = origin)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.8) +
  scale_colour_manual(values = origin_cols, labels = c("Native", "Cultivar")) +
  scale_fill_manual(values = origin_cols, guide = "none") +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Proportion alive") +
  theme_ms

ggsave(file.path(paths$figures, "fig_survival.png"), fig_survival,
       width = 6.5, height = 4, dpi = 300)
ggsave(file.path(paths$figures, "fig_survival.pdf"), fig_survival, width = 6.5, height = 4)

# ---- Fig_BLUPs ---------------------------------------------------------

blups <- read_csv(file.path(paths$derived, "blups_cultivar_survival.csv"), show_col_types = FALSE) %>%
  mutate(family = factor(family, levels = family[order(blup_logit)]))

fig_blups <- ggplot(blups, aes(blup_logit, family)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
  geom_errorbar(aes(xmin = blup_logit - se, xmax = blup_logit + se),
                orientation = "y", width = 0, colour = origin_cols["cultivar"]) +
  geom_point(colour = origin_cols["cultivar"], size = 1.8) +
  labs(x = "Family conditional mode (logit scale)", y = "Cultivar-origin family") +
  theme_ms + theme(legend.position = "none")

ggsave(file.path(paths$figures, "fig_blups.png"), fig_blups,
       width = 5, height = 4.5, dpi = 300)
ggsave(file.path(paths$figures, "fig_blups.pdf"), fig_blups, width = 5, height = 4.5)

message("wrote fig_survival and fig_blups to ", paths$figures)
