# 09_diagnostics.R
# -----------------------------------------------------------------------------
# Model diagnostics for the linear mixed models, supporting Appendix S1.
#
# For the top-supported model of each single-time trait (origin + family) and
# each repeated-measures trait (plant + family), the script produces:
#   * residual vs fitted plot          (homoscedasticity)
#   * scale-location plot               (homoscedasticity)
#   * normal Q-Q plot of residuals      (normality)
#   * per-family influence (leave-one-out Cook's distance on the fixed effects)
# and a summary table with a Shapiro-Wilk normality test and a
# residual-variance-vs-fitted correlation as a heteroscedasticity check.
#
# The survival GLMM diagnostics (DHARMa) are produced by 05_survival.R.
#
# Inputs : data/derived/traits_harvest.csv, data/derived/traits_repeated.csv,
#          data/admixture.csv
# Outputs: media/diagnostics/<model>/{resid_fitted,scale_location,qq}.png
#          data/results/lmm_diagnostics.csv
#          data/results/lmm_influential_families.csv
# -----------------------------------------------------------------------------

source("R/00_setup.R")
source("R/helpers.R")
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(purrr)
  library(lme4); library(ggplot2)
})

LMM_CTRL <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
diag_dir <- file.path(paths$figures, "diagnostics")
dir.create(diag_dir, showWarnings = FALSE, recursive = TRUE)

theme_d <- theme_classic(base_size = 10)

# ---- Data ---------------------------------------------------------------

adm <- read_csv(file.path(paths$data, "admixture.csv"), show_col_types = FALSE) %>%
  mutate(id = as.character(id))

harvest <- read_csv(file.path(paths$derived, "traits_harvest.csv"), show_col_types = FALSE) %>%
  mutate(id = as.character(id)) %>% left_join(adm, by = c("mom", "id")) %>%
  mutate(origin = factor(origin, levels = c("native", "cultivar")), family = factor(mom))

repdat <- read_csv(file.path(paths$derived, "traits_repeated.csv"), show_col_types = FALSE) %>%
  mutate(id = as.character(id)) %>% left_join(adm, by = c("mom", "id")) %>%
  filter(!is.na(value), value > 0) %>%
  mutate(origin = factor(origin, levels = c("native", "cultivar")),
         family = factor(mom), plant = factor(paste(mom, id, sep = "_")), census = factor(census))

# ---- Diagnostic plots + summary for one fitted model -----------------

diagnose <- function(model, name) {
  d <- tibble(fitted = fitted(model), resid = residuals(model)) %>%
    mutate(std_resid = resid / sd(resid), sqrt_abs = sqrt(abs(std_resid)))

  out_dir <- file.path(diag_dir, name)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  ggsave(file.path(out_dir, "resid_fitted.png"),
    ggplot(d, aes(fitted, resid)) + geom_point(alpha = 0.4, size = 0.9) +
      geom_hline(yintercept = 0, linetype = 2) + geom_smooth(se = FALSE, colour = "firebrick") +
      labs(title = name, x = "Fitted", y = "Residual") + theme_d,
    width = 4, height = 3.2, dpi = 200)

  ggsave(file.path(out_dir, "scale_location.png"),
    ggplot(d, aes(fitted, sqrt_abs)) + geom_point(alpha = 0.4, size = 0.9) +
      geom_smooth(se = FALSE, colour = "firebrick") +
      labs(title = name, x = "Fitted", y = expression(sqrt(abs("std. residual")))) + theme_d,
    width = 4, height = 3.2, dpi = 200)

  qq <- qqnorm(d$std_resid, plot.it = FALSE)
  ggsave(file.path(out_dir, "qq.png"),
    ggplot(tibble(theoretical = qq$x, sample = qq$y), aes(theoretical, sample)) +
      geom_abline(slope = 1, intercept = 0, linetype = 2) + geom_point(alpha = 0.4, size = 0.9) +
      labs(title = name, x = "Theoretical quantiles", y = "Standardised residual") + theme_d,
    width = 4, height = 3.2, dpi = 200)

  sw <- tryCatch(shapiro.test(sample(d$std_resid, min(5000, nrow(d))))$p.value,
                 error = function(e) NA_real_)
  het <- suppressWarnings(cor(d$fitted, d$sqrt_abs, use = "complete.obs"))

  tibble(model = name, n = nrow(d),
         shapiro_p = sw, hetero_cor = het,
         resid_sd = sd(d$resid))
}

# ---- Single-time trait models --------------------------------------

single_specs <- c(height = "height_cm_ln", stem_diam = "stemd_mm_ln",
                  leaf_number = "leaves_n_ln", AGB = "agb_g_ln", BGB = "bgb_g_ln",
                  leaf_biomass = "leaf_biomass_g_ln", stem_biomass = "stem_biomass_g_ln")

single_diag <- imap_dfr(single_specs, function(y, nm) {
  d <- harvest %>% filter(!is.na(.data[[y]]))
  m <- lmer(reformulate(c("days_z", "seed_weight_z", "(1|origin)", "(1|family)"), y),
            data = d, REML = TRUE, control = LMM_CTRL)
  diagnose(m, paste0("single_", nm))
})

# ---- Repeated-measures trait models -------------------------------

repeated_diag <- imap_dfr(c(height = "PH", stem_diam = "SD", leaf_number = "NL"), function(tr, nm) {
  d <- repdat %>% filter(trait == tr)
  m <- lmer(reformulate(c("days_z", "seed_weight_z", "census", "(1|plant)", "(1|family)"), "value_ln"),
            data = d, REML = TRUE, control = LMM_CTRL)
  diagnose(m, paste0("repeated_", nm))
})

lmm_diagnostics <- bind_rows(single_diag, repeated_diag) %>%
  mutate(across(c(shapiro_p, hetero_cor, resid_sd), ~ round(.x, 4)))
write_csv(lmm_diagnostics, file.path(paths$results, "lmm_diagnostics.csv"))

# ---- Per-family influence for the single-time family models -------

infl <- imap_dfr(single_specs, function(y, nm) {
  d  <- harvest %>% filter(!is.na(.data[[y]])) %>% mutate(family = droplevels(family))
  m0 <- lmer(reformulate(c("days_z", "seed_weight_z", "(1|family)"), y),
             data = d, REML = TRUE, control = LMM_CTRL)
  refit_fn <- function(dd) lmer(stats::formula(m0), data = dd, REML = TRUE, control = LMM_CTRL)
  cd <- family_cooks_d(m0, d, refit_fn = refit_fn)
  thr <- 4 / nlevels(d$family)
  tibble(trait = nm, family = names(cd), cooks_d = round(as.numeric(cd), 4),
         above_4_over_n = as.numeric(cd) > thr)

})
write_csv(infl, file.path(paths$results, "lmm_influential_families.csv"))

message("\n--- LMM diagnostics ---")
print(as.data.frame(lmm_diagnostics), row.names = FALSE)
message("\nshapiro_p < 0.05 flags non-normal residuals; |hetero_cor| large flags heteroscedasticity")
message("\n--- influential families (Cook's D > 4/n) ---")
print(as.data.frame(infl), row.names = FALSE)
