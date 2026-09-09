# 05_survival.R
# -----------------------------------------------------------------------------
# Seedling survival to the end of the second growing season, analysed with
# binomial generalized linear mixed models (logit link).
#
# Response : alive_final (1 = alive at the 2019-09-18 census, 0 = dead)
# Fixed    : days_z (scaled days to germinate), seed_weight_z (scaled maternal
#            seed weight, origin-mean imputed); AP_z (scaled Admixture
#            Percentage) in the admixture models only
# Random   : maternal family; origin (native/cultivar) in the origin models
#
# For each data subset the script compares
#   fixed-only  <  + (1|family)  <  + (1|origin/family)
# by AIC / AICc, and for the mixed models estimates latent-scale variance
# components and half-sib heritability with 1000-replicate parametric
# bootstrap confidence intervals. Family conditional modes (BLUPs) from the
# cultivar-subset family model give Fig_BLUPs.
#
# Inputs : data/derived/survival.csv, data/admixture.csv
# Outputs: data/results/survival_aic.csv
#          data/results/survival_quantgen.csv
#          data/results/survival_fixed_effects.csv
#          data/results/survival_lrt.csv
#          data/results/survival_diagnostics.csv
#          data/derived/blups_cultivar_survival.csv   (Fig_BLUPs plotting data)
# -----------------------------------------------------------------------------

source("R/00_setup.R")
source("R/helpers.R")
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(lme4); library(DHARMa)
})

check_packages()
NSIM <- as.integer(Sys.getenv("NSIM", "1000"))
GLMM_CTRL <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))

# ---- 1. Data --------------------------------------------------------------

surv <- read_csv(file.path(paths$derived, "survival.csv"), show_col_types = FALSE) %>%
  mutate(id = as.character(id)) %>%
  left_join(read_csv(file.path(paths$data, "admixture.csv"), show_col_types = FALSE) %>%
              mutate(id = as.character(id)),
            by = c("mom", "id")) %>%
  mutate(
    origin = factor(origin, levels = c("native", "cultivar")),
    family = factor(mom),
    AP_z   = as.numeric(scale(AP))
  )

subsets <- list(
  full     = surv,
  native   = surv %>% filter(origin == "native")   %>% mutate(family = droplevels(family)),
  cultivar = surv %>% filter(origin == "cultivar") %>% mutate(family = droplevels(family)),
  genotyped = surv %>% filter(!is.na(AP))          %>% mutate(family = droplevels(family))
)

message("subset sizes: ",
        paste(sprintf("%s=%d", names(subsets), vapply(subsets, nrow, integer(1))), collapse = "  "))

# ---- 2. Fit the model sets ---------------------------------------------

fit_fixed  <- function(d, extra = NULL)
  glm(reformulate(c("days_z", "seed_weight_z", extra), "alive_final"),
      data = d, family = binomial())

fit_family <- function(d, extra = NULL)
  glmer(reformulate(c("days_z", "seed_weight_z", extra, "(1|family)"), "alive_final"),
        data = d, family = binomial(), control = GLMM_CTRL)

fit_origin <- function(d)
  glmer(alive_final ~ days_z + seed_weight_z + (1 | origin/family),
        data = d, family = binomial(), control = GLMM_CTRL)

models <- list(
  full_fixed         = fit_fixed(subsets$full),
  full_family        = fit_family(subsets$full),
  full_origin_family = fit_origin(subsets$full),
  native_fixed       = fit_fixed(subsets$native),
  native_family      = fit_family(subsets$native),
  cultivar_fixed     = fit_fixed(subsets$cultivar),
  cultivar_family    = fit_family(subsets$cultivar),
  geno_family        = fit_family(subsets$genotyped),
  geno_AP_family     = fit_family(subsets$genotyped, extra = "AP_z")
)

# ---- 3. AIC / AICc comparison ----------------------------------------

aic_sets <- list(
  Total     = c("fixed only" = "full_fixed", "family" = "full_family",
                "origin + family" = "full_origin_family"),
  Native    = c("fixed only" = "native_fixed", "family" = "native_family"),
  Cultivar  = c("fixed only" = "cultivar_fixed", "family" = "cultivar_family"),
  AP        = c("family" = "geno_family", "AP + family" = "geno_AP_family")
)

survival_aic <- imap_dfr(aic_sets, function(mset, subset_name) {
  aic_table(models[mset]) %>%
    mutate(model = names(mset)[match(model, mset)], subset = subset_name, .before = 1)
})
write_csv(survival_aic, file.path(paths$results, "survival_aic.csv"))

# ---- 4. Variance components, heritability, bootstrap CIs ------------

boot_targets <- c(
  Total_family        = "full_family",
  Total_origin_family = "full_origin_family",
  Native_family       = "native_family",
  Cultivar_family     = "cultivar_family",
  AP_family           = "geno_family"
)

survival_quantgen <- imap_dfr(boot_targets, function(mod_name, label) {
  message("bootstrapping ", label, " (", NSIM, " sims) ...")
  bc <- boot_quantgen(models[[mod_name]], resid_var = LATENT_RESID_VAR, nsim = NSIM, seed = SEED)
  bc %>% mutate(set = label, model = mod_name, .before = 1)
})
write_csv(survival_quantgen, file.path(paths$results, "survival_quantgen.csv"))

# ---- 5. Fixed effects (odds ratios) --------------------------------

survival_fixed_effects <- imap_dfr(models, function(m, nm) {
  co <- if (inherits(m, "glmerMod")) summary(m)$coefficients else summary(m)$coefficients
  as_tibble(co, rownames = "term") %>%
    rename(estimate = Estimate, se = `Std. Error`, z = `z value`, p = `Pr(>|z|)`) %>%
    mutate(model = nm, odds_ratio = exp(estimate), .before = 1)
})
write_csv(survival_fixed_effects, file.path(paths$results, "survival_fixed_effects.csv"))

# ---- 6. Likelihood-ratio tests for the random effects -------------

survival_lrt <- bind_rows(
  lrt(models$full_fixed,   models$full_family)        %>% mutate(subset = "Total",    effect = "family"),
  lrt(models$full_family,  models$full_origin_family) %>% mutate(subset = "Total",    effect = "origin"),
  lrt(models$native_fixed, models$native_family)      %>% mutate(subset = "Native",   effect = "family"),
  lrt(models$cultivar_fixed, models$cultivar_family)  %>% mutate(subset = "Cultivar", effect = "family")
) %>% select(subset, effect, chisq, df, p)
write_csv(survival_lrt, file.path(paths$results, "survival_lrt.csv"))

# ---- 7. GLMM diagnostics (DHARMa) ---------------------------------

survival_diagnostics <- imap_dfr(
  list(Total = models$full_family, Cultivar = models$cultivar_family,
       Native = models$native_family),
  function(m, nm) {
    sr <- simulateResiduals(m, n = 1000, seed = SEED)
    tibble(
      subset       = nm,
      uniformity_p = testUniformity(sr, plot = FALSE)$p.value,
      dispersion_p = testDispersion(sr, plot = FALSE)$p.value,
      outlier_p    = testOutliers(sr, plot = FALSE)$p.value
    )
  })
write_csv(survival_diagnostics, file.path(paths$results, "survival_diagnostics.csv"))

# ---- 8. Family BLUPs for the cultivar survival model (Fig_BLUPs) --

blups <- ranef(models$cultivar_family, condVar = TRUE)$family
blup_sd <- sqrt(attr(blups, "postVar")[1, 1, ])
blups_cultivar <- tibble(
  family     = rownames(blups),
  blup_logit = blups[, 1],
  se         = blup_sd
) %>% arrange(blup_logit)
write_csv(blups_cultivar, file.path(paths$derived, "blups_cultivar_survival.csv"))

# ---- 9. Console summary ------------------------------------------

message("\n--- AIC ---")
print(as.data.frame(survival_aic), row.names = FALSE)
message("\n--- variance components / h2 (latent scale) ---")
print(as.data.frame(survival_quantgen %>%
        filter(parameter %in% c("V_origin", "V_family", "h2")) %>%
        mutate(across(c(estimate, ci_low, ci_high), ~ round(.x, 3)))), row.names = FALSE)
message("\n--- LRT ---")
print(as.data.frame(survival_lrt), row.names = FALSE)
message("\n--- diagnostics (p > 0.05 = no detectable violation) ---")
print(as.data.frame(survival_diagnostics), row.names = FALSE)
