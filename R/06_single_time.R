# 06_single_time.R
# -----------------------------------------------------------------------------
# End-of-experiment growth traits, analysed one measurement per seedling with
# linear mixed models on the ln scale.
#
# Traits : height, stem diameter, leaf number, above-ground biomass (AGB),
#          below-ground biomass (BGB), leaf biomass, stem biomass
# Fixed  : days_z, seed_weight_z; AP_z in the admixture models
# Random : maternal family; origin (native/cultivar) in the origin models
#
# For each trait:
#   * model comparison  fixed-only < +(1|family) < +(1|origin)+(1|family)
#     by AIC/AICc                                     -> Tab_GrowthAIC
#   * an admixture arm   +(1|family) vs +AP_z+(1|family) on genotyped seedlings
#   * variance components, half-sib h2, parametric-bootstrap CIs, ranova
#     for the origin+family and family models              -> Tab_GrowthHeritability
#   * the same within each origin (family model only)        -> Tab_OriginHeritability
#   * a sensitivity pass: refit after dropping individual outliers
#     (|studentised residual| > 3) and after dropping the single most
#     influential family (largest leave-one-out Cook's distance)
#   * fixed-effect tests for the covariates                  -> Tab_CovariateEffects
#
# Inputs : data/derived/traits_harvest.csv, data/admixture.csv
# Outputs: data/results/growth_aic.csv
#          data/results/growth_quantgen.csv          (full + origin models)
#          data/results/growth_origin_quantgen.csv   (within-origin models)
#          data/results/growth_covariates.csv
# -----------------------------------------------------------------------------

source("R/00_setup.R")
source("R/helpers.R")
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(purrr)
  library(lme4); library(lmerTest)
})

check_packages()
NSIM     <- as.integer(Sys.getenv("NSIM", "1000"))
LMM_CTRL <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))

TRAITS <- c(
  height       = "height_cm_ln",
  stem_diam    = "stemd_mm_ln",
  leaf_number  = "leaves_n_ln",
  AGB          = "agb_g_ln",
  BGB          = "bgb_g_ln",
  leaf_biomass = "leaf_biomass_g_ln",
  stem_biomass = "stem_biomass_g_ln"
)

# ---- 1. Data ------------------------------------------------------------

dat <- read_csv(file.path(paths$derived, "traits_harvest.csv"), show_col_types = FALSE) %>%
  mutate(id = as.character(id)) %>%
  left_join(read_csv(file.path(paths$data, "admixture.csv"), show_col_types = FALSE) %>%
              mutate(id = as.character(id)),
            by = c("mom", "id")) %>%
  mutate(
    origin = factor(origin, levels = c("native", "cultivar")),
    family = factor(mom),
    AP_z   = as.numeric(scale(AP))
  )

# ---- 2. Model fitting for one trait ----------------------------------

fit_lm   <- function(d, y, extra = NULL)
  lm(reformulate(c("days_z", "seed_weight_z", extra), y), data = d)

fit_lmer <- function(d, y, rand, extra = NULL)
  lmer(reformulate(c("days_z", "seed_weight_z", extra, rand), y),
       data = d, REML = TRUE, control = LMM_CTRL)

trait_models <- function(y) {
  d_all  <- dat %>% filter(!is.na(.data[[y]]))
  d_gen  <- d_all %>% filter(!is.na(AP)) %>% mutate(family = droplevels(family))
  d_nat  <- d_all %>% filter(origin == "native")   %>% mutate(family = droplevels(family))
  d_cul  <- d_all %>% filter(origin == "cultivar") %>% mutate(family = droplevels(family))
  list(
    data = list(all = d_all, gen = d_gen, native = d_nat, cultivar = d_cul),
    models = list(
      Total_fixed         = fit_lm(d_all, y),
      Total_family        = fit_lmer(d_all, y, "(1|family)"),
      Total_origin_family = fit_lmer(d_all, y, c("(1|origin)", "(1|family)")),
      AP_family           = fit_lmer(d_gen, y, "(1|family)"),
      AP_AP_family        = fit_lmer(d_gen, y, "(1|family)", extra = "AP_z"),
      Native_fixed        = fit_lm(d_nat, y),
      Native_family       = fit_lmer(d_nat, y, "(1|family)"),
      Cultivar_fixed      = fit_lm(d_cul, y),
      Cultivar_family     = fit_lmer(d_cul, y, "(1|family)")
    )
  )
}

fitted_all <- map(TRAITS, trait_models)

# ---- 3. AIC / AICc  -> Tab_GrowthAIC --------------------------------

aic_sets <- list(
  Total    = c("fixed only" = "Total_fixed", "family" = "Total_family",
               "origin + family" = "Total_origin_family"),
  AP       = c("family" = "AP_family", "AP + family" = "AP_AP_family"),
  Native   = c("fixed only" = "Native_fixed", "family" = "Native_family"),
  Cultivar = c("fixed only" = "Cultivar_fixed", "family" = "Cultivar_family")
)

growth_aic <- imap_dfr(TRAITS, function(y, trait_name) {
  mods <- fitted_all[[trait_name]]$models
  imap_dfr(aic_sets, function(mset, subset_name) {
    aic_table(mods[mset]) %>%
      mutate(model = names(mset)[match(model, mset)],
             trait = trait_name, subset = subset_name, .before = 1)
  })
})
write_csv(growth_aic, file.path(paths$results, "growth_aic.csv"))

# ---- 4. Variance / heritability with sensitivity  -> Tab_GrowthHeritability ----

# one variance-component row set for a fitted model, with bootstrap CIs
vc_row <- function(model, label, sensitivity, removed) {
  qg <- quantgen_summary(model)                       # Gaussian: residual from model
  bc <- boot_quantgen(model, nsim = NSIM, seed = SEED)
  rv <- tryCatch(as.data.frame(lmerTest::ranova(model)), error = function(e) NULL)
  p_origin <- if (!is.null(rv) && "(1 | origin)" %in% rownames(rv)) rv["(1 | origin)", "Pr(>Chisq)"] else NA_real_
  p_family <- if (!is.null(rv) && "(1 | family)" %in% rownames(rv)) rv["(1 | family)", "Pr(>Chisq)"] else NA_real_
  bc %>%
    filter(parameter %in% c("V_origin", "V_family", "V_residual",
                            "P_origin", "P_family", "P_residual", "h2")) %>%
    select(parameter, estimate, ci_low, ci_high) %>%
    pivot_wider(names_from = parameter, values_from = c(estimate, ci_low, ci_high)) %>%
    mutate(model = label, sensitivity = sensitivity, removed = removed,
           p_origin = p_origin, p_family = p_family, .before = 1)
}

refit <- function(model, data) lmer(stats::formula(model), data = droplevels(data),
                                    REML = TRUE, control = LMM_CTRL)

sensitivity_rows <- function(y, mod_name, label, data_key) {
  d   <- fitted_all[[y]]$data[[data_key]]
  m0  <- fitted_all[[y]]$models[[mod_name]]

  # original
  out <- vc_row(m0, label, "original", NA_character_)

  # drop individual outliers
  bad_obs <- outlier_obs(m0, d, "id")
  if (length(bad_obs)) {
    m_o <- refit(m0, d %>% filter(!id %in% bad_obs))
    out <- bind_rows(out, vc_row(m_o, label, "obs_removed", as.character(length(bad_obs))))
  }

  # drop the single most influential family (largest leave-one-out Cook's D)
  mi <- most_influential_family(m0, d, refit_fn = function(dd) refit(m0, dd))
  if (!is.null(mi)) {
    m_f <- refit(m0, d %>% filter(as.character(family) != mi$family))
    out <- bind_rows(out, vc_row(m_f, label, "family_removed", mi$family))
  }
  out
}

growth_quantgen <- imap_dfr(TRAITS, function(y, trait_name) {
  message("quantgen: ", trait_name)
  bind_rows(
    sensitivity_rows(trait_name, "Total_origin_family", "origin + family", "all"),
    sensitivity_rows(trait_name, "Total_family",        "family",          "all")
  ) %>% mutate(trait = trait_name, .before = 1)
})
write_csv(growth_quantgen, file.path(paths$results, "growth_quantgen.csv"))

# ---- 5. Within-origin models  -> Tab_OriginHeritability ------------

growth_origin_quantgen <- imap_dfr(TRAITS, function(y, trait_name) {
  message("within-origin quantgen: ", trait_name)
  bind_rows(
    sensitivity_rows(trait_name, "Native_family",   "family", "native")   %>% mutate(subset = "Native"),
    sensitivity_rows(trait_name, "Cultivar_family", "family", "cultivar") %>% mutate(subset = "Cultivar")
  ) %>% mutate(trait = trait_name, .before = 1)
})
write_csv(growth_origin_quantgen, file.path(paths$results, "growth_origin_quantgen.csv"))

# ---- 6. Covariate fixed effects  -> Tab_CovariateEffects ----------

growth_covariates <- imap_dfr(TRAITS, function(y, trait_name) {
  imap_dfr(list(Total = "Total_family", Native = "Native_family", Cultivar = "Cultivar_family"),
           function(mod_name, subset_name) {
    m <- fitted_all[[trait_name]]$models[[mod_name]]
    co <- as.data.frame(summary(m)$coefficients)
    tibble(
      trait = trait_name, subset = subset_name, term = rownames(co),
      estimate = co[, "Estimate"], se = co[, "Std. Error"],
      df = co[, "df"], t = co[, "t value"], p = co[, "Pr(>|t|)"]
    ) %>% filter(term != "(Intercept)")
  })
})
write_csv(growth_covariates, file.path(paths$results, "growth_covariates.csv"))

# ---- 6b. Germination timing (GD)  -------------------------------
# Does origin predict the timing of germination itself? Response is
# developmental time since germination (dev_time; higher = earlier). No
# germination-timing covariate here because it is the response.

gd <- dat %>% filter(!is.na(dev_time))
gd_models <- list(
  fixed         = lm(dev_time ~ seed_weight_z, data = gd),
  family        = lmer(dev_time ~ seed_weight_z + (1 | family), data = gd, REML = TRUE, control = LMM_CTRL),
  origin_family = lmer(dev_time ~ seed_weight_z + (1 | origin) + (1 | family), data = gd, REML = TRUE, control = LMM_CTRL)
)
gd_aic <- aic_table(gd_models) %>%
  mutate(model = recode(model, fixed = "fixed only", family = "family",
                        origin_family = "origin + family"))
gd_qg <- boot_quantgen(gd_models$origin_family, nsim = NSIM, seed = SEED) %>%
  filter(parameter %in% c("V_origin", "V_family", "V_residual", "h2"))
gd_ranova <- as.data.frame(lmerTest::ranova(gd_models$origin_family))

germination_timing <- bind_rows(
  gd_aic %>% transmute(what = "AIC", model, value = sprintf("AICc=%.1f dAICc=%.2f wAICc=%.2f", AICc, dAICc, wAICc)),
  gd_qg  %>% transmute(what = "quantgen", model = parameter, value = sprintf("%.3f [%.3f, %.3f]", estimate, ci_low, ci_high)),
  tibble(what = "ranova_p",
         model = rownames(gd_ranova)[-1],
         value = sprintf("%.4f", gd_ranova[["Pr(>Chisq)"]][-1]))
)
write_csv(germination_timing, file.path(paths$results, "germination_timing.csv"))
message("\n--- germination timing (GD) ---")
print(as.data.frame(germination_timing), row.names = FALSE)

# ---- 6c. Origin conditional modes (BLUPs) from the origin + family models ----
# Supports the Discussion point that model-predicted growth is higher for
# cultivar-origin than native-origin seedlings.

origin_blups <- imap_dfr(TRAITS, function(y, trait_name) {
  re <- lme4::ranef(fitted_all[[trait_name]]$models[["Total_origin_family"]])$origin
  tibble(trait = trait_name,
         native   = re["native", 1],
         cultivar = re["cultivar", 1],
         cultivar_minus_native = re["cultivar", 1] - re["native", 1])
})
write_csv(origin_blups, file.path(paths$results, "origin_blups.csv"))
message("\n--- origin conditional modes (ln scale) ---")
print(as.data.frame(origin_blups %>% mutate(across(where(is.numeric), ~ round(.x, 4)))), row.names = FALSE)

# ---- 7. Console summary -----------------------------------------

message("\n--- Tab_GrowthAIC (top model per trait/subset) ---")
print(as.data.frame(growth_aic %>% group_by(trait, subset) %>% slice_min(AICc, n = 1) %>%
        ungroup() %>% select(trait, subset, model, wAICc)), row.names = FALSE)

message("\n--- h2 (origin+family and family, original) ---")
print(as.data.frame(growth_quantgen %>% filter(sensitivity == "original") %>%
        transmute(trait, model,
                  h2 = sprintf("%.3f [%.3f, %.3f]", estimate_h2, ci_low_h2, ci_high_h2),
                  p_family = round(p_family, 4))), row.names = FALSE)
