# 07_repeated_measures.R
# -----------------------------------------------------------------------------
# Height, stem diameter and leaf number analysed as repeated measures across
# four censuses (start and end of each growing season), with linear mixed
# models on the ln scale.
#
# Response : value_ln at each census
# Fixed    : days_z, seed_weight_z, census (factor)
# Random   : individual plant; maternal family; origin (native/cultivar)
#
# For each trait:
#   * model comparison over
#       fixed-only,  +(1|plant),  +(1|plant)+(1|family),
#       +(1|plant)+(1|origin),  +(1|plant)+(1|origin)+(1|family)
#     by AIC / AICc                                    -> Tab_RepeatedAIC
#   * variance components (origin, family, plant, residual), half-sib h2,
#     bootstrap CIs and ranova for the origin+family and family models,
#     with an outlier / influential-family sensitivity pass
#                                                       -> Tab_RepeatedHeritability
#
# Input  : data/derived/traits_repeated.csv, data/admixture.csv
# Outputs: data/results/repeated_aic.csv
#          data/results/repeated_quantgen.csv
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

TRAITS <- c(height = "PH", stem_diam = "SD", leaf_number = "NL")

# ---- 1. Data -----------------------------------------------------------

rep_dat <- read_csv(file.path(paths$derived, "traits_repeated.csv"), show_col_types = FALSE) %>%
  mutate(id = as.character(id)) %>%
  left_join(read_csv(file.path(paths$data, "admixture.csv"), show_col_types = FALSE) %>%
              mutate(id = as.character(id)),
            by = c("mom", "id")) %>%
  filter(!is.na(value), value > 0) %>%
  mutate(
    origin = factor(origin, levels = c("native", "cultivar")),
    family = factor(mom),
    plant  = factor(paste(mom, id, sep = "_")),
    census = factor(census),
    AP_z   = as.numeric(scale(AP))
  )

# ---- 2. Model set for one trait ------------------------------------

fit_lm   <- function(d, extra = NULL)
  lm(reformulate(c("days_z", "seed_weight_z", "census", extra), "value_ln"), data = d)

fit_lmer <- function(d, rand, extra = NULL)
  lmer(reformulate(c("days_z", "seed_weight_z", "census", extra, rand), "value_ln"),
       data = d, REML = TRUE, control = LMM_CTRL)

trait_models <- function(tr) {
  d <- rep_dat %>% filter(trait == tr)
  list(
    data = d,
    models = list(
      fixed               = fit_lm(d),
      plant               = fit_lmer(d, "(1|plant)"),
      plant_family        = fit_lmer(d, c("(1|plant)", "(1|family)")),
      plant_origin        = fit_lmer(d, c("(1|plant)", "(1|origin)")),
      plant_origin_family = fit_lmer(d, c("(1|plant)", "(1|origin)", "(1|family)"))
    )
  )
}

fitted_all <- map(TRAITS, trait_models)

# ---- 3. AIC / AICc  -> Tab_RepeatedAIC ---------------------------

model_labels <- c(
  fixed = "fixed only", plant = "plant", plant_family = "plant + family",
  plant_origin = "plant + origin", plant_origin_family = "plant + family + origin"
)

repeated_aic <- imap_dfr(TRAITS, function(tr, trait_name) {
  aic_table(fitted_all[[trait_name]]$models) %>%
    mutate(model = model_labels[model], trait = trait_name, .before = 1)
})
write_csv(repeated_aic, file.path(paths$results, "repeated_aic.csv"))

# ---- 4. Variance / heritability with sensitivity  -> Tab_RepeatedHeritability ----

refit <- function(model, data) lmer(stats::formula(model), data = droplevels(data),
                                    REML = TRUE, control = LMM_CTRL)

vc_row <- function(model, label, sensitivity, removed) {
  bc <- boot_quantgen(model, nsim = NSIM, seed = SEED)
  rv <- tryCatch(as.data.frame(lmerTest::ranova(model)), error = function(e) NULL)
  getp <- function(term) if (!is.null(rv) && term %in% rownames(rv)) rv[term, "Pr(>Chisq)"] else NA_real_
  bc %>%
    filter(parameter %in% c("V_origin", "V_family", "V_plant", "V_residual",
                            "P_origin", "P_family", "P_plant", "P_residual", "h2")) %>%
    select(parameter, estimate, ci_low, ci_high) %>%
    pivot_wider(names_from = parameter, values_from = c(estimate, ci_low, ci_high)) %>%
    mutate(model = label, sensitivity = sensitivity, removed = removed,
           p_origin = getp("(1 | origin)"), p_family = getp("(1 | family)"),
           p_plant  = getp("(1 | plant)"), .before = 1)
}

sensitivity_rows <- function(trait_name, mod_name, label) {
  d  <- fitted_all[[trait_name]]$data
  m0 <- fitted_all[[trait_name]]$models[[mod_name]]

  out <- vc_row(m0, label, "original", NA_character_)

  # individual measurement outliers: |studentised residual| > 3
  bad_idx <- which(abs(as.numeric(scale(residuals(m0)))) > 3)
  if (length(bad_idx)) {
    m_o <- refit(m0, d[-bad_idx, ])
    out <- bind_rows(out, vc_row(m_o, label, "obs_removed", as.character(length(bad_idx))))
  }

  mi <- most_influential_family(m0, d, refit_fn = function(dd) refit(m0, dd))
  if (!is.null(mi)) {
    m_f <- refit(m0, d %>% filter(as.character(family) != mi$family))
    out <- bind_rows(out, vc_row(m_f, label, "family_removed", mi$family))
  }
  out
}

repeated_quantgen <- imap_dfr(TRAITS, function(tr, trait_name) {
  message("quantgen: ", trait_name)
  bind_rows(
    sensitivity_rows(trait_name, "plant_origin_family", "origin + family"),
    sensitivity_rows(trait_name, "plant_family",        "family")
  ) %>% mutate(trait = trait_name, .before = 1)
})
write_csv(repeated_quantgen, file.path(paths$results, "repeated_quantgen.csv"))

# ---- 5. Console summary ------------------------------------------

message("\n--- Tab_RepeatedAIC ---")
print(as.data.frame(repeated_aic %>% select(trait, model, AICc, dAICc, wAICc) %>%
        mutate(across(c(AICc, dAICc, wAICc), ~ round(.x, 2)))), row.names = FALSE)

message("\n--- h2 (original) ---")
print(as.data.frame(repeated_quantgen %>% filter(sensitivity == "original") %>%
        transmute(trait, model,
                  h2 = sprintf("%.3f [%.3f, %.3f]", estimate_h2, ci_low_h2, ci_high_h2),
                  p_family = round(p_family, 4), p_origin = round(p_origin, 4))),
      row.names = FALSE)
