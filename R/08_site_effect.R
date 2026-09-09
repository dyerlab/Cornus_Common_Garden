# 08_site_effect.R
# -----------------------------------------------------------------------------
# Native maternal trees were collected at two locations: Pocahontas State Park
# and the VCU Rice Rivers Center. This script tests whether collection site
# explains variation in survival or any growth trait among native-origin
# seedlings, to justify pooling the two into a single "native" class.
#
# For each response, on native-origin seedlings only:
#   reduced :  ~ days_z + seed_weight_z + (1|family)
#   full    :  reduced + site
# tested by a likelihood-ratio test of the site term.
#
# Inputs : data/derived/traits_harvest.csv, data/derived/survival.csv
# Output : data/results/site_effect.csv   (Tab_SiteEffect)
# -----------------------------------------------------------------------------

source("R/00_setup.R")
source("R/helpers.R")
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(purrr); library(lme4); library(lmerTest)
})

LMM_CTRL  <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
GLMM_CTRL <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))

harvest <- read_csv(file.path(paths$derived, "traits_harvest.csv"), show_col_types = FALSE)
surv    <- read_csv(file.path(paths$derived, "survival.csv"),       show_col_types = FALSE)

TRAITS <- c(
  height       = "height_cm_ln",  stem_diam    = "stemd_mm_ln",
  leaf_number  = "leaves_n_ln",   AGB          = "agb_g_ln",
  BGB          = "bgb_g_ln",      leaf_biomass = "leaf_biomass_g_ln",
  stem_biomass = "stem_biomass_g_ln"
)

native_traits <- harvest %>%
  filter(origin == "native", !is.na(site)) %>%
  mutate(family = factor(mom), site = factor(site))

native_surv <- surv %>%
  filter(origin == "native", !is.na(site)) %>%
  mutate(family = factor(mom), site = factor(site))

message("native seedlings with a site: ", nrow(native_surv),
        "  (", paste(names(table(native_surv$site)), table(native_surv$site),
                     sep = " = ", collapse = ", "), ")")
message("native families per site: ",
        paste(names(table(distinct(native_surv, mom, site)$site)),
              table(distinct(native_surv, mom, site)$site), sep = " = ", collapse = ", "))

# ---- growth traits (LMM) ----------------------------------------------

lmm_site <- imap_dfr(TRAITS, function(y, trait_name) {
  d <- native_traits %>% filter(!is.na(.data[[y]]))
  reduced <- lmer(reformulate(c("days_z", "seed_weight_z", "(1|family)"), y),
                  data = d, REML = FALSE, control = LMM_CTRL)
  full    <- lmer(reformulate(c("days_z", "seed_weight_z", "site", "(1|family)"), y),
                  data = d, REML = FALSE, control = LMM_CTRL)
  co <- summary(full)$coefficients
  site_row <- grep("^site", rownames(co))
  test <- lrt(reduced, full)
  tibble(
    response = trait_name, n = nrow(d),
    site_estimate = co[site_row, "Estimate"],
    site_se       = co[site_row, "Std. Error"],
    chisq = test$chisq, df = test$df, p = 2 * test$p   # two-sided: undo the RE halving
  )
})

# ---- survival (GLMM) ------------------------------------------------

surv_reduced <- glmer(alive_final ~ days_z + seed_weight_z + (1 | family),
                      data = native_surv, family = binomial(), control = GLMM_CTRL)
surv_full    <- glmer(alive_final ~ days_z + seed_weight_z + site + (1 | family),
                      data = native_surv, family = binomial(), control = GLMM_CTRL)
co <- summary(surv_full)$coefficients
test <- lrt(surv_reduced, surv_full)

site_effect <- bind_rows(
  tibble(
    response = "survival", n = nrow(native_surv),
    site_estimate = co[grep("^site", rownames(co)), "Estimate"],
    site_se       = co[grep("^site", rownames(co)), "Std. Error"],
    chisq = test$chisq, df = test$df, p = 2 * test$p
  ),
  lmm_site
) %>%
  mutate(across(c(site_estimate, site_se, chisq), ~ round(.x, 3)),
         p = round(p, 3))

write_csv(site_effect, file.path(paths$results, "site_effect.csv"))

message("\n--- Tab_SiteEffect (Pocahontas vs Rice Rivers, native seedlings) ---")
print(as.data.frame(site_effect), row.names = FALSE)
message("\nall p > 0.05 supports pooling the two native collection sites")
