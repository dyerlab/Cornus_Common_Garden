# 12_within_family_variance.R
# -----------------------------------------------------------------------------
# Post-hoc (Discussion) analysis: is the *spread* of offspring phenotypes
# within a maternal family larger for cultivar-origin arrays than for
# native-origin arrays?
#
# Motivation (not an a priori hypothesis): in the main analysis the maternal
# family random effect improved model fit for the cultivar subset but not the
# native subset for growth traits and for survival (Tab_OriginHeritability,
# Tab_Survival). Every other test in the manuscript compares *means*. This
# script compares within-family *variance*, as a lead-in to the follow-up
# study of pollen-pool diversity and structure along the urban -> forest
# gradient: urban maternal environments span a wide range of light, water and
# nutrient provisioning and (potentially) a more heterogeneous pollen pool,
# whereas the second-growth forest understory is more buffered and homogeneous.
#
# Method
#   For each harvest trait (ln scale):
#     1. Remove the origin mean, each mother's mean, and the germination-timing
#        covariate:  y_ln ~ days_z + mom  (mom as a fixed factor), fit
#        SEPARATELY within each origin class. Keep the residuals e_ij. Because
#        the fit is per origin and mom is a fixed factor, e_ij carries no
#        origin- or family-level mean signal, only within-family spread.
#     2. Per mother summarise that spread:
#          s2_i  = unbiased variance of e_ij   (denominator n_i - 1)
#          mad_i = median absolute deviation of e_ij from the family median
#     3. Compare cultivar vs native:
#        (A) family-level Levene / Brown-Forsythe -- the MOTHER is the unit of
#            replication (n = up to 20 cultivar, 16 native). Mann-Whitney and a
#            family-label permutation test on log(s2_i) and on mad_i.
#        (B) pooled-residual variance ratio  R = Var(e|cultivar) / Var(e|native),
#            with the null from permuting origin across mothers (sib clusters
#            kept intact) -- uses every residual but respects the clustering.
#        (C) individual-level Levene test, shown ONLY for contrast: it treats
#            sibs as independent, so its p-value is anti-conservative. A
#            mother-clustered version, lmer(|e| ~ origin + (1|mom)), is the
#            honest individual-level analogue.
#   Robustness
#     * restrict to mothers with n_i >= 8 (small native arrays lose more
#       within-family variance to the n_i - 1 correction, which would bias
#       toward the hypothesis);
#     * rarefy every retained family to a common offspring count;
#     * check for a mother-level mean-variance relationship.
#
# Inputs : data/derived/traits_harvest.csv
# Outputs: data/results/within_family_variance.csv          (one row per trait/test)
#          data/results/within_family_variance_family.csv    (per-mother dispersion)
#          media/fig_within_family_variance.png / .pdf
# -----------------------------------------------------------------------------

source("R/00_setup.R")
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(purrr); library(tibble)
  library(lme4); library(lmerTest); library(nlme); library(ggplot2)
})
check_packages()

NPERM <- as.integer(Sys.getenv("NPERM", "9999"))
set.seed(SEED)

TRAITS <- c(
  height       = "height_cm_ln",
  stem_diam    = "stemd_mm_ln",
  leaf_number  = "leaves_n_ln",
  AGB          = "agb_g_ln",
  BGB          = "bgb_g_ln",
  leaf_biomass = "leaf_biomass_g_ln",
  stem_biomass = "stem_biomass_g_ln"
)

MIN_N_ROBUST <- 8     # stricter mother cut for the robustness pass
RAREFY_K     <- 6     # common offspring count for the rarefaction pass
RAREFY_B     <- 999

# ---- 1. Data ----------------------------------------------------------------

dat <- read_csv(file.path(paths$derived, "traits_harvest.csv"), show_col_types = FALSE) %>%
  mutate(origin = factor(origin, levels = c("native", "cultivar")),
         mom    = factor(mom))

# ---- 2. Within-family residuals -------------------------------------------
# y_ln ~ days_z + mom, fit within each origin class. Mothers with < 2
# usable observations for a trait contribute no within-family residual.

family_residuals <- function(d, y) {
  d <- d %>%
    transmute(origin, mom, yval = .data[[y]], days_z) %>%
    filter(!is.na(yval), !is.na(days_z)) %>%
    droplevels()
  d %>%
    group_split(origin) %>%
    map_dfr(function(dd) {
      dd <- dd %>% group_by(mom) %>% filter(n() >= 2) %>% ungroup() %>% droplevels()
      if (nlevels(dd$mom) < 2) return(NULL)
      dd$e <- residuals(lm(yval ~ days_z + mom, data = dd))
      dd
    })
}

# ---- 3. Per-mother dispersion -------------------------------------------

family_dispersion <- function(res) {
  res %>%
    group_by(origin, mom) %>%
    summarise(
      n          = n(),
      mom_mean_y = mean(yval),
      s2         = var(e),
      log_s2     = log(var(e)),
      mad_dev    = median(abs(e - median(e))),
      .groups = "drop"
    ) %>%
    filter(is.finite(log_s2), s2 > 0)
}

# ---- 4. Tests -----------------------------------------------------------

# (A) family-level: mother is the unit
famlevel_test <- function(fd, value, nperm = NPERM) {
  x_c <- fd[[value]][fd$origin == "cultivar"]
  x_n <- fd[[value]][fd$origin == "native"]
  if (length(x_c) < 2 || length(x_n) < 2)
    return(tibble(test = paste0("family_", value), ratio = NA_real_,
                  mw_p = NA_real_, perm_p_greater = NA_real_, perm_p_two = NA_real_,
                  n_cultivar_moms = length(x_c), n_native_moms = length(x_n)))
  obs   <- mean(x_c) - mean(x_n)
  ratio <- if (value == "log_s2") exp(obs) else mean(x_c) / mean(x_n)   # geometric-mean variance ratio for log_s2
  mw    <- suppressWarnings(wilcox.test(x_c, x_n, alternative = "greater"))
  pool  <- c(x_c, x_n)
  g0    <- rep(c("cultivar", "native"), c(length(x_c), length(x_n)))
  perm  <- replicate(nperm, {
    gg <- sample(g0)
    mean(pool[gg == "cultivar"]) - mean(pool[gg == "native"])
  })
  tibble(
    test            = paste0("family_", value),
    ratio           = ratio,
    mw_p            = mw$p.value,
    perm_p_greater  = (1 + sum(perm >= obs))       / (nperm + 1),
    perm_p_two      = (1 + sum(abs(perm) >= abs(obs))) / (nperm + 1),
    n_cultivar_moms = length(x_c),
    n_native_moms   = length(x_n)
  )
}

# (B) pooled-residual variance ratio, null = permute origin across mothers
pooled_ratio_test <- function(res, nperm = NPERM) {
  mo      <- res %>% distinct(mom, origin)
  origins <- as.character(mo$origin)
  names(origins) <- as.character(mo$mom)
  ratio_for <- function(lbl) {
    o <- lbl[as.character(res$mom)]
    var(res$e[o == "cultivar"]) / var(res$e[o == "native"])
  }
  obs  <- ratio_for(origins)
  perm <- replicate(nperm, ratio_for(setNames(sample(origins), names(origins))))
  tibble(
    test           = "pooled_variance_ratio",
    ratio          = obs,
    mw_p           = NA_real_,
    perm_p_greater = (1 + sum(perm >= obs)) / (nperm + 1),
    perm_p_two     = (1 + sum(abs(log(perm)) >= abs(log(obs)))) / (nperm + 1),
    n_cultivar_moms = sum(origins == "cultivar"),
    n_native_moms   = sum(origins == "native")
  )
}

# (C) individual-level Levene (anti-conservative) + mother-clustered analogue
levene_individual <- function(res) {
  r <- res %>% group_by(origin, mom) %>%
    mutate(absdev = abs(e - median(e))) %>% ungroup()
  p_classic <- anova(lm(absdev ~ origin, data = r))[["Pr(>F)"]][1]
  mt <- lmerTest::lmer(absdev ~ origin + (1 | mom), data = r)
  co <- summary(mt)$coefficients
  tibble(
    test            = c("levene_individual_naive", "levene_mother_clustered"),
    ratio           = NA_real_,
    mw_p            = NA_real_,
    perm_p_greater  = c(p_classic / 2, co["origincultivar", "Pr(>|t|)"] / 2),
    perm_p_two      = c(p_classic,     co["origincultivar", "Pr(>|t|)"]),
    n_cultivar_moms = NA_integer_,
    n_native_moms   = NA_integer_
  )
}

# (D) location-scale model: the principled arbiter. lme(y ~ days_z + origin,
# random = ~1|mom) with a residual variance that is allowed to differ by origin
# (varIdent). The likelihood-ratio test of that extra parameter is a direct
# test of "within-family variance differs by origin", using every observation,
# weighting families correctly, and holding the among-family structure fixed.
location_scale_test <- function(res) {
  d <- res %>% mutate(mom = factor(mom), origin = factor(origin, levels = c("native", "cultivar")))
  ctrl <- nlme::lmeControl(opt = "optim", maxIter = 300, msMaxIter = 300, returnObject = TRUE)
  out <- tryCatch({
    m0 <- nlme::lme(yval ~ days_z + origin, random = ~ 1 | mom, data = d,
                    method = "ML", control = ctrl)
    m1 <- update(m0, weights = nlme::varIdent(form = ~ 1 | origin))
    lr <- anova(m0, m1)
    sd_ratio <- as.numeric(coef(m1$modelStruct$varStruct, unconstrained = FALSE)["cultivar"])
    p_two <- lr$`p-value`[2]
    tibble(test = "location_scale_lme",
           ratio = sd_ratio^2,                       # cultivar / native residual variance
           mw_p = NA_real_,
           perm_p_greater = if (sd_ratio > 1) p_two / 2 else 1 - p_two / 2,
           perm_p_two = p_two,
           n_cultivar_moms = NA_integer_, n_native_moms = NA_integer_)
  }, error = function(e) tibble(test = "location_scale_lme", ratio = NA_real_,
           mw_p = NA_real_, perm_p_greater = NA_real_, perm_p_two = NA_real_,
           n_cultivar_moms = NA_integer_, n_native_moms = NA_integer_))
  out
}

# rarefaction: subsample every retained family to k, integrate a family-label
# permutation test over the subsampling
rarefy_test <- function(res, k = RAREFY_K, B = RAREFY_B) {
  mo      <- res %>% distinct(mom, origin)
  origins <- as.character(mo$origin); names(origins) <- as.character(mo$mom)
  moms_ok <- res %>% count(mom) %>% filter(n >= k) %>% pull(mom) %>% as.character()
  if (length(unique(origins[moms_ok])) < 2 ||
      min(table(origins[moms_ok])) < 2)
    return(tibble(test = "rarefied_variance_ratio", ratio = NA_real_,
                  mw_p = NA_real_, perm_p_greater = NA_real_, perm_p_two = NA_real_,
                  n_cultivar_moms = sum(origins[moms_ok] == "cultivar"),
                  n_native_moms   = sum(origins[moms_ok] == "native")))
  r <- res %>% filter(as.character(mom) %in% moms_ok)
  obs <- nul <- numeric(B)
  for (b in seq_len(B)) {
    sub <- r %>% group_by(mom) %>% slice_sample(n = k) %>% ungroup()
    o   <- origins[as.character(sub$mom)]
    obs[b] <- var(sub$e[o == "cultivar"]) / var(sub$e[o == "native"])
    po  <- setNames(sample(origins[moms_ok]), moms_ok)[as.character(sub$mom)]
    nul[b] <- var(sub$e[po == "cultivar"]) / var(sub$e[po == "native"])
  }
  tibble(
    test            = "rarefied_variance_ratio",
    ratio           = median(obs),
    mw_p            = NA_real_,
    perm_p_greater  = (1 + sum(nul >= median(obs))) / (B + 1),
    perm_p_two      = (1 + sum(abs(log(nul)) >= abs(log(median(obs))))) / (B + 1),
    n_cultivar_moms = sum(origins[moms_ok] == "cultivar"),
    n_native_moms   = sum(origins[moms_ok] == "native")
  )
}

# ---- 5. Run for every trait ------------------------------------------

disp_all   <- list()
result_all <- list()

for (nm in names(TRAITS)) {
  y   <- TRAITS[[nm]]
  res <- family_residuals(dat, y)
  fd  <- family_dispersion(res)
  disp_all[[nm]] <- fd %>% mutate(trait = nm, .before = 1)

  fd_rob  <- fd  %>% filter(n >= MIN_N_ROBUST)
  res_rob <- res %>% semi_join(fd_rob, by = "mom")

  # mother-level diagnostics: mean-variance coupling, and (crucially) whether
  # the per-family variance estimate scales with family size -- a strong n vs
  # log(s2) correlation flags the small-native-family artefact.
  diag <- fd %>% group_by(origin) %>%
    summarise(mv_rho  = suppressWarnings(cor(mom_mean_y, log_s2, method = "spearman")),
              n_rho   = suppressWarnings(cor(n, log_s2, method = "spearman")),
              .groups = "drop") %>%
    tidyr::pivot_wider(names_from = origin, values_from = c(mv_rho, n_rho))

  rows <- bind_rows(
    famlevel_test(fd, "log_s2")                      %>% mutate(pass = "main"),
    famlevel_test(fd, "mad_dev")                     %>% mutate(pass = "main"),
    pooled_ratio_test(res)                           %>% mutate(pass = "main"),
    location_scale_test(res)                         %>% mutate(pass = "main"),
    levene_individual(res)                           %>% mutate(pass = "main"),
    rarefy_test(res)                                 %>% mutate(pass = "main"),
    famlevel_test(fd_rob, "log_s2")                  %>% mutate(pass = "n>=8"),
    pooled_ratio_test(res_rob)                       %>% mutate(pass = "n>=8"),
    location_scale_test(res_rob)                     %>% mutate(pass = "n>=8")
  ) %>%
    mutate(trait = nm, .before = 1) %>%
    bind_cols(diag[rep(1, nrow(.)), ])

  result_all[[nm]] <- rows
}

family_dispersion_tbl <- bind_rows(disp_all)
within_family_variance <- bind_rows(result_all) %>%
  relocate(pass, .after = trait)

write_csv(family_dispersion_tbl, file.path(paths$results, "within_family_variance_family.csv"))
write_csv(within_family_variance, file.path(paths$results, "within_family_variance.csv"))

# ---- 6. Cross-trait summary -----------------------------------------
# directional consistency of the pooled variance ratio across the 7 traits

pooled <- within_family_variance %>% filter(test == "pooled_variance_ratio", pass == "main")
n_gt1  <- sum(pooled$ratio > 1, na.rm = TRUE)
sign_p <- binom.test(n_gt1, sum(is.finite(pooled$ratio)), alternative = "greater")$p.value

message("\n=== Within-family variance: cultivar vs native ===")
message(sprintf("Pooled variance ratio (cultivar / native) > 1 in %d of %d traits (sign test p = %.3f)",
                n_gt1, sum(is.finite(pooled$ratio)), sign_p))
message("\n-- variance ratio (cultivar / native) and one-sided p, per trait --")
print(as.data.frame(
  within_family_variance %>%
    filter(pass == "main", test %in% c("family_log_s2", "family_mad_dev",
                                       "pooled_variance_ratio", "location_scale_lme",
                                       "rarefied_variance_ratio",
                                       "levene_mother_clustered")) %>%
    transmute(trait, test,
              ratio = round(ratio, 2),
              p_greater = round(perm_p_greater, 3))
), row.names = FALSE)

message("\n-- robustness (mothers with n >= 8) --")
print(as.data.frame(
  within_family_variance %>% filter(pass == "n>=8") %>%
    transmute(trait, test, ratio = round(ratio, 2),
              p_greater = round(perm_p_greater, 3),
              n_cultivar_moms, n_native_moms)
), row.names = FALSE)

message("\n-- diagnostics: Spearman rho of per-family log(s2) with family mean")
message("   (mean-variance coupling) and with family size n (small-family artefact) --")
print(as.data.frame(
  within_family_variance %>% filter(pass == "main", test == "pooled_variance_ratio") %>%
    transmute(trait,
              meanvar_native   = round(mv_rho_native, 2),
              meanvar_cultivar = round(mv_rho_cultivar, 2),
              size_native      = round(n_rho_native, 2),
              size_cultivar    = round(n_rho_cultivar, 2))
), row.names = FALSE)

# ---- 7. Figure -------------------------------------------------------

origin_cols <- c(native = "#2c7fb8", cultivar = "#d95f0e")
theme_ms <- theme_classic(base_size = 11) +
  theme(legend.position = "top", legend.title = element_blank(),
        axis.text = element_text(colour = "black"),
        strip.background = element_blank())

trait_labels <- c(height = "Height", stem_diam = "Stem diameter",
                  leaf_number = "Leaf number", AGB = "Above-ground biomass",
                  BGB = "Below-ground biomass", leaf_biomass = "Leaf biomass",
                  stem_biomass = "Stem biomass")

fig_dat <- family_dispersion_tbl %>%
  mutate(trait = factor(trait, levels = names(TRAITS), labels = trait_labels[names(TRAITS)]))

fig_wfv <- ggplot(fig_dat, aes(origin, log_s2, colour = origin, fill = origin)) +
  geom_boxplot(width = 0.5, alpha = 0.15, outlier.shape = NA) +
  geom_jitter(aes(size = n), width = 0.12, height = 0, alpha = 0.75) +
  facet_wrap(~ trait, scales = "free_y", nrow = 2) +
  scale_colour_manual(values = origin_cols, labels = c("Native", "Cultivar")) +
  scale_fill_manual(values = origin_cols, guide = "none") +
  scale_size_area(name = "offspring in family", max_size = 4,
                  breaks = c(3, 10, 20, 28)) +
  scale_x_discrete(labels = c(native = "Native", cultivar = "Cultivar")) +
  labs(x = NULL, y = expression("Within-family residual variance,  " * log(s[i]^2)),
       caption = paste("Each point is one maternal family. The low native outliers are the",
                       "smallest families (n = 2-4), whose\nvariance is estimated with 1-3 df;",
                       "they drive the equal-weight statistics but not the size-weighted ones.")) +
  theme_ms +
  theme(plot.caption = element_text(hjust = 0, colour = "grey30"))

ggsave(file.path(paths$figures, "fig_within_family_variance.png"), fig_wfv,
       width = 9, height = 5.5, dpi = 300)
ggsave(file.path(paths$figures, "fig_within_family_variance.pdf"), fig_wfv,
       width = 9, height = 5.5)

message("\nwrote within_family_variance.csv, within_family_variance_family.csv, ",
        "and fig_within_family_variance")
