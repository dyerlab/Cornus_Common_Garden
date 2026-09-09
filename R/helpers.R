# helpers.R
# -----------------------------------------------------------------------------
# Shared machinery for the mixed-model analyses (survival, single-time, and
# repeated-measures). Sourced after 00_setup.R.
#
# The quantitative-genetic quantities follow a half-sib design:
#   V_A  (additive genetic variance)      = 4 * V_family
#   V_P  (total phenotypic variance)      = V_origin + V_family + V_residual
#   h2   (narrow-sense heritability)      = V_A / V_P
# For a binomial GLMM the residual variance on the latent (logit) scale is
# fixed at pi^2 / 3 (Gilmour et al. 1985).
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({ library(lme4); library(dplyr); library(tibble); library(purrr) })

LATENT_RESID_VAR <- pi^2 / 3

# ---- Variance-component extraction --------------------------------------

# Return the variance for a named grouping factor, matching either order of a
# nested interaction name (e.g. "family:origin" vs "origin:family").
vc_value <- function(model, group) {
  vc <- as.data.frame(lme4::VarCorr(model))
  vc <- vc[is.na(vc$var2), ]                       # variances, not covariances
  alt <- if (grepl(":", group)) paste(rev(strsplit(group, ":")[[1]]), collapse = ":") else NULL
  hit <- vc$vcov[vc$grp %in% c(group, alt)]
  if (length(hit)) hit[1] else 0
}

# family variance, whether the model used (1|family) or (1|origin/family)
family_variance <- function(model) {
  v <- vc_value(model, "family:origin")
  if (v == 0) v <- vc_value(model, "family")
  v
}

origin_variance <- function(model) vc_value(model, "origin")
plant_variance  <- function(model) vc_value(model, "plant")

# ---- Heritability and variance partition ------------------------------

# `resid_var`: supply LATENT_RESID_VAR for a binomial GLMM; for a Gaussian LMM
# leave NULL and the residual variance is read from the model.
# A repeated-measures model with (1|plant) contributes V_plant to V_P.
quantgen_summary <- function(model, resid_var = NULL) {
  v_fam <- family_variance(model)
  v_org <- origin_variance(model)
  v_plt <- plant_variance(model)
  v_res <- if (is.null(resid_var)) {
    vc <- as.data.frame(lme4::VarCorr(model))
    vc$vcov[vc$grp == "Residual"]
  } else resid_var

  v_p  <- v_org + v_fam + v_plt + v_res
  c(
    V_origin   = v_org,
    V_family   = v_fam,
    V_plant    = v_plt,
    V_residual = v_res,
    P_origin   = if (v_p > 0) v_org / v_p else NA_real_,
    P_family   = if (v_p > 0) v_fam / v_p else NA_real_,
    P_plant    = if (v_p > 0) v_plt / v_p else NA_real_,
    P_residual = if (v_p > 0) v_res / v_p else NA_real_,
    h2         = if (v_p > 0) 4 * v_fam / v_p else NA_real_
  )
}

# ---- Parametric bootstrap CIs (bootMer) -------------------------------

boot_quantgen <- function(model, resid_var = NULL, nsim = 1000, seed = SEED) {
  fn <- function(m) quantgen_summary(m, resid_var = resid_var)
  bb <- lme4::bootMer(model, FUN = fn, nsim = nsim, type = "parametric",
                      use.u = FALSE, seed = seed,
                      .progress = "none")
  t0 <- bb$t0
  ci <- apply(bb$t, 2, function(col) stats::quantile(col, c(0.025, 0.975), na.rm = TRUE))
  tibble(
    parameter = names(t0),
    estimate  = as.numeric(t0),
    ci_low    = ci[1, ],
    ci_high   = ci[2, ],
    n_bad     = apply(bb$t, 2, function(col) sum(!is.finite(col)))
  )
}

# ---- AIC / AICc model comparison -------------------------------------

# `models`: named list. Mixed models and plain (g)lm are handled together.
# REML-fitted LMMs are refit by ML first, because a REML likelihood is not
# comparable with the ML likelihood of a plain lm (or of an LMM with different
# fixed effects). GLMMs are already ML.
aic_table <- function(models) {
  models <- lapply(models, function(m) {
    if (inherits(m, "lmerMod") && lme4::isREML(m)) update(m, REML = FALSE) else m
  })
  aic  <- vapply(models, function(m) AIC(m), numeric(1))
  aicc <- vapply(models, function(m) MuMIn::AICc(m), numeric(1))
  k    <- vapply(models, function(m) attr(logLik(m), "df"), numeric(1))
  tibble(model = names(models), k = k, AIC = aic, AICc = aicc) %>%
    mutate(
      dAIC  = AIC  - min(AIC),
      dAICc = AICc - min(AICc),
      wAIC  = exp(-0.5 * dAIC)  / sum(exp(-0.5 * dAIC)),
      wAICc = exp(-0.5 * dAICc) / sum(exp(-0.5 * dAICc))
    ) %>%
    arrange(AICc)
}

# ---- Likelihood-ratio test for a single random effect ---------------

# Compares two nested models by likelihood ratio. Both are refit with ML when
# needed; a mix of (g)lm and (g)lmer is handled by comparing log-likelihoods
# directly, halving the p-value because the null puts a variance on the
# boundary (Self & Liang 1987).
lrt <- function(reduced, full) {
  refit_ml <- function(m) if (inherits(m, "merMod") && lme4::isREML(m)) update(m, REML = FALSE) else m
  reduced <- refit_ml(reduced); full <- refit_ml(full)
  ll_r <- as.numeric(logLik(reduced)); ll_f <- as.numeric(logLik(full))
  df   <- attr(logLik(full), "df") - attr(logLik(reduced), "df")
  chisq <- max(0, 2 * (ll_f - ll_r))
  tibble(chisq = chisq, df = df, p = 0.5 * pchisq(chisq, df = df, lower.tail = FALSE))
}

# ---- Influence / sensitivity ---------------------------------------

# Individual observations with |studentised residual| > 3.
outlier_obs <- function(model, data, id_col) {
  r <- residuals(model, type = if (inherits(model, "glmerMod")) "pearson" else "response")
  s <- r / sd(r, na.rm = TRUE)
  data[[id_col]][which(abs(s) > 3)]
}

# Group-level Cook's distance from leave-one-group-out refits:
#   D_g = (b - b_(g))' V(b)^-1 (b - b_(g)) / p
# where V(b) is the fixed-effect covariance of the full model (HLMdiag's
# definition). Returns every group's D, sorted descending.
# `refit_fn(data)` must refit the model on a supplied data frame.
family_cooks_d <- function(model, data, refit_fn, group_col = "family") {
  groups <- unique(as.character(data[[group_col]]))
  b   <- lme4::fixef(model)
  Vi  <- tryCatch(solve(as.matrix(vcov(model))), error = function(e) NULL)
  p   <- length(b)
  cd <- vapply(groups, function(g) {
    d2 <- data[as.character(data[[group_col]]) != g, , drop = FALSE]
    m2 <- tryCatch(refit_fn(droplevels(d2)), error = function(e) NULL)
    if (is.null(m2)) return(NA_real_)
    db <- lme4::fixef(m2)[names(b)] - b
    if (is.null(Vi)) return(sum(db^2) / p)
    as.numeric(t(db) %*% Vi %*% db) / p
  }, numeric(1))
  sort(cd[is.finite(cd)], decreasing = TRUE)
}

# The single most influential group (largest Cook's distance).
most_influential_family <- function(model, data, refit_fn, group_col = "family") {
  cd <- family_cooks_d(model, data, refit_fn, group_col)
  if (!length(cd)) return(NULL)
  list(family = names(cd)[1], cooks_d = unname(cd[1]), all = cd)
}

fmt_ci <- function(est, lo, hi, digits = 3) {
  ifelse(is.na(est), "—",
         sprintf(paste0("%.", digits, "f [%.", digits, "f, %.", digits, "f]"), est, lo, hi))
}
