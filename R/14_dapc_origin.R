# 14_dapc_origin.R
# -----------------------------------------------------------------------------
# First-pass test for whether a multivariate combination of the nine loci
# separates putative-native from putative-cultivar maternal trees, as a third
# complement to 02_admixture.R (individual allele sharing; AP) and
# 13_multilocus_ld.R (within-group multilocus association; index of
# association). This asks a different question from either: not "do the
# groups share fewer/more alleles on average" (AP) and not "is there hidden
# structure within a group" (Ia), but "is there ANY linear combination of
# the nine loci, taken together, that separates the two groups better than
# chance."
#
# Method: Discriminant Analysis of Principal Components (DAPC; Jombart,
# Devillard & Balloux 2010), restricted to the 58 unrelated maternal trees
# for the same reason as 13_multilocus_ld.R -- pooling in the offspring would
# let shared maternal-family membership (not origin) drive apparent
# separation, since siblings correlate at every locus together.
#
# DAPC on its own is not a significance test: with enough retained PCA axes
# relative to sample size, a DAPC will show visual separation between almost
# any two labeled groups, including random ones (a well-documented overfitting
# risk -- see the adegenet DAPC vignette). This script guards against that in
# two ways:
#   1. The number of PCA axes retained is chosen by cross-validation
#      (xvalDapc) on the real labels, not by hand.
#   2. Significance is assessed by permutation: with that same number of PCA
#      axes held fixed, the origin labels are shuffled 999 times (preserving
#      the 26/32 group sizes) and the same repeated-holdout classification
#      accuracy is recomputed each time. The reported p-value is the fraction
#      of permutations whose accuracy meets or exceeds the observed accuracy.
#
# Input : data/genotypes.csv
# Output: data/results/dapc_origin.csv   (one row: observed vs. permutation
#         summary) and data/results/dapc_origin_permutations.csv (the full
#         null distribution, for a histogram/sanity check if wanted)
# -----------------------------------------------------------------------------

source("R/00_setup.R")
suppressPackageStartupMessages({ library(readr); library(dplyr) })

if (!requireNamespace("adegenet", quietly = TRUE)) {
  stop("Missing required package: adegenet\nInstall with install.packages(\"adegenet\")")
}
suppressPackageStartupMessages({ library(adegenet) })

LOCI <- c("cf020", "cf125", "cf213", "cf273", "cf581", "cf585", "cf597", "cf634", "cf701")

geno <- read_csv(file.path(paths$data, "genotypes.csv"), show_col_types = FALSE) %>%
  filter(OffID == 0)   # maternal trees only -- see header

gid <- df2genind(
  geno[LOCI],
  sep       = ":",
  ind.names = geno$ID,
  pop       = geno$Origin,
  ploidy    = 2,
  NA.char   = ""
)

message(nInd(gid), " maternal trees (", table(pop(gid))["native"], " native, ",
        table(pop(gid))["cultivar"], " cultivar)")

# Allele-count table, one column per allele; missing genotypes replaced by
# the mean allele frequency (adegenet's recommended handling before PCA/DAPC
# -- PCA cannot itself accommodate missing values).
mat <- tab(gid, NA.method = "mean")
grp <- pop(gid)

# ---- 1. Choose the number of PCA axes to retain, by cross-validation ------
# Kept deliberately small relative to N = 58 (n.pca.max = 20) -- retaining
# PCs approaching N-1 is exactly the overfitting regime DAPC is criticized
# for; n.da = 1 is the only possible value for a 2-group discriminant.

set.seed(SEED)
n_pca <- tryCatch({
  xval <- xvalDapc(mat, grp, n.pca.max = 20, n.da = 1, training.set = 0.9,
                    n.rep = 100, xval.plot = FALSE)
  xval$DAPC$n.pca
}, error = function(e) {
  message("xvalDapc failed (", conditionMessage(e), "); falling back to n.pca = 5")
  5L
})
message("PCA axes retained (cross-validated): ", n_pca)

# ---- 2. Repeated-holdout classification accuracy at that fixed n.pca -----
# 80% train / 20% test, stratified by group, repeated n.rep times; reports
# the mean out-of-sample assignment accuracy.

cv_success <- function(mat, grp, n.pca, n.rep = 30, training.frac = 0.8) {
  n <- nrow(mat)
  successes <- numeric(n.rep)
  for (i in seq_len(n.rep)) {
    train_idx <- unlist(lapply(levels(grp), function(g) {
      idx_g <- which(grp == g)
      sample(idx_g, size = max(1, round(length(idx_g) * training.frac)))
    }))
    test_idx <- setdiff(seq_len(n), train_idx)
    if (length(test_idx) == 0) { successes[i] <- NA; next }
    dp   <- dapc(mat[train_idx, , drop = FALSE], grp[train_idx], n.pca = n.pca, n.da = 1)
    pred <- predict(dp, newdata = mat[test_idx, , drop = FALSE])
    successes[i] <- mean(pred$assign == grp[test_idx])
  }
  mean(successes, na.rm = TRUE)
}

set.seed(SEED)
obs_accuracy <- cv_success(mat, grp, n.pca = n_pca, n.rep = 30)

no_info_rate <- max(table(grp)) / length(grp)
message(sprintf("observed repeated-holdout accuracy: %.3f  (no-information rate: %.3f)",
                 obs_accuracy, no_info_rate))

# ---- 3. Permutation test: shuffle origin labels, same n.pca, same procedure ----

N_PERM <- 999
message("running ", N_PERM, " label permutations (this may take a minute)...")
set.seed(SEED)
perm_accuracy <- vapply(seq_len(N_PERM), function(i) {
  cv_success(mat, sample(grp), n.pca = n_pca, n.rep = 30)
}, numeric(1))

p_value <- (sum(perm_accuracy >= obs_accuracy) + 1) / (N_PERM + 1)

message(sprintf("\npermutation null: mean accuracy = %.3f, 95th pctile = %.3f",
                 mean(perm_accuracy), quantile(perm_accuracy, 0.95)))
message(sprintf("empirical p-value (observed vs. %d permutations): %.4f", N_PERM, p_value))

# ---- 4. Write results ------------------------------------------------------

result <- tibble(
  n_maternal        = nInd(gid),
  n_native          = as.integer(table(grp)["native"]),
  n_cultivar        = as.integer(table(grp)["cultivar"]),
  n_pca_retained    = n_pca,
  observed_accuracy = obs_accuracy,
  no_info_rate      = no_info_rate,
  perm_mean         = mean(perm_accuracy),
  perm_95th_pctile  = quantile(perm_accuracy, 0.95),
  n_permutations    = N_PERM,
  p_value           = p_value
)
write_csv(result, file.path(paths$results, "dapc_origin.csv"))
write_csv(tibble(permutation = seq_len(N_PERM), accuracy = perm_accuracy),
          file.path(paths$results, "dapc_origin_permutations.csv"))

message("\np_value < 0.05 would indicate the two groups are more separable, on the full",
        " multivariate genotype, than a random same-sized split of these 58 trees.")
