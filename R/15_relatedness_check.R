# 15_relatedness_check.R
# -----------------------------------------------------------------------------
# Side-check (referenced in the manuscript Discussion as Tab_Relatedness): is
# the asymmetry in family-level variance we see between the putative-native and
# putative-cultivar common-garden subsets (Tab_OriginHeritability) explained,
# at least in part, by the two groups' maternal trees being differently
# related to each other -- rather than by anything intrinsic to "cultivar"
# vs. "native" ancestry itself? Native-origin mothers come from two nearby
# forest sites and *C. florida* shows isolation-by-distance (Dyer et al.
# 2012), so they may be more closely related to one another than the
# cultivar-origin mothers, which are 20 independently-bred named accessions.
# If so, more within-group relatedness among native mothers -> less
# *independent* genetic variation among native families -> less power to
# detect a family variance component there, confounded with (not caused by)
# the origin label.
#
# Restricted to the maternal trees whose families were planted in THIS
# common garden (study_genotypes_file() in 00_setup.R): 28 of the 36 garden
# mothers were genotyped (9 native, 19 cultivar; nat18, nat24, nat28,
# nat30-33 and urb10 were not). An earlier version used all 58 genotyped
# adults in the larger collection, but 30 of those have no family in the
# garden, so they cannot bear on the garden-based family-variance asymmetry
# this check is meant to address. 02-04 and 13-16 all use the same
# study-only set.
#
# Method: pairwise relatedness from Nason's coancestry coefficient Fij
# (gstudio::rel_nason, averaged across loci, reported on the relatedness
# scale r = 2*Fij, following rel_nason()'s own as.relatedness=TRUE convention)
# -- the same coefficient already underlying the isolation-by-distance
# framing (Dyer et al. 2012) used elsewhere in this manuscript, computed
# once from the fixed genotypes for all C(28,2) = 378 pairs among the 28
# genotyped garden mothers.
#
# (Wang's (2002) dyadic-moment estimator has the lowest sampling variance
# among the available pairwise estimators and was the first choice here, but
# the R package that implements it, `related`, requires compiling Fortran
# source, isn't on CRAN, and its GitHub repo isn't laid out as an installable
# package -- so this uses gstudio's Nason coefficient instead, which is
# already part of this project's toolchain and needs no extra install.)
#
# Rather than treating the 36 within-native or 171 within-cultivar pairs as
# independent samples (they are not -- every individual appears in multiple
# pairs), the origin label is permuted across the same 28 individuals,
# holding the true group sizes (9/19) fixed, and the same pairwise values
# are relabeled "within-native" / "within-cultivar" under each permutation.
# This gives a null distribution for (a) the difference in mean within-group
# relatedness and (b) a two-sample Kolmogorov-Smirnov statistic between the
# two within-group distributions, without ever assuming the pairs are
# independent and without needing a separate rarefaction step -- permuting
# real individuals into groups of the true sizes already reproduces the
# correct sampling variability for those sizes under the null.
#
# Input : data/genotypes.csv via study_genotypes_file() (garden mothers, OffID == 0)
# Output: data/results/relatedness_check.csv          (observed stats + p-values)
#         data/results/relatedness_permutations.csv   (null distributions)
# -----------------------------------------------------------------------------

source("R/00_setup.R")
suppressPackageStartupMessages({ library(readr); library(dplyr); library(purrr); library(tibble); library(gstudio) })

set.seed(SEED)
NPERM <- 9999
LOCI  <- c("cf020", "cf125", "cf213", "cf273", "cf581", "cf585", "cf597", "cf634", "cf701")

# ---- 1. Load genotyped garden mothers -----------------------------------

pop <- read_population(study_genotypes_file(),
                        type = "separated", locus.columns = 6:14)

geno <- pop %>%
  filter(OffID == 0) %>%
  arrange(ID)

message(nrow(geno), " genotyped adult trees (",
        sum(geno$Origin == "native"), " native, ",
        sum(geno$Origin == "cultivar"), " cultivar)")

origin <- setNames(geno$Origin, geno$ID)
ids    <- geno$ID

# ---- 2. Nason pairwise relatedness for all pairs ---------------------
# genetic_relatedness(..., mode="Nason") returns the multilocus-averaged
# coancestry Fij; multiply by 2 to put it on the relatedness scale r (what
# rel_nason()'s own as.relatedness=TRUE argument does for a single locus).

# NB: genetic_relatedness()'s `loci` arg is buggy for multi-locus vectors --
# `if (is.na(loci))` errors when loci has length > 1. LOCI matches the loci
# gstudio would auto-detect via column_class(geno, "locus") anyway, so we
# just omit the argument and let it default.
fij_mat <- genetic_relatedness(geno, mode = "Nason")
r_mat   <- 2 * fij_mat

pair_idx <- which(upper.tri(r_mat), arr.ind = TRUE)
rel <- tibble(
  ind1    = ids[pair_idx[, 1]],
  ind2    = ids[pair_idx[, 2]],
  r_nason = r_mat[pair_idx]
)

# true within/cross group label for each pair, from the real origin classification
pair_type <- function(o1, o2) ifelse(o1 == o2, o1, "cross")
rel <- rel %>%
  mutate(o1 = origin[ind1], o2 = origin[ind2],
         type = pair_type(o1, o2))

message("\n--- observed pair counts ---")
print(table(rel$type))

# ---- 3. Observed statistics ----------------------------------------------

ks_stat <- function(x, y) suppressWarnings(ks.test(x, y)$statistic)

obs_native   <- rel$r_nason[rel$type == "native"]
obs_cultivar <- rel$r_nason[rel$type == "cultivar"]

obs_mean_diff <- mean(obs_native) - mean(obs_cultivar)
obs_ks        <- ks_stat(obs_native, obs_cultivar)

# ---- 4. Permutation: shuffle origin labels across the 28 mothers, holding ---
#         group sizes (9/19) fixed; relabel the SAME fixed pairwise values --

n_native <- sum(origin == "native")

perm_stats <- map_dfr(seq_len(NPERM), function(i) {
  perm_origin <- setNames(rep("cultivar", length(ids)), ids)
  perm_origin[sample(ids, n_native)] <- "native"

  o1 <- perm_origin[rel$ind1]; o2 <- perm_origin[rel$ind2]
  ptype <- pair_type(o1, o2)

  x <- rel$r_nason[ptype == "native"]
  y <- rel$r_nason[ptype == "cultivar"]
  tibble(mean_diff = mean(x) - mean(y), ks = ks_stat(x, y))
})

p_mean_diff <- mean(abs(perm_stats$mean_diff) >= abs(obs_mean_diff))
p_ks        <- mean(perm_stats$ks >= obs_ks)

write_csv(perm_stats, file.path(paths$results, "relatedness_permutations.csv"))

# ---- 5. Summary  -----------------------------------------------------------

summary_tbl <- tibble(
  what = c("n_native_adults", "n_cultivar_adults", "n_pairs_native", "n_pairs_cultivar",
           "n_pairs_cross", "mean_r_native", "mean_r_cultivar", "mean_r_cross",
           "obs_mean_diff", "perm_p_mean_diff", "obs_ks", "perm_p_ks", "n_permutations"),
  value = c(sum(origin == "native"), sum(origin == "cultivar"),
            length(obs_native), length(obs_cultivar), sum(rel$type == "cross"),
            mean(obs_native), mean(obs_cultivar), mean(rel$r_nason[rel$type == "cross"]),
            obs_mean_diff, p_mean_diff, as.numeric(obs_ks), p_ks, NPERM)
)
write_csv(summary_tbl, file.path(paths$results, "relatedness_check.csv"))

message("\n--- relatedness_check summary ---")
print(as.data.frame(summary_tbl), row.names = FALSE)
