# 13_multilocus_ld.R
# -----------------------------------------------------------------------------
# First-pass test for non-random multilocus allelic association (linkage
# disequilibrium) between the native and cultivar groups, as a complement to
# the individual-allele-sharing approach in 02_admixture.R.
#
# Rationale: AP (02_admixture.R) and the full-vs-Wadl-panel comparison
# (Tab_AdmixturePanelComparison) both look at individual alleles one locus at
# a time, and both found a diffuse pattern with no sharp separation between
# groups. That leaves open a different question: do the nine loci co-vary
# with one another more than expected under random mating within a group,
# even though no single locus (or simple sum across loci) separates the
# groups? Elevated multilocus linkage disequilibrium can appear from cryptic
# structure of exactly that kind. It is also the standard signature used to
# detect clonal or partially-clonal reproduction (cultivars are clonally
# propagated), which is why this uses the same tool developed for that
# purpose (poppr's index of association).
#
# Method: standardized index of association (Ia, rbarD; Agapow & Burt 2001),
# with a 999-permutation test of the null of no multilocus association
# (linkage equilibrium), via poppr (Kamvar, Tabima & Grunwald 2014), computed
# for two datasets:
#   (1) maternal trees only (58 unrelated individuals) -- the clean test:
#       no two of these individuals are siblings, so any Ia signal here
#       cannot be a family/kinship artifact.
#   (2) maternal trees + their offspring, pooled (886 individuals) -- included
#       for completeness only. This dataset contains many maternal half- and
#       full-sib families, and sibling groups are *expected* to show elevated
#       multilocus association purely from shared parentage (a Wahlund-type
#       effect from family structure), regardless of any cultivar signal. A
#       positive result here, without a matching positive result in (1), is
#       not evidence of a cultivar/native effect.
# Both are computed separately within each Origin (native, cultivar) and for
# the pooled/"Total" sample, via poppr()'s built-in per-population summary.
#
# This is an exploratory first pass, not yet wired into run_all.R or
# 11_tables.R. If it turns up a real, unconfounded (native vs. cultivar,
# maternal-only) signal, promote it into the numbered pipeline and add a
# Results paragraph / Tab_MultilocusLD. If not, the intended manuscript
# treatment is a one-line parenthetical noting the test was done and returned
# no evidence of multilocus disequilibrium (data not shown).
#
# Input : data/genotypes.csv
# Output: data/results/multilocus_ld.csv
# -----------------------------------------------------------------------------

source("R/00_setup.R")
suppressPackageStartupMessages({ library(readr); library(dplyr) })

for (pkg in c("adegenet", "poppr")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Missing required package: ", pkg, "\nInstall with install.packages(\"", pkg, "\")")
  }
}
suppressPackageStartupMessages({ library(adegenet); library(poppr) })

LOCI <- c("cf020", "cf125", "cf213", "cf273", "cf581", "cf585", "cf597", "cf634", "cf701")

geno <- read_csv(file.path(paths$data, "genotypes.csv"), show_col_types = FALSE)

# genotypes.csv stores missing loci as an empty string, not the literal "NA".
build_genind <- function(df) {
  gid <- df2genind(
    df[LOCI],
    sep       = ":",
    ind.names = paste(df$ID, df$OffID, sep = "_"),
    pop       = df$Origin,
    ploidy    = 2,
    NA.char   = ""
  )
  # Drop individuals missing more than 20% of loci (>1 of 9) rather than
  # impute; poppr's own recommended step before ia()/poppr()
  # (type = "geno" filters individuals, as opposed to "loci" which filters loci).
  missingno(gid, type = "geno", cutoff = 0.2, quiet = TRUE)
}

gid_maternal <- build_genind(geno %>% filter(OffID == 0))
gid_pooled   <- build_genind(geno)

message("maternal-only dataset: ", nInd(gid_maternal), " individuals (",
        table(pop(gid_maternal))["native"], " native, ",
        table(pop(gid_maternal))["cultivar"], " cultivar) after missing-data filter")
message("pooled (maternal + offspring) dataset: ", nInd(gid_pooled), " individuals (",
        table(pop(gid_pooled))["native"], " native, ",
        table(pop(gid_pooled))["cultivar"], " cultivar) after missing-data filter")

set.seed(SEED)
ld_maternal <- poppr(gid_maternal, sample = 999, quiet = TRUE) %>%
  mutate(dataset = "maternal (unrelated individuals)")

set.seed(SEED)
ld_pooled <- poppr(gid_pooled, sample = 999, quiet = TRUE) %>%
  mutate(dataset = "pooled maternal + offspring (family-structured; confounded, see header)")

ld_result <- bind_rows(ld_maternal, ld_pooled) %>%
  select(dataset, Pop, N, Ia, p.Ia, rbarD, p.rD)

write_csv(ld_result, file.path(paths$results, "multilocus_ld.csv"))

message("\nMultilocus linkage disequilibrium (index of association; 999 permutations):")
print(as.data.frame(ld_result), row.names = FALSE)
message("\np.Ia / p.rD < 0.05 rejects the null of linkage equilibrium (no multilocus association)",
        " for that population. Judge the native/cultivar question from the 'maternal' rows only --",
        " the 'pooled' rows are expected to show significant LD from family structure alone.")
