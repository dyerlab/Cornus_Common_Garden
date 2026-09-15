# 16_offspring_relatedness.R
# -----------------------------------------------------------------------------
# Side-check / supplementary analysis: every heritability estimate in the
# manuscript (Tab_GrowthHeritability, Tab_OriginHeritability, Tab_
# RepeatedHeritability) assumes maternal half-sibs share the expected
# relatedness of a half-sib design (r = 0.25, so V_A = V_family / r =
# 4 * V_family). Open pollination means some full-sib pairs (shared father,
# r = 0.5) are almost certainly mixed in with half-sibs within a mother's
# offspring array, which would push the true multiplier below 4 (down to a
# floor of 1/0.5 = 2 if every pair were full sibs). This script asks how much
# that matters in practice, using realized marker-based relatedness among
# genotyped siblings in place of the assumed constant.
#
# Method, in three steps:
#   1. Restrict to offspring with a COMPLETE genotype (no missing locus).
#      gstudio::genetic_relatedness(mode = "Nason") imputes a missing
#      genotype with population allele frequencies, and that imputation is
#      both statistically undesirable here (it shrinks a missing-data
#      individual's Fij with everyone toward the population mean,
#      mechanically biasing pairwise r downward) and, worse, buggy in the
#      installed gstudio: internally it assigns `loci[is.na(x)] <-
#      freq$Frequency`, and because `is.na(x)` is a length-N logical vector
#      recycled across an N x k matrix, only the FIRST missing individual at
#      a locus gets the correct frequency vector -- every subsequent missing
#      individual at that locus gets a rotated/shifted version of it
#      (verified directly against gstudio::rel_nason source). This is not a
#      rare edge case in this dataset: 71 of 522 mother x locus combinations
#      (among mothers with >= 2 genotyped offspring) have more than one
#      offspring missing at the same locus. Dropping incomplete genotypes
#      sidesteps both problems; we are not interested in per-mother point
#      estimates (many native families are small), so losing some offspring
#      per mother is immaterial to the question asked here.
#   2. For each genotyped mother (>= MIN_OFFSPRING complete-genotype
#      offspring), compute pairwise Nason's coancestry (r = 2*Fij; same
#      method and package as Tab_Relatedness in 15_relatedness_check.R)
#      among her own offspring -- full- and half-sib pairs are not
#      distinguished, this is realized relatedness across the whole
#      sibship, not a paternity call. All pairs are then POOLED within
#      origin (native / cultivar), not averaged per mother first, since the
#      question is the shape of the relatedness distribution across each
#      origin's half-sib arrays, weighted by how many pairs actually
#      support it.
#   3. For each trait, take the already-fitted P_family = V_family / V_P
#      point estimate from the origin-specific family model
#      (growth_origin_quantgen.csv, from 06_single_time.R -- the same
#      quantity behind Tab_OriginHeritability's h2 = 4 * P_family) and
#      re-express heritability under the pooled realized relatedness for
#      that origin instead of the assumed r = 0.25:
#          h2 = P_family / r
#      This holds V_family/V_P fixed (it is not re-estimated here) and only
#      varies the assumed-relatedness multiplier. We report this as a RANGE
#      (h2 at the full-sib floor r = 0.5, to h2 at the half-sib ceiling
#      r = 0.25, i.e. the value already published) bracketing a single point
#      estimate -- h2 at each origin's pooled MEDIAN r (not the mean; a
#      fraction of pairwise Fij values are outliers well outside [0, 1] with
#      only 9 loci, and those drag the mean off the half-sib/full-sib band --
#      see the Section 5 comment in the script) -- rather than a per-mother
#      distribution. This is meant only to give a reasonable expectation of
#      where heritability likely falls given realized sibship structure, not
#      a formal re-estimate; pollen-donor structure is left to a follow-up
#      manuscript.
#
# Input : data/genotypes.csv (all genotyped offspring, OffID != 0)
#         data/results/growth_origin_quantgen.csv (P_family per trait/origin,
#           from 06_single_time.R, "family" model, "original" sensitivity)
# Output: data/results/offspring_relatedness_pairs.csv     (every within-mother sib pair, complete genotypes only)
#         data/results/offspring_relatedness_by_mom.csv    (per-mother diagnostic summary, not used downstream)
#         data/results/offspring_relatedness_by_origin.csv (pooled pairwise relatedness summary by origin)
#         data/results/offspring_heritability_summary.csv  (h2 range + mean/median-r point estimates, per trait x origin)
#         media/fig_offspring_heritability.png / .pdf
#         media/fig_offspring_relatedness.png / .pdf
# -----------------------------------------------------------------------------

source("R/00_setup.R")
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(purrr)
  library(tibble)
  library(gstudio)
  library(ggplot2)
})

set.seed(SEED)
MIN_OFFSPRING <- 2 # a mother needs >= 2 complete-genotype offspring for a pairwise comparison

# ---- 1. Load genotyped offspring, drop incomplete genotypes ---------------

pop <- read_population(
  file.path(paths$data, "genotypes.csv"),
  type = "separated",
  locus.columns = 6:14
)
off_raw <- pop %>%
  filter(OffID != 0) %>%
  arrange(ID)

loci <- column_class(off_raw, "locus")
n_missing_loci <- Reduce(
  `+`,
  lapply(loci, function(l) as.integer(is.na(off_raw[[l]])))
)
off <- off_raw[n_missing_loci == 0, ]
message(
  nrow(off_raw) - nrow(off),
  " of ",
  nrow(off_raw),
  " genotyped offspring dropped for having ",
  "a missing genotype at >= 1 of ",
  length(loci),
  " loci (see header comment: gstudio's ",
  "missing-data imputation is both biased and, for >1 missing individual per locus, buggy). ",
  nrow(off),
  " offspring with complete genotypes retained across ",
  n_distinct(off$ID),
  " mothers (",
  n_distinct(off$ID[off$Origin == "native"]),
  " native, ",
  n_distinct(off$ID[off$Origin == "cultivar"]),
  " cultivar)."
)

moms <- off %>%
  count(ID, Origin, name = "n_offspring") %>%
  filter(n_offspring >= MIN_OFFSPRING)
message(
  nrow(moms),
  " mothers have >= ",
  MIN_OFFSPRING,
  " complete-genotype offspring (",
  sum(moms$Origin == "native"),
  " native, ",
  sum(moms$Origin == "cultivar"),
  " cultivar)"
)

# ---- 2. Within-mother pairwise Nason relatedness ---------------------------
# Same method as Tab_Relatedness (15_relatedness_check.R): Fij averaged
# across loci, reported on the relatedness scale r = 2*Fij. genetic_
# relatedness()'s `loci` argument is omitted for the same reason documented
# there (buggy for multi-locus vectors) -- loci are auto-detected via
# column_class(geno, "locus"). No missing genotypes reach this call (step 1),
# so gstudio's imputation branch never fires.
#
# Allele frequencies (pbar in Nason's Fij) are computed ONCE per origin, over
# ALL of that origin's complete-genotype offspring, rather than per mother.
# genetic_relatedness()/rel_nason() uses frequencies(x) computed on whatever
# data.frame x is passed in, so calling it separately per mother (as an
# earlier version of this script did) makes pbar a family-local statistic
# from as few as 2-3 offspring -- with that few individuals, loci are
# frequently fixed within the family (pbar = 1), and Fij's denominator
# pbar*(1-pbar) divides by zero, producing NA (confirmed: 20 mothers, 642 of
# 3231 pairs). Neither origin-level pool has a fixed locus (verified: 0 of 9
# loci monomorphic in either the native, n=245, or cultivar, n=316,
# complete-genotype pool), and using the origin as the reference population
# for its own families' kinship coefficients is the more standard choice
# for Nason's Fij besides.

relatedness_for_origin <- function(this_origin) {
  d <- off %>% filter(Origin == this_origin) %>% arrange(ID, OffID)
  fij_mat <- genetic_relatedness(d, mode = "Nason")
  r_mat <- 2 * fij_mat
  rownames(r_mat) <- colnames(r_mat) <- as.character(d$OffID)
  list(d = d, r_mat = r_mat)
}

origin_rel <- setNames(
  lapply(unique(moms$Origin), relatedness_for_origin),
  unique(moms$Origin)
)

sib_pairs_for_mom <- function(mom_id, this_origin) {
  info <- origin_rel[[this_origin]]
  d <- info$d %>% filter(ID == mom_id)
  idx <- match(as.character(d$OffID), colnames(info$r_mat))
  r_sub <- info$r_mat[idx, idx, drop = FALSE]
  pair_idx <- which(upper.tri(r_sub), arr.ind = TRUE)
  tibble(
    mom = mom_id,
    origin = this_origin,
    off1 = d$OffID[pair_idx[, 1]],
    off2 = d$OffID[pair_idx[, 2]],
    r_nason = r_sub[pair_idx]
  )
}

sib_pairs <- map2_dfr(moms$ID, moms$Origin, sib_pairs_for_mom)
message(
  sum(is.na(sib_pairs$r_nason)),
  " of ",
  nrow(sib_pairs),
  " pairs came back NA (should be 0 -- see comment above on origin-level allele frequencies)."
)
write_csv(
  sib_pairs,
  file.path(paths$results, "offspring_relatedness_pairs.csv")
)

# ---- 3. Per-mother summary (diagnostic only; not used downstream) ---------

by_mom <- sib_pairs %>%
  group_by(mom, origin) %>%
  summarise(
    n_pairs = n(),
    mean_r = mean(r_nason),
    sd_r = sd(r_nason),
    min_r = min(r_nason),
    max_r = max(r_nason),
    .groups = "drop"
  ) %>%
  left_join(moms %>% select(ID, n_offspring), by = c("mom" = "ID")) %>%
  mutate(
    band = case_when(
      mean_r < 0.25 ~ "below half-sib expectation (0.25)",
      mean_r <= 0.5 ~ "within half-sib/full-sib range [0.25, 0.5]",
      TRUE ~ "above full-sib expectation (0.5)"
    )
  )
write_csv(by_mom, file.path(paths$results, "offspring_relatedness_by_mom.csv"))

# ---- 4. Pooled relatedness distribution by origin --------------------------
# All within-mother pairs pooled within origin (not averaged per mother
# first) -- this is the distribution of interest: the shape of realized
# sibling relatedness across each origin's half-sib arrays.

by_origin <- sib_pairs %>%
  group_by(origin) %>%
  summarise(
    n_pairs = n(),
    n_mothers = n_distinct(mom),
    mean_r = mean(r_nason),
    median_r = median(r_nason),
    sd_r = sd(r_nason),
    min_r = min(r_nason),
    max_r = max(r_nason),
    n_below_0.25 = sum(r_nason < 0.25),
    n_above_0.5 = sum(r_nason > 0.5),
    .groups = "drop"
  )
write_csv(
  by_origin,
  file.path(paths$results, "offspring_relatedness_by_origin.csv")
)
message("\n--- pooled within-mother pairwise relatedness, by origin ---")
print(
  as.data.frame(
    by_origin %>% mutate(across(where(is.numeric), ~ round(.x, 4)))
  ),
  row.names = FALSE
)

# ---- 5. Ex post facto heritability under realized relatedness -------------
# h2 = P_family / r, for each trait's origin-specific family model
# (P_family fixed; only the assumed-relatedness multiplier 1/r varies).
# Reported as a range (r = 0.5 full-sib floor, to r = 0.25 half-sib ceiling
# -- the value already published) bracketing a single point estimate at
# each origin's pooled MEDIAN r. Median, not mean: with only 9 loci, Nason's
# pairwise Fij is noisy enough that a meaningful fraction of pairs fall
# outside the biologically sensible [0, 1] range (11% of cultivar pairs
# exceed r = 1, up to r = 7.2 for one pair), and those outliers drag the
# mean well off the half-sib/full-sib band while the median stays inside
# it. The full pooled distribution (mean included) is still in
# offspring_relatedness_by_origin.csv for reference.

p_family <- read_csv(
  file.path(paths$results, "growth_origin_quantgen.csv"),
  show_col_types = FALSE
) %>%
  filter(model == "family", sensitivity == "original") %>%
  transmute(
    trait,
    origin_cap = subset,
    p_family = estimate_P_family,
    h2_ceiling_halfsib = estimate_h2, # r = 0.25, already published
    h2_floor_fullsib = 2 * estimate_P_family
  ) # r = 0.5, for reference only

by_origin_cap <- by_origin %>%
  mutate(origin_cap = tools::toTitleCase(origin)) %>%
  select(origin_cap, median_r)

offspring_heritability_summary <- p_family %>%
  inner_join(by_origin_cap, by = "origin_cap") %>%
  rename(origin = origin_cap) %>%
  mutate(h2_at_median_r = p_family / median_r) %>%
  select(
    trait,
    origin,
    p_family,
    median_r,
    h2_at_median_r,
    h2_ceiling_halfsib,
    h2_floor_fullsib
  ) %>%
  arrange(origin, trait)
write_csv(
  offspring_heritability_summary,
  file.path(paths$results, "offspring_heritability_summary.csv")
)

message(
  "\n--- h2 range (full-sib r=0.5 to half-sib r=0.25) and pooled-median-r point estimate, by trait x origin ---"
)
print(
  as.data.frame(
    offspring_heritability_summary %>%
      mutate(across(where(is.numeric), ~ round(.x, 3)))
  ),
  row.names = FALSE
)

# Bootstrap uncertainty around each origin's pooled median r, for the
# fig_offspring_heritability violin. Resampling raw pairwise r_nason values
# directly is not usable here (h2 = p_family/r blows up near r = 0: the 99th
# percentile of a direct pointwise transform is already ~3.7x the published
# half-sib h2, max ~170x) -- that shows estimator instability, not biology.
# Instead this bootstraps the MEDIAN r itself: mothers are resampled with
# replacement within origin (preserving each drawn mother's full set of
# pairs, since pairs sharing a mother are not independent), the pooled
# median r is recomputed per replicate, and converted to h2 = p_family /
# median_r_boot. The resulting violin shows sampling uncertainty in the
# point estimate given the mothers actually sampled, not spread across
# sibships.
N_BOOT <- 1000
boot_median_r <- function(this_origin) {
  r_by_mom <- split(
    sib_pairs$r_nason[sib_pairs$origin == this_origin],
    sib_pairs$mom[sib_pairs$origin == this_origin]
  )
  moms_here <- names(r_by_mom)
  medians <- vapply(
    seq_len(N_BOOT),
    function(i) {
      samp <- sample(moms_here, length(moms_here), replace = TRUE)
      median(unlist(r_by_mom[samp], use.names = FALSE))
    },
    numeric(1)
  )
  tibble(origin_cap = tools::toTitleCase(this_origin), median_r_boot = medians)
}
r_boot <- map_dfr(unique(moms$Origin), boot_median_r)

h2_boot <- p_family %>%
  inner_join(r_boot, by = "origin_cap", relationship = "many-to-many") %>%
  rename(origin = origin_cap) %>%
  mutate(h2_boot = p_family / median_r_boot) %>%
  select(trait, origin, h2_boot)

# ---- 6. Figures -------------------------------------------------------------

theme_ms <- theme_classic(base_size = 11) +
  theme(
    legend.position = "none",
    axis.text = element_text(colour = "black"),
    strip.background = element_rect(fill = "grey92", colour = NA)
  )
origin_cols <- c(Native = "#2c7fb8", Cultivar = "#d95f0e")
trait_labels <- c(
  height = "Height",
  stem_diam = "Stem diam.",
  leaf_number = "Leaf no.",
  AGB = "AGB",
  BGB = "BGB",
  leaf_biomass = "Leaf biomass",
  stem_biomass = "Stem biomass"
)

# Fig_OffspringRelatedness: pooled pairwise relatedness within each origin --
# the Supplementary Materials figure showing the shape of the distribution
# between the half-sib (0.25) and full-sib (0.5) expectations. The density
# is estimated from the full pooled distribution (all pairs, outliers
# included); the x-axis is then zoomed via coord_cartesian so the 0.25/0.5
# reference lines and the bulk of the distribution stay legible -- a small
# tail of noisy pairwise estimates (see Section 5 comment) extends further
# than the plotted range but is not excluded from the calculation.
relatedness_plot_dat <- sib_pairs %>%
  mutate(
    origin_cap = factor(
      tools::toTitleCase(origin),
      levels = c("Native", "Cultivar")
    )
  )

fig_offspring_relatedness <- ggplot(
  relatedness_plot_dat,
  aes(r_nason, fill = origin_cap)
) +
  geom_density(alpha = 0.7, colour = NA) +
  geom_vline(xintercept = 0.25, linetype = 2, colour = "grey30") +
  geom_vline(xintercept = 0.5, linetype = 3, colour = "grey30") +
  scale_fill_manual(values = origin_cols) +
  coord_cartesian(xlim = c(-0.25, 1)) +
  facet_grid(origin_cap ~ .) +
  labs(
    x = "Pairwise within-mother offspring relatedness (Nason's r)",
    y = "Density"
  ) +
  theme_ms

ggsave(
  file.path(paths$figures, "fig_offspring_relatedness.png"),
  fig_offspring_relatedness,
  width = 5,
  height = 4.5,
  dpi = 300
)
ggsave(
  file.path(paths$figures, "fig_offspring_relatedness.pdf"),
  fig_offspring_relatedness,
  width = 5,
  height = 4.5
)

# Fig_OffspringHeritability: for each trait x origin, a violin of bootstrap
# uncertainty in h2 at the pooled median r (see comment above), the h2 range
# implied by the full-sib/half-sib bracket (thin bar), and the point
# estimate itself. Mean r is not used here -- see Section 5 comment on
# outlier-driven pairwise Fij with only 9 loci.
h2_plot_dat <- offspring_heritability_summary %>%
  mutate(
    origin = factor(origin, levels = c("Native", "Cultivar")),
    trait = factor(trait, levels = names(trait_labels), labels = trait_labels)
  )
h2_boot_plot <- h2_boot %>%
  mutate(
    origin = factor(origin, levels = c("Native", "Cultivar")),
    trait = factor(trait, levels = names(trait_labels), labels = trait_labels)
  )

fig_offspring_h2 <- ggplot(
  h2_plot_dat,
  aes(x = trait, colour = origin, fill = origin)
) +
  geom_violin(
    data = h2_boot_plot,
    aes(y = h2_boot),
    position = position_dodge(width = 0.6),
    trim = TRUE,
    scale = "width",
    width = 0.6,
    colour = NA,
    alpha = 0.35
  ) +
  geom_linerange(
    aes(ymin = h2_floor_fullsib, ymax = h2_ceiling_halfsib),
    position = position_dodge(width = 0.6),
    linewidth = 0.6,
    alpha = 0.7
  ) +
  geom_point(
    aes(y = h2_at_median_r),
    position = position_dodge(width = 0.6),
    size = 2,
    shape = 16
  ) +
  scale_colour_manual(values = origin_cols) +
  scale_fill_manual(values = origin_cols) +
  labs(x = NULL, y = expression(h^2)) +
  theme_ms +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

ggsave(
  file.path(paths$figures, "fig_offspring_heritability.png"),
  fig_offspring_h2,
  width = 7,
  height = 4,
  dpi = 300
)
ggsave(
  file.path(paths$figures, "fig_offspring_heritability.pdf"),
  fig_offspring_h2,
  width = 7,
  height = 4
)

message(
  "\nwrote offspring relatedness / heritability outputs to ",
  paths$results,
  " and ",
  paths$figures
)
