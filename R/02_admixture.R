# 02_admixture.R
# -----------------------------------------------------------------------------
# Compute Admixture Percentage (AP) for each genotyped seedling.
#
# AP is the percentage of a seedling's called microsatellite alleles that are
# shared with the single best-matching known cultivar voucher, across the nine
# diagnostic loci (cf020, cf125, cf213, cf273, cf581, cf585, cf597, cf634, cf701).
# It is used as a continuous, genetics-based measure of cultivar ancestry, in
# parallel with the binary native/cultivar classification from maternal location.
#
# Method (following Wadl et al. 2008, as adapted for this study):
#   For a seedling and a cultivar voucher, each locus scores
#       2  genotypes identical (same unordered allele pair)
#       1  at least one allele shared
#       0  no allele shared
#   The seedling is matched to the voucher with the highest total score.
#   AP = 100 * (total score) / (2 * number of loci compared with data on both sides)
#
# Inputs
#   data/genotypes.csv                                   seedling + maternal genotypes (anonymized)
#   SENSITIVE/data_raw/plate7_group1 - plate7_group1.csv  cultivar voucher panel, loci group 1
#   SENSITIVE/data_raw/plate7_group2 - plate7_group2.csv  cultivar voucher panel, loci group 2
#     (the plate7 files carry original sample names and are not part of the
#      public repository; only the derived AP values are published)
#
# Output
#   data/admixture.csv   mom, id, AP, n_loci   (one row per genotyped seedling)
#
# data/admixture.csv is the one generated file that IS committed to the
# repository: it is derived from the non-public voucher panel, so a clone
# without SENSITIVE/ cannot rebuild it. When SENSITIVE/data_raw/ is absent this
# script skips and the committed file is used as-is by steps 05-12.
#
# Validation (skipped automatically if SENSITIVE/ is absent)
#   Compares AP/100 against the legacy pipeline's per-voucher score grid
#   (SENSITIVE/legacy/common_garden_hd/data/processed/offspring_cultivar_comparison_1.csv),
#   taking each seedling's best voucher by total shared-allele score. The
#   legacy "match" file (offspring_cultivar_match_1.csv) is stale relative to
#   that grid and is not used here.
# -----------------------------------------------------------------------------

source("R/00_setup.R")
suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr); library(stringr); library(purrr) })

LOCI <- c("cf020", "cf125", "cf213", "cf273", "cf581", "cf585", "cf597", "cf634", "cf701")

# ---- 0. Guard: needs the non-public voucher panel -------------------------
# The plate7 voucher genotypes live only in SENSITIVE/. Their product,
# data/admixture.csv, is committed, so a clone without SENSITIVE/ skips this
# step and the downstream scripts use the committed file.
plate7 <- file.path(paths$sensitive, "data_raw",
                    c("plate7_group1 - plate7_group1.csv",
                      "plate7_group2 - plate7_group2.csv"))
ap_file <- file.path(paths$data, "admixture.csv")
if (!all(file.exists(plate7))) {
  if (!file.exists(ap_file))
    stop("SENSITIVE/data_raw/ absent and ", ap_file, " missing — cannot continue.")
  message("SENSITIVE/data_raw/ not found — skipping AP recomputation; ",
          "keeping committed ", ap_file)
  quit(save = "no", status = 0)
}

# ---- 1. Seedling genotypes ------------------------------------------------
# genotypes.csv stores each locus as "allele1:allele2". Split into two integer
# columns per locus. OffID 0 marks the maternal tree; keep only offspring.

geno_raw <- read_csv(file.path(paths$data, "genotypes.csv"), show_col_types = FALSE)

split_alleles <- function(df, loci) {
  for (l in loci) {
    parts <- str_split_fixed(df[[l]], ":", 2)
    a1 <- suppressWarnings(as.integer(parts[, 1]))
    a2 <- suppressWarnings(as.integer(parts[, 2]))
    a2[is.na(a2) & !is.na(a1)] <- a1[is.na(a2) & !is.na(a1)]   # homozygous shorthand
    df[[paste0(l, "_1")]] <- pmin(a1, a2)
    df[[paste0(l, "_2")]] <- pmax(a1, a2)
  }
  df
}

offspring <- geno_raw %>%
  filter(OffID != 0) %>%
  transmute(mom = ID, id = OffID) %>%
  bind_cols(split_alleles(geno_raw %>% filter(OffID != 0), LOCI) %>% select(ends_with("_1"), ends_with("_2")))

# ---- 2. Cultivar voucher panel ------------------------------------------
# The plate7 files hold a mix of named cultivar vouchers, maternal-tree
# re-runs, and controls. Keep only the named cultivar lines: the sample name
# is in "<name>_multi<n>_<well>_Mar-4-2021" format and starts with a letter,
# and is not a negative control or the "tree2" re-run.

read_plate7 <- function(file, loci_cols) {
  x <- read_csv(file, show_col_types = FALSE, name_repair = "minimal")
  names(x) <- c("Name", loci_cols)
  x %>%
    mutate(
      tag   = str_extract(Name, "^[A-Za-z][A-Za-z0-9]*_multi\\d+_[A-Z]\\d+_Mar-4-2021"),
      voucher = str_extract(tag, "^[A-Za-z][A-Za-z0-9]*"),
      well    = str_extract(tag, "(?<=_)[A-Z]\\d+(?=_)")
    ) %>%
    filter(!is.na(tag),
           !str_detect(tolower(voucher), "negcont|tree2|common")) %>%
    mutate(across(all_of(loci_cols), ~ suppressWarnings(as.integer(.x))))
}

g1 <- read_plate7(file.path(paths$sensitive, "data_raw", "plate7_group1 - plate7_group1.csv"),
                  c("cf213_1","cf213_2","cf585_1","cf585_2","cf634_1","cf634_2",
                    "cf597_1","cf597_2","cf581_1","cf581_2"))
g2 <- read_plate7(file.path(paths$sensitive, "data_raw", "plate7_group2 - plate7_group2.csv"),
                  c("cf020_1","cf020_2","cf701_1","cf701_2","cf125_1","cf125_2",
                    "cf273_1","cf273_2"))

# plate7 loci group 2 was run on a separate plate with a 1 bp calibration
# offset at cf701 relative to the family genotypes; align it here.
g2 <- g2 %>% mutate(cf701_1 = cf701_1 - 1L, cf701_2 = cf701_2 - 1L)

vouchers <- inner_join(
  g1 %>% select(voucher, well, starts_with("cf")),
  g2 %>% select(voucher, well, starts_with("cf")),
  by = c("voucher", "well")
)

# order the two alleles per locus (min, max) so genotype comparison is unordered
for (l in LOCI) {
  lo <- pmin(vouchers[[paste0(l, "_1")]], vouchers[[paste0(l, "_2")]])
  hi <- pmax(vouchers[[paste0(l, "_1")]], vouchers[[paste0(l, "_2")]])
  vouchers[[paste0(l, "_1")]] <- lo
  vouchers[[paste0(l, "_2")]] <- hi
}

message("cultivar voucher genotypes: ", nrow(vouchers),
        " (", n_distinct(vouchers$voucher), " distinct lines)")

# ---- 3. Score each seedling against each voucher -----------------------

# per-locus score for one seedling allele pair vs one voucher allele pair
locus_score <- function(s1, s2, v1, v2) {
  if (is.na(s1) || is.na(v1)) return(NA_integer_)
  s <- c(s1, s2); v <- c(v1, v2)
  if (setequal(s, v)) return(2L)
  if (length(intersect(s, v)) > 0) return(1L)
  0L
}

voucher_mat <- as.matrix(vouchers[, paste0(rep(LOCI, each = 2), c("_1", "_2"))])

# The seedling is matched to the voucher with the highest total shared-allele
# score. Ties are broken toward the voucher compared on more loci (the more
# conservative choice: a larger denominator), so the result does not depend on
# the order of the voucher table.
ap_one_seedling <- function(row) {
  s <- unlist(row[paste0(rep(LOCI, each = 2), c("_1", "_2"))])
  best_score <- -1; best_loci <- -1
  for (j in seq_len(nrow(voucher_mat))) {
    v <- voucher_mat[j, ]
    per_locus <- vapply(seq_along(LOCI), function(k) {
      idx <- 2 * k - 1
      locus_score(s[idx], s[idx + 1], v[idx], v[idx + 1])
    }, integer(1))
    tot  <- sum(per_locus, na.rm = TRUE)
    ncmp <- sum(!is.na(per_locus))
    if (ncmp > 0 && (tot > best_score || (tot == best_score && ncmp > best_loci))) {
      best_score <- tot; best_loci <- ncmp
    }
  }
  tibble(AP = 100 * best_score / (2 * best_loci), n_loci = best_loci)
}

admixture <- offspring %>%
  mutate(res = purrr::pmap(pick(everything()), ~ ap_one_seedling(list(...)))) %>%
  select(mom, id, res) %>%
  unnest(res) %>%
  filter(!is.na(AP))

write_csv(admixture, file.path(paths$data, "admixture.csv"))

message("AP computed for ", nrow(admixture), " seedlings; ",
        "mean ", round(mean(admixture$AP), 1), "%, range ",
        round(min(admixture$AP), 1), "-", round(max(admixture$AP), 1), "%")

# ---- 4. Validation against the legacy pipeline ------------------------

legacy_grid <- file.path(paths$sensitive, "legacy/common_garden_hd/data/processed/offspring_cultivar_comparison_1.csv")
key_file    <- file.path(paths$sensitive, "mom_anonymization_key.csv")

if (file.exists(legacy_grid) && file.exists(key_file)) {
  key <- read_csv(key_file, show_col_types = FALSE)            # mom_original, origin, mom_ID

  legacy <- read_csv(legacy_grid, show_col_types = FALSE) %>%
    group_by(unique_ID) %>%
    slice_max(order_by = total_match_sa, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(
      mom_original = str_split_fixed(sample_ID, "\\.", 2)[, 1],
      id           = str_split_fixed(sample_ID, "\\.", 2)[, 2]
    ) %>%
    left_join(key, by = "mom_original") %>%
    transmute(mom = mom_ID, id = as.character(id),
              prop_sa_legacy = prop_sa, n_loci_legacy = loci_compared)

  cmp <- admixture %>%
    mutate(id = as.character(id)) %>%
    inner_join(legacy, by = c("mom", "id")) %>%
    mutate(diff = AP / 100 - prop_sa_legacy)

  message("\nvalidation vs legacy score grid: ", nrow(cmp), " seedlings in common")
  message("  AP identical (|diff| < 0.001): ", sum(abs(cmp$diff) < 0.001, na.rm = TRUE))
  message("  loci compared identical      : ", sum(cmp$n_loci == cmp$n_loci_legacy, na.rm = TRUE))
  message("  mean |diff| in AP/100        : ", round(mean(abs(cmp$diff), na.rm = TRUE), 4))
  write_csv(cmp, file.path(paths$derived, "admixture_vs_legacy.csv"))
} else {
  message("\nvalidation skipped: legacy grid or anonymization key not found")
}
