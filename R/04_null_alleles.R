# 04_null_alleles.R
# -----------------------------------------------------------------------------
# Two related marker checks:
#
#  (1) Null-allele frequency per locus and origin class, estimated by maximum
#      likelihood (Genepop's default estimator, which assumes apparent null
#      homozygotes are true nulls), with a 95% CI. Loci whose CI excludes zero
#      in both origin classes are dropped to form a "reduced" marker set.
#
#  (2) Combined multilocus exclusion probability (gstudio, Dyer 2009) for the
#      full nine-locus set and for the reduced set.
#
# Input : data/genotypes.csv
# Output: data/results/null_alleles.csv        (Tab_NullAlleles)
#         data/results/exclusion_probability.csv
#         data/derived/reduced_marker_set.txt   (loci kept in the reduced set)
#
# Uses genepop::nulls (Rousset 2008) and gstudio::genetic_diversity(mode = "Pe").
# -----------------------------------------------------------------------------

source("R/00_setup.R")
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr); library(stringr); library(gstudio) })

LOCI <- c("cf020", "cf125", "cf213", "cf273", "cf581", "cf585", "cf597", "cf634", "cf701")

geno <- read_csv(file.path(paths$data, "genotypes.csv"), show_col_types = FALSE) %>%
  mutate(group = if_else(OffID == 0, "maternal", "offspring"),
         ind   = paste(ID, OffID, sep = "_"))

# ---- 1. Write a Genepop-format file --------------------------------------
# Alleles are padded to three digits; a genotype is the two alleles
# concatenated (six digits); missing data is "000000". Individuals are
# grouped into populations by origin.

allele_to_3digit <- function(x) {
  a <- suppressWarnings(as.integer(x))
  ifelse(is.na(a), "000", sprintf("%03d", a))
}

encode_locus <- function(col) {
  parts <- str_split_fixed(col, ":", 2)
  a1 <- allele_to_3digit(parts[, 1])
  a2 <- allele_to_3digit(parts[, 2])
  a2 <- ifelse(a2 == "000" & a1 != "000", a1, a2)   # homozygous shorthand
  paste0(a1, a2)
}

genepop_file <- file.path(paths$derived, "genotypes_genepop.txt")

gp_body <- geno %>%
  mutate(across(all_of(LOCI), encode_locus)) %>%
  arrange(Origin, ind)

con <- file(genepop_file, "w")
writeLines("Cornus florida common garden: nine cultivar-diagnostic microsatellites", con)
writeLines(LOCI, con)
for (o in unique(gp_body$Origin)) {
  writeLines("POP", con)
  sub <- gp_body %>% filter(Origin == o)
  lines <- paste0(sub$ind, " , ", apply(sub[, LOCI], 1, paste, collapse = " "))
  writeLines(lines, con)
}
close(con)
message("wrote ", genepop_file, " (", nrow(gp_body), " individuals, ",
        n_distinct(gp_body$Origin), " populations)")

# ---- 2. Run the null-allele estimator ----------------------------------

null_out <- file.path(paths$derived, "null_alleles_genepop.txt")
genepop::nulls(genepop_file, outputFile = null_out, verbose = FALSE)

# genepop::nulls() scribbles cmdline.txt and fichier.in into the working
# directory; remove them so they don't pollute the repository (run_all.R does
# the same sweep at the end, this covers running step 04 on its own).
file.remove(intersect(c("cmdline.txt", "fichier.in"), list.files()))

# Parse the confidence-interval summary table that Genepop prints at the end:
#
#   Locus      Population   estimate   bound    bound
#   cf020      <pop1>       0.0429     0.0264   0.0624
#              <pop2>       0.0226     0.0063   0.0448
#   ...
#
# Rows are: a locus name (optional) + population label + three numbers, or the
# text "(No info for CI)" in place of the two bounds.
raw <- readLines(null_out)

start <- grep("Frequency\\s+0\\.0?250\\s+0\\.9750|estimate\\s+bound", raw)
tail_block <- if (length(start)) raw[(min(start)):length(raw)] else raw

parse_ci_table <- function(lines, loci) {
  cur_locus <- NA_character_
  rows <- list()
  for (ln in lines) {
    lc <- str_extract(ln, "cf\\d+")
    if (!is.na(lc)) cur_locus <- lc
    nums <- suppressWarnings(as.numeric(
      unlist(regmatches(ln, gregexpr("\\d+\\.\\d+", ln)))))
    if (length(nums) >= 1 && !is.na(cur_locus)) {
      no_ci <- grepl("No info for CI", ln)
      rows[[length(rows) + 1]] <- tibble(
        Locus    = cur_locus,
        null_est = nums[1],
        ci_low   = if (no_ci || length(nums) < 2) NA_real_ else nums[2],
        ci_high  = if (no_ci || length(nums) < 3) NA_real_ else nums[3]
      )
    }
  }
  bind_rows(rows) %>% filter(Locus %in% loci)
}

null_long <- parse_ci_table(tail_block, LOCI)

# label the two populations by origin (Genepop names a population after its
# first individual; the file was written cultivar block first, native second)
pop_labels <- unique(null_long$Locus)  # placeholder to keep order
null_long <- null_long %>%
  group_by(Locus) %>%
  mutate(population = if_else(row_number() == 1, "cultivar", "native")) %>%
  ungroup()

null_tbl <- null_long %>%
  mutate(Locus = factor(Locus, levels = LOCI)) %>%
  arrange(Locus, population) %>%
  select(Locus, population, null_est, ci_low, ci_high)

message("\nnull-allele frequency estimates (ML, 95% CI):")
print(as.data.frame(null_tbl), row.names = FALSE)

write_csv(null_tbl, file.path(paths$results, "null_alleles.csv"))

# ---- 3. Reduced marker set -------------------------------------------
# A locus is dropped from the reduced set when its null-allele estimate has a
# 95% CI excluding zero in BOTH origin classes (consistent evidence of null
# alleles). Loci flagged in only one class, or in neither, are kept.

flagged_by_pop <- null_tbl %>%
  mutate(flagged = !is.na(ci_low) & ci_low > 0) %>%
  group_by(Locus) %>%
  summarise(n_flagged = sum(flagged), .groups = "drop")

flagged <- flagged_by_pop %>% filter(n_flagged == 2) %>% pull(Locus) %>% as.character()
reduced <- setdiff(LOCI, flagged)

writeLines(reduced, file.path(paths$derived, "reduced_marker_set.txt"))
message("\nloci with null-allele evidence in both origin classes: ", paste(flagged, collapse = ", "))
message("reduced marker set (", length(reduced), " loci): ", paste(reduced, collapse = ", "))

# ---- 4. Exclusion probabilities --------------------------------------

pop <- read_population(file.path(paths$data, "genotypes.csv"),
                       type = "separated", locus.columns = 6:14)

pe <- genetic_diversity(pop, mode = "Pe")
names(pe)[2] <- "Pe"

combined_pe <- function(loci) 1 - prod(1 - pe$Pe[pe$Locus %in% loci])

excl <- tibble(
  set   = c("full (9 loci)", paste0("reduced (", length(reduced), " loci)")),
  loci  = c(paste(LOCI, collapse = ", "), paste(reduced, collapse = ", ")),
  P_excl_multilocus = c(combined_pe(LOCI), combined_pe(reduced))
)

write_csv(bind_rows(
  pe %>% mutate(set = "per-locus") %>% select(set, Locus, Pe),
  excl %>% transmute(set, Locus = "multilocus", Pe = P_excl_multilocus)
), file.path(paths$results, "exclusion_probability.csv"))

message("\nper-locus exclusion probability (gstudio mode = Pe):")
print(as.data.frame(pe %>% mutate(Pe = round(Pe, 4))), row.names = FALSE)
message("\ncombined multilocus exclusion probability:")
print(as.data.frame(excl %>% mutate(P_excl_multilocus = signif(P_excl_multilocus, 6))), row.names = FALSE)
