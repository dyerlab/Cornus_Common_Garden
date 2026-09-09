# 03_marker_diversity.R
# -----------------------------------------------------------------------------
# Per-locus genetic diversity of the nine cultivar-diagnostic microsatellite
# loci, computed separately for the maternal trees and their offspring arrays.
#
# Reports, per locus:
#   A   number of alleles
#   Ae  effective number of alleles = 1 / sum(p_i^2)
#   Ho  observed heterozygosity
#   He  expected heterozygosity (Nei's gene diversity)
#   Fis inbreeding coefficient (1 - Ho/He)
# and the multilocus mean / sum across loci.
#
# Input : data/genotypes.csv   (ID, OffID, Origin, coords, cf* genotypes)
#          OffID == 0 marks the maternal tree.
# Output: data/results/marker_diversity.csv        (Tab_MarkerDiversity)
#         data/results/marker_diversity_by_origin.csv  (supplementary split)
#
# Uses gstudio::genetic_diversity (Dyer 2009).
# -----------------------------------------------------------------------------

source("R/00_setup.R")
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr); library(gstudio) })

LOCI <- c("cf020", "cf125", "cf213", "cf273", "cf581", "cf585", "cf597", "cf634", "cf701")

# ---- 1. Load genotypes as gstudio locus objects -------------------------

pop <- read_population(file.path(paths$data, "genotypes.csv"),
                       type = "separated", locus.columns = 6:14)

pop <- pop %>%
  mutate(group = if_else(OffID == 0, "Maternal", "Offspring"))

message("maternal genotypes: ", sum(pop$group == "Maternal"),
        "   offspring genotypes: ", sum(pop$group == "Offspring"))

# ---- 2. Per-locus diversity for one set of individuals ----------------

diversity_table <- function(df, label) {
  modes <- c("A", "Ae", "Ho", "He", "Fis")
  out <- Reduce(function(acc, m) {
    d <- genetic_diversity(df, mode = m)
    names(d)[2] <- m
    if (is.null(acc)) d else left_join(acc, d, by = "Locus")
  }, modes, init = NULL)

  out <- out %>%
    mutate(Locus = factor(Locus, levels = LOCI)) %>%
    arrange(Locus)

  # multilocus summary row: A and Ae summed, heterozygosities averaged
  summary_row <- tibble(
    Locus = "Multilocus",
    A   = sum(out$A),
    Ae  = sum(out$Ae),
    Ho  = mean(out$Ho, na.rm = TRUE),
    He  = mean(out$He, na.rm = TRUE),
    Fis = mean(out$Fis, na.rm = TRUE)
  )

  bind_rows(mutate(out, Locus = as.character(Locus)), summary_row) %>%
    mutate(group = label)
}

maternal_div  <- diversity_table(filter(pop, group == "Maternal"),  "Maternal")
offspring_div <- diversity_table(filter(pop, group == "Offspring"), "Offspring")

marker_diversity <- bind_rows(maternal_div, offspring_div) %>%
  mutate(across(c(Ae, Ho, He, Fis), ~ round(.x, 3))) %>%
  select(group, Locus, A, Ae, Ho, He, Fis)

write_csv(marker_diversity, file.path(paths$results, "marker_diversity.csv"))

# ---- 3. Supplementary: same split additionally by origin -------------

by_origin <- pop %>%
  distinct(Origin) %>%
  pull(Origin) %>%
  purrr::set_names() %>%
  purrr::map_dfr(function(o) {
    bind_rows(
      diversity_table(filter(pop, group == "Maternal",  Origin == o), "Maternal"),
      diversity_table(filter(pop, group == "Offspring", Origin == o), "Offspring")
    ) %>% mutate(origin = o)
  }) %>%
  mutate(across(c(Ae, Ho, He, Fis), ~ round(.x, 3))) %>%
  select(origin, group, Locus, A, Ae, Ho, He, Fis)

write_csv(by_origin, file.path(paths$results, "marker_diversity_by_origin.csv"))

# ---- 4. Console summary ------------------------------------------------

am <- marker_diversity %>% filter(group == "Maternal", Locus != "Multilocus")
message("\nmaternal alleles per locus: ", min(am$A), " - ", max(am$A),
        "   (manuscript states 5-22)")
print(as.data.frame(marker_diversity), row.names = FALSE)
