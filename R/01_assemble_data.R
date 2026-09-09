# 01_assemble_data.R
# -----------------------------------------------------------------------------
# Build the analysis tables used by every downstream model script.
#
# Inputs  (data/, anonymized; one row per seedling-measurement or per mom):
#   germination.csv      mom, id, origin, cold_stratification, germination,
#                        days (developmental time: days from germination to a fixed
#                        reference near the end of the germination window; a higher
#                        value means the seed germinated earlier and had more growing
#                        time before measurement), cold_days, bilocular
#   phenotypes.csv       mom, id, origin, date, phenotype (PH/SD/NL/...), value, note
#   biomass.csv          mom, id, origin, phenotype (AGB/BGB/LB/SB), value
#   seed_weight.csv      mom, origin, mean_seed_weight, n_seeds
#   native_mom_sites.csv mom, site   (native mothers only)
#
# Outputs (data/derived/):
#   seedlings.csv          one row per seedling: identifiers, covariates
#   traits_harvest.csv     one row per seedling: end-of-experiment trait values
#   traits_repeated.csv    long: seedling x census x trait, for repeated-measures models
#   survival.csv           one row per seedling: alive at the final census
#   survival_by_census.csv proportion alive per census and origin, for the survival figure
#
# Conventions
#   * "family" = maternal tree (mom). Half-sib design: offspring of one mother.
#   * "origin" = native vs cultivar, from the maternal tree's location.
#   * Covariate `days_z`      = z-scored developmental time (higher = earlier germination).
#   * Covariate `seed_weight_z` = z-scored mean maternal seed weight, missing
#     values imputed with the origin-specific mean before scaling.
#   * Trait responses are ln-transformed for modelling; AGB and stem biomass
#     values of 0 are treated as missing (a 0 means the tissue was not recovered).
# -----------------------------------------------------------------------------

source("R/00_setup.R")
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr); library(purrr)
})

# ---- 1. Load inputs --------------------------------------------------------

germ    <- read_csv(file.path(paths$data, "germination.csv"),      show_col_types = FALSE)
pheno   <- read_csv(file.path(paths$data, "phenotypes.csv"),       show_col_types = FALSE)
biomass <- read_csv(file.path(paths$data, "biomass.csv"),          show_col_types = FALSE)
seedwt  <- read_csv(file.path(paths$data, "seed_weight.csv"),      show_col_types = FALSE)
sites   <- read_csv(file.path(paths$data, "native_mom_sites.csv"), show_col_types = FALSE)

# The germination table is the seedling registry: every seedling that was
# stratified, germinated, and planted into the garden has exactly one row here.
germ <- germ %>%
  filter(!is.na(id)) %>%
  distinct(mom, id, .keep_all = TRUE)

n_seedlings <- nrow(germ)
message("seedlings in registry: ", n_seedlings,
        "  (", sum(germ$origin == "cultivar"), " cultivar / ",
        sum(germ$origin == "native"), " native)")
message("maternal families: ", dplyr::n_distinct(germ$mom),
        "  (", dplyr::n_distinct(germ$mom[germ$origin == "cultivar"]), " cultivar / ",
        dplyr::n_distinct(germ$mom[germ$origin == "native"]), " native)")

# ---- 2. Seed weight: impute missing families, then z-score ----------------
# seed_weight.csv has one row per mother that was weighed. Families with no
# record get the mean of the other families of the same origin. Scaling is done
# on the seedling-level vector so the covariate is centred for the analysis set.

origin_mean_wt <- seedwt %>%
  group_by(origin) %>%
  summarise(origin_mean = mean(mean_seed_weight, na.rm = TRUE), .groups = "drop")

mom_wt <- germ %>%
  distinct(mom, origin) %>%
  left_join(seedwt %>% select(mom, mean_seed_weight), by = "mom") %>%
  left_join(origin_mean_wt, by = "origin") %>%
  mutate(
    seed_weight     = mean_seed_weight,
    seed_weight_imp = coalesce(mean_seed_weight, origin_mean),
    imputed_wt      = is.na(mean_seed_weight)
  ) %>%
  select(mom, seed_weight, seed_weight_imp, imputed_wt)

message("families with imputed seed weight: ", sum(mom_wt$imputed_wt),
        " (", paste(sort(mom_wt$mom[mom_wt$imputed_wt]), collapse = ", "), ")")

# ---- 3. Seedling registry -------------------------------------------------

seedlings <- germ %>%
  transmute(
    mom, id, origin,
    bilocular  = as.integer(bilocular %||% 0),
    dev_time   = days,     # developmental time; higher = earlier germination
    germination_date = germination
  ) %>%
  left_join(sites, by = "mom") %>%                     # site: native mothers only
  left_join(mom_wt, by = "mom") %>%
  mutate(
    days_z         = zscore(dev_time),
    seed_weight_z  = zscore(seed_weight_imp)
  ) %>%
  arrange(origin, mom, id)

write_csv(seedlings, file.path(paths$derived, "seedlings.csv"))

# ---- 4. End-of-experiment (single-measurement) traits --------------------
# Height, stem diameter, and leaf number are taken at their last census
# (2019-09-18, 2019-09-17, 2019-09-08 respectively). Biomass fractions come
# from the harvest at the end of the second season.

last_pheno <- function(trait_code, on_date) {
  pheno %>%
    filter(phenotype == trait_code, date == as.Date(on_date)) %>%
    select(mom, id, value)
}

harvest_traits <- seedlings %>%
  select(mom, id, origin) %>%
  left_join(last_pheno("PH", "2019-09-18") %>% rename(height_cm = value),  by = c("mom", "id")) %>%
  left_join(last_pheno("SD", "2019-09-17") %>% rename(stemd_mm  = value),  by = c("mom", "id")) %>%
  left_join(last_pheno("NL", "2019-09-08") %>% rename(leaves_n  = value),  by = c("mom", "id")) %>%
  left_join(
    biomass %>%
      mutate(phenotype = recode(phenotype, LB = "leaf_biomass_g", SB = "stem_biomass_g",
                                AGB = "agb_g", BGB = "bgb_g")) %>%
      pivot_wider(names_from = phenotype, values_from = value),
    by = c("mom", "id", "origin")
  ) %>%
  mutate(
    # a recorded 0 for above-ground or stem biomass means tissue was not
    # recovered, not a true zero mass
    agb_g          = na_if(agb_g, 0),
    stem_biomass_g = na_if(stem_biomass_g, 0),
    across(c(height_cm, stemd_mm, leaves_n, agb_g, bgb_g, leaf_biomass_g, stem_biomass_g),
           list(ln = ~ log(.x)), .names = "{.col}_ln")
  )

traits_harvest <- seedlings %>%
  left_join(harvest_traits %>% select(-origin), by = c("mom", "id"))

write_csv(traits_harvest, file.path(paths$derived, "traits_harvest.csv"))

# ---- 5. Repeated-measures long table ------------------------------------
# Height, stem diameter and leaf number were recorded at the start and end of
# each of the two growing seasons. These are the censuses used in the
# repeated-measures models (one row per seedling x census x trait).
#
# The first stem-diameter census (2018-05-30) is the diameter measured just
# below the cotyledons (stored as phenotype "SD_cotyledon"); the three later
# censuses are the standard "SD" measurement at the cotyledon scar.

pheno_sd <- pheno %>%
  mutate(phenotype = if_else(phenotype == "SD_cotyledon" & date == as.Date("2018-05-30"),
                             "SD", phenotype))

repeated_censuses <- tribble(
  ~trait,      ~date,
  "PH", "2018-06-24", "PH", "2018-10-20", "PH", "2019-05-02", "PH", "2019-09-18",
  "SD", "2018-05-30", "SD", "2019-02-18", "SD", "2019-05-02", "SD", "2019-09-17",
  "NL", "2018-07-03", "NL", "2018-10-20", "NL", "2019-07-09", "NL", "2019-09-08"
) %>% mutate(date = as.Date(date))

traits_repeated <- pheno_sd %>%
  filter(phenotype %in% c("PH", "SD", "NL")) %>%
  inner_join(repeated_censuses, by = c("phenotype" = "trait", "date")) %>%
  transmute(mom, id, origin, trait = phenotype, date, value) %>%
  group_by(trait) %>%
  mutate(census = factor(as.integer(factor(date)))) %>%   # 1..4 within each trait
  ungroup() %>%
  left_join(seedlings %>% select(mom, id, days_z, seed_weight_z, site), by = c("mom", "id")) %>%
  mutate(value_ln = log(value)) %>%
  arrange(trait, date, origin, mom, id)

write_csv(traits_repeated, file.path(paths$derived, "traits_repeated.csv"))

# ---- 6. Survival ------------------------------------------------------------
# A seedling counts as alive at a census if a height was recorded for it at
# that census. The survival endpoint is being alive at the final census
# (2019-09-18). The per-census proportions feed the survival figure.

height_censuses <- pheno %>%
  filter(phenotype == "PH") %>%
  distinct(date) %>%
  arrange(date) %>%
  pull(date)

alive_grid <- tidyr::expand_grid(
  seedlings %>% select(mom, id, origin),
  date = height_censuses
) %>%
  left_join(
    pheno %>% filter(phenotype == "PH") %>% transmute(mom, id, date, measured = 1L),
    by = c("mom", "id", "date")
  ) %>%
  mutate(alive = coalesce(measured, 0L))

final_date <- max(height_censuses)

survival <- seedlings %>%
  left_join(
    alive_grid %>% filter(date == final_date) %>% select(mom, id, alive_final = alive),
    by = c("mom", "id")
  )

write_csv(survival, file.path(paths$derived, "survival.csv"))

survival_by_census <- alive_grid %>%
  group_by(origin, date) %>%
  summarise(n = n(), alive = sum(alive), .groups = "drop") %>%
  mutate(
    prop_alive = alive / n,
    se         = sqrt(prop_alive * (1 - prop_alive) / n)
  )

write_csv(survival_by_census, file.path(paths$derived, "survival_by_census.csv"))

# ---- 7. Summary to console ----------------------------------------------

message("\n--- assembled tables ---")
message("seedlings.csv           : ", nrow(seedlings), " rows")
message("traits_harvest.csv      : ", nrow(traits_harvest), " rows, ",
        sum(!is.na(traits_harvest$agb_g)), " with AGB")
message("traits_repeated.csv     : ", nrow(traits_repeated), " rows (",
        paste(sprintf("%s=%d", c("PH","SD","NL"),
                      c(sum(traits_repeated$trait == "PH"),
                        sum(traits_repeated$trait == "SD"),
                        sum(traits_repeated$trait == "NL"))), collapse = ", "), ")")
message("survival.csv            : ", nrow(survival), " rows, ",
        sum(survival$alive_final), " alive at ", format(final_date),
        " (", round(100 * mean(survival$alive_final), 1), "%)")
message("  by origin: ",
        paste(sprintf("%s %.1f%%",
                      tapply(survival$alive_final, survival$origin, length) |> names(),
                      100 * tapply(survival$alive_final, survival$origin, mean)),
              collapse = "  "))
