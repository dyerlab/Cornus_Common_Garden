# 00_setup.R
# Shared configuration for the common-garden analysis pipeline.
# Sourced at the top of every analysis script (01_*.R ... ).
#
# It defines: project paths, the random-number seed, and a list of the
# CRAN packages each downstream script expects. It does not fit any models.

# ---- Paths -------------------------------------------------------------------
# All scripts are run with the repository root as the working directory
# (e.g. `Rscript R/01_assemble_data.R` from the project folder).

paths <- list(
  data      = "data",              # inputs: 6 anonymized CSVs + admixture.csv (committed AP values)
  derived   = "data/derived",      # pipeline intermediates the model scripts consume — git-ignored
  sensitive = "SENSITIVE",         # non-anonymized source, keys, notes, legacy repos — git-ignored, local only
  results   = "data/results",      # generated result tables (CSV) + manuscript_tables.md — git-ignored
  figures   = "media"              # generated manuscript figures (+ diagnostics/ subtree) — committed
)

for (p in c(paths$derived, paths$results, paths$figures)) {
  if (!dir.exists(p)) dir.create(p, recursive = TRUE)
}

# ---- Reproducibility --------------------------------------------------------
# Used for parametric bootstraps and any simulation-based diagnostics.
SEED <- 20240907

# ---- Packages --------------------------------------------------------------
# Checked here so a missing package fails early with a clear message rather
# than midway through a long model fit.

required_packages <- c(
  "readr", "dplyr", "tidyr", "stringr", "purrr", "tibble",  # data handling
  "lme4", "lmerTest", "MuMIn", "boot", "nlme",              # mixed models
  "DHARMa",                                                 # GLMM diagnostics
  "gstudio", "hierfstat", "genepop",                        # population genetics
  "ggplot2"                                                 # figures
)

check_packages <- function(pkgs = required_packages) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Missing required packages: ", paste(missing, collapse = ", "),
         "\nInstall them with install.packages(c(",
         paste(sprintf('\"%s\"', missing), collapse = ", "), "))")
  }
  invisible(TRUE)
}

# ---- Study genotypes ------------------------------------------------------
# data/genotypes.csv holds every genotyped tree in the larger collection (58
# maternal trees + 828 offspring), but only some of those mothers contributed
# families to THIS common garden. Every genetic analysis in the manuscript is
# restricted to the study material:
#   * maternal trees (OffID == 0) whose family was planted in the garden
#     (mom present in data/germination.csv: 36 families, 28 of them genotyped)
#   * offspring that are garden seedlings (mom + id present in germination.csv)
# The filtered rows are written verbatim (same text, same columns) to
# data/derived/genotypes_study.csv so gstudio::read_population() and
# read_csv() parse it identically to the original. Returns that path.
study_genotypes_file <- function() {
  out  <- file.path(paths$derived, "genotypes_study.csv")
  germ <- utils::read.csv(file.path(paths$data, "germination.csv"),
                          colClasses = "character")
  garden_moms <- unique(germ$mom)
  garden_kids <- paste(germ$mom, germ$id)

  lines  <- readLines(file.path(paths$data, "genotypes.csv"))
  fields <- strsplit(lines[-1], ",", fixed = TRUE)
  id     <- vapply(fields, `[`, "", 1)
  offid  <- vapply(fields, `[`, "", 2)
  keep   <- ifelse(offid == "0", id %in% garden_moms,
                   paste(id, offid) %in% garden_kids)

  writeLines(c(lines[1], lines[-1][keep]), out)
  message("study genotypes: ", sum(keep & offid == "0"), " maternal trees, ",
          sum(keep & offid != "0"), " garden seedlings (of ", length(keep),
          " genotyped individuals in data/genotypes.csv)")
  out
}

# ---- Small utilities ------------------------------------------------------
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

# z-score that ignores NA, returned as a plain numeric vector
zscore <- function(x) as.vector(scale(x, center = TRUE, scale = TRUE))

message("setup loaded: project paths, seed = ", SEED)
