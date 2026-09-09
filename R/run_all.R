# run_all.R
# -----------------------------------------------------------------------------
# Run the full analysis pipeline in order. From the repository root:
#
#     Rscript R/run_all.R
#
# Set the bootstrap replicate count with the NSIM environment variable
# (default 1000; use a small value for a fast dry run):
#
#     NSIM=50 Rscript R/run_all.R
#
# Each step reads from data/ or data/derived/ and writes to data/derived/,
# data/results/ or media/. Steps are independent given their inputs,
# but 02-11 all depend on 01, and 05-07 depend on 02.
# -----------------------------------------------------------------------------

steps <- c(
  "R/01_assemble_data.R",       # data/ -> data/derived/ analysis tables
  "R/02_admixture.R",           # Admixture Percentage per genotyped seedling
  "R/03_marker_diversity.R",    # Tab_MarkerDiversity
  "R/04_null_alleles.R",        # Tab_NullAlleles, exclusion probabilities
  "R/05_survival.R",            # Tab_Survival, Fig_BLUPs data, survival diagnostics
  "R/06_single_time.R",         # Tab_GrowthAIC, Tab_GrowthHeritability, Tab_OriginHeritability, Tab_CovariateEffects
  "R/07_repeated_measures.R",   # Tab_RepeatedAIC, Tab_RepeatedHeritability
  "R/08_site_effect.R",         # Tab_SiteEffect
  "R/09_diagnostics.R",         # Appendix S1 LMM diagnostics
  "R/10_figures.R",             # Fig_Survival, Fig_BLUPs
  "R/11_tables.R",              # assemble data/results/manuscript_tables.md
  "R/12_within_family_variance.R"  # Discussion: within-family variance, cultivar vs native
)

t0 <- Sys.time()
for (s in steps) {
  message("\n========== ", s, "  (", format(Sys.time(), "%H:%M:%S"), ") ==========")
  st <- Sys.time()
  res <- system2("Rscript", s, stdout = "", stderr = "")
  if (!identical(res, 0L)) stop("step failed: ", s)
  message("  ", s, " done in ", round(difftime(Sys.time(), st, units = "mins"), 1), " min")
}

# ---- cleanup ---------------------------------------------------------------
# Some tools scribble scratch files into the working directory that are not
# outputs and should not end up in the repository:
#   cmdline.txt, fichier.in  <- genepop::nulls() (step 04)
#   Rplots.pdf               <- a plotting call with no explicit device
# Their real outputs live under data/derived/, data/results/ and media/.
scratch  <- c("cmdline.txt", "fichier.in", "Rplots.pdf")
leftover <- scratch[file.exists(scratch)]
if (length(leftover)) {
  file.remove(leftover)
  message("cleaned up scratch files: ", paste(leftover, collapse = ", "))
}

message("\npipeline complete in ",
        round(difftime(Sys.time(), t0, units = "mins"), 1), " min")
