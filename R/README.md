# Analysis pipeline

Reproduces every table and figure in the manuscript from the anonymized data in
`../data/`. Run from the repository root:

```sh
Rscript R/run_all.R           # full run (bootstrap NSIM = 1000, ~15-20 min)
NSIM=50 Rscript R/run_all.R    # fast dry run
```

`run_all.R` runs steps 01–12 in order, then deletes any scratch files
(`cmdline.txt`, `fichier.in`, `Rplots.pdf`) tools leave in the working directory.

## Layout

| Path | Contents |
|---|---|
| `../data/*.csv` | inputs — 6 anonymized collected tables (genotypes, germination, phenotypes, biomass, seed weight, native mom sites) plus `admixture.csv` (see below) |
| `../data/derived/` | pipeline intermediates the model scripts consume — **git-ignored**, rebuilt every run |
| `../data/results/` | one CSV per result set, plus `manuscript_tables.md` — **git-ignored**, rebuilt every run (also rendered into `manuscript.md`) |
| `../media/` | `fig_survival`, `fig_blups`, `fig_within_family_variance`, and per-model diagnostic panels under `diagnostics/` — **committed** so `manuscript.md` renders without R |
| `../SENSITIVE/` | non-anonymized source (`data_raw/`), the anonymization key, notes, legacy repos — **git-ignored, local only, never published** |

`data/derived/` and `data/results/` are both pipeline output and a fresh checkout
rebuilds them from `data/*.csv` with `Rscript R/run_all.R`. The one committed
generated file is **`data/admixture.csv`** (per-seedling Admixture Percentage):
`02_admixture.R` builds it from the non-public plate7 voucher panel in
`SENSITIVE/`, so a clone cannot regenerate it — step 02 detects a missing
`SENSITIVE/` and skips, and steps 05–12 use the committed file.

## Scripts

| Script | Produces |
|---|---|
| `00_setup.R` | paths, seed, package check (sourced by the others) |
| `helpers.R` | shared model machinery: variance components, half-sib h², bootstrap CIs, AIC tables, LRT, influence measures |
| `01_assemble_data.R` | `seedlings`, `traits_harvest`, `traits_repeated`, `survival`, `survival_by_census` |
| `02_admixture.R` | `data/admixture.csv` — Admixture Percentage (AP) per genotyped seedling (committed; needs `SENSITIVE/`, skips without it) |
| `03_marker_diversity.R` | Tab_MarkerDiversity |
| `04_null_alleles.R` | Tab_NullAlleles, exclusion probabilities |
| `05_survival.R` | Tab_Survival, Fig_BLUPs data, survival GLMM diagnostics |
| `06_single_time.R` | Tab_GrowthAIC, Tab_GrowthHeritability, Tab_OriginHeritability, Tab_CovariateEffects |
| `07_repeated_measures.R` | Tab_RepeatedAIC, Tab_RepeatedHeritability |
| `08_site_effect.R` | Tab_SiteEffect |
| `09_diagnostics.R` | Appendix S1 LMM diagnostic panels and summary |
| `10_figures.R` | Fig_Survival, Fig_BLUPs |
| `11_tables.R` | `manuscript_tables.md` — all Tab_* blocks in Markdown |
| `12_within_family_variance.R` | Discussion analysis: is within-family phenotypic variance larger for cultivar than native arrays? (`NPERM` permutation reps, default 9999) |

## Analysis choices

- **Origin** is the binary native/cultivar classification from maternal-tree location.
- **Admixture Percentage (AP)** is a continuous fixed covariate (percentage of a
  seedling's alleles shared with its best-matching cultivar voucher). It replaces
  the earlier `prop_sa` / `prop_sa_cat` / `prop_sa_g1` variants.
- Initial rearing location is **not** a model term (tested previously, no effect).
- Growth-trait responses are ln-transformed; AGB and stem biomass values of 0 are
  treated as missing.
- Half-sib design: V_A = 4 · V_family; h² = V_A / V_P. Survival h² is on the latent
  (logit) scale with residual variance fixed at π²/3.
- Confidence intervals are 1000-replicate parametric bootstraps (`lme4::bootMer`).
- Random-effect significance is a boundary-corrected likelihood-ratio test
  (`lmerTest::ranova` for LMMs).

## Requirements

R (>= 4.2) with: `readr dplyr tidyr stringr purrr tibble lme4 lmerTest MuMIn boot
nlme DHARMa gstudio hierfstat genepop ggplot2`. `00_setup.R` checks these on load
and stops with an install line if any are missing. `nlme` ships with R;
`gstudio` is on GitHub (`dyerlab/gstudio`), the rest are on CRAN.
`package_versions.txt` records the versions used for the published run.
