# Analysis pipeline

Reproduces every table and figure in the manuscript from the anonymized data in
`../data/`. Run from the repository root:

```sh
Rscript R/run_all.R           # full run (bootstrap NSIM = 1000, ~15-20 min)
NSIM=50 Rscript R/run_all.R    # fast dry run
```

`run_all.R` runs steps 01–16 in order, then deletes any scratch files
(`cmdline.txt`, `fichier.in`, `Rplots.pdf`) tools leave in the working directory.

## Layout

| Path | Contents |
|---|---|
| `../data/*.csv` | inputs — 6 anonymized collected tables (genotypes, germination, phenotypes, biomass, seed weight, native mom sites) plus `admixture.csv` and `admixture_maternal.csv` (see below) |
| `../data/derived/` | pipeline intermediates the model scripts consume — **git-ignored**, rebuilt every run |
| `../data/results/` | one CSV per result set, plus `manuscript_tables.md` — **git-ignored**, rebuilt every run (also rendered into `manuscript.md`) |
| `../media/` | `fig_survival`, `fig_blups`, `fig_within_family_variance`, `fig_lmm_qq`, and per-model diagnostic panels under `diagnostics/` — **committed** so `manuscript.md` renders without R |
| `../SENSITIVE/` | non-anonymized source (`data_raw/`), the anonymization key, notes, legacy repos — **git-ignored, local only, never published** |

`data/derived/` and `data/results/` are both pipeline output and a fresh checkout
rebuilds them from `data/*.csv` with `Rscript R/run_all.R`. The two committed
generated files are **`data/admixture.csv`** (per-seedling Admixture
Percentage) and **`data/admixture_maternal.csv`** (per-maternal-tree Admixture
Percentage): `02_admixture.R` builds both from the non-public plate7 voucher
panel in `SENSITIVE/`, so a clone cannot regenerate them — step 02 detects a
missing `SENSITIVE/` and skips, and steps 05–12 use the committed files.

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
| `09_diagnostics.R` | `Fig_LMM_QQ` (composite Q-Q for all 10 LMMs) + per-model diagnostic panels and the Shapiro/heteroscedasticity summary — Supplementary Materials |
| `10_figures.R` | Fig_Survival, Fig_BLUPs |
| `11_tables.R` | `manuscript_tables.md` — all Tab_* blocks in Markdown |
| `12_within_family_variance.R` | Discussion analysis: is within-family phenotypic variance larger for cultivar than native arrays? (`NPERM` permutation reps, default 9999) |
| `13_multilocus_ld.R` | Tab_LinkageDisequilibrium — index of association (Ia, rbarD) within each origin group, on the 58 unrelated maternal trees |
| `14_dapc_origin.R` | Tab_DAPC — discriminant analysis of principal components, testing whether native/cultivar maternal trees separate on full multilocus genotype |
| `15_relatedness_check.R` | Tab_Relatedness — pairwise relatedness among maternal trees, checking whether native mothers are more closely related to one another than cultivar mothers |
| `16_offspring_relatedness.R` | Tab_HeritabilityRelatedness — within-mother offspring relatedness (same Nason/`gstudio` method as Tab_Relatedness, applied to each mother's own genotyped offspring, pooled by origin) and heritability re-expressed at the pooled median realized relatedness, bracketed by the half-sib (r = 0.25) and full-sib (r = 0.5) assumptions |

## Analysis choices

- **Origin** is the binary native/cultivar classification from maternal-tree location.  
- Admixture Percentage (AP)** is a continuous fixed covariate (percentage of a seedling's alleles shared with its best-matching cultivar voucher). It replaces the earlier `prop_sa` / `prop_sa_cat` / `prop_sa_g1` variants.
- Initial rearing location is **not** a model term (tested previously, no effect). Growth-trait responses are ln-transformed; AGB and stem biomass values of 0 are treated as missing.
- Half-sib design: V_A = 4 · V_family; h² = V_A / V_P. Survival h² is on the latent (logit) scale with residual variance fixed at π²/3. `16_offspring_relatedness.R` checks this half-sib assumption against realized marker-based relatedness among
  each mother's own genotyped offspring.
- Confidence intervals are 1000-replicate parametric bootstraps (`lme4::bootMer`).
- Random-effect significance is a boundary-corrected likelihood-ratio test (`lmerTest::ranova` for LMMs).

## Requirements

R (>= 4.2) with: `readr dplyr tidyr stringr purrr tibble lme4 lmerTest MuMIn boot nlme DHARMa gstudio hierfstat genepop ggplot2`. `00_setup.R` checks these on load and stops with an install line if any are missing. `nlme` ships with R; `gstudio` is on GitHub (`dyerlab/gstudio`), the rest are on CRAN. `package_versions.txt` records the versions used for the published run.



## ToDo

The following items need to be addressed:

1. In the materials and methods, blocks are mentioned at one point, but this is after the description of the garden, so were there formal blocks, like a randomized block design, or did you just randomize seedlings within the entire garden? 
2. Typically, I am asked to give the equations for the linear mixed models I fit in the main text. This helps everyone kind of figure out what you fit exactly, although your writing is quite clear. Not sure if you want to do that, but it’s a common ask from reviewers.
3. Although you point out that more description for some of the methods is located in the supplement, unless I missed it, it does not state in the main manuscript that you used bootstrapping to get confidence in intervals. That might be helpful since there are numerous ways to get these intervals. 
4. A common interpretation of low additive genetic variance in the context of what this manuscript is about often includes selection producing well-fitted populations to their local environments. Given that you argue the cultivars are sort of babied by humans, no such babying occurs for the trees in the natural stands, and they could’ve been under intense selection pressures for these early life-hood traits to meet a fitness optimum, and therefore have lower level levels of diversity. Maybe it’s the selectionist in me, but that might need to be mentioned.
5. Paragraph 1: (Pierce et al. 2008; Suchecki and Gibson 2008) These are about forest successional structure rather than stiltgrass specifically What I originally had but maybe should have placed the citations earlier in the sentence: Changes in forest successional composition toward more shade tolerant species increases disease severity and competition from exotic invasive species like *Microstegium vimineum* (Trin) A. Camus has limited seedling recruitment (Pierce, Bromer, and Rabenold 2008; Suchecki and Gibson 2008).
6. annual sales of $31 million came from USDA. 2020. “2019 Census of Horticultural Specialties.”
7. Paragraph 2: I'm not sure the Herrera 1981 citation fits. It's about endozoochory, but not in *C. florida* specifically
8. Methods: seeds were not weighed in 2017
9. I believe the molecular key in Wadl et al. 2008 is based on only 4 loci: CF213, CF581, CF585, CF597. CF634 and CF273 are in the paper but not in the key
10. Results: survival was higher in this experiment than in Redwine 2013



