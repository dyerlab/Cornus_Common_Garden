# Early fitness consequences of cultivar gene escape in flowering dogwood (*Cornus florida*, Cornaceae L.)

Jane Remfert, Andrew Eckert, Rodney Dyer 

Open-pollinated dogwood seedlings from managed-landscape ("cultivar") and second-growth-forest ("native") mothers were grown together for two years in a deeply shaded forest common garden and scored for survival, growth, and biomass. Cultivar-origin seedlings were **not** disadvantaged — survival was high and did not vary with origin or genetic admixture, and origin explained little trait variation — so early viability selection is unlikely to stop escaped cultivar alleles from establishing, and longer-term introgression into wild populations deserves attention.

This repository holds the manuscript (`manuscript.md`, `supplemental.md`), the data, and the `R` code that reproduces every table and figure. 

An online version of the manuscript is [here](https://dyerlab.github.io/Cornus_Common_Garden/manuscript.html) and supplemental materials are [here](https://dyerlab.github.io/Cornus_Common_Garden/supplemental.html)

The rest of this file is a standalone summary of the study; repository and reproduction details are at the bottom.

---

## Competing hypotheses

| | Prediction for cultivar-origin seedlings in a wild understory |
| :-- | :-- |
| **Maladaptation** | Selection under cultivation (showy bloom, vigor in open high light) relaxed selection for shade tolerance / understory competitiveness / early growth → reduced survival, growth, and/or biomass vs. native-origin seedlings. Early viability selection would then weed out escaped alleles. |
| **Cultivar persistence** | Cultivation-selected vigor carries over as a general performance advantage → cultivar-origin seedlings do at least as well as native. Early viability selection does *not* restrain gene escape; sets the stage for introgression. |

**Result: supports cultivar persistence** — no evidence for maladaptation.

![](https://dyerlab.github.io/Cornus_Common_Garden/media/fig_offspring_heritability.png)
*Heritability re-expressed under realized within-mother offspring relatedness, for native and cultivar family arrays across each measured trait. The violin shows bootstrap uncertainty in heritability at each origin’s median pairwise sibling relatedness (maternal families resampled with replacement within origin); the point marks that median-relatedness estimate, and the vertical line spans the range implied by a strict half-sib to strict full-sib assumption.*


## Study design & data

- **569 seedlings** from **36 open-pollinated maternal families**, seeds collected fall 2017 around Richmond, VA.
  - **20 cultivar-origin families / 395 seedlings** — mothers in residential/commercial managed landscapes.
  - **16 native-origin families / 174 seedlings** — mothers in second-growth forest (Pocahontas State Park + VCU Rice Rivers Center).
- **Common garden**: understory of a second-growth hardwood forest at the VCU Rice Rivers Center, deer exclosure, seedlings in pots placed at random. Two growing seasons; survival censuses at peak growth and end of season each year; final harvest/census **2019-09-18**. 468 alive at the end.
- **Traits**: germination timing (developmental "days"), survival, plant height (PH), stem diameter (SD), number of leaves (NL), above- and below-ground biomass (AGB, BGB), AGB partitioned into leaf (LB) and stem (SB) biomass. PH/SD/NL measured at four censuses (repeated measures); biomass at harvest.
- **Genetics**: 9 cultivar-diagnostic microsatellite loci (Wadl et al. 2008; Wang et al. 2009) on maternal trees, offspring, and ~19 known cultivar voucher lines. Used to build **Admixture Percentage (AP)** = proportion of a seedling's alleles shared with the best-matching cultivar voucher (all cultivars pooled).
- **Two ancestry classifications, analyzed in parallel**: (1) *origin* = maternal-tree location (cultivar vs. native); (2) *AP* = continuous genetic admixture.
- **Covariates**: scaled developmental time since germination (proxy for germination timing); scaled maternal-family mean seed weight (provisioning).

## Analysis & key assumptions

- **Mixed models** (`lme4` / `lmerTest`, R 4.5.3). Gaussian LMMs for ln-transformed growth traits; **binomial GLMM (logit)** for survival.
  - Single-time: `y ~ days + seed_weight + (1|family) [+ (1|origin)]`
  - Repeated-measures: adds a fixed `census` term and `(1|plant)`.
  - Admixture models: replace the `origin` random effect with **AP as a continuous fixed covariate**.
- **Model selection** by AIC/AICc + Akaike weights: fixed-only vs. +family vs. +family+origin (or family ± AP). Models within ΔAIC ≤ 2 treated as comparable; both decompositions reported when there is no clear winner. AIC fits by ML, variance components by REML.
- **Heritability**: half-sib design, **h² = 4·V_family / V_P**. Survival h² on the latent liability scale with residual variance fixed at π²/3 ≈ 3.29 (Gilmour et al. 1985).
- **Uncertainty**: 1,000-replicate parametric bootstrap CIs (`bootMer`); random-effect significance by LRT (`ranova`, boundary-corrected for the GLMM). An h² is called "significant" when its bootstrap CI excludes zero.
- **Diagnostics / sensitivity**: residual-vs-fitted, Q–Q, Cook's D; `DHARMa` simulated residuals for the GLMM (no meaningful violations). Every model refit (a) dropping individual outliers |resid z| > 3 and (b) dropping the single most influential maternal family.
- **Key assumptions / caveats**: random effects independent and normally distributed; families treated as pure half-sibs (paternity unknown; open-pollinated); missing seed weights (1 cultivar mom / 14 seedlings; 8 native moms / 61 seedlings) imputed with the origin mean → conservative for detecting seed-weight effects; cultivar voucher panel incomplete (~19 lines) and diagnostic alleles not proven exclusive to cultivars; two years only — no data on reproduction or later-generation hybrid breakdown; rearing/nursery location not modeled.

## Findings

**1. Native collection-site effect (can the two native sources be pooled?)** LRT of a `site` fixed term added to the family model, per trait, native-origin seedlings only (Pocahontas 11 families/147 vs. Rice Rivers Center 5/27). No site effect on survival or any trait (all *p* = 0.15–0.68; Tab_SiteEffect) → native sources pooled into one *native* category.

**2. Marker diversity & power.** Highly polymorphic — 5–26 alleles/locus in mothers, similar in offspring (Tab_MarkerDiversity). Null-allele frequencies low everywhere (all point estimates < 0.07) but CIs excluded zero at four loci (cf020, cf585, cf597, cf701) in both origin classes (Tab_NullAlleles). Combined exclusion probability high for the full nine-locus panel (*P*excl > 0.9999) and a reduced five-locus panel (*P*excl = 0.999) → the panel keeps strong discriminating power.

**3. Does origin/admixture affect survival?** Binomial GLMM; model selection ± origin, ± family; separate AP models; native/cultivar subset fits.
- Overall survival **~82 %** (native 79.3 %, cultivar 83.5 %); mortality roughly constant through time (Fig_Survival).
- Adding origin did **not** improve fit over family alone (family ΔAICc = 0.0, *w*AICc = 0.54; origin+family ΔAICc = 2.0); origin variance ≈ 0 (LRT χ² ≈ 0, *p* ≈ 0.5).
- AP not a significant survival predictor (*p* = 0.88).
- **Family split**: the family random effect improved fit in the **cultivar** subset (ΔAICc = 1.5; LRT χ² = 3.5, df = 1, *p* = 0.030) but **not** the native subset. Cultivar-family survival h² = 0.51, very imprecise (95 % CI [0, 1.07]); ~13 % of latent-scale variance among cultivar families. Family conditional modes (BLUPs) show heterogeneous survival among cultivar families (Fig_BLUPs).

**4. Does origin/admixture affect germination timing?** Origin fit as a random effect (the germination-timing covariate is the response, so omitted). Origin did not improve fit (ΔAICc = 2.0; LRT *p* = 0.48). But germination timing itself is **heritable** among families: h² = 0.42 [0.15, 0.72] (Tab_GerminationTiming).

**5. Does origin/admixture affect growth traits (PH, SD, NL, AGB, BGB, LB, SB)?**
- `family` and `family + origin` within ΔAICc ≤ 2 for PH, SD, NL, AGB, BGB, SB → origin adds little beyond family (Tab_GrowthAIC). The bootstrap CI for V_origin includes zero for **every** trait (Tab_GrowthHeritability).
- **Leaf biomass** is the one exception where the origin-inclusive model is clearly favored (ΔAICc = 4.9, *w*AICc = 0.92) — but its V_origin CI still includes zero.
- Repeated-measures models show the same pattern (Tab_RepeatedAIC / Tab_RepeatedHeritability).
- AP never improved fit over family alone; the AP coefficient is not significant for any trait (smallest *p* = 0.11, plant height).
- Origin conditional modes (from origin+family models) are consistently **positive for cultivar, negative for native** for PH, AGB, BGB, SB → if anything, cultivar-origin seedlings grow slightly *better*, not worse.

**6. Family-level variation & heritability in growth (within-origin).** Data split by origin; family model fit within each; same sensitivity passes.
- **Native subset**: fixed-effects-only model best for every growth trait; family-variance CIs include zero throughout.
- **Cultivar subset**: family model best-supported for every growth trait except LB; significant, modest family variance for **PH (11.9 % of V_P), AGB (8.1 %), SB (11.0 %)**.
- **Significant narrow-sense h² (cultivar subset)**: PH 0.47 [0.11, 0.93], SB ≈ 0.44 [0.07, 0.89], AGB 0.32 [0.03, 0.72] — moderate, comparable to other tree species.
- **Sensitivity**: robust to removing outlier individuals (the family effect if anything strengthens; biggest shift BGB h² 0.23 → 0.36). Removing the single most influential family (**urb09** for most traits) matters more — e.g. SB: V_family 0.046 → 0.016, h² 0.51 → 0.19, LRT still *p* = 0.022. The heritable signal is real but partly leveraged by one family.

**7. Covariate effects.**
- **Earlier germination** is a strong positive predictor of survival and of every growth trait except NL, in both origin classes (Tab_CovariateEffects) → a developmental head start is a general early advantage.
- **Maternal seed weight** is weaker and inconsistent: negatively associated with survival overall and in the native subset; among cultivar seedlings negatively associated with AGB, BGB, LB (smaller-seeded families accumulated more biomass); not significant for growth in native seedlings.

**8. Marker-panel validation — do native and cultivar maternal trees form distinguishable genetic groups?** Three complementary checks on the 58 unrelated maternal trees, run because AP (above) looks at each locus independently and found only a diffuse difference between groups. No evidence of non-random multilocus allelic association within either origin group (index of association; native *I*_A = −0.071, *p* = 0.717; cultivar *I*_A = 0.063, *p* = 0.235; Tab_LinkageDisequilibrium). Discriminant analysis of principal components (DAPC) on the full nine-locus genotype does not separate the groups better than chance (56.061% cross-validated classification accuracy vs. 55.172% no-information rate; permutation *p* = 0.237; Tab_DAPC). And native mothers are not more closely related to one another than cultivar mothers are (Tab_Relatedness), so the origin-dependent family-variance asymmetry (Finding 6) is unlikely to simply reflect unequal kinship structure. Together these suggest the cultivar-diagnostic loci may be less exclusively informative in this sample than assumed, consistent with the broad AP overlap noted in Finding 3.

## Bottom line

- No evidence that cultivar ancestry (by maternal location *or* genetic admixture) hurts early survival or growth in a wild understory; cultivar-origin seedlings trend slightly better for growth.
- Early viability selection is therefore **unlikely to be self-limiting** for cultivar gene escape → potential for longer-term introgression into wild *C. florida*.
- Meaningful heritable variation *among cultivar families* for height and biomass (and heterogeneous cultivar-family survival) → *which* cultivar families escape may matter; there is variation for selection to act on.
- Caveats bounding the claim: two growing seasons only (no reproduction, no later-generation hybrid breakdown), an incomplete cultivar voucher panel, and no disease challenge — cultivar disease-resistance alleles (e.g. to *Discula destructiva*) could even benefit wild populations.

---

## Repository

### Contents

| Path | Contents |
|---|---|
| `manuscript.md` | Full manuscript: main text, Literature Cited, tables, figure legends |
| `supplemental.md` | Supplementary Materials — statistical model structure, variance decomposition, marker validation (linkage disequilibrium, DAPC, relatedness), diagnostics |
| `data/*.csv` | Input data (see **Data files** below) |
| `media/` | Manuscript figures (`fig_survival`, `fig_blups`, `fig_lmm_qq`, `fig_within_family_variance`) and per-model diagnostic panels under `diagnostics/` |
| `R/` | Analysis pipeline — one script per step, plus `run_all.R`; see `R/README.md` |

`data/derived/` (analysis intermediates) and `data/results/` (result tables + `manuscript_tables.md`) are both regenerated by the pipeline from `data/*.csv` and are not tracked; the result tables also appear, rendered, in `manuscript.md`. Non-anonymized source files, the de-anonymization key, and project notes are kept locally in `SENSITIVE/` and are never published.

### Reproducing the analysis

From the repository root:

```sh
Rscript R/run_all.R           # full run, 1000-replicate bootstraps (~15-20 min)
NSIM=50 Rscript R/run_all.R    # fast dry run
```

This rebuilds `data/derived/`, `data/results/`, and `media/` from the input data. Only `02_admixture.R` needs the non-anonymized cultivar-voucher panel (kept in `SENSITIVE/`, not published); it detects a missing `SENSITIVE/` and skips, and its one output — `data/admixture.csv` — is committed so every downstream step still runs. Requirements and a per-script guide are in `R/README.md`.

### Data files

All maternal trees are identified by anonymized code: `urb**` = cultivar origin (managed landscape), `nat**` = native origin (second-growth forest).

| File | One row per | Key columns |
|---|---|---|
| `germination.csv` | seedling | `mom`, `id`, `origin`, germination date, developmental time, seed weight |
| `phenotypes.csv` | seedling × census | height, stem diameter, leaf number, stem-base/cotyledon measures, date |
| `biomass.csv` | seedling | above- and below-ground biomass, leaf and stem biomass (harvest) |
| `seed_weight.csv` | maternal family | mean weighed seed mass |
| `native_mom_sites.csv` | native maternal family | collection site (Pocahontas State Park / VCU Rice Rivers Center) |
| `genotypes.csv` | individual | 9 cultivar-diagnostic microsatellite loci for mothers (`OffID = 0`) and offspring |
| `admixture.csv` | genotyped seedling | Admixture Percentage — % of alleles shared with the best-matching cultivar voucher. Generated by `02_admixture.R` from a non-public voucher panel; committed because it cannot be rebuilt from a clone. |
| `admixture_maternal.csv` | genotyped maternal tree | Admixture Percentage per mother (same method as `admixture.csv`, aggregated to the 58 maternal trees). Also generated by `02_admixture.R`; committed for the same reason. |

### Citation

If you use the data or code, please cite the paper:

- Remfert J, Eckert A, Dyer R. Early fitness consequences of cultivar gene escape in flowering dogwood (*Cornus florida* L.). A DOI will be minted from the tagged release.

---

Any mistakes in these data are the sole responsibility of [R. Dyer](mailto:rjdyer@vcu.edu).
