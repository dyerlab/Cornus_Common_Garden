# Appendix S1: Statistical Models, Marker Analyses, and Supplementary Results

This appendix provides full methodological detail and results supporting analyses described in the Methods and Discussion. The first part (Model Structure through Sensitivity Analyses) supports the *Cultivar Impact Analysis* subsection of the Methods, including model equations, variance-component and heritability derivations, diagnostic procedures, and sensitivity analyses. The second part (Multilocus Linkage Disequilibrium through Summary of Marker Analyses) provides full detail for two supplementary tests referenced in the Discussion: a test for non-random multilocus allelic association (linkage disequilibrium) between the native and cultivar groups, and a discriminant analysis of principal components (DAPC) testing whether native and cultivar maternal trees separate on their full multilocus genotype. Both marker-based tests complement the individual-allele-sharing approach (Admixture Percentage, AP) used in the main text, which considers each of the nine microsatellite loci independently. Annotated analysis code and the underlying data are available in the online supplementary materials.

## Model Structure

Two general mixed-effects model forms were used, depending on whether a trait was measured once per individual or repeatedly across the study.

Germination date (GD), seedling survival (SS), and the traits measured at harvest (AGB, BGB, leaf biomass, stem biomass) were modeled as:

$$
y_{ijk} = \mu + Days + Sw + Family_j + Origin_k + \varepsilon_{ijk}
$$

where $y_{ijk}$ is the trait value of the $i^{th}$ individual from the $j^{th}$ maternal family (*Family*) and the $k^{th}$ origin class (*Origin*). Fixed effects included the scaled developmental time since germination (*Days*, a proxy for germination timing that accounts for differences in developmental stage at measurement; omitted when the response variable was GD itself) and the scaled mean seed weight of the $j^{th}$ maternal family ($S_w$, to account for maternal provisioning effects). Random effects were maternal family (*Family*) and origin (*Origin*), with residual error $\varepsilon_{ijk}$. Survival was fit under this same general structure using a binomial error distribution and logit link rather than a Gaussian residual.

Height, stem diameter, and number of leaves were measured at the beginning and end of each growing season and were modeled with a repeated-measures extension:
$$
y_{ijkm} = \mu + Days + Sw + Time_m + Plant_i + Family_j + Origin_k + \varepsilon_{ijkm}
$$
which adds a fixed effect of census occasion (*Time*) and a random effect of individual plant (*Plant*) to account for repeated sampling of the same seedlings.

Random effects were assumed independent and normally distributed:
$$
Family_j \sim N(0,\sigma^2_{family}), \quad Origin_k \sim N(0,\sigma^2_{origin}), \quad Plant_i \sim N(0,\sigma^2_{plant}), \quad \varepsilon \sim N(0,\sigma^2)
$$
Cultivar ancestry was represented two different ways — the binary, maternal-location *Origin* classification (a random effect) and the continuous Admixture Percentage (AP; a fixed covariate) — and each trait was modeled separately under both, allowing us to compare inference across classification schemes.

Missing maternal seed weights ($n = 1$ maternal tree, 14 cultivar-origin seedlings; $n = 8$ maternal trees, 61 native-origin seedlings) were imputed using the origin-specific mean seed weight, which reduces variance in seed weight among affected families and should be considered a conservative choice with respect to detecting seed-weight effects.

## Model Comparison

For each trait we compared a fixed-effects-only model, a model adding *Family* as a random effect, and a model adding both *Family* and *Origin*, using Akaike Information Criterion (AIC and, for the small samples, AICc), and calculated Akaike weights across the candidate set. For the admixture-based analysis the candidate set was the *Family* model with and without AP as a fixed covariate. Models within $\Delta AIC \leq 2$ were treated as having comparable empirical support; where no single model was clearly preferred, we report variance and heritability estimates under both variance decompositions rather than selecting a single "best" model (see Results). Models compared by AIC were fit by maximum likelihood; final variance-component and heritability estimates were taken from models fit by restricted maximum likelihood (REML).

## Variance Partitioning and Heritability Estimation

Variance components for origin, maternal family, and residual variation were extracted from each fitted model using `VarCorr` (*lme4*; Bates et al. 2015). We assumed a half-sib design, under which additive genetic variance is estimated as four times the maternal-family variance component ($V_A = 4V_{family}$). Total phenotypic variance was estimated as $V_P = V_{origin} + V_{family} + V_{residual}$ across origin classes, and as $V_P(s) = V_{family}(s) + V_{residual}(s)$ within each origin class $s$. The proportion of phenotypic variance attributed to each random effect was calculated by dividing that component by the corresponding total.

Narrow-sense heritability was estimated as:

$$
h^2 = \frac{4V_{family}}{V_P}, \qquad h^2(s) = \frac{4V_{family}(s)}{V_P(s)}
$$

Survival heritability was estimated on the latent liability scale, with residual variance fixed at the distribution-specific value $\pi^2/3 = 3.29$ (Gilmour et al. 1985), following the same general form.

For the repeated-measures models, $V_P$ additionally included the among-plant variance component $V_{plant}$. Ninety-five percent confidence intervals for variance components and heritability estimates were obtained from 1,000 parametric bootstrap replicates using `bootMer` (*lme4*; Bates et al. 2015). The statistical significance of each random effect was evaluated using likelihood-ratio tests (`ranova`, *lmerTest*; Kuznetsova et al. 2017, for the linear mixed models; the equivalent boundary-corrected test for the survival GLMM).

## Model Diagnostics

For the linear mixed models, residual normality was assessed with normal quantile–quantile plots (Fig_LMM_QQ), homoscedasticity with residual-vs-fitted and scale–location plots, and group-level influence with leave-one-out Cook's distance on the fixed effects. Response variables were log-transformed where these plots indicated right-skew or non-constant variance. Formal Shapiro–Wilk tests rejected residual normality for most traits. At these sample sizes ($n$ = 465–2021) the test detects small departures, and the Q–Q plots show that the central distribution of residuals is well behaved for every model; the departures are confined to the lower tail and reflect mild left-skew — a small number of unusually low residuals — most pronounced for below-ground and leaf biomass. The outlier-removal sensitivity analysis below was run specifically to check that these tail observations do not drive the variance-component estimates: the estimates are qualitatively unchanged and, where they shift, the maternal-family signal strengthens rather than weakens. The correlation between residual scale and fitted value was small for every model ($|r| \leq 0.17$), indicating no meaningful heteroscedasticity. Residual-vs-fitted, scale–location, and Q–Q panels for every model are available in the code repository.

For the generalized linear mixed models (survival), fit was assessed using simulation-based residual diagnostics (1,000 simulated residuals per model via `simulateResiduals`, *DHARMa*; Hartig 2016). Tests for residual uniformity, dispersion, and outliers were non-significant for the full sample and both origin subsets (all $p > 0.2$), with the single exception of a marginal outlier test in the cultivar subset ($p = 0.02$); we detected no meaningful violations of GLMM assumptions. Family-specific conditional modes (BLUPs) were extracted from the cultivar-origin survival model using `ranef` (*lme4*; Bates et al. 2015) to visualize among-family variation in survival (Fig_BLUPs).

![](media/fig_lmm_qq.png)

**Fig_LMM_QQ.** Normal quantile–quantile plots of the standardised residuals from the top-supported model of each trait: the seven end-of-experiment single-measurement traits and the three repeated-measures traits (RM). Points on the dashed line indicate normally distributed residuals. The central distribution is well behaved for every model; departures are confined to the lower tail (mild left-skew), most visible for below-ground and leaf biomass. Shapiro–Wilk tests reject normality for most models (see text).

## Sensitivity Analyses

To assess the robustness of fixed- and random-effects inference, we refit each model after temporarily removing (a) individual-level outliers, identified by standardized residuals with $|z| > 3$, and (b) influential maternal families, identified by large Cook's distance values.

Removing outlier individual-level observations did not substantially change heritability estimates for most traits, and the significance of family-level random effects generally increased slightly; the largest shift was in below-ground biomass ($h^2 = 0.23 \rightarrow 0.36$).

Removing the single most influential maternal family (identified by the largest leave-one-out Cook's distance on the fixed effects) had a larger effect on variance-component and heritability estimates than removing outlier individuals. For most growth traits this family was urb09. For example, in the stem biomass model, removing urb09 reduced the family variance component ($V_{family} = 0.046 \rightarrow 0.016$) and heritability ($h^2 = 0.51 \rightarrow 0.19$), and weakened but did not eliminate support for the family random effect ($\chi^2 = 5.3$, df $= 1$, $p = 0.022$).

## Multilocus Linkage Disequilibrium (Index of Association)

AP and the comparison of the full nine-locus panel against the reduced six-locus panel of Wadl *et al.* (2008) (Results, Genetic Analysis) both examine individual alleles one locus at a time, and both found a diffuse pattern of allele sharing with no sharp separation between the native and cultivar groups. This leaves open a different question: do the nine loci co-vary with one another more than expected under random mating within a group, even though no single locus, or a simple sum across loci, separates the groups? Elevated multilocus linkage disequilibrium can arise from cryptic structure of exactly this kind, and is also the standard signature used to detect clonal or partially clonal reproduction — relevant here because cultivars are clonally propagated.

We used the standardized index of association ($\bar{r}_D$, and the related $I_A$; Agapow and Burt 2001), with a 999-permutation test of the null hypothesis of linkage equilibrium (no multilocus association), implemented in `poppr` (Kamvar, Tabima, and Grünwald 2014). This was computed for two datasets:

1. **Maternal trees only** (58 unrelated individuals) — the clean test. No two of these individuals are siblings, so any $I_A$ signal here cannot be a family or kinship artifact.
2. **Maternal trees pooled with their offspring** (886 individuals) — included for reference only. This dataset contains many maternal half- and full-sib families, and sibling groups are *expected* to show elevated multilocus association purely from shared parentage (a Wahlund-type effect of family structure), regardless of any cultivar signal. A positive result here, without a matching positive result in the maternal-only dataset, is not evidence of a native/cultivar effect.

Both datasets were analyzed separately within each origin group (native, cultivar) and for the pooled/total sample.

**Table S1.** Index of association ($I_A$) and standardized index of association ($\bar{r}_D$; Agapow and Burt 2001), with permutation-based *p*-values (999 permutations; `poppr`, Kamvar, Tabima, and Grünwald 2014).

| Dataset | Group | N | $I_A$ | *p* ($I_A$) | $\bar{r}_D$ | *p* ($\bar{r}_D$) |
|---|---|---|---|---|---|---|
| Maternal trees only (unrelated individuals) | Native | 26 | -0.071 | 0.717 | -0.009 | 0.718 |
| Maternal trees only (unrelated individuals) | Cultivar | 32 | 0.063 | 0.235 | 0.008 | 0.235 |
| Maternal trees only (unrelated individuals) | Total | 58 | -0.012 | 0.542 | -0.002 | 0.542 |
| Pooled maternal + offspring (family-structured) | Native | 330 | 0.063 | 0.004 | 0.008 | 0.004 |
| Pooled maternal + offspring (family-structured) | Cultivar | 471 | 0.167 | 0.001 | 0.021 | 0.001 |
| Pooled maternal + offspring (family-structured) | Total | 801 | 0.061 | 0.001 | 0.008 | 0.001 |

In the maternal-only dataset — the test unconfounded by family structure — neither the native nor the cultivar group showed evidence of multilocus linkage disequilibrium ($p > 0.2$ in both groups). The pooled dataset showed strongly significant linkage disequilibrium in both groups ($p \leq 0.004$), exactly as expected from the family (sibship) structure it contains rather than from any native/cultivar effect, and is reported here only to confirm that the test correctly detects the association it is known to be sensitive to. Judged from the maternal-only rows, we find no evidence of non-random multilocus allelic association distinguishing the native and cultivar groups.

## Discriminant Analysis of Principal Components (DAPC)

As a complementary multivariate test, we asked whether native and cultivar maternal trees separate on their full nine-locus genotype, considered jointly rather than one locus or one linkage relationship at a time, using discriminant analysis of principal components (DAPC; Jombart, Devillard, and Balloux 2010). This analysis was restricted to the 58 unrelated maternal trees for the same reason given above: no two of these individuals are siblings, so the result cannot reflect family structure.

DAPC is known to be prone to manufacturing apparent group separation when too many principal-component axes are retained relative to sample size, since additional axes give the discriminant step more ways to fit idiosyncratic differences between arbitrarily labeled groups. We guarded against this in two ways. First, rather than fixing the number of retained axes by inspection, we selected it by cross-validation (`xvalDapc`, capped at a maximum of 20 axes given $N = 58$), which chooses the number of axes that maximizes out-of-sample classification accuracy. Second, we assessed statistical significance by permutation rather than by treating any visual separation as evidence on its own: repeated (30 replicates) 80/20 stratified holdout classification accuracy was computed for the true native/cultivar labels, and this entire procedure was then repeated on 999 random relabelings of the same 58 trees (holding the number of retained axes fixed) to build an empirical null distribution.

**Table S2.** DAPC classification of native versus cultivar maternal trees on the full nine-locus genotype.

| $N$ (native, cultivar) | PCA axes retained (cross-validated) | Observed accuracy | No-information rate | Permutation mean | Permutation 95th percentile | *p*-value (999 permutations) |
|---|---|---|---|---|---|---|
| 58 (26, 32) | 18 | 0.561 | 0.552 | 0.509 | 0.615 | 0.237 |

Cross-validation retained 18 of the 20 allowed principal-component axes — a large number relative to $N = 58$, and the regime in which DAPC is most prone to overfitting. Even so, observed classification accuracy (56.1%) was close to the no-information rate set by the class imbalance (55.2% native/cultivar split) and indistinguishable from the permutation null (permutation mean 50.9%, 95th percentile 61.5%, $p = 0.237$). That the test had ample room to overfit and still did not exceed chance strengthens, rather than weakens, the conclusion that native and cultivar maternal trees are not distinguishable on overall multilocus genotype at these nine loci.

## Summary of Marker Analyses

Together with the admixture-percentage comparisons in the main text (Results, Genetic Analysis), these two analyses — one examining correlation structure among loci, the other examining multivariate separability of the full genotype — find no evidence that native and cultivar maternal trees form distinguishable genetic groups at these nine microsatellite loci.

## Software

All statistical analyses were performed in `R` (version 4.5.3; R Core Team 2026). Mixed-effects models used `lme4` (version 2.0.6; Bates et al. 2015), `lmerTest` (version 3.2.1; Kuznetsova et al. 2017), `MuMIn` (version 1.48.19) for AICc, and `DHARMa` (version 0.5.0; Hartig 2016) for GLMM diagnostics. Marker-based analyses used `gstudio` (version 1.14; Dyer 2009) and `genepop` (version 1.2.17; Rousset 2008) for exclusion probabilities and null-allele frequency estimation, `poppr` (version 2.9.8; Kamvar, Tabima, and Grünwald 2014) for the multilocus linkage-disequilibrium analysis, and `adegenet` (version 2.1.11; Jombart 2008) for the DAPC analysis. Annotated code and data are available in the online supplementary materials.
