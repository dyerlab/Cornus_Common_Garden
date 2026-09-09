# Appendix S1: Statistical Model Structure, Variance Decomposition, and Model Diagnostics

This appendix provides the full statistical detail supporting the *Cultivar Impact Analysis* subsection of the Methods, including model equations, variance-component and heritability derivations, diagnostic procedures, and sensitivity analyses. Annotated analysis code and the underlying data are available in the online supplementary materials.

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

For linear mixed models, assumptions were evaluated using residual-vs-fitted plots (homoscedasticity), Q–Q plots of residuals (normality), and Cook's distance (influential groups). Response variables were log-transformed where diagnostic plots indicated heteroscedasticity or non-normality. For generalized linear mixed models (survival), fit was assessed using simulation-based residual diagnostics (1,000 simulated residuals per model via `simulateResiduals`, *DHARMa*; Hartig 2016), with visual inspection and formal tests for residual uniformity, dispersion, and outliers; no meaningful violations of GLMM assumptions were detected. Family-specific conditional modes (BLUPs) were extracted from the cultivar-origin survival model using `ranef` (*lme4*; Bates et al. 2015) to visualize among-family variation in survival (Fig_BLUPs).
## Sensitivity Analyses
To assess the robustness of fixed- and random-effects inference, we refit each model after temporarily removing (a) individual-level outliers, identified by standardized residuals with $|z| > 3$, and (b) influential maternal families, identified by large Cook's distance values.

Removing outlier individual-level observations did not substantially change heritability estimates for most traits, and the significance of family-level random effects generally increased slightly; the largest shift was in below-ground biomass ($h^2 = 0.23 \rightarrow 0.36$).

Removing the single most influential maternal family (identified by the largest leave-one-out Cook's distance on the fixed effects) had a larger effect on variance-component and heritability estimates than removing outlier individuals. For most growth traits this family was urb09. For example, in the stem biomass model, removing urb09 reduced the family variance component ($V_{family} = 0.046 \rightarrow 0.016$) and heritability ($h^2 = 0.51 \rightarrow 0.19$), and weakened but did not eliminate support for the family random effect ($\chi^2 = 5.3$, df $= 1$, $p = 0.022$).

## Software

All statistical analyses were performed in `R` (version 4.5.3; R Core Team 2026) using `lme4` (version 2.0.6; Bates et al. 2015), `lmerTest` (version 3.2.1; Kuznetsova et al. 2017), `MuMIn` (version 1.48.19) for AICc, `DHARMa` (version 0.5.0; Hartig 2016), and `gstudio` (version 1.14; Dyer 2009) and `genepop` (version 1.2.17; Rousset 2008) for the marker analyses. Annotated code and data are available in the online supplementary materials.
