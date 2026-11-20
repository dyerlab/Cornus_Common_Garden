# Text changes: restricting genetic analyses to the garden mothers

All genetic analyses have been rerun on the study material only: the maternal
trees whose families were planted in the common garden, plus their genotyped
garden seedlings. Previously they used the full genotyped collection (58
mothers, 828 offspring), and 30 of those mothers had no family in the garden.

**What's in the restricted set**
- Garden families: 36 (16 native, 20 cultivar). Trait and survival models are unchanged.
- Garden families with a genotyped mother: **28 (9 native, 19 cultivar)**.
- Families with no genetic data at all (no genotyped mother or seedlings): nat18, nat24, nat28, nat30, nat31, nat32, nat33, urb10.
- Genotyped garden seedlings: **397** (103 native, 294 cultivar). These are the same 397 seedlings already used in the AP-subset models.
- 15 genotyped offspring of garden mothers do not appear among the 569 planted seedlings, so they are excluded.

**Unchanged:** every trait, survival, heritability, site and diagnostic result
(Tab_Samples, Tab_Survival, Tab_Growth\*, Tab_OriginHeritability, Tab_Repeated\*,
Tab_CovariateEffects, Tab_SiteEffect, Fig_Survival, Fig_BLUPs). Their pipeline
output is byte-identical to the previous run.

Locations below are given as **section → paragraph**, followed by a short
string to search for (Ctrl-F) in the Google Doc.

---

## 0. Methods: say what was genotyped (new text)

**Manuscript → Methods → Genetic Analyses, end of paragraph 1**
(find: `individuals in commercial nurseries in the greater Richmond`)

Add:

> Maternal trees and seedlings were genotyped for 28 of the 36 garden families
> (9 native, 19 cultivar; 397 seedlings). The maternal trees of the remaining
> eight families (seven native, one cultivar) were not genotyped, and none of
> their seedlings were. All marker-based analyses are therefore restricted to
> these 28 families.

*Optional:* in the Tab_Survival and Tab_GrowthAIC captions, after "AP
(seedlings with a genotype)", add "(n = 397 seedlings from 28 families)".

---

## 1. Marker diversity (Tab_MarkerDiversity)

**Manuscript → Results → Genetic Analysis, paragraph 1** (find: `5–26 alleles per locus`)
- `5–26 alleles per locus` → **`5–21 alleles per locus`**

**Supplement → Marker Panel Diversity and Validation → Tab_MarkerDiversity**: replace the whole table.

| Sample | Locus | *A* | *A*ₑ | *H*ₒ | *H*ₑ | *F*IS | *P*ₑ |
| :-- | :-- | --: | --: | --: | --: | --: | --: |
| Maternal | cf020 | 10 | 7.839 | 0.815 | 0.872 | 0.066 | 0.877 |
| Maternal | cf125 | 5 | 3.424 | 0.714 | 0.708 | −0.009 | 0.681 |
| Maternal | cf213 | 21 | 15.525 | 0.964 | 0.936 | −0.031 | 0.931 |
| Maternal | cf273 | 9 | 2.330 | 0.643 | 0.571 | −0.126 | 0.582 |
| Maternal | cf581 | 5 | 1.298 | 0.143 | 0.230 | 0.378 | 0.215 |
| Maternal | cf585 | 15 | 9.278 | 0.929 | 0.892 | −0.041 | 0.881 |
| Maternal | cf597 | 11 | 7.127 | 0.857 | 0.860 | 0.003 | 0.875 |
| Maternal | cf634 | 11 | 6.222 | 0.750 | 0.839 | 0.106 | 0.841 |
| Maternal | cf701 | 7 | 4.052 | 0.643 | 0.753 | 0.146 | 0.753 |
| Maternal | **all loci** | 94 | 57.094 | 0.718 | 0.740 | 0.055 | — |
| Offspring | cf020 | 14 | 8.143 | 0.795 | 0.877 | 0.094 | 0.877 |
| Offspring | cf125 | 6 | 3.105 | 0.696 | 0.678 | −0.026 | 0.681 |
| Offspring | cf213 | 27 | 14.342 | 0.902 | 0.930 | 0.031 | 0.931 |
| Offspring | cf273 | 10 | 2.397 | 0.566 | 0.583 | 0.028 | 0.582 |
| Offspring | cf581 | 5 | 1.271 | 0.196 | 0.213 | 0.082 | 0.215 |
| Offspring | cf585 | 20 | 8.305 | 0.772 | 0.880 | 0.122 | 0.881 |
| Offspring | cf597 | 15 | 8.032 | 0.787 | 0.875 | 0.101 | 0.875 |
| Offspring | cf634 | 14 | 6.277 | 0.813 | 0.841 | 0.032 | 0.841 |
| Offspring | cf701 | 12 | 4.045 | 0.697 | 0.753 | 0.074 | 0.753 |
| Offspring | **all loci** | 123 | 55.917 | 0.692 | 0.737 | 0.060 | — |

---

## 2. Null alleles and exclusion probability (Tab_NullAlleles) — the result changes

The reduced panel changes. Three loci are now flagged in both origins
(cf020, cf597, cf634), where previously four were (cf020, cf585, cf597, cf701).
That leaves a **six-locus** reduced panel (cf125, cf213, cf273, cf581, cf585,
cf701), not five. Native point estimates are now as high as 0.12.

**Manuscript → Results → Genetic Analysis, paragraph 1** (find: `all point estimates < 0.07`)

Replace the sentence beginning "Maximum-likelihood estimates of null-allele frequency…" with:

> Maximum-likelihood estimates of null-allele frequency were generally low
> (point estimates ≤ 0.12, highest in the native sample at cf020 and cf634),
> and their confidence intervals excluded zero in both the native and
> cultivar samples at three loci (cf020, cf597, cf634; Tab_NullAlleles,
> Supplementary Materials).

Same paragraph (find: `P_{excl,multi} = 0.9990`, or `0.9990` in the rendered doc):
- `0.9990` → **`0.9998`**

**Supplement → Tab_NullAlleles**: replace the table.

| Locus | Cultivar | Native | In reduced panel |
| :-- | :-- | :-- | :--: |
| cf020 | 0.039 [0.017, 0.065] | 0.119 [0.073, 0.181] | no |
| cf125 | 0.015 [0.000, 0.050] | 0.000 (no CI) | yes |
| cf213 | 0.021 [0.005, 0.042] | 0.000 (no CI) | yes |
| cf273 | 0.029 [0.008, 0.060] | 0.018 [0.000, 0.071] | yes |
| cf581 | 0.013 [0.000, 0.060] | 0.077 [0.022, 0.138] | yes |
| cf585 | 0.068 [0.044, 0.096] | 0.000 (no CI) | yes |
| cf597 | 0.036 [0.016, 0.061] | 0.050 [0.011, 0.102] | no |
| cf634 | 0.022 [0.003, 0.047] | 0.120 [0.071, 0.187] | no |
| cf701 | 0.028 [0.004, 0.056] | 0.009 (no CI) | yes |

**Supplement, line below that table** (find: `reduced five-locus panel`)
- → "full nine-locus panel *P*excl = 1.00000; **reduced six-locus panel (cf125, cf213, cf273, cf581, cf585, cf701) *P*excl = 0.99979**."

*Note:* Genepop returned no CI for four native loci. The native sample is now
small (9 mothers + 103 seedlings), so the native estimates are noisier.

---

## 3. Admixture Percentage (Tab_AdmixturePanelComparison)

**Manuscript → Results → Genetic Analysis, paragraph 3** (find: `Among the 58 genotyped maternal trees`)
- `Among the 58 genotyped maternal trees, AP averaged 49.210% (range 33.333–88.889%), and across the 828 genotyped offspring it averaged 49.263% (range 27.778–75.000%)`
  → **`Among the 28 genotyped maternal trees, AP averaged 49.752% (range 33.333–88.889%), and across the 397 genotyped offspring it averaged 50.363% (range 27.778–75.000%)`**
- `AP averaged 51.667% (range 33.333–91.667%) among maternal trees and 52.104% (range 25.000–100.000%) among offspring`
  → **`AP averaged 53.274% (range 33.333–91.667%) among maternal trees and 54.039% (range 25.000–100.000%) among offspring`**
- The surrounding interpretation ("slightly higher on average… broad overlap") still holds.

**Manuscript → Tables → Tab_AdmixturePanelComparison**: replace the table body.

| Group | Marker set | *n* | Mean AP (%) | SD | Range (%) |
| :-- | :-- | --: | --: | --: | :-- |
| Maternal | full (9 loci) | 28 | 49.752 | 11.881 | 33.333–88.889 |
| Offspring | full (9 loci) | 397 | 50.363 | 8.858 | 27.778–75.000 |
| Maternal | Wadl (6 loci) | 28 | 53.274 | 10.957 | 33.333–91.667 |
| Offspring | Wadl (6 loci) | 397 | 54.039 | 10.569 | 25.000–100.000 |

---

## 4. Multilocus linkage disequilibrium (Tab_LinkageDisequilibrium)

**Manuscript → Discussion, paragraph beginning "A surprising result in this work"**
(find: `restricted to the 58 unrelated maternal trees`)
- `58` → **`28`**
- `native ($I_A$ = -0.071, *p* = 0.717) or cultivar ($I_A$ = 0.063, *p* = 0.235)`
  → **`native ($I_A$ = -0.116, *p* = 0.607) or cultivar ($I_A$ = 0.079, *p* = 0.262)`**

**Supplement → Multilocus Linkage Disequilibrium**
- List item 1 (find: `58 unrelated individuals`) → **`28 unrelated individuals`**
- List item 2 (find: `886 individuals`) → **`425 individuals`**
- Replace the table:

| Dataset | Group | N | $I_A$ | *p* ($I_A$) | $\bar{r}_D$ | *p* ($\bar{r}_D$) |
|---|---|---|---|---|---|---|
| Maternal trees only (unrelated individuals) | Native | 9 | -0.116 | 0.607 | -0.016 | 0.613 |
| Maternal trees only (unrelated individuals) | Cultivar | 19 | 0.079 | 0.262 | 0.010 | 0.261 |
| Maternal trees only (unrelated individuals) | Total | 28 | -0.033 | 0.621 | -0.004 | 0.622 |
| Pooled maternal + offspring (family-structured) | Native | 95 | 0.215 | 0.001 | 0.027 | 0.001 |
| Pooled maternal + offspring (family-structured) | Cultivar | 254 | 0.197 | 0.001 | 0.025 | 0.001 |
| Pooled maternal + offspring (family-structured) | Total | 349 | 0.121 | 0.001 | 0.015 | 0.001 |

- Paragraph after the table (find: `$p \leq 0.004$`) → **`$p = 0.001$`**. The `$p > 0.2$ in both groups` wording still holds.

(Pooled N is lower than 425 because poppr drops individuals missing more than 1 of 9 loci. The same was true before: 801 of 886.)

---

## 5. DAPC (Tab_DAPC) — wording must change

Observed accuracy now falls **below** the no-information rate. The conclusion
(no separation) is the same or stronger, but "close to the no-information
rate" no longer describes it.

**Manuscript → Discussion, same "A surprising result" paragraph** (find: `18 of a maximum 20`)
- `(18 of a maximum 20, chosen by cross-validation)` → **`(14 of a maximum 20, chosen by cross-validation)`**
- `cross-validated classification accuracy (56.061%) was indistinguishable from chance given the sample's class balance (no-information rate 55.172%)`
  → **`cross-validated classification accuracy (64.444%) did not exceed chance given the sample's class balance (no-information rate 67.857%)`**
- `999 relabelings of the same 58 trees (*p* = 0.237)` → **`999 relabelings of the same 28 trees (*p* = 0.273)`**

**Supplement → DAPC**
- Paragraph 1 (find: `restricted to the 58 unrelated maternal trees`) → **`28`**
- Paragraph 2 (find: `given $N = 58$`) → **`given $N = 28$`**
- Paragraph 2 (find: `relabelings of the same 58 trees`) → **`28`**
- Table row → **`| 28 (9, 19) | 14 | 0.644 | 0.679 | 0.600 | 0.722 | 0.273 |`**
- Replace the paragraph after the table (find: `Cross-validation retained 18`) with:

> Cross-validation retained 14 of the 20 allowed principal-component axes — a
> large number relative to $N = 28$, and the regime in which DAPC is most prone
> to overfitting. Even so, observed classification accuracy (64.4%) fell below
> the no-information rate set by the class imbalance (67.9%; 19 of 28 trees are
> cultivar-origin) and was indistinguishable from the permutation null
> (permutation mean 60.0%, 95th percentile 72.2%, $p = 0.273$). That the test
> had ample room to overfit and still did not exceed chance strengthens, rather
> than weakens, the conclusion that native and cultivar maternal trees are not
> distinguishable on overall multilocus genotype at these nine loci.

*Note:* with N = 28, a cap of 20 axes is generous. Reviewers may ask about it.
I can rerun with a lower cap (for example 10) if you'd like.

---

## 6. Pairwise relatedness among maternal trees (Tab_Relatedness)

**Manuscript → Results, paragraph beginning "Native-origin mothers were sampled from only two"**
The conclusion is unchanged, so no number edits are needed. Consider adding a
power caveat, for example after "We found no evidence of such a difference":
"(although with only nine genotyped native mothers this test has limited power)".

**Supplement → Pairwise Relatedness Among Maternal Trees, paragraph 2**
- (find: `the same 58 genotyped maternal trees`) → **`the same 28 genotyped maternal trees`**
- (find: `across the same 58 individuals`) → **`28`**
- (find: `(26 native, 32 cultivar)`) → **`(9 native, 19 cultivar)`**

**Tab_Relatedness row** → **`| 28 (9, 19) | -0.025 | 0.002 | -0.027 | 0.672 | 0.133 | 0.887 |`**

**Paragraph after the table**: replace the first sentence (find: `smaller than expected under random relabeling`) with:

> Neither statistic showed evidence of a difference in relatedness structure
> between the two groups: neither the observed difference in mean pairwise
> relatedness (*p* = 0.672) nor the difference between the two within-group
> relatedness distributions (*p* = 0.887) was larger than expected under random
> relabeling of the same 28 individuals.

(The old claim that the two groups were "if anything, more similar" relied on
*p* = 0.995. That no longer applies, so drop it.) The rest of the paragraph
stands.

---

## 7. Realized offspring relatedness and heritability (Tab_HeritabilityRelatedness, Fig_OffspringRelatedness, Fig_Heritability)

**Supplement → Tab_HeritabilityRelatedness caption** and **Fig_OffspringRelatedness caption**: both contain the same parenthetical.
- `native median $r$ = 0.259, $n$ = 1421 pairs among 24 mothers; cultivar median $r$ = 0.291, $n$ = 1810 pairs among 31 mothers`
  → **`native median $r$ = 0.268, $n$ = 256 pairs among 8 mothers; cultivar median $r$ = 0.274, $n$ = 878 pairs among 18 mothers`**
- In Fig_OffspringRelatedness the form is "Native median $r$ = 0.259 ($n$ = 1421 pairs among 24 mothers); cultivar median $r$ = 0.291 ($n$ = 1810 pairs among 31 mothers)". Change it the same way.
- The 8 and 18 mothers are the garden mothers with at least 2 complete-genotype seedlings.

**Tab_HeritabilityRelatedness**: only the "Relatedness" column changes. Half- and full-sib columns are unchanged.

| Trait                | Origin   | Half Sibling | Relatedness | Full Sibling |
|----------------------|----------|:------------:|:-----------:|:------------:|
| Plant height         | Native   | 0.238        | 0.222       | 0.119        |
| Plant height         | Cultivar | 0.474        | 0.432       | 0.237        |
| Stem diameter        | Native   | 0.008        | 0.007       | 0.004        |
| Stem diameter        | Cultivar | 0.247        | 0.226       | 0.124        |
| Number of leaves     | Native   | 0.000        | 0.000       | 0.000        |
| Number of leaves     | Cultivar | 0.263        | 0.239       | 0.131        |
| Above-ground biomass | Native   | 0.000        | 0.000       | 0.000        |
| Above-ground biomass | Cultivar | 0.322        | 0.294       | 0.161        |
| Below-ground biomass | Native   | 0.103        | 0.096       | 0.051        |
| Below-ground biomass | Cultivar | 0.233        | 0.213       | 0.117        |
| Leaf biomass         | Native   | 0.000        | 0.000       | 0.000        |
| Leaf biomass         | Cultivar | 0.130        | 0.119       | 0.065        |
| Stem biomass         | Native   | 0.090        | 0.084       | 0.045        |
| Stem biomass         | Cultivar | 0.438        | 0.399       | 0.219        |

The surrounding supplement text and the main-text Discussion paragraph (find:
`As a first, more limited step toward that verification`) remain accurate as
written: estimates still shift modestly downward and the pattern is intact.

**Figures:** replace `media/fig_offspring_relatedness.png` (Fig_OffspringRelatedness)
and `media/fig_offspring_heritability.png` (Fig_Heritability) in the Google Doc
with the regenerated versions. The Fig_Heritability caption is unchanged.

---

## 8. Nothing to change

- Abstract (36 families, 20/16, 569 seedlings): all garden counts, still correct.
- Tab_CultivarVoucherPanel, "21 named lines": unchanged.
- Discussion, "$n = 174$ across 16 families vs. $n = 395$ across 20": garden counts, still correct.
- Supplement, "Summary of Marker Analyses": no numbers, and the conclusions still hold.
