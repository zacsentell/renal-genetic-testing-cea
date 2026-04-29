# Methods

Cost-effectiveness analysis (CEA) of genetic testing strategies for adults with chronic kidney disease (CKD) referred for suspected hereditary kidney disease. Perspective: Quebec public payer (RAMQ). Time horizon: 12-month diagnostic episode. Reporting follows CHEERS 2022.

---

## 1. Testing strategies

Three strategies are compared:

| # | Strategy | Description |
|---|---|---|
| 1 | **Panel** | Upfront targeted renal gene panel, matched to phenotype |
| 2 | **ES** | Upfront clinical exome sequencing |
| 3 | **Reflex (Panel→ES)** | Targeted panel first; unsolved cases reflex to clinical exome |

**Genome sequencing (GS)** is excluded from the base case because standard short-read pipelines lack validated detection of *MUC1* VNTR expansions (~3–5% of the renal genetics population). GS is evaluated in a scenario analysis (Section 8) modelling enhanced bioinformatics detection.

---

## 2. Target population

Adults with CKD referred to a specialized renal genetics clinic (MUHC, Montreal). Five phenotype categories are modelled: Cystic, Glomerular, Tubulointerstitial, Tubulopathies, and CKD of unknown etiology (CKDu). All strategies model singleton proband testing; cascade testing of first-degree relatives is included as a cost component for diagnosed probands. Costs are in 2024 CAD. Long-term outcomes (CKD progression, dialysis, QALYs) are not modelled.

---

## 3. Model structure

Individual-level probabilistic microsimulation with second-order Monte Carlo simulation. Two nested loops:

- **Outer loop (PSA):** 1,000 iterations; each draws uncertain parameters from prior distributions.
- **Inner loop:** 1,000 probands per iteration, fixed random-seed sequence for cross-strategy comparability.

All counts are configurable in `config.yml`.

### Per-patient pathway

Each simulated proband goes through the following steps:

1. **Phenotype assignment** — Dirichlet-multinomial draw from pooled phenotype proportions (literature synthesis).
2. **Monogenic status** — Beta-distributed probability conditional on phenotype.
3. **Causal architecture** — If monogenic, draw a (Gene, Inheritance mode, Variant class) triple from phenotype-specific empirical joint proportions derived from the curated variant table.
4. **Strategy application** — check gene coverage, apply analytic detection probability by modality and variant class. If detected: diagnose, accumulate costs, stop. If negative and strategy has reflex: draw a Bernoulli trial with probability equal to the reflex uptake parameter. Probands who accept proceed to the ES step; those who decline are counted as panel-only outcomes with no further testing costs.
5. **VUS outcome** — independent Bernoulli draw per testing step using modality-specific Beta-distributed probabilities.
6. **Cost accumulation** — sum clinical visits, laboratory testing, VUS follow-up, and cascade testing for diagnosed probands.

### Variant classes

Three mutually exclusive categories:

| Class | Definition |
|---|---|
| **SNV/Indel** | Small variants outside difficult loci |
| **CNV/SV** | Exon-level copy-number variants |
| **Difficult Locus** | Analytically challenging regions: *PKD1* (pseudogene/segmental duplication) and *MUC1* (VNTR) |

Phenotype-specific panels include bundled difficult-locus assays where clinically relevant:
- Cystic panel → targeted *PKD1* long-range PCR assay
- Tubulointerstitial panel → targeted *MUC1* VNTR assay
- ES and GS rely on intrinsic pipeline detection only (no bundled assays)

---

## 4. Parameterisation

### 4.1 Cohort phenotype composition and diagnostic yields

Parameterised from a structured MEDLINE search of adult renal genetics cohorts (2020–September 2025) meeting pre-specified criteria (≥150 probands, adult-predominant, broad CKD phenotype spectrum, first-line NGS, high-income country). Six cohorts were retained:

| Study | PMID | N | Modality |
|---|---|---|---|
| Vaisitti 2021 | 33226606 | 86 | Panel |
| Dahl 2023 | 37794564 | 1,619 | Panel |
| Jayasinghe 2021 | 32939031 | 165 | Exome |
| Blasco 2024 | 38972501 | 748 | Panel |
| Elhassan 2022 | 35099770 | 428 | Panel |
| Giovanella 2026 | 40533884 | 647 | Panel |

**Phenotype proportions:** Dirichlet parameters = pooled proband counts per phenotype across all six studies. Curation decisions documented in `data/audit/01_cohort_imported.csv`.

**Diagnostic yields:** Phenotype-specific yields pooled by random-effects meta-analysis (`scripts/02_meta_analysis.R`; R package `meta`). Beta distribution parameters (α, β) fit by method-of-moments to the pooled estimate and 95% CI:

| Phenotype | Studies | N | Pooled Yield [95% CI] | I² | α | β |
|---|---|---|---|---|---|---|
| Glomerular | 6 | 1,059 | 30.0% [21.4–40.2%] | 82% | 26.6 | 62.2 |
| Cystic | 6 | 931 | 62.1% [49.3–73.5%] | 91% | 36.2 | 22.1 |
| Tubulopathies | 5 | 194 | 39.5% [16.8–67.9%] | 76% | 4.2 | 6.4 |
| CKDu | 5 | 1,367 | 19.6% [6.9–44.4%] | 94% | 3.2 | 13.1 |
| Tubulointerstitial | 3 | 187 | 27.7% [15.6–44.3%] | 80% | 9.7 | 25.3 |

### 4.2 Causal architecture

Phenotype-specific empirical joint proportions for (Gene, Inheritance mode, Variant class) derived from the curated variant table (`data/raw/variant_curation.xlsx`). One row per P/LP variant: gene, HGVS notation, variant class, zygosity, inheritance mode, test modality. Audit trail in `data/audit/01_variant_imported.csv`, `02_variant_difficult.csv`, `02_variant_standard.csv`.

### 4.3 Analytic detection

Beta distribution parameters by modality and variant class, derived from laboratory validation datasets and published reference cohorts (`data/raw/variant_performance_curation.xlsx`). Posterior mean sensitivities:

| Variant Class | Panel (cystic) | Panel (glom.) | Panel (tub.int.) | Panel (tub.path.) | Panel (CKDu) | ES | GS |
|---|---|---|---|---|---|---|---|
| SNV/Indel | 99.6% | 99.6% | 99.6% | 99.6% | 99.6% | 99.6% | 99.9% |
| CNV/SV | 97.4% | 97.4% | 97.4% | 97.4% | 97.4% | 97.4% | 99.6% |
| PKD1 (difficult) | 98.9% | 12.5% | 12.5% | 12.5% | 12.5% | 12.5% | 92.9% |
| MUC1 (difficult) | 0% | 0% | 95.2% | 0% | 0% | 0% | 0% |

0% = no validated detection capability modelled for that modality–variant class combination.

### 4.4 VUS probabilities

Binary per-proband outcome (≥1 reportable VUS). Beta priors by modality and gene-count bin, fit from published VUS rates. Panel strategies are assigned a bin based on panel gene count; ES and GS use exome/genome-wide parameters. Parameters in `data/params/vus_betas.rds`.

| Modality | Gene Panel Size | Mean Pr(≥1 VUS) |
|---|---|---|
| Panel (tubulointerstitial) | 11–25 genes | 18.1% |
| Panel (cystic, glomerular, tubulopathies) | 51–100 genes | 47.1% |
| Panel (CKDu/comprehensive) | >200 genes | 76.2% |
| ES | Exome-wide | 27.6% |
| GS | Genome-wide | 27.6% |

### 4.5 Uptake probabilities

Two probabilities govern whether a proband completes an offered testing step. Both are loaded from `data/raw/uptake_parameters.csv`.

| Parameter | Base Case | PSA Prior | DSA Range | Source |
|---|---|---|---|---|
| Reflex ES uptake | 50% | Beta(5,5); mean 0.50 [0.19–0.81] | 20%–100% | Bindhu et al. |
| Cascade uptake | 100% | Fixed (structural assumption) | — | — |

Reflex ES uptake is the probability that a panel-negative proband proceeds to ES when offered. A Bernoulli trial is drawn independently per proband at each PSA iteration. Cascade uptake is held at 100% as a structural upper-bound assumption; this creates a conservative cost-effectiveness bias against higher-yield strategies because high-yield strategies trigger more cascade testing at 100% uptake.

### 4.6 Unit costs (2024 CAD)

All cost parameters from `data/raw/genetic_test_costs.csv`. Probabilistic uncertainty applied via Gamma distribution (CV = 0.25):

| Item | Mean (CAD) |
|---|---|
| Targeted multigene panel | $1,400 |
| Clinical exome | $2,500 |
| Clinical genome | $4,500 |
| Targeted single-variant test | $300 |
| Pre-test genetics consultation (RAMQ 09001) | $350 |
| Post-test genetics consultation (RAMQ 09001/09003) | $200 |

Cascade testing per eligible relative: 1 pre-test consultation + 1 targeted familial variant test. Base case: 2 eligible first-degree relatives per diagnosed proband (one-way DSA range: 0–4). VUS follow-up: 1 post-test visit + 1 targeted familial variant test per VUS-positive undiagnosed proband.

---

## 5. Outcomes

**Primary:** diagnostic yield (proportion with P/LP result), total cost per proband, incremental cost per additional diagnosis (ICPD).  
**Secondary:** VUS burden (proportion with ≥1 reportable VUS), cost composition by component.

Incremental analysis uses the efficiency frontier: strategies are ordered by cost; strictly and extendedly dominated strategies are removed; ICPD is computed between adjacent non-dominated strategies.

### Clinical impact classification

A secondary outcome quantifies the proportion of diagnosed probands whose molecular result drives a defined clinical action beyond genetic counselling and reproductive planning. The framework is adapted from Knoers et al. (NDT 2022), which identifies domains of clinical benefit from genetic testing in CKD. We use three non-exclusive primary categories that each correspond to a distinct downstream clinical action:

1. **Therapeutic.** Diagnosis enables a gene-specific treatment decision. Partitioned into three non-exclusive sub-categories: *targeted drug* (e.g., tolvaptan for *PKD1/PKD2*, eculizumab for complement genes, ERT for *GLA*, RAS inhibition initiated on the basis of an Alport diagnosis for *COL4A3/4/5*); *avoidance of ineffective therapy* (e.g., immunosuppression in hereditary SRNS or Alport, parathyroidectomy in FHH, AL-amyloid chemotherapy in hereditary amyloidosis); and *supportive or dietary management* (e.g., diet for hyperoxaluria, electrolyte replacement in tubulopathies).

*Caveat on RAS inhibition for Alport.* RAS inhibition is a drug class used broadly for proteinuric kidney disease, not a gene-specific pharmacotherapy in the narrow sense of tolvaptan or eculizumab. It is classified as *targeted drug* for *COL4A3/4/5* because initiation is triggered by the molecular Alport diagnosis (including in pre-proteinuric carriers), which is the clinical action the category is meant to capture. Readers who prefer a strict gene-specific-agent definition can subtract the Alport genes from the *targeted drug* count using `gene_coverage_audit.csv`.
2. **Extrarenal surveillance.** Diagnosis triggers a defined screening protocol for an extrarenal complication such as hearing loss, ocular abnormalities, diabetes, or tumour predisposition.
3. **Transplant and donor management.** Diagnosis alters living-donor evaluation, post-transplant recurrence risk, or transplant timing.

Genetic counselling and reproductive planning was excluded because it applies universally to any monogenic diagnosis and therefore cannot discriminate between strategies. Diagnostic reclassification was considered as a fourth category but excluded because it is clinically consequential only when it changes management, and that change is already captured by one of the three action-oriented categories above.

Gene-level assignments were curated from the "Diagnosis utility" column of the Knoers supplemental tables (S1 through S5) as the primary annotation source, with the paper body "Clinical benefit of genetic testing" subsections used as secondary evidence. Generic phrases ("adequation of treatment and follow-up", "renal protection strategies" alone, genetic counselling variants) did not qualify a gene for any category. Twenty-seven simulation genes absent from Knoers were annotated from independent clinical literature. Every assignment was manually reviewed against current clinical evidence. Full curation methodology and data sources are documented at `data/raw/clinical_utility/curation_methodology.md`; the curated table is at `data/raw/clinical_utility/gene_clinical_utility.csv`.

A gene was classified as high clinical impact if it mapped to at least one primary category. For each PSA iteration and testing strategy, the proportion of diagnosed probands whose causal gene falls into each category was computed (diagnosed-proband denominator). Cost per high-impact diagnosis was derived as total cost per proband divided by the proportion of all probands receiving a high-impact diagnosis (all-proband denominator), because testing costs are borne by all probands regardless of diagnostic outcome. No weights were applied across categories and no composite utility score was constructed.

**Phenotype attribution.** Several simulation genes contribute diagnoses to more than one phenotype category (for example *COL4A3*, *COL4A5*, *PKD1*, *UMOD*, *HNF1B*, *PAX2*). To preserve exact phenotype attribution, the simulation engine tabulates diagnosed probands at the (gene, phenotype) level at the end of each iteration, rather than summing over phenotype and allocating post hoc. This allows phenotype-stratified clinical impact proportions to reflect the same cohort realisation as the iteration-level cost and yield outcomes.

**Clinical impact figure.** One composite three-panel figure is produced at `clinical_impact_composite.png`. Panel A shows, for all four strategies (Panel, ES, Reflex, GS), the proportion of diagnosed probands in each primary impact category (High Impact, Therapeutic, Extrarenal Surveillance, Transplant/Donor), with 95% uncertainty intervals. Panel B shows, for the Reflex reference strategy, the proportion of all probands in each phenotype who fall into each impact category (all-proband denominator, so bar height absorbs phenotype-specific yield). Panel C is a condensed gene-by-category dot matrix of the top ten genes by High Impact contribution under Reflex; dot size encodes mean cohort proportion diagnosed. Impact categories are non-exclusive.

---

## 6. Uncertainty and sensitivity analyses

| Analysis | Method | Script |
|---|---|---|
| Probabilistic sensitivity (PSA) | Outer Monte Carlo loop; robustness = fraction of iterations where each strategy is non-dominated | `06a`, `08a` |
| Cost-effectiveness acceptability (CEAC) | Net Monetary Benefit at WTP $0–$100,000 CAD (step $1,000) | `08_wtp_analysis.R` |
| Global sensitivity (PRCC) | Partial rank correlation between 31 sampled inputs and incremental cost/yield outcomes; two comparisons: Reflex vs Panel, GS scenario vs Reflex | `07_global_sensitivity_analysis.R` |

---

## 7. Subgroup and scenario analyses

- **Cystic subgroup:** base-case CEA restricted to cystic phenotype (`06c_cystic_subgroup.R`). Clinically relevant because the cystic panel bundles a targeted PKD1 assay, giving Panel a structural detection advantage over standalone ES.
- **Phenotype stratification:** Reflex (Panel→ES) outcomes by phenotype category (`06e_phenotype_figure.R`).
- **Phenotype-specific WTP winners:** most probable cost-effective strategy by phenotype at selected WTP anchors, reported as a matrix figure and findings tables (`06e_phenotype_figure.R`).
- **GS scenario:** +10% relative yield uplift applied to GS arm to model enhanced bioinformatics detection of difficult loci. GS scenario compared against Reflex as the incumbent frontier strategy (`09a_gs_scenario_analysis.R`). The uplift is a deterministic scenario assumption; uncertainty in its magnitude is not propagated through the PSA.
- **PKD1 full detection scenario:** PKD1 difficult-locus detection overridden to 100% for all panels and ES, representing a near-future state where optimised capture and pseudogene-aware bioinformatics eliminate the current PKD1 detection barrier (`09d_pkd1_full_detection_scenario.R`). A full PSA re-run (1,000 iterations × 1,000 probands) is performed. MUC1 VNTR detection is unchanged. Test costs are unchanged.
- **Cascade eligible relatives DSA:** one-way sensitivity varying the number of eligible first-degree relatives from 0 to 4 (`09c_cascade_dsa.R`). Examines whether frontier ranking is preserved across cascade assumptions.
- **Reflex uptake DSA:** one-way sensitivity varying reflex ES uptake from 20% to 100% (`09e_reflex_uptake_dsa.R`). Examines incremental cost per additional diagnosis as a function of the proportion of panel-negative probands who proceed to ES.

### 7.1 Phenotype-specific WTP winner analysis

This analysis identifies the most probable cost-effective strategy for each phenotype at selected willingness-to-pay (WTP) anchors. Input data are iteration-level phenotype outcomes from the base case.

For each phenotype, strategy, and WTP anchor, net monetary benefit (NMB) is computed as:

NMB = (WTP × diagnoses per proband) − total cost per proband.

Within each iteration, the strategy with the highest NMB is selected as the winner. If strategies tie, deterministic ordering by strategy ID is used to ensure reproducible winner assignment. Winner probability is the proportion of iterations in which each strategy is optimal for that phenotype and WTP anchor.

WTP anchors are derived deterministically from the configured WTP range and step in `config.yml`, including the minimum and maximum values. Two output tables are produced in `outputs/results/supplement/phenotype_stratified/`: a long-format table (`phenotype_wtp_winner_long.csv`) with one row per phenotype and anchor, and a wide findings table (`phenotype_wtp_winner_wide.csv`) for direct report display. The corresponding matrix figure is exported as `phenotype_wtp_winner_matrix.png` and `.svg`.

---

## 8. Limitations

1. **Time horizon.** 12-month episode only; no long-term outcomes, QALYs, or downstream treatment effects.
2. **Gene coverage.** ES and GS are modelled with 100% causal-gene coverage. Coverage gaps in practice may overestimate ES/GS yield.
3. **GS MUC1 detection.** Modelled as 0% for standard short-read GS. Long-read and bespoke bioinformatics pipelines are not parameterised.
4. **VUS modelling.** Binary per-proband outcome; VUS count, gene-specific rates, and correlation with diagnosis are not modelled.
5. **Cost sources.** Send-out pricing and RAMQ list prices; in-house and negotiated institutional pricing may differ.
6. **Referral population.** Synthetic cohort derived from published case series; phenotype proportions may not generalise to all renal genetics referral populations.
7. **Uptake assumptions.** Cascade testing is modelled at 100% uptake (structural upper-bound assumption), which creates a conservative bias against higher-yield strategies. Reflex ES uptake is set to 50% in the base case, based on published cascade and reflex acceptance estimates; the impact of varying this from 20% to 100% is characterised in the reflex uptake DSA.

