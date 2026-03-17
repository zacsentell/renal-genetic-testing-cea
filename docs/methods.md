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
4. **Strategy application** — check gene coverage, apply analytic detection probability by modality and variant class. If detected: diagnose, accumulate costs, stop. If negative and strategy has reflex: execute ES step.
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
| Vaisitti 2021 | 33226590 | 131 | Panel |
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

### 4.5 Unit costs (2024 CAD)

All cost parameters from `data/raw/genetic_test_costs.csv`. Probabilistic uncertainty applied via Gamma distribution (CV = 0.25):

| Item | Mean (CAD) |
|---|---|
| Targeted multigene panel | $1,400 |
| Clinical exome | $2,500 |
| Clinical genome | $4,500 |
| Targeted single-variant test | $300 |
| Pre-test genetics consultation (RAMQ 09001) | $350 |
| Post-test genetics consultation (RAMQ 09001/09003) | $200 |

Cascade testing per eligible relative: 1 pre-test consultation + 1 targeted familial variant test. Base case: 2 eligible first-degree relatives per diagnosed proband (sensitivity range: 0–4). VUS follow-up: 1 post-test visit + 1 targeted familial variant test.

---

## 5. Outcomes

**Primary:** diagnostic yield (proportion with P/LP result), total cost per proband, incremental cost per additional diagnosis (ICPD).  
**Secondary:** VUS burden (proportion with ≥1 reportable VUS), cost composition by component.

Incremental analysis uses the efficiency frontier: strategies are ordered by cost; strictly and extendedly dominated strategies are removed; ICPD is computed between adjacent non-dominated strategies.

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
- **GS scenario:** +10% relative yield uplift applied to GS arm to model enhanced bioinformatics detection of difficult loci. GS scenario compared against Reflex as the incumbent frontier strategy (`09a_gs_scenario_analysis.R`). The uplift is a deterministic scenario assumption; uncertainty in its magnitude is not propagated through the PSA.

---

## 8. Limitations

1. **Time horizon.** 12-month episode only; no long-term outcomes, QALYs, or downstream treatment effects.
2. **Gene coverage.** ES and GS are modelled with 100% causal-gene coverage. Coverage gaps in practice may overestimate ES/GS yield.
3. **GS MUC1 detection.** Modelled as 0% for standard short-read GS. Long-read and bespoke bioinformatics pipelines are not parameterised.
4. **VUS modelling.** Binary per-proband outcome; VUS count, gene-specific rates, and correlation with diagnosis are not modelled.
5. **Cost sources.** Send-out pricing and RAMQ list prices; in-house and negotiated institutional pricing may differ.
6. **Referral population.** Synthetic cohort derived from published case series; phenotype proportions may not generalise to all renal genetics referral populations.

