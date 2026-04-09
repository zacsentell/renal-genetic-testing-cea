# Supplementary Methods and Results
## Cost-Effectiveness of Genetic Evaluation in Chronic Kidney Disease
### Renal Genetics Clinic, MUHC

---

**Purpose.** This document provides a self-contained account of the cost-effectiveness analysis (CEA) supporting the manuscript *Genetic Evaluation in Chronic Kidney Disease: A Canadian Single-Centre Experience.* It follows CHEERS 2022 reporting guidelines (see Appendix A). All numerical results are generated directly from simulation outputs; the CHEERS checklist and key-finding callouts are reviewed and confirmed after each pipeline run.

**Run configuration.** n_iter = 1000; n_probands = 1000; seed_start = 2025; cost CV = 0.25; GS yield uplift (scenario) = +10%.
**Generated:** 2026-03-29 (post reflex-uptake parameterization; base-case reflex uptake = 50%). See `docs/report/report.md` for the canonical version including §7.3 (PKD1 full detection scenario) and Appendix E (reflex uptake DSA).

---

## §1  Decision Problem and Study Context

### 1.1  Clinical question

We evaluate the cost-effectiveness of alternative first-line genetic testing strategies for adults with chronic kidney disease (CKD) referred to a specialized Renal Genetics Clinic for suspected hereditary kidney disease. The primary decision is which sequencing modality — or which sequencing sequence — maximises diagnostic yield at acceptable incremental cost. The analysis is conducted from the perspective of the Quebec public payer (RAMQ) over a 12-month diagnostic episode.

**Primary outcomes**

1. **Diagnostic yield:** proportion of probands receiving a pathogenic or likely pathogenic (P/LP) variant result.
2. **Total cost per proband:** mean diagnostic-episode cost.
3. **Incremental cost per additional diagnosis (ICPD):** the marginal cost of achieving one additional diagnosis by moving to the next strategy on the efficiency frontier.

**Secondary outcomes:** VUS burden (proportion with ≥1 reportable variant of uncertain significance), cost composition by category.

### 1.2  Testing strategies

Three strategies reflecting current MUHC practice in adult renal genetics are compared:

| # | Strategy | Description |
| --- | --- | --- |
| 1 | **Panel** | Upfront targeted renal gene panel, matched to phenotype category |
| 2 | **ES** | Upfront clinical exome sequencing |
| 3 | **Reflex (Panel→ES)** | Targeted panel first; unsolved cases reflex to clinical exome sequencing |

> **Note on genome sequencing (GS).** Standard short-read GS pipelines lack validated detection of *MUC1* VNTR expansions, a key cause of ADTKD present in approximately 3–5% of the renal genetics population. Including GS in the base case would misrepresent its diagnostic performance relative to phenotype-specific panels that bundle targeted *MUC1* assays. GS is therefore evaluated in a separate scenario analysis (§7.2) modelling near-future enhanced bioinformatics detection.

### 1.3  Target population and setting

Adults with CKD referred to the MUHC Renal Genetics Clinic for evaluation of suspected hereditary kidney disease. Five phenotype categories are modelled: Cystic, Glomerular, Tubulointerstitial (ADTKD), Tubulopathies, and CKD of unknown etiology (CKDu). All strategies are singleton proband testing; cascade testing of first-degree relatives is included as a cost component for probands with a P/LP result. Costs are valued in 2024 Canadian dollars (CAD). Long-term outcomes (CKD progression, dialysis, transplantation, QALYs) are not modelled.

---

## §2  Model Structure

### 2.1  Overview

The CEA uses an individual-level probabilistic cohort model with two nested Monte Carlo loops. The **outer loop** (PSA) draws uncertain parameters from their prior distributions for each of 1000 iterations. The **inner loop** simulates 1000 unique probands per iteration under each strategy, using a fixed random seed sequence for cross-strategy comparability. Aggregating inner-loop probands produces one cost and yield estimate per strategy per outer iteration; the outer-loop distribution defines uncertainty intervals.

### 2.2  Per-patient diagnostic pathway

For each simulated proband and strategy:

1. **Assign phenotype** — draw from a Dirichlet-multinomial distribution parameterised by pooled phenotype proportions from the literature synthesis.
2. **Assign monogenic status** — conditional on phenotype, draw from a Beta distribution representing the phenotype-specific detectable monogenic probability under a high-sensitivity reference strategy.
3. **Assign causal architecture** — if monogenic, draw a (Gene, Inheritance mode, Variant class) triple from phenotype-specific empirical joint proportions derived from the curated variant table.
4. **Apply strategy** — traverse the planned test sequence: apply gene-coverage check, then analytic detection probability for the modality and variant class. If detected, record diagnosis and accumulate costs; if negative and the strategy includes reflex, draw a Bernoulli trial with probability equal to the reflex uptake parameter to determine whether the proband proceeds to ES.
5. **VUS outcome** — independently draw whether ≥1 reportable VUS is returned at each testing step, using modality-specific Beta-distributed probabilities.
6. **Accumulate costs** — sum clinical visits, laboratory testing, VUS follow-up (if applicable), and cascade testing for diagnosed probands. Cascade cost is computed as Nₜᵣₗ × (pre-test consultation + targeted familial variant test), applied only to probands with a P/LP result. Three simplifying assumptions are made: (a) **100% cascade uptake**; (b) **no health benefit to relatives is modelled**; (c) **50% reflex uptake (base case)** — of panel-negative probands offered ES, half accept. Parameterized as a Bernoulli trial per proband; Beta(5,5) prior in the PSA (mean 0.50, 95% interval 0.19–0.81). The base case assumes 2 eligible first-degree relatives; one-way sensitivity analyses vary Nₜᵣₗ from 0 to 4 (Appendix D) and reflex uptake from 20% to 100% (Appendix E in the canonical report).

### 2.3  Variant classes and detection logic

Three mutually exclusive variant-class categories are modelled: **SNV/Indel** (small variants outside difficult loci), **CNV/SV** (exon-level copy-number variants), and **Difficult Locus** (variants in analytically challenging regions, currently *PKD1* pseudogene region and *MUC1* VNTR). Phenotype-specific panels include bundled difficult-locus assays where clinically relevant: the cystic panel includes a targeted *PKD1* assay; the tubulointerstitial panel includes a targeted *MUC1* assay. ES and GS rely on intrinsic pipeline detection for all variant classes.

---

## §3  Parameterization

### 3.1  Cohort phenotype composition

Phenotype proportions are derived from a systematic literature synthesis of contemporary adult or adult-predominant renal genetics cohorts (2020 onward). Dirichlet parameters are estimated from pooled proportions across eligible cohorts.

### 3.2  Diagnostic yield estimates

Phenotype-specific diagnostic yields were pooled using a random-effects meta-analysis of phenotype-stratified yields from eligible cohorts. Beta distribution parameters (Alpha, Beta) are derived from the pooled estimate and its 95% credible interval using method-of-moments fitting.

| Phenotype | Studies | Probands | Pooled Yield (95% CI) | I² (%) | α | β |
| --- | --- | --- | --- | --- | --- | --- |
| Glomerular | 6 | 1059 | 30.0% [21.4-40.2] | 82% | 26.6 | 62.2 |
| Cystic | 6 | 931 | 62.1% [49.3-73.5] | 91% | 36.2 | 22.1 |
| Tubulopathies | 5 | 194 | 39.5% [16.8-67.9] | 76% | 4.2 | 6.4 |
| CKDu | 5 | 1367 | 19.6% [6.9-44.4] | 94% | 3.2 | 13.1 |
| Tubulointerstitial | 3 | 187 | 27.7% [15.6-44.3] | 80% | 9.7 | 25.3 |

_k = number of studies; N = total probands; I² = between-study heterogeneity; α, β = Beta distribution parameters (method-of-moments fit to pooled yield and 95% CI)._

![Genetic landscape: top 20 causal genes by phenotype and variant class](../outputs/figures/gene_architecture_by_phenotype.png)

_Figure A. Genetic landscape of the top 20 causal genes ranked by weighted global prevalence across phenotypes. The central panel shows total weighted diagnoses per gene; the left panel shows the distribution of variant classes (SNV/Indel, CNV/SV, Difficult Locus); the right panel shows phenotype distribution. Asterisks (*) denote difficult loci modelled separately from the standard analytic detection layer._

### 3.3  Analytic detection by modality and variant class

Detection probabilities are derived from laboratory validation datasets and published reference cohorts. Beta distribution parameters are estimated per modality and variant class. The table below shows mean analytic sensitivity across all five panel variants, ES, and GS. Phenotype-specific panels include bundled difficult-locus assays where clinically relevant: the cystic panel includes a targeted *PKD1* assay; the tubulointerstitial panel includes a targeted *MUC1* assay. Panels without a bundled assay for a difficult locus show 0% sensitivity for that variant class.

| Variant Class | Panel (cystic) | Panel (glom.) | Panel (tub.int.) | Panel (tub.path.) | Panel (CKDu) | ES | GS |
| --- | --- | --- | --- | --- | --- | --- | --- |
| SNV / Indel | 99.6% | 99.6% | 99.6% | 99.6% | 99.6% | 99.6% | 99.9% |
| CNV / SV | 97.4% | 97.4% | 97.4% | 97.4% | 97.4% | 97.4% | 99.6% |
| Difficult locus — PKD1 | 98.9% | 12.5% | 12.5% | 12.5% | 12.5% | 12.5% | 92.9% |
| Difficult locus — MUC1 | 0% | 0% | 95.2% | 0% | 0% | 0% | 0% |

_Sensitivity values are posterior means of Beta-distributed priors. A 0% entry indicates no detection capability is modelled for that modality–variant class combination in the base case (e.g., standard short-read GS and ES for MUC1 VNTR; panels without bundled difficult-locus assay for PKD1 or MUC1)._

### 3.4  VUS probabilities

The probability of receiving ≥1 reportable VUS per proband is modelled as a modality-specific binary outcome. Beta priors are loaded from `data/params/vus_betas.rds`, which was constructed from a curated VUS rate table stratified by panel gene-count bin, ES, and GS. The VUS outcome is drawn independently from the diagnostic outcome.

**Table A.** VUS Beta prior parameters by modality and gene-count bin.

| Modality | Gene Panel Size | α | β | Mean Pr(≥1 VUS) |
| --- | --- | --- | --- | --- |
| Panel | 11–25 genes | 44,900 | 203,168 | 18.1% |
| Panel | 26–50 genes | 115,542 | 258,379 | 30.9% |
| Panel | 51–100 genes | 2,450 | 2,750 | 47.1% |
| Panel | 101–200 genes | 49,230 | 24,139 | 67.1% |
| Panel | >200 genes | 64,250 | 20,068 | 76.2% |
| ES | Exome-wide | 5,569 | 14,603 | 27.6% |
| GS | Genome-wide | 5,569 | 14,603 | 27.6% |

_α, β = Beta distribution shape parameters (method-of-moments fit). Mean Pr(≥1 VUS) = α / (α + β). Panel bins are assigned to testing strategies based on gene-panel size: tubulointerstitial panel → 11–25 gene bin; cystic, glomerular, and tubulopathies panels → 51–100 gene bin; CKDu/comprehensive panel → >200 gene bin._

### 3.5  Unit costs

All costs are in 2024 CAD. Probabilistic uncertainty is applied using a Gamma distribution centred on the point estimate with coefficient of variation (CV) = 0.25.

| Cost Item | Applied When | Mean (CAD) | Range (CAD) | Source |
| --- | --- | --- | --- | --- |
| Targeted multigene panel | Panel-first strategy | $1,400 | $700–$2,100 | Cost Parameters.docx refs 16,17 |
| Clinical exome (or equivalent) | Exome- or broad-test-first strategy | $2,500 | $2,000–$3,000 | Cost Parameters.docx refs 16,17 |
| Clinical genome | Genome-first strategy | $4,500 | $4,000–$5,000 | Cost Parameters.docx refs 16,17 |
| Targeted single-variant test | Familial variant testing plus VUS segregation | $300 | $200–$400 | Cost Parameters.docx |
| Pre-test genetics consultation | Initial genetics visit for proband | $350 | $250–$450 | RAMQ code 09001 |
| Post-test genetics consultation | Follow-up visit(s) | $200 | $150–$250 | RAMQ codes 09001/09003 |

_Cascade testing cost per eligible relative: one pre-test consultation + one targeted familial variant test. VUS follow-up cost per VUS-positive proband: one post-test visit + one targeted familial variant test (single relative)._

### 3.6  Strategy comparison summary

The table below summarises the key analytic parameters for each testing modality variant, allowing direct comparison of detection performance, gene scope, cost, and VUS burden.

**Table B.** Strategy comparison: analytic detection, gene scope, unit cost, and VUS probability.

| Strategy | Phenotype | Panel Version | Genes (n) | SNV/Indel (%) | CNV/SV (%) | PKD1 (%) | MUC1 (%) | Test Cost (CAD) | VUS Pr(≥1) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Panel (cystic) | Cystic | v4.3 | 77 | 99.6% | 97.4% | 98.9% | 0% | $1,400 | 47.1% |
| Panel (glomerular) | Glomerular | v5.2 | 57 | 99.6% | 97.4% | 12.5% | 0% | $1,400 | 47.1% |
| Panel (tubulointerstitial) | Tubulointerstitial | v3.6 | 20 | 99.6% | 97.4% | 12.5% | 95.2% | $1,400 | 18.1% |
| Panel (tubulopathies) | Tubulopathies | v5.6 + nephr. v5.1 | 66 | 99.6% | 97.4% | 12.5% | 0% | $1,400 | 47.1% |
| Panel (CKDu) | CKDu / all | Union all panels | 267 | 99.6% | 97.4% | 12.5% | 0% | $1,400 | 76.2% |
| ES | Any | — | Exome-wide | 99.6% | 97.4% | 12.5% | 0% | $2,500 | 27.6% |
| GS | Any | — | Genome-wide | 99.9% | 99.6% | 92.9% | 0% | $4,500 | 27.6% |

_Detection % values are posterior mean analytic sensitivities. PKD1 and MUC1 columns show difficult-locus detection: 0% indicates the modality lacks a bundled assay or validated pipeline for that locus. VUS Pr(≥1) is the mean of the Beta prior for that modality's gene-count bin._

---

## §4  Base-Case Results

### 4.1  Cost and diagnostic yield by strategy

**Table 1.** Base-case cost and diagnostic outcomes by strategy (mean across 1,000 PSA iterations; 95% uncertainty intervals in brackets).

| Strategy | Total Cost per Proband, Mean [95% UI] | Diagnostic Yield, Mean [95% UI] | Cost per Diagnosis, Mean [95% UI] | Pr(≥1 VUS), Mean [95% UI] |
| --- | --- | --- | --- | --- |
| Panel | $2,554 [$1,863, $3,394] | 30.7% [23.3%, 39.4%] | $8,447 [$5,800, $12,069] | 56.4% [53.3%, 59.4%] |
| ES | $3,527 [$2,473, $4,847] | 31.1% [23.8%, 40.6%] | $11,540 [$7,610, $17,021] | 27.6% [24.7%, 30.5%] |
| Reflex (Panel→ES) | $3,522 [$2,532, $4,694] | 32.5% [24.6%, 41.9%] | $11,021 [$7,392, $15,748] | 60.4% [56.4%, 64.5%] |

> **Key finding:** Panel is the least-cost strategy (mean $2,554/proband) and achieves 30.7% diagnostic yield. At 50% reflex uptake, Reflex (Panel→ES) costs $3,522/proband — slightly cheaper than ES ($3,527) — while achieving higher yield (32.5% vs 31.1%). ES is therefore strictly dominated by Reflex and is excluded from the choice set entirely.

### 4.2  Cost composition

**Table 2.** Mean cost by component, by strategy (base case).

| Strategy | Component | Mean Cost (CAD) | Share of Total |
| --- | --- | --- | --- |
| Panel | Cascade Testing | $398 | 15.6% |
| Panel | Genetics Visits | $547 | 21.4% |
| Panel | Laboratory Testing | $1,406 | 55.1% |
| Panel | VUS Follow-up | $203 | 7.9% |
| ES | Cascade Testing | $403 | 11.4% |
| ES | Genetics Visits | $547 | 15.5% |
| ES | Laboratory Testing | $2,482 | 70.4% |
| ES | VUS Follow-up | $96 | 2.7% |
| Reflex (Panel→ES) | Cascade Testing | $421 | 12.0% |
| Reflex (Panel→ES) | Genetics Visits | $616 | 17.5% |
| Reflex (Panel→ES) | Laboratory Testing | $2,267 | 64.4% |
| Reflex (Panel→ES) | VUS Follow-up | $217 | 6.2% |

![Cost Composition by Strategy](../outputs/results/base_case/cost_composition/mean_cost_composition_by_strategy_base_case.png)

_Figure 1. Stacked bars show mean total diagnostic-episode cost per proband partitioned into cost components. Laboratory testing dominates in all strategies (55–70% of total cost). VUS follow-up is largest under Reflex due to the higher VUS burden associated with combined panel and exome reporting._

### 4.3  Incremental analysis and efficiency frontier

**Table 3.** Dominance analysis and incremental cost per additional diagnosis (ICPD).

| Strategy | Dominance Status | Incremental Cost | Incremental Yield | ICPD (CAD) |
| --- | --- | --- | --- | --- |
| Panel | Non Dominated | — | — | — |
| ES | Strictly Dominated | — | — | — |
| Reflex (Panel→ES) | Non Dominated | $968 | 0.0181 | $53,619 |

_Frontier: Reflex (Panel→ES) vs Panel: +18.1 diagnoses per 1,000 probands at +$968/proband; ICPD $53,619 per additional diagnosis._

![Efficiency Frontier](../outputs/results/incremental_analysis/base_case/efficiency_frontier_base_case.png)

_Figure 2. Efficiency frontier on the cost–diagnosis plane. Each point is a strategy positioned by mean diagnostic yield (x-axis) and mean total cost per proband (y-axis). The frontier line connects non-dominated strategies. ES is strictly dominated by Reflex (lower yield and higher cost). The slope between Panel and Reflex equals the ICPD._

> **Key finding:** The efficiency frontier runs from Panel (lowest cost, lowest yield) to Reflex (highest yield). ES is strictly dominated — Reflex achieves higher yield at lower cost than ES — and is excluded from the frontier. ICPD for Reflex vs Panel: $53,619 per additional diagnosis.

---

## §5  VUS Burden

The model tracks whether each proband receives ≥1 reportable VUS as a secondary outcome. VUS are treated as non-diagnostic in the base case but generate additional follow-up costs. Outcome proportions are approximated from iteration-level estimates using the independence assumption: Pr(VUS only) = Pr(≥1 VUS) × (1 − Pr(diagnosed)); this slightly underestimates the true VUS-only proportion.

**Table 4.** Mean outcome proportions by strategy.

| Strategy | Diagnosed (P/LP) | VUS Only | Negative |
| --- | --- | --- | --- |
| Panel | 30.7% | 39.1% | 30.3% |
| ES | 31.1% | 19.1% | 49.9% |
| Reflex (Panel→ES) | 32.5% | 40.8% | 26.8% |

![Outcome Proportions by Strategy](../outputs/results/base_case/vus_burden/outcome_proportions_stacked.png)

_Figure 3. Stacked horizontal bars showing the proportion of probands in each outcome category: Diagnosed (P/LP), VUS only, and Negative._

> **Key finding:** ES has the lowest VUS burden (27.6% Pr[≥1 VUS]), compared to Panel (56.4%) and Reflex (60.4%). At 50% reflex uptake, the Reflex VUS-only proportion is 40.8% — lower than at full uptake (42.3%) but still the highest of the three strategies.

---

## §6  Uncertainty and Sensitivity Analyses

### 6.1  Probabilistic sensitivity analysis (PSA)

PSA is fully represented by the outer Monte Carlo loop. For each of the 1,000 iterations, we construct the efficiency frontier and record which strategies are non-dominated. The fraction of iterations in which each strategy is non-dominated is the PSA robustness metric.

**Table 5.** Proportion of iterations in which each strategy is non-dominated.

| Strategy | Pr(Non-dominated) |
| --- | --- |
| Panel | 94.7% |
| ES | 22.8% |
| Reflex (Panel→ES) | 95.5% |

![PSA Cost–Diagnosis Scatter](../outputs/results/uncertainty_sensitivity/probabilistic_sensitivity_analysis/psa_cost_diagnosis_scatter.png)

_Figure 4. PSA scatter on the cost–diagnosis plane. Each small point is one Monte Carlo iteration. Large solid circles mark the mean. The cloud spread reflects joint uncertainty in costs and diagnostic yields._

> **Key finding:** Reflex is non-dominated in 95.5% of PSA iterations; Panel in 94.7%. ES is non-dominated in 22.8% of iterations. The ~4.5% of iterations where Reflex is dominated occur when sampled reflex uptake is low, making Reflex yield close to Panel yield at higher cost.

### 6.2  Cost-effectiveness acceptability analysis (CEAC)

The CEAC quantifies the probability that each strategy is cost-effective (highest Net Monetary Benefit) at each willingness-to-pay (WTP) threshold for one additional diagnosis, ranging from $0 to $100,000 CAD.

**NMB** = (WTP × diagnoses per proband) − total cost per proband.

**Table 6.** Probability cost-effective at selected WTP thresholds (CAD per additional diagnosis).

| WTP Threshold (CAD) | Pr Cost-effective: Panel | Pr Cost-effective: ES | Pr Cost-effective: Reflex (Panel→ES) |
| --- | --- | --- | --- |
| $0 | 92.6% | 7.4% | 0.0% |
| $10,000 | 91.5% | 8.5% | 0.0% |
| $20,000 | 89.3% | 10.5% | 0.2% |
| $30,000 | 84.7% | 12.3% | 3.0% |
| $40,000 | 73.1% | 12.3% | 14.6% |
| $50,000 | 54.8% | 12.0% | 33.2% |
| $60,000 | 37.4% | 10.9% | 51.7% |
| $70,000 | 23.8% | 10.4% | 65.8% |
| $80,000 | 16.2% | 10.1% | 73.7% |
| $90,000 | 10.0% | 9.6% | 80.4% |
| $100,000 | 7.0% | 9.1% | 83.9% |

![CEAC](../outputs/results/uncertainty_sensitivity/willingness_to_pay/ceac_plot.png)

_Figure 5. Cost-effectiveness acceptability curve. Lines show the probability each strategy is cost-effective across WTP thresholds from $0 to $100,000 CAD per additional diagnosis._

> **Key finding:** Panel is the most probable cost-effective strategy at WTP below approximately the Reflex ICPD ($53,619/diagnosis). Reflex surpasses Panel near this threshold and reaches 83.9% probability at $100,000/diagnosis. ES retains residual cost-effectiveness probability (up to ~12%) because strict dominance by Reflex is probabilistic, not absolute across all PSA iterations.

---

## §7  Subgroup and Scenario Analyses

### 7.1  Phenotype-dependent cost-effectiveness

Cost-effectiveness of genetic testing is highly sensitive to the phenotypic composition of the referral population, because the a priori probability of a monogenic diagnosis — and therefore the denominator of cost per diagnosis — varies dramatically by phenotype category. The cystic subgroup illustrates this most clearly: the cystic panel bundles a targeted *PKD1* difficult-locus assay, giving panel testing a structural detection advantage over ES in this phenotype.

**Table 7.** Cystic subgroup outcomes by strategy.

| Strategy | Total Cost per Proband, Mean [95% UI] | Diagnostic Yield, Mean [95% UI] | Cost per Diagnosis, Mean [95% UI] |
| --- | --- | --- | --- |
| Panel | $2,828 [$2,088, $3,723] | 60.3% [45.6%, 73.2%] | $4,740 [$3,361, $6,654] |
| ES | $3,752 [$2,684, $5,083] | 50.6% [38.0%, 62.1%] | $7,525 [$5,064, $10,798] |
| Reflex (Panel→ES) | $3,382 [$2,525, $4,346] | 61.3% [46.3%, 74.3%] | $5,603 [$3,940, $7,837] |

In the cystic subgroup, Panel achieves 60.3% diagnostic yield at the lowest cost — ES is strictly dominated (lower yield and higher cost). At 50% reflex uptake, the Reflex cystic cost falls to $3,382 (from $3,942 at full uptake), because only half of cystic panel-negatives receive the ES step. Reflex remains on the frontier with ICPD $57,691 vs the all-phenotype base case ICPD of $53,619. Under Reflex, diagnostic yield ranges from 17.1% (CKDu) to 61.3% (Cystic), while cost per proband varies modestly ($3,382–$3,636). Cost per diagnosis tracks almost entirely with yield, ranging from $5,603 (Cystic) to $29,441 (CKDu).

![Cost vs Yield by Phenotype](../outputs/results/supplement/phenotype_stratified/reflex_cost_yield_by_phenotype.png)

_Figure 6. Mean total cost per proband (y-axis) vs diagnostic yield (x-axis) for the Reflex (Panel→ES) strategy, stratified by phenotype. Points are labeled with phenotype name and mean yield._

> **Key finding:** Cost-effectiveness varies up to 5-fold across phenotype categories. Centres with cystic-enriched referral populations will achieve structurally lower cost per diagnosis; those with predominantly undifferentiated CKD referrals face the highest ICPD. Phenotype mix is a modifiable driver of cost-effectiveness through referral criteria.

### 7.2  Genome sequencing scenario (+10% yield uplift)

Genome sequencing (GS) is excluded from the base case for two reasons. First, diagnostic performance of GS in hereditary kidney disease is insufficiently characterized: few studies have benchmarked GS against gene panels in renal cohorts, and variant interpretation pipelines for difficult loci (e.g., *MUC1* VNTR expansions) are still maturing. Modelling GS alongside well-characterized modalities would require assumptions unsupported by the current evidence base. Second, GS is not currently offered as a clinical diagnostic pathway in Quebec, and the base case reflects testing strategies that are both empirically grounded and operationally available within the system under evaluation. This scenario analysis models a near-future state where GS bioinformatics achieves improved resolution of difficult loci, implemented as a deterministic +10% relative increase in GS diagnostic yield across all phenotypes. Scenario GS is compared against Reflex (Panel→ES) as the incumbent frontier strategy.

**Table 8.** GS scenario vs Reflex: summary outcomes.

| Scenario | Strategy | Diagnostic Yield | Total Cost per Proband | Cost per Diagnosis |
| --- | --- | --- | --- | --- |
| Base Case | Reflex (Panel→ES) | 32.5% | $3,522 | $11,021 |
| GS Scenario (+10% yield) | GS | 37.4% | $5,576 | $15,175 |

_Incremental: GS (scenario) vs Reflex (Panel→ES): delta yield +4.9% (+48.8 per 1,000 probands); delta cost +$2,054/proband; ICPD $42,065 per additional diagnosis._

**Threshold analysis.** Without any yield uplift, GS base yield (34.0%) now exceeds the 50% uptake Reflex yield (32.5%), so GS without uplift is non-dominated but has an impractically high ICPD vs Reflex (~$137,000). A yield improvement of approximately +7% relative to GS base yield (reaching ~36.3% yield) is required for GS ICPD vs Reflex to fall below the Reflex vs Panel ICPD of $53,619. At +10% uplift, the ICPD of $42,065 crosses this threshold comfortably.

![GS Scenario Efficiency Frontier](../outputs/results/scenario_analysis/gs_uplift/efficiency_frontier_gs_scenario.png)

_Figure 7. Efficiency frontier including all base-case strategies and the GS scenario arm. Bubble size represents VUS burden._

![CEAC: GS Scenario](../outputs/results/scenario_analysis/gs_uplift/ceac_plot_gs_scenario.png)

_Figure 8. CEAC including the GS scenario arm alongside base-case strategies (WTP range $0–$100,000 CAD)._

> **Key finding:** Under the +10% yield uplift, GS achieves an ICPD of $42,065 vs Reflex, compared to $53,619 for Reflex vs Panel. A minimum yield improvement of approximately +7% is still required for GS ICPD to fall below the Panel-to-Reflex threshold. GS unit cost (PRCC 0.978) and reflex uptake (PRCC -0.751) are the dominant cost uncertainty drivers.

---

## §8  Limitations and Key Assumptions

1. **Time horizon.** The analysis covers a 12-month diagnostic episode. Downstream clinical consequences (treatment changes, QALYs) and reanalysis value are not captured. The GS yield uplift is a deterministic multiplier applied to observed diagnostic yield, not to detection probability; it is not propagated through the PSA and is applied uniformly across phenotypes.

2. **Evidence base.** Yields and variant frequencies derive from 6 published cohorts (4,738 probands, 2020–2026) with high heterogeneity (I² 76–94%). Published series over-represent specialty centres and solved cases. Only five adult CKD phenotype categories are modelled; CAKUT and aHUS are excluded.

3. **Phenotype composition.** Probands classified as "other/mixed" in source cohorts, including complement-mediated disease and aHUS, are remapped to CKDu. Nephrolithiasis is remapped to tubulopathies. The CKDu panel is defined as the union of all phenotype-specific gene panels, which may not correspond to any single commercial product.

4. **Variant class model.** Causal variants are classified into three classes (SNV/indel, CNV/SV, difficult locus) from 1,978 curated variants across 6 cohorts. Detection is at the modality-by-class level; gene-level sequencing difficulty, CNV size, and variant-specific factors are collapsed. Panel and ES use identical SNV/indel and CNV/SV detection parameters despite panels typically achieving higher per-gene sequencing depth. ES and GS assume 100% gene coverage.

5. **Difficult-locus detection.** The finding that Panel outperforms ES in cystic disease depends on two specific values: 98.9% sensitivity for the bundled PKD1 long-range PCR assay and 12.5% for ES, each from a small number of validation studies and applied uniformly across PKD1 variant types. MUC1 VNTR detection is 0% for ES and GS. Bundled assay costs are assumed zero above the panel unit price.

6. **VUS model.** VUS probabilities derive from a single published dataset, stratified by modality and panel gene count. ES and GS use identical VUS parameters. VUS is binary and drawn independently of diagnostic outcome for all probands. Follow-up costs are applied only to undiagnosed probands with at least one VUS. Gene-specific VUS rates are not captured.

7. **Cost model.** All panels cost $1,400 CAD regardless of gene count (11 to 200+ genes). ES ($2,500) and GS ($4,500) are commercial send-out list prices; in-house testing is substantially cheaper. Consultation costs use Quebec RAMQ fee codes. Cascade testing assumes 100% uptake with 2 eligible relatives and does not distinguish by inheritance pattern; autosomal recessive diagnoses trigger the same cascade cost as autosomal dominant. VUS follow-up is costed as a single segregation test; actual management ranges from no action to extensive functional workup. The impact of varying eligible relatives from 0 to 4 is characterised in Appendix D; the frontier ranking is preserved across this range.

8. **Model structure.** All testing is singleton; trio sequencing is excluded. Diagnosis is binary (P/LP or undiagnosed). Secondary and incidental findings from ES/GS outside the renal indication are not captured. The reflex strategy models 50% uptake of the ES step (base case), parameterized as a Bernoulli trial per proband using a Beta(5,5) prior. Uptake uncertainty and one-way sensitivity are characterized in §6.1 (PSA), Appendix C (PRCC), and Appendix E (one-way DSA; see canonical report).

9. **Generalisability.** Parameters are calibrated to a Quebec tertiary referral setting. Phenotype mix, test pricing, and reimbursement structures vary by jurisdiction.

---

## Appendix A: CHEERS 2022 Reporting Checklist

_Husereau D et al. CHEERS 2022 Statement. BMJ. 2022;376:e067975._

| Item | Element | Where Reported |
| --- | --- | --- |
| **TITLE** | | |
| 1 | Title | Manuscript |
| **ABSTRACT** | | |
| 2 | Abstract | Manuscript |
| **INTRODUCTION** | | |
| 3 | Background and objectives | §1.1 |
| **METHODS** | | |
| 4 | Health economic analysis plan | No pre-registered analysis plan. This is an exploratory model-based economic evaluation. |
| 5 | Study population | §1.3 |
| 6 | Setting and location | §1.3 |
| 7 | Comparators | §1.2 |
| 8 | Perspective | §1.1 |
| 9 | Time horizon | §1.1 |
| 10 | Discount rate | Not applicable — 12-month episode; no discounting required |
| 11 | Selection of outcomes | §1.1 |
| 12 | Measurement of outcomes | §2.2 |
| 13 | Valuation of outcomes | §1.1 — diagnoses counted as a physical outcome unit; not converted to utility. ICPD is the primary measure |
| 14 | Measurement and valuation of resources and costs | §3.5 |
| 15 | Currency, price date, and conversion | §3.5 — 2024 CAD; no currency conversion required |
| 16 | Rationale and description of model | §2.1, §2.2; model code publicly available (Appendix B) |
| 17 | Analytics and assumptions | §2.3, §3, Appendix D, §8 |
| 18 | Characterizing heterogeneity | §7.1 — phenotype-stratified subgroup analysis |
| 19 | Characterizing distributional effects | Not modelled — aggregate public payer perspective adopted; no distributional equity analysis |
| 20 | Characterizing uncertainty | §6.1 (PSA), §6.2 (CEAC), Appendix C (PRCC), §7.2 (GS scenario), Appendix D (cascade DSA) |
| 21 | Approach to engagement with patients and others affected by the study | Nephrologists at the MUHC Renal Genetics Clinic were engaged in parameter selection and model validation throughout the study. No formal patient engagement was undertaken. |
| **RESULTS** | | |
| 22 | Study parameters | §3 — all inputs with distributional parameters (Tables A–B, detection matrix, unit costs) |
| 23 | Summary of main results | §4–§5 (Tables 1–4, Figures 1–3) |
| 24 | Effect of uncertainty | §6.1–§6.2, Appendix C, §7.2, Appendix D |
| 25 | Effect of engagement with patients and others affected by the study | Clinician input shaped the selection of testing strategies, cost parameters, and the rationale for excluding GS from the base case. No patient engagement was conducted and no adjustment to findings was made on that basis. |
| **DISCUSSION** | | |
| 26 | Study findings, limitations, generalizability, and current knowledge | §8; Manuscript discussion |
| **OTHER RELEVANT INFORMATION** | | |
| 27 | Source of funding | Manuscript |
| 28 | Conflicts of interest | Manuscript |

---

## Appendix B: Code Availability

All analysis code, raw data inputs, and parameter files required to reproduce this analysis are publicly available at: **https://github.com/zacsentell/renal-genetic-testing-cea**

The pipeline was executed using **R version 4.4.0** (Puppy Cup). 

The full pipeline can be reproduced by running:

```
Rscript scripts/run_pipeline.R
```

from the repository root. All outputs (figures, tables, and this report) are regenerated deterministically from the raw data inputs in `data/raw/` and the configuration in `config.yml`. Required R packages are listed in `README.md`.

---

## Appendix C: Global Sensitivity Analysis (PRCC)

Partial Rank Correlation Coefficients (PRCC) quantify the contribution of each sampled input parameter to variation in the two key incremental decisions: (A) Reflex vs Panel (primary frontier transition) and (B) GS scenario vs Reflex (secondary frontier transition, §7.2). Table C.1 reports material values across all four analyses. Figure C.1 shows tornado plots for the primary comparison (A). Parameters with |PRCC| < 0.1 or 95% CI crossing zero are omitted.

**Table C.1.** Material PRCC values across all four analyses.

| Domain | Parameter | A: Cost | A: Yield | B: Cost | B: Yield |
| --- | --- | --- | --- | --- | --- |
| Uptake | Reflex ES uptake probability | 0.947 | 0.719 | -0.751 | -0.664 |
| Cost | Clinical exome unit cost | 0.897 | — | -0.621 | — |
| Cost | Clinical genome unit cost | — | — | 0.978 | — |
| Cost | Targeted multigene panel unit cost | — | — | -0.815 | — |
| Yield | Monogenic probability: CKDu | -0.322 | 0.436 | 0.153 | 0.715 |
| Yield | Monogenic probability: Cystic | -0.227 | — | 0.082 | 0.259 |
| Cost | Post-test genetics consultation cost | 0.230 | — | -0.138 | — |
| Yield | Monogenic probability: Glomerular | -0.114 | 0.286 | — | 0.375 |
| Yield | Monogenic probability: Tubulopathies | — | — | 0.088 | 0.136 |

_A = Reflex (Panel→ES) vs Panel; B = GS scenario vs Reflex. Values shown only for parameters with |PRCC| ≥ 0.1 and 95% CI excluding zero in at least one analysis. — = not material. Reflex ES uptake probability is new in this analysis (base case = 50%, Beta(5,5) PSA prior)._

![PRCC Tornado: composite](../outputs/results/uncertainty_sensitivity/global_sensitivity_analysis_prcc/prcc_tornado_composite.png)

_Figure C.1. PRCC tornado plots for the primary frontier transition (Reflex vs Panel). Panel A: incremental cost. Panel B: incremental diagnoses. Bar length = |PRCC|; direction = sign of association. Faded bars have 95% CI crossing zero._

> **Key finding:** Reflex ES uptake probability is the dominant driver of both Reflex vs Panel incremental cost (PRCC 0.947) and incremental yield (PRCC 0.719), displacing clinical exome unit cost (PRCC 0.897). Monogenic probability: CKDu is the second-ranked yield driver (PRCC 0.436). For GS vs Reflex, uptake is also material (PRCC -0.751 on cost; -0.664 on yield): higher reflex uptake makes Reflex a stronger comparator. Clinical genome unit cost (PRCC 0.978) dominates GS scenario cost uncertainty.

---

## Appendix D: Cascade Testing Family-Size Sensitivity

The base case assumes 2 eligible first-degree relatives per diagnosed proband. This one-way DSA varies that count from 0 to 4 by analytical rescaling of the cascade cost component (linear in N). The table below reports the frontier-relevant ICPD for Reflex (Panel→ES) vs Panel with 95% UI from the PSA.

**Table D.1.** ICPD for Reflex (Panel→ES) vs Panel (CAD) [95% UI] by cascade family-size assumption (N = 0–4).

| Eligible relatives (n) | Reflex (Panel→ES) vs Panel |
| --- | --- |
| 0 | $59,070 [$23,982–$130,292] |
| 1 | $59,719 [$24,530–$130,870] |
| **2 (base)** | $60,367 [$25,073–$131,448] |
| 3 | $61,016 [$25,621–$132,166] |
| 4 | $61,665 [$26,171–$133,021] |

_Across the full range of N = 0–4, mean ICPD shifts by $2,595 for Reflex (Panel→ES) vs Panel. The efficiency frontier ranking is unchanged across all values of N. Note: values here are mean of per-iteration ICPD ratios, which differs from the main table estimate ($53,619); frontier ranking conclusion is identical under both approaches._
