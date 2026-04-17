# Gene Clinical Utility Table — Curation Methodology

This document describes how `gene_clinical_utility.csv` was constructed.
The CSV is the sole clinical-utility input consumed by the pipeline
(`scripts/10_clinical_impact_analysis.R`). It maps each gene modelled
in the simulation to a binary clinical-impact category framework.

## 1. Purpose

Quantify the proportion of diagnosed probands whose molecular result
drives a defined clinical action beyond genetic counselling and
reproductive planning, so that testing strategies can be compared on
clinical value rather than raw diagnostic yield alone.

## 2. Source data

Primary source: **Knoers N, et al. The genetic testing odyssey of adult
patients with chronic kidney disease. Nephrology Dialysis Transplantation
2022.** Five disease-grouped supplemental tables located at
`data/raw/clinical_utility/knoers_2022/Table S{1..5}.xlsx`.

| Table | Disease group | Columns used |
|---|---|---|
| S1 | Glomerulopathies | gene symbol, "Diagnosis utility" |
| S2 | Tubulopathies | gene symbol, "Diagnosis utility" |
| S3 | Complement-mediated | gene symbol, "Diagnosis utility" |
| S4 | CAKUT | gene symbol, "Diagnosis utility" |
| S5 | Ciliopathies | gene symbol, "Diagnosis utility" |

After deduplication across tables, 273 unique Knoers genes were extracted.
Gene symbols were normalized: trailing asterisks stripped, parenthetical
annotations removed (`CASR (inactivating mutations)` → `CASR`),
slash-delimited multi-gene entries split, and non-standard structural-variant
entries excluded.

Secondary source: the paper body, particularly the per-disease
"Clinical benefit of genetic testing" subsections, was used to resolve
ambiguous utility-column entries and to confirm gene-directed
interventions that are not explicit in the tables.

## 3. Category framework

Three non-exclusive primary categories capture distinct clinical actions
triggered by the molecular diagnosis:

1. **Therapeutic** — diagnosis enables a gene-specific treatment decision.
   Partitioned into three non-exclusive sub-categories:
   - **Targeted drug** — initiation of a gene-directed pharmacotherapy
     (e.g., tolvaptan for *PKD1/PKD2*, eculizumab for complement genes,
     CoQ10 for *COQ* pathway, ERT for *GLA*).
   - **Avoid ineffective therapy** — molecular diagnosis rules out a
     treatment that would otherwise be indicated (e.g., immunosuppression
     in hereditary SRNS or Alport, parathyroidectomy in FHH, AL-amyloid
     chemotherapy in hereditary amyloidosis).
   - **Supportive or dietary** — diagnosis drives supportive or dietary
     management (e.g., hyperoxaluria diet, electrolyte replacement in
     tubulopathies, RAS inhibition for Alport).
2. **Extrarenal surveillance** — diagnosis triggers a defined screening
   protocol for an extrarenal complication (hearing, eye, diabetes,
   tumour predisposition, etc.).
3. **Transplant and donor management** — diagnosis alters living-donor
   evaluation, post-transplant recurrence risk, or transplant timing.

The Knoers fifth domain, genetic counselling and reproductive planning,
was excluded because it applies universally to any monogenic diagnosis
and cannot discriminate between strategies.

Diagnostic reclassification was considered but excluded as a standalone
category. Reclassification is only clinically consequential when it
changes management, and that change is already captured by one of the
three action-oriented categories above. Treating reclassification as a
separate category introduced double-counting and vague assignments.

## 4. High-impact rule

`high_impact = 1` if and only if the gene maps to at least one of the
three primary categories. `therapeutic = 1` if and only if at least one
of the three therapeutic sub-categories is 1. Both invariants are
enforced at load time by `scripts/10_clinical_impact_analysis.R`.

## 5. Curation procedure

The committed table is the product of one manual review pass over every
gene, with assignments guided by the following rules.

**Starting cues from the "Diagnosis utility" column.** Specific
action-oriented phrases were taken as evidence for the corresponding
category. Examples: "CoQ10 supplementation" → therapeutic / targeted
drug; "Adapted dietary management" or "Prevention of recurrent
nephrolithiasis" → therapeutic / supportive; "Early detection and
treatment of sensorineural hearing loss", "Ophthalmologic screening",
"diabetes screening" → extrarenal surveillance; "Predict risk of
recurrence in the graft" → transplant and donor management;
"Differential diagnosis with hyperparathyroidism (avoid unnecessary
surgery)" → avoid ineffective therapy.

**Excluded phrases.** Generic phrases were treated as insufficient
evidence and did not, on their own, qualify a gene for any category:
"adequation of treatment and follow-up", any variant of "genetic
counselling", "identification of female carriers", and the phrase "renal
protection strategies" when used without a specific intervention (this
phrase appears for nearly every glomerular and tubulopathy gene in the
Knoers tables).

**Paper body as secondary evidence.** Where the utility column was
ambiguous or generic, the relevant "Clinical benefit of genetic testing"
subsection of the paper body was used to confirm whether a specific
gene-directed intervention exists. Examples include RAS inhibition for
*COL4A3/4/5*, tolvaptan for *PKD1/PKD2*, and eculizumab for *CFH* and
*C3*.

**Manual review against current clinical evidence.** Every Knoers-derived
assignment was manually reviewed against current clinical literature to
confirm it reflects present-day management, not the 2022 snapshot alone.
Where present-day evidence conflicted with the Knoers utility text, the
current evidence was taken as authoritative.

## 6. Extension genes

Twenty-seven genes present in the simulation gene list but absent from
the Knoers tables were annotated from independent clinical literature
and disease-specific management guidelines. High-impact assignments in
this group include *GLA* (Fabry, ERT + surveillance), *TSC1/TSC2*
(everolimus + surveillance), *VHL* and *FLCN* (tumour surveillance),
*WFS1* (Wolfram surveillance), *COL4A1* (cerebrovascular surveillance),
*LDLR* (cardiovascular), *FGA* and *APOA1* (hereditary amyloidosis,
avoiding AL-directed chemotherapy), and *HNF1B*-related genes. Genes
with insufficient direct-action evidence (e.g., *APOL1*, *CUBN*,
*SLC5A2*, *COL4A6* in isolation) were assigned `high_impact = 0`.

## 7. Schema

Columns of `gene_clinical_utility.csv` (15 columns):

| Column | Type | Description |
|---|---|---|
| `gene_symbol` | string | HGNC-approved symbol |
| `high_impact` | 0/1 | Any primary category = 1 |
| `therapeutic` | 0/1 | Any therapeutic sub-category = 1 |
| `targeted_drug` | 0/1 | Therapeutic sub-category |
| `avoid_ineffective_therapy` | 0/1 | Therapeutic sub-category |
| `supportive_dietary` | 0/1 | Therapeutic sub-category |
| `extrarenal_surveillance` | 0/1 | Primary category |
| `transplant_donor` | 0/1 | Primary category |
| `rationale` | string | Free-text justification |
| `knoers_table` | string | Source table(s) in Knoers 2022, or empty |
| `knoers_utility_text` | string | Verbatim utility column text, or empty |
| `notes` | string | Curator comments |
| `in_simulation` | TRUE/FALSE | Gene is modelled in the simulation |
| `needs_review` | TRUE/FALSE | Flagged for follow-up review |
| `disease_groups` | string | Knoers disease group(s), or empty |

**Derived at analysis time (not stored in the CSV).** The `source` field
(`knoers_2022` vs. `extension`) is derived inside
`scripts/10_clinical_impact_analysis.R` from whether `knoers_table` is
non-empty. Do not add a `source` column to the CSV.

**Load-time invariants enforced by `scripts/10_clinical_impact_analysis.R`.**
Any edit to the CSV that breaks these will fail validation:

1. Every category column must contain only `0` or `1`.
2. `therapeutic = 1` requires at least one of `targeted_drug`,
   `avoid_ineffective_therapy`, or `supportive_dietary` to be `1`.
3. `high_impact = 1` if and only if at least one primary category
   (`therapeutic`, `extrarenal_surveillance`, `transplant_donor`) is `1`.

## 8. Provenance of the initial draft

An earlier version of this project used a keyword-matching R script
(`scripts/archive/10a_build_gene_utility_draft.R`) to generate a
first-pass draft of this table directly from the Knoers supplemental
tables. That script is frozen in `scripts/archive/` for provenance only.
It is **not** part of the pipeline and will not produce a table
compatible with the current schema. The committed
`gene_clinical_utility.csv` is authoritative and is the product of
manual curation.
