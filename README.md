# Renal Genetic CEA

Individual-level microsimulation model for the cost-effectiveness of genetic testing strategies in adults with chronic kidney disease (CKD). Compares targeted gene panels, clinical exome sequencing, and a reflex (panel to exome) strategy from a Canadian payer perspective.

## Requirements

- R ≥ 4.3
- [renv](https://rstudio.github.io/renv/) (dependency management)

## Setup

```r
# In R, from project root:
renv::restore()
```

## Running the pipeline

```bash
Rscript scripts/run_pipeline.R
```

This runs all pipeline steps in sequence (defined in `config.yml`). Default configuration: 1,000 PSA iterations × 1,000 probands.


## Repository structure

```
data/raw/                       Source inputs (read-only): costs, gene panels, curated literature
data/params/                    Derived parameter distributions (.rds)
data/audit/                     Import and validation audit trail (CSV)
scripts/                        Pipeline scripts (00–20) and shared utilities
config.yml                      Pipeline scripts and required outputs
docs/methods.md                 CEA methods (CHEERS 2022)
outputs/figures/                Manual figures (PRISMA, workflow)
outputs/results/parameters/     Tracked parameter reference tables and summaries
outputs/results/base_case/      Iteration-level traces, summary tables, cost composition, VUS burden
outputs/results/incremental_analysis/   Incremental tables and efficiency frontiers
outputs/results/uncertainty_sensitivity/ PRCC, PSA, CEAC, cascade DSA, reflex-uptake DSA
outputs/results/scenario_analysis/      GS uplift, PKD1 full detection, ES augmented
outputs/results/supplement/     Phenotype-stratified and cystic subgroup outputs
outputs/results/clinical_impact/        Gene-level clinical impact analyses
```

## License

MIT — see [LICENSE](LICENSE)
