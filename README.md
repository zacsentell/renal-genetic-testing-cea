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
data/raw/          Source inputs: costs, gene panels, curated literature data
data/params/       Derived parameter distributions (.rds)
data/audit/        Import and validation audit trail (CSV)
scripts/           All pipeline scripts and utilities
docs/methods.md    CEA methods (CHEERS 2022)
config.yml         Pipeline configuration
outputs/parameters/ Parameter summaries (meta-analysis, detection matrix)
```

## License

MIT — see [LICENSE](LICENSE)
