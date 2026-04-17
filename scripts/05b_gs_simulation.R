# 05b_gs_simulation.R
# Purpose: GS baseline simulation for scenario analysis (§8.5)
# Produces gs_baseline_iteration_outcomes.csv consumed by 09a_gs_scenario_analysis.R
#
# Reproducibility note: iter_seeds are derived identically to 05_simulation_engine.R
# (set.seed(SEED_START); sample.int(1e9, N_ITER)), and parameters are sampled in the
# same order per iteration, ensuring cohorts and outcomes are identical to what was
# previously computed inside the main simulation loop.

library(dplyr)
library(readr)

# Source main engine to load all helper functions and parameters.
# The if (sys.nframe() == 0) guard in 05 ensures the main loop does not execute.
source("scripts/05_simulation_engine.R")

# ==============================================================================
# IO Paths
# ==============================================================================
GS_STRATEGY_ID <- 4L
OUTPUT_DIR <- "outputs/results/scenario_analysis/gs_uplift"
OUTPUT_GS_BASELINE_ITER_CSV <- file.path(OUTPUT_DIR, "gs_baseline_iteration_outcomes.csv")
OUTPUT_GS_GENE_COUNTS_CSV <- "outputs/results/clinical_impact/gs_iteration_gene_counts.csv"

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
if (!dir.exists(dirname(OUTPUT_GS_GENE_COUNTS_CSV))) {
    dir.create(dirname(OUTPUT_GS_GENE_COUNTS_CSV), recursive = TRUE)
}

# ==============================================================================
# Configuration — must match 05_simulation_engine.R exactly
# ==============================================================================
N_ITER     <- cfg$simulation$n_iter
N_PROBANDS <- cfg$simulation$n_probands
SEED_START <- cfg$simulation$seed_start
COST_CV    <- cfg$simulation$cost_cv

# Derive iter_seeds identically to 05_simulation_engine.R
set.seed(SEED_START)
iter_seeds <- sample.int(1e9, N_ITER)

PHENOTYPE_CATEGORIES <- names(params$phenotype_alphas)

cat(sprintf(
    "Starting GS simulation: %d iterations, %d probands\n",
    N_ITER, N_PROBANDS
))

# ==============================================================================
# GS Simulation Loop
# ==============================================================================
gs_results_list <- vector("list", N_ITER)
gs_gene_count_list <- list()
gs_gene_count_idx <- 0

for (i in 1:N_ITER) {
    if (i %% 50 == 0) cat(sprintf("  Processing Iteration %d / %d...\n", i, N_ITER))

    curr_seed <- iter_seeds[i]
    set.seed(curr_seed)

    # Sample parameters in exactly the same order as 05_simulation_engine.R
    s_det <- sample_detection_matrix(params$detection_betas$strategy_detection_matrix)
    s_vus <- sample_vus_probs(params$vus_betas)

    curr_costs <- rapply(params$costs, function(x) {
        if (is.numeric(x) && length(x) == 1) {
            sample_gamma_cost(x, cv = COST_CV)
        } else {
            x
        }
    }, how = "replace")

    sampled_yields <- sapply(PHENOTYPE_CATEGORIES, function(p) {
        bp <- params$yield_betas[[p]]$params
        if (is.null(bp)) return(0)
        rbeta(1, shape1 = bp["shape1"], shape2 = bp["shape2"])
    })

    # Sample uptake probabilities (maintains RNG alignment with 05_simulation_engine.R)
    uptake_reflex_i <- if (!is.na(params$uptake$reflex$beta_shape1)) {
        rbeta(1, params$uptake$reflex$beta_shape1, params$uptake$reflex$beta_shape2)
    } else {
        params$uptake$reflex$base_case
    }
    uptake_cascade_i <- if (!is.na(params$uptake$cascade$beta_shape1)) {
        rbeta(1, params$uptake$cascade$beta_shape1, params$uptake$cascade$beta_shape2)
    } else {
        params$uptake$cascade$base_case
    }

    cohort_res <- generate_cohort(
        n_probands = N_PROBANDS,
        params = params,
        sampled_yields = sampled_yields,
        seed = curr_seed
    )
    cohort <- cohort_res$cohort

    curr_pheno_costs <- build_phenotype_cost_map(params$panels, curr_costs$tests$panel)

    # Run GS strategy
    res <- run_simulation_step(
        cohort, "GS", s_det, s_vus, curr_costs, params$panels, curr_pheno_costs,
        cascade_n_relatives = params$cascade$eligible_relatives_base,
        uptake_reflex = uptake_reflex_i,
        uptake_cascade = uptake_cascade_i
    )

    cost_per_diagnosis <- if (res$diagnoses_per_proband > 0) {
        res$total_cost_per_proband_cad / res$diagnoses_per_proband
    } else {
        NA_real_
    }

    gs_results_list[[i]] <- data.frame(
        iteration_id = i,
        strategy_id = GS_STRATEGY_ID,
        strategy_label = "GS",
        scenario_label = "baseline",
        total_cost_per_proband_cad = res$total_cost_per_proband_cad,
        diagnoses_per_proband = res$diagnoses_per_proband,
        cost_per_diagnosis_cad = cost_per_diagnosis,
        pr_at_least_one_vus = res$pr_at_least_one_vus,
        stringsAsFactors = FALSE
    )

    # Gene-level diagnosed counts for clinical impact analysis
    # (pre-tabulated inside run_simulation_step).
    gene_counts <- res$diagnosed_gene_counts
    gene_counts$iteration_id <- i
    gene_counts$strategy_id <- GS_STRATEGY_ID
    gene_counts$strategy_label <- "GS"
    gene_counts$n_probands <- N_PROBANDS
    gs_gene_count_idx <- gs_gene_count_idx + 1
    gs_gene_count_list[[gs_gene_count_idx]] <- gene_counts
}

gs_baseline_df <- do.call(rbind, gs_results_list)
gs_gene_count_df <- do.call(rbind, gs_gene_count_list)

# ==============================================================================
# Validate and Write
# ==============================================================================
gs_required <- c(
    "iteration_id", "strategy_id", "strategy_label", "scenario_label",
    "total_cost_per_proband_cad", "diagnoses_per_proband", "cost_per_diagnosis_cad",
    "pr_at_least_one_vus"
)
assert_required_columns(gs_baseline_df, gs_required, "gs_baseline_iteration_outcomes", exact = TRUE)
assert_no_na(
    gs_baseline_df,
    c("iteration_id", "strategy_id", "strategy_label", "scenario_label",
      "total_cost_per_proband_cad", "diagnoses_per_proband"),
    "gs_baseline_iteration_outcomes"
)
write_csv_validated(gs_baseline_df, OUTPUT_GS_BASELINE_ITER_CSV, "gs_baseline_iteration_outcomes")

# Gene counts for clinical impact analysis
gs_gene_count_required <- c("gene_symbol", "phenotype", "n_diagnosed",
                            "iteration_id", "strategy_id", "strategy_label",
                            "n_probands")
assert_required_columns(gs_gene_count_df, gs_gene_count_required,
                        "gs_iteration_gene_counts", exact = TRUE)
write_csv_validated(gs_gene_count_df, OUTPUT_GS_GENE_COUNTS_CSV, "gs_iteration_gene_counts")

cat("\nGS simulation complete.\n")
cat("Results saved to:", OUTPUT_GS_BASELINE_ITER_CSV, "\n")
cat("Gene counts saved to:", OUTPUT_GS_GENE_COUNTS_CSV, "\n")
