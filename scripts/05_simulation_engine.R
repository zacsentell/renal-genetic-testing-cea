# 05_simulation_engine.R
# Purpose: Run probabilistic CEA simulation (Monte Carlo)
# Author: Renal Genetics CEA Team
# Date: 2025-12-29

if (!require("dplyr")) install.packages("dplyr")
if (!require("tidyr")) install.packages("tidyr")
if (!require("readr")) install.packages("readr")
if (!require("purrr")) install.packages("purrr")
if (!requireNamespace("config", quietly = TRUE)) {
    stop("Package 'config' is required. Install it with install.packages('config').")
}

library(dplyr)
library(tidyr)
library(readr)
library(purrr)

# Source Cohort Generator (to get generate_cohort and sampling functions)
source("scripts/04_cohort_generator.R")
source("scripts/utils_schema.R")

# ==============================================================================
# IO Paths & Configuration
# ==============================================================================
PARAMS_DIR <- "data/params"
OUTPUT_DIR <- "data/intermediate"
OUTPUT_RESULTS <- file.path(OUTPUT_DIR, "05_cea_results.rds")

# New output paths per methods.md §7.1
OUTPUT_ITERATION_CSV <- "outputs/results/base_case/iteration_level/strategy_iteration_outcomes.csv"
OUTPUT_COMPONENT_CSV <- "outputs/results/base_case/cost_composition/iteration_cost_components.csv"
OUTPUT_PHENOTYPE_CSV <- "outputs/results/supplement/phenotype_stratified/base_case/phenotype_iteration_outcomes.csv"
OUTPUT_PARAM_TRACE_CSV <- "outputs/results/uncertainty_sensitivity/global_sensitivity_analysis_prcc/iteration_parameter_trace.csv"

# Ensure output directories exist
for (out_path in c(
    OUTPUT_ITERATION_CSV, OUTPUT_COMPONENT_CSV, OUTPUT_PHENOTYPE_CSV,
    OUTPUT_PARAM_TRACE_CSV
)) {
    if (!dir.exists(dirname(out_path))) dir.create(dirname(out_path), recursive = TRUE)
}

# Strategy ID mapping (consistent across all outputs)
# NOTE: GS excluded from base case (see methods.md for rationale)
STRATEGY_MAP <- c(
    "Panel" = 1L,
    "ES" = 2L,
    "Panel_Reflex_ES" = 3L
)

# Load runtime configuration from config.yml
cfg <- config::get()

# ==============================================================================
# 1. Load Parameters
# ==============================================================================
cat("Loading parameters from", PARAMS_DIR, "...\n")

read_param <- function(name) {
    path <- file.path(PARAMS_DIR, paste0(name, ".rds"))
    if (!file.exists(path)) stop("Parameter file not found: ", path)
    readRDS(path)
}

params <- list(
    phenotype_alphas = read_param("phenotype_alphas"),
    yield_betas = read_param("yield_betas"),
    gene_variant_joint = read_param("gene_variant_joint"),
    difficult_loci = read_param("difficult_loci"),
    difficult_by_phenotype = read_param("difficult_by_phenotype"),
    detection_betas = read_param("detection_betas"),
    vus_betas = read_param("vus_betas"),
    costs = read_param("costs"),
    panels = read_param("panels"),
    # Cascade params are loaded separately from config.yml to ensure they are
    # never passed through sample_gamma_cost() in the PSA rapply block.
    cascade = read_param("cascade_params"),
    # Uptake probabilities (reflex ES, cascade testing) loaded from
    # data/raw/uptake_parameters.csv via 01b_load_parameters.R.
    uptake = read_param("uptake_params")
)

cat("Parameters loaded.\n")

# ==============================================================================
# 2. Helper Functions
# ==============================================================================

# 2.1 Detection Sampling
# --------------------
sample_detection_matrix <- function(detection_betas) {
    # Expects input: detection_betas$strategy_detection_matrix

    sampled_matrix <- list()

    for (mod in names(detection_betas)) {
        sampled_matrix[[mod]] <- list()

        # If Modality is Panel or ES_Augmented, next level is Phenotype
        if (mod %in% c("Panel", "ES_Augmented")) {
            for (pheno in names(detection_betas[[mod]])) {
                sampled_matrix[[mod]][[pheno]] <- list()
                pheno_obj <- detection_betas[[mod]][[pheno]]

                sampled_matrix[[mod]][[pheno]]$bundled_assays <- pheno_obj$bundled_assays

                for (vc in names(pheno_obj)) {
                    if (vc == "bundled_assays") next
                    vc_obj <- pheno_obj[[vc]]

                    if (vc %in% c("SNV_Indel", "CNV_SV")) {
                        if (isTRUE(vc_obj$applicable)) {
                            prob <- rbeta(1, vc_obj$alpha, vc_obj$beta)
                            sampled_matrix[[mod]][[pheno]][[vc]] <- prob
                        } else {
                            sampled_matrix[[mod]][[pheno]][[vc]] <- 0
                        }
                    } else if (vc == "Difficult_Locus") {
                        sampled_matrix[[mod]][[pheno]][[vc]] <- list()
                        for (gene in names(vc_obj)) {
                            g_obj <- vc_obj[[gene]]
                            if (isTRUE(g_obj$applicable)) {
                                prob <- rbeta(1, g_obj$alpha, g_obj$beta)
                                sampled_matrix[[mod]][[pheno]][[vc]][[gene]] <- prob
                            } else {
                                sampled_matrix[[mod]][[pheno]][[vc]][[gene]] <- 0
                            }
                        }
                    }
                }
            }
        } else {
            # ES / GS (Universal)
            pheno <- "Universal"
            sampled_matrix[[mod]][[pheno]] <- list()
            mod_obj <- detection_betas[[mod]]

            sampled_matrix[[mod]][[pheno]]$bundled_assays <- mod_obj$bundled_assays

            for (vc in names(mod_obj)) {
                if (vc == "bundled_assays") next
                vc_obj <- mod_obj[[vc]]

                if (vc %in% c("SNV_Indel", "CNV_SV")) {
                    if (isTRUE(vc_obj$applicable)) {
                        prob <- rbeta(1, vc_obj$alpha, vc_obj$beta)
                        sampled_matrix[[mod]][[pheno]][[vc]] <- prob
                    } else {
                        sampled_matrix[[mod]][[pheno]][[vc]] <- 0
                    }
                } else if (vc == "Difficult_Locus") {
                    sampled_matrix[[mod]][[pheno]][[vc]] <- list()
                    for (gene in names(vc_obj)) {
                        g_obj <- vc_obj[[gene]]
                        if (isTRUE(g_obj$applicable)) {
                            prob <- rbeta(1, g_obj$alpha, g_obj$beta)
                            sampled_matrix[[mod]][[pheno]][[vc]][[gene]] <- prob
                        } else {
                            sampled_matrix[[mod]][[pheno]][[vc]][[gene]] <- 0
                        }
                    }
                }
            }
        }
    }
    return(sampled_matrix)
}

# 2.2 VUS Sampling
# ---------------
sample_vus_probs <- function(vus_betas) {
    probs <- list()
    probs$panels <- lapply(vus_betas$panels, function(x) {
        rbeta(1, x$shape1, x$shape2)
    })
    probs$es <- rbeta(1, vus_betas$es$shape1, vus_betas$es$shape2)
    probs$gs <- rbeta(1, vus_betas$gs$shape1, vus_betas$gs$shape2)
    return(probs)
}

# 2.3 Cost Sampling
# -----------------
sample_gamma_cost <- function(mean_cost, cv = 0.2) {
    if (is.na(mean_cost) || mean_cost <= 0) {
        return(mean_cost)
    }
    shape <- 1 / (cv^2)
    scale <- mean_cost * (cv^2)
    rgamma(1, shape = shape, scale = scale)
}

# Build per-phenotype panel costs using a single panel unit cost.
build_phenotype_cost_map <- function(panels_list, panel_cost) {
    purrr::imap(panels_list, function(genes, name) {
        list(panel_cost = panel_cost)
    })
}


# 2.3 Vectorized Evaluation
# ------------------------
# Helper to get detection probability vector for the cohort
get_detection_probs <- function(cohort, modality, sampled_detection) {
    # Pre-allocate
    n <- nrow(cohort)
    probs <- numeric(n)

    # If modality is Panel, we need to respect phenotype.
    # If ES/GS, phenotype is Universal.

    # We loop over the unique combinations present in the cohort to minimize lookups
    # Group by Phenotype (if Panel), Variant Class

    # Define mapping strategy
    if (modality %in% c("Panel", "ES_Augmented")) {
        # Loop over phenotypes (Panel and ES_Augmented both use phenotype-specific detection)
        phenos <- unique(cohort$phenotype)
        for (p in phenos) {
            # Get Detection Object for this phenotype
            det_obj <- sampled_detection[[modality]][[p]]
            if (is.null(det_obj)) next

            # Sub-indices
            idx_pheno <- which(cohort$phenotype == p)

            # Simple Classes
            idx_snv <- idx_pheno[cohort$variant_class[idx_pheno] == "SNV_Indel"]
            if (length(idx_snv) > 0) probs[idx_snv] <- det_obj$SNV_Indel

            idx_cnv <- idx_pheno[cohort$variant_class[idx_pheno] == "CNV_SV"]
            if (length(idx_cnv) > 0) probs[idx_cnv] <- det_obj$CNV_SV

            # Difficult Loci - Gene specific
            idx_diff <- idx_pheno[cohort$variant_class[idx_pheno] == "Difficult_Locus"]
            if (length(idx_diff) > 0) {
                # Need to check gene
                sub_cohort <- cohort[idx_diff, ]
                genes <- unique(sub_cohort$difficult_locus_gene)
                for (g in genes) {
                    if (is.na(g)) next

                    p_gene <- 0
                    if (!is.null(det_obj$Difficult_Locus[[g]])) {
                        p_gene <- det_obj$Difficult_Locus[[g]]
                    }

                    idx_g <- idx_diff[sub_cohort$difficult_locus_gene == g]
                    probs[idx_g] <- p_gene
                }
            }
        }
    } else {
        # ES / GS - Universal
        det_obj <- sampled_detection[[modality]][["Universal"]]

        # Simple Classes
        idx_snv <- which(cohort$variant_class == "SNV_Indel")
        probs[idx_snv] <- det_obj$SNV_Indel

        idx_cnv <- which(cohort$variant_class == "CNV_SV")
        probs[idx_cnv] <- det_obj$CNV_SV

        # Difficult Loci
        idx_diff <- which(cohort$variant_class == "Difficult_Locus")
        if (length(idx_diff) > 0) {
            genes <- unique(cohort$difficult_locus_gene[idx_diff])
            for (g in genes) {
                if (is.na(g)) next

                p_gene <- 0
                if (!is.null(det_obj$Difficult_Locus[[g]])) {
                    p_gene <- det_obj$Difficult_Locus[[g]]
                }
                idx_g <- idx_diff[cohort$difficult_locus_gene[idx_diff] == g]
                probs[idx_g] <- p_gene
            }
        }
    }

    return(probs)
}

# Helper to get VUS probability vector
get_vus_probs <- function(cohort, modality, sampled_vus, panel_rec) {
    n <- nrow(cohort)
    if (modality %in% c("ES", "ES_Augmented")) {
        # ES_Augmented uses ES VUS probability (assay co-ordering does not change VUS rate)
        return(rep(sampled_vus$es, n))
    }
    if (modality == "GS") {
        return(rep(sampled_vus$gs, n))
    }

    # Panel - depends on phenotype -> gene list size
    # We pre-calculate map: Phenotype -> VUS Prob
    phenos <- names(panel_rec)
    pheno_vus_map <- numeric(length(phenos))
    names(pheno_vus_map) <- phenos

    for (p in phenos) {
        n_genes <- length(panel_rec[[p]])
        bin <- "bin_gt_200"
        if (n_genes <= 25) {
            bin <- "bin_11_25"
        } else if (n_genes <= 50) {
            bin <- "bin_26_50"
        } else if (n_genes <= 100) {
            bin <- "bin_51_100"
        } else if (n_genes <= 200) bin <- "bin_101_200"

        val <- sampled_vus$panels[[bin]]
        if (is.null(val)) val <- sampled_vus$panels$bin_gt_200
        pheno_vus_map[p] <- val
    }

    # Map to vector
    return(pheno_vus_map[cohort$phenotype])
}

# Helper to check coverage
get_coverage <- function(cohort, modality, panel_rec) {
    if (modality != "Panel") {
        return(rep(TRUE, nrow(cohort)))
    }

    # Vectorized check: is causal_gene in panel_genes_list[phenotype]
    # This is a bit slow row-wise, can map.
    # But since gene lists are static, we can loop phenotypes.

    covered <- logical(nrow(cohort))
    # Default is FALSE.

    phenos <- unique(cohort$phenotype)
    for (p in phenos) {
        genes <- panel_rec[[p]]
        idx <- which(cohort$phenotype == p)
        # Check membership
        # Note: causal_gene can be NA
        covered[idx] <- cohort$causal_gene[idx] %in% genes
    }
    return(covered)
}

# Main Vectorized Strategy Runner
# cascade_n_relatives: number of eligible first-degree relatives per diagnosed
#   proband for cascade testing costing. Kept as an explicit argument (not read
#   from current_costs) so the DSA can sweep this independently of the PSA.
# uptake_reflex: probability that a panel-negative proband proceeds to ES
#   in reflex strategies. Drawn from Beta(5,5) in PSA; swept in DSA.
# uptake_cascade: probability that a diagnosed proband undergoes cascade
#   testing. Fixed at 1.0 in the current analysis; parameterized for future use.
run_simulation_step <- function(cohort, strategy,
                                sampled_detection, sampled_vus, current_costs,
                                panel_genes_list, phenotype_costs_map,
                                cascade_n_relatives = 2,
                                uptake_reflex = 1.0,
                                uptake_cascade = 1.0) {
    n_probands <- nrow(cohort)

    # 1. Base Costs (Consultations) in a vector
    # ----------------------------------------
    total_costs <- rep(current_costs$consultations$pretest + current_costs$consultations$posttest, n_probands)
    is_diagnosed <- logical(n_probands) # FALSE
    has_vus <- logical(n_probands) # FALSE

    # Helper to apply a test
    # Returns a list of vectors: diagnosed, vus, cost_added
    apply_test <- function(modality, mask_eligible) {
        # mask_eligible: boolean vector of who gets this test

        n_test <- sum(mask_eligible)
        if (n_test == 0) {
            return(list(diag = logical(n_probands), vus = logical(n_probands), cost = numeric(n_probands)))
        }

        # Sub-cohort for calcs
        # NOTE: to keep alignment, we calculate for all, but only use masked

        # A. Cost
        cost_vec <- numeric(n_probands)
        if (modality == "Panel") {
            # Map phenotype to cost
            # Optimization: Create named vector
            p_cost_map <- sapply(phenotype_costs_map, function(x) x$panel_cost)
            cost_vec <- p_cost_map[cohort$phenotype]
        } else if (modality %in% c("ES", "ES_Augmented")) {
            # ES_Augmented uses same unit cost as ES; bundled assay costs absorbed
            cost_vec <- rep(current_costs$tests$es, n_probands)
        } else if (modality == "GS") {
            cost_vec <- rep(current_costs$tests$gs, n_probands)
        }
        cost_vec[!mask_eligible] <- 0

        # B. Diagnosis
        # Logic: Monogenic & Covered & Detected
        covered <- get_coverage(cohort, modality, panel_genes_list)
        probs <- get_detection_probs(cohort, modality, sampled_detection)

        rng <- runif(n_probands)
        detected <- cohort$is_monogenic & covered & (rng < probs)
        detected[!mask_eligible] <- FALSE

        # C. VUS
        # Logic: runif < vus_prob
        vus_p_vec <- get_vus_probs(cohort, modality, sampled_vus, panel_genes_list)
        rng_vus <- runif(n_probands)
        vus_res <- rng_vus < vus_p_vec
        vus_res[!mask_eligible] <- FALSE

        list(diag = detected, vus = vus_res, cost = cost_vec)
    }

    # 2. Execute Strategy
    # -------------------
    reflex_mask <- rep(FALSE, n_probands) # tracks who actually received reflex ES

    if (strategy %in% c("Panel", "ES", "GS", "ES_Augmented")) {
        res <- apply_test(strategy, rep(TRUE, n_probands))

        is_diagnosed <- res$diag
        has_vus <- res$vus
        total_costs <- total_costs + res$cost
        cost_testing <- res$cost
    } else if (strategy == "Panel_Reflex_ES") {
        # Step 1: Panel (All)
        res1 <- apply_test("Panel", rep(TRUE, n_probands))

        # Update State
        is_diagnosed <- res1$diag
        has_vus <- res1$vus # Tracks if VUS occurred in Step 1
        total_costs <- total_costs + res1$cost
        cost_testing <- res1$cost

        # Step 2: Reflex (Only Undiagnosed who accept reflex testing)
        mask_step2 <- !is_diagnosed & (runif(n_probands) < uptake_reflex)
        reflex_mask <- mask_step2
        if (any(mask_step2)) {
            res2 <- apply_test("ES", mask_step2)

            # Update diagnosis
            is_diagnosed[mask_step2] <- res2$diag[mask_step2]

            # Update VUS (Cumulative: if VUS in either step)
            has_vus[mask_step2] <- has_vus[mask_step2] | res2$vus[mask_step2]

            # Reflex uses full exome unit cost for panel-negative probands.
            total_costs <- total_costs + res2$cost
            cost_testing <- cost_testing + res2$cost

            # Add extra disclosure visit for reflex step
            # Probands who proceed to reflex get: pretest + panel disclosure + ES disclosure
            total_costs[mask_step2] <- total_costs[mask_step2] + current_costs$consultations$posttest
        }
    } else if (strategy == "Reflex_Augmented") {
        # Step 1: Panel (all probands) — unchanged, bundled assays via panel_by_phenotype
        res1 <- apply_test("Panel", rep(TRUE, n_probands))
        is_diagnosed <- res1$diag
        has_vus <- res1$vus
        total_costs <- total_costs + res1$cost
        cost_testing <- res1$cost

        # Step 2: ES_Augmented (panel-negatives who accept reflex testing)
        # PKD1 assay universal; MUC1 assay for tubulointerstitial.
        # Note: for cystic/tubulointerstitial probands the relevant assay is applied at
        # both steps; <2% of cystic probands reach step 2 (panel PKD1 sensitivity ~99%),
        # so marginal yield inflation is negligible (minor approximation).
        mask_step2 <- !is_diagnosed & (runif(n_probands) < uptake_reflex)
        reflex_mask <- mask_step2
        if (any(mask_step2)) {
            res2 <- apply_test("ES_Augmented", mask_step2)
            is_diagnosed[mask_step2] <- is_diagnosed[mask_step2] | res2$diag[mask_step2]
            has_vus[mask_step2] <- has_vus[mask_step2] | res2$vus[mask_step2]
            total_costs <- total_costs + res2$cost
            cost_testing <- cost_testing + res2$cost
            # Extra post-test disclosure for reflex probands
            total_costs[mask_step2] <- total_costs[mask_step2] + current_costs$consultations$posttest
        }
    }

    # Track which probands had reflex (for consultation cost tracking).
    # Uses reflex_mask which reflects both panel-negative status AND uptake acceptance.
    had_reflex <- reflex_mask

    # 3. Downstream Outcomes (with per-component cost tracking)
    # ---------------------------------------------------------

    # Per-component cost tracking (§7.4)

    # Consultation costs: pretest + posttest for all probands
    # Plus extra posttest for reflex probands (already added to total_costs above)
    cost_consultation <- rep(current_costs$consultations$pretest +
        current_costs$consultations$posttest, n_probands)
    cost_consultation[had_reflex] <- cost_consultation[had_reflex] +
        current_costs$consultations$posttest

    # A. Cascade Testing (If Diagnosed AND accepts cascade testing)
    # cascade_n_relatives * (pretest + familial_test)
    # n_relatives is passed explicitly as cascade_n_relatives to keep it out of
    # the Gamma-sampling rapply block; it is an integer count, not a unit cost.
    cost_cascade_unit <- current_costs$consultations$pretest + current_costs$tests$familial
    cascade_accepted <- is_diagnosed & (runif(n_probands) < uptake_cascade)
    cost_cascade <- ifelse(cascade_accepted, cascade_n_relatives * cost_cascade_unit, 0)
    total_costs <- total_costs + cost_cascade

    # B. VUS Follow-up (If VUS AND NOT Diagnosed)
    # (posttest + familial_test)
    cost_vus_unit <- current_costs$consultations$posttest + current_costs$tests$familial

    cost_vus_followup <- ifelse(has_vus & !is_diagnosed, cost_vus_unit, 0)
    total_costs <- total_costs + cost_vus_followup

    # Return enhanced output with per-component costs
    list(
        # Primary outcomes
        diagnoses_per_proband = mean(is_diagnosed),
        total_cost_per_proband_cad = mean(total_costs),
        pr_at_least_one_vus = mean(has_vus),

        # Per-component costs (§7.4)
        cost_consultation_mean = mean(cost_consultation),
        cost_testing_mean = mean(cost_testing),
        cost_cascade_mean = mean(cost_cascade),
        cost_vus_followup_mean = mean(cost_vus_followup),

        # Proband-level vectors for detailed analysis
        diagnosed_vec = is_diagnosed,
        cost_vec = total_costs,
        vus_vec = has_vus,
        phenotype_vec = cohort$phenotype, # For phenotype-stratified analysis

        # Component cost vectors (added for phenotype-stratified cost composition)
        cost_consultation_vec = cost_consultation,
        cost_testing_vec = cost_testing,
        cost_cascade_vec = cost_cascade,
        cost_vus_followup_vec = cost_vus_followup
    )
}

if (sys.nframe() == 0) {
    # ==============================================================================
    # 3. Main Simulation Loop
    # ==============================================================================

    # Configuration
    N_ITER <- cfg$simulation$n_iter
    N_PROBANDS <- cfg$simulation$n_probands
    SEED_START <- cfg$simulation$seed_start
    COST_CV <- cfg$simulation$cost_cv

    # Strategies to run (GS excluded from base case - see scenario analysis)
    STRATEGIES <- c("Panel", "ES", "Panel_Reflex_ES")

    cat(sprintf(
        "Starting Simulation: %d Iterations, %d Probands, %d Strategies\n",
        N_ITER, N_PROBANDS, length(STRATEGIES)
    ))

    # Storage for results - using new column names per §7.1
    iter_results_list <- vector("list", N_ITER * length(STRATEGIES))
    cost_components_list <- vector("list", N_ITER * length(STRATEGIES))
    phenotype_results_list <- list() # 5.5.5: Phenotype-stratified outcomes
    param_trace_list <- vector("list", N_ITER) # 5.5.6: Parameter trace for PRCC
    list_idx <- 0
    pheno_idx <- 0

    set.seed(SEED_START)
    iter_seeds <- sample.int(1e9, N_ITER)

    # Get phenotype categories
    PHENOTYPE_CATEGORIES <- names(params$phenotype_alphas)

    for (i in 1:N_ITER) {
        if (i %% 50 == 0) cat(sprintf("  Processing Iteration %d / %d...\n", i, N_ITER))

        curr_seed <- iter_seeds[i]
        set.seed(curr_seed)

        # 1. Sample parameters for this iteration
        s_det <- sample_detection_matrix(params$detection_betas$strategy_detection_matrix)
        s_vus <- sample_vus_probs(params$vus_betas)

        # Sample costs using Gamma distribution.
        # curr_costs will mirror params$costs structure
        curr_costs <- rapply(params$costs, function(x) {
            if (is.numeric(x) && length(x) == 1) {
                sample_gamma_cost(x, cv = COST_CV)
            } else {
                x
            }
        }, how = "replace")

        sampled_yields <- sapply(PHENOTYPE_CATEGORIES, function(p) {
            bp <- params$yield_betas[[p]]$params
            if (is.null(bp)) {
                return(0)
            }
            rbeta(1, shape1 = bp["shape1"], shape2 = bp["shape2"])
        })

        # Sample uptake probabilities (Beta if parameterized, else fixed)
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

        # 2. Cohort using the exact sampled yields for this iteration.
        cohort_res <- generate_cohort(
            n_probands = N_PROBANDS,
            params = params,
            sampled_yields = sampled_yields,
            seed = curr_seed
        )
        cohort <- cohort_res$cohort
        sampled_pheno_props <- cohort_res$sampled_proportions

        # Build phenotype panel costs for this iteration using one panel unit cost.
        curr_pheno_costs <- build_phenotype_cost_map(params$panels, curr_costs$tests$panel)

        iteration_truth <- list(
            sampled_pheno_props = sampled_pheno_props,
            sampled_yields = cohort_res$sampled_yields,
            sampled_detection = s_det,
            sampled_vus = s_vus,
            sampled_costs = curr_costs
        )

        # Extract VUS probabilities
        vus_panel_vals <- unlist(s_vus$panels)
        prob_vus_panel_mean <- if (length(vus_panel_vals) > 0) mean(vus_panel_vals, na.rm = TRUE) else NA_real_

        # Extract detection sensitivities for ES and GS
        det_snv_es <- s_det$ES$Universal$SNV_Indel
        det_cnv_es <- s_det$ES$Universal$CNV_SV
        det_snv_gs <- s_det$GS$Universal$SNV_Indel
        det_cnv_gs <- s_det$GS$Universal$CNV_SV

        # Extract difficult locus detection for ES and GS
        det_pkd1_es <- if (!is.null(s_det$ES$Universal$Difficult_Locus$PKD1)) s_det$ES$Universal$Difficult_Locus$PKD1 else NA_real_
        det_pkd1_gs <- if (!is.null(s_det$GS$Universal$Difficult_Locus$PKD1)) s_det$GS$Universal$Difficult_Locus$PKD1 else NA_real_
        det_muc1_es <- if (!is.null(s_det$ES$Universal$Difficult_Locus$MUC1)) s_det$ES$Universal$Difficult_Locus$MUC1 else NA_real_
        det_muc1_gs <- if (!is.null(s_det$GS$Universal$Difficult_Locus$MUC1)) s_det$GS$Universal$Difficult_Locus$MUC1 else NA_real_

        # Extract Panel difficult locus detection (bundled assays)
        det_pkd1_panel <- if (!is.null(s_det$Panel$cystic$Difficult_Locus$PKD1)) s_det$Panel$cystic$Difficult_Locus$PKD1 else NA_real_
        det_muc1_panel <- if (!is.null(s_det$Panel$tubulointerstitial$Difficult_Locus$MUC1)) s_det$Panel$tubulointerstitial$Difficult_Locus$MUC1 else NA_real_

        # Panel detection aggregated as phenotype-weighted means
        panel_phenos <- names(s_det$Panel)
        panel_snv_vals <- sapply(panel_phenos, function(p) {
            val <- s_det$Panel[[p]][["SNV_Indel"]]
            if (is.null(val)) NA_real_ else val
        })
        panel_cnv_vals <- sapply(panel_phenos, function(p) {
            val <- s_det$Panel[[p]][["CNV_SV"]]
            if (is.null(val)) NA_real_ else val
        })

        # Use phenotype proportions as weights (only for phenotypes present in both)
        common_phenos <- intersect(names(sampled_pheno_props), panel_phenos)
        if (length(common_phenos) > 0) {
            weights <- sampled_pheno_props[common_phenos]
            weights <- weights / sum(weights, na.rm = TRUE) # Renormalize
            det_snv_panel_weighted <- sum(panel_snv_vals[common_phenos] * weights, na.rm = TRUE)
            det_cnv_panel_weighted <- sum(panel_cnv_vals[common_phenos] * weights, na.rm = TRUE)
        } else {
            det_snv_panel_weighted <- mean(panel_snv_vals, na.rm = TRUE)
            det_cnv_panel_weighted <- mean(panel_cnv_vals, na.rm = TRUE)
        }

        # Combine all traced parameters (expanded set)
        trace_params <- c(
            # Phenotype prevalence (Dirichlet-sampled)
            setNames(sampled_pheno_props, paste0("prev_", names(sampled_pheno_props))),
            # Diagnostic yields
            setNames(iteration_truth$sampled_yields[PHENOTYPE_CATEGORIES], paste0("yield_", PHENOTYPE_CATEGORIES)),
            # VUS probabilities
            prob_vus_panel_mean = prob_vus_panel_mean,
            prob_vus_es = iteration_truth$sampled_vus$es,
            prob_vus_gs = iteration_truth$sampled_vus$gs,
            # ES/GS detection (standard variant classes)
            det_snv_es = det_snv_es,
            det_cnv_es = det_cnv_es,
            det_snv_gs = det_snv_gs,
            det_cnv_gs = det_cnv_gs,
            # Panel detection (phenotype-weighted)
            det_snv_panel_weighted = det_snv_panel_weighted,
            det_cnv_panel_weighted = det_cnv_panel_weighted,
            # Difficult locus detection
            det_pkd1_es = det_pkd1_es,
            det_pkd1_gs = det_pkd1_gs,
            det_muc1_es = det_muc1_es,
            det_muc1_gs = det_muc1_gs,
            det_pkd1_panel = det_pkd1_panel,
            det_muc1_panel = det_muc1_panel,
            # Unit Costs (Gamma Sampled)
            cost_panel = iteration_truth$sampled_costs$tests$panel,
            cost_es = iteration_truth$sampled_costs$tests$es,
            cost_gs = iteration_truth$sampled_costs$tests$gs,
            cost_familial = iteration_truth$sampled_costs$tests$familial,
            cost_consult_pre = iteration_truth$sampled_costs$consultations$pretest,
            cost_consult_post = iteration_truth$sampled_costs$consultations$posttest,
            # Uptake probabilities
            uptake_reflex = uptake_reflex_i,
            uptake_cascade = uptake_cascade_i
        )

        # Store parameter trace (5.5.6)
        param_trace_list[[i]] <- data.frame(
            iteration_id = i,
            parameter_name = names(trace_params),
            parameter_value = unname(trace_params),
            stringsAsFactors = FALSE
        )

        for (strat in STRATEGIES) {
            list_idx <- list_idx + 1

            # Run Vectorized Step
            res <- run_simulation_step(
                cohort, strat, iteration_truth$sampled_detection, iteration_truth$sampled_vus, iteration_truth$sampled_costs,
                params$panels, curr_pheno_costs,
                cascade_n_relatives = params$cascade$eligible_relatives_base,
                uptake_reflex = uptake_reflex_i,
                uptake_cascade = uptake_cascade_i
            )

            # Compute derived columns
            cost_per_diagnosis <- if (res$diagnoses_per_proband > 0) {
                res$total_cost_per_proband_cad / res$diagnoses_per_proband
            } else {
                NA_real_
            }

            # Store iteration outcomes (§7.1 spec)
            iter_results_list[[list_idx]] <- data.frame(
                iteration_id = i,
                strategy_id = STRATEGY_MAP[[strat]],
                strategy_label = strat,
                total_cost_per_proband_cad = res$total_cost_per_proband_cad,
                diagnoses_per_proband = res$diagnoses_per_proband,
                cost_per_diagnosis_cad = cost_per_diagnosis,
                pr_at_least_one_vus = res$pr_at_least_one_vus,
                stringsAsFactors = FALSE
            )

            # Store cost components (§7.4)
            cost_components_list[[list_idx]] <- data.frame(
                iteration_id = i,
                strategy_id = STRATEGY_MAP[[strat]],
                strategy_label = strat,
                cost_consultation_cad = res$cost_consultation_mean,
                cost_testing_cad = res$cost_testing_mean,
                cost_cascade_cad = res$cost_cascade_mean,
                cost_vus_followup_cad = res$cost_vus_followup_mean,
                total_cost_per_proband_cad = res$total_cost_per_proband_cad,
                stringsAsFactors = FALSE
            )

            # 5.5.5: Phenotype-stratified aggregation (§7.7)
            for (pheno in PHENOTYPE_CATEGORIES) {
                pheno_idx <- pheno_idx + 1
                pheno_mask <- res$phenotype_vec == pheno
                n_pheno <- sum(pheno_mask)

                if (n_pheno > 0) {
                    pheno_yield <- mean(res$diagnosed_vec[pheno_mask])
                    pheno_cost <- mean(res$cost_vec[pheno_mask])
                    pheno_cpd <- if (pheno_yield > 0) pheno_cost / pheno_yield else NA_real_

                    # Component costs
                    pheno_consult <- mean(res$cost_consultation_vec[pheno_mask])
                    pheno_testing <- mean(res$cost_testing_vec[pheno_mask])
                    pheno_cascade <- mean(res$cost_cascade_vec[pheno_mask])
                    pheno_vus <- mean(res$cost_vus_followup_vec[pheno_mask])
                } else {
                    pheno_yield <- NA_real_
                    pheno_cost <- NA_real_
                    pheno_cpd <- NA_real_

                    pheno_consult <- NA_real_
                    pheno_testing <- NA_real_
                    pheno_cascade <- NA_real_
                    pheno_vus <- NA_real_
                }

                phenotype_results_list[[pheno_idx]] <- data.frame(
                    iteration_id = i,
                    phenotype_category = pheno,
                    strategy_id = STRATEGY_MAP[[strat]],
                    strategy_label = strat,
                    n_probands = n_pheno,

                    # Outcomes
                    total_cost_per_proband_cad = pheno_cost,
                    diagnoses_per_proband = pheno_yield,
                    cost_per_diagnosis_cad = pheno_cpd,

                    # Component Costs
                    cost_consultation_cad = pheno_consult,
                    cost_testing_cad = pheno_testing,
                    cost_cascade_cad = pheno_cascade,
                    cost_vus_followup_cad = pheno_vus,
                    stringsAsFactors = FALSE
                )
            }
        }
    }

    # Combine results
    iter_results_df <- do.call(rbind, iter_results_list)
    cost_components_df <- do.call(rbind, cost_components_list)
    phenotype_results_df <- do.call(rbind, phenotype_results_list)
    param_trace_df <- do.call(rbind, param_trace_list)

    # ==============================================================================
    # 4. Export Results
    # ==============================================================================
    iter_required <- c(
        "iteration_id", "strategy_id", "strategy_label",
        "total_cost_per_proband_cad", "diagnoses_per_proband",
        "cost_per_diagnosis_cad", "pr_at_least_one_vus"
    )
    component_required <- c(
        "iteration_id", "strategy_id", "strategy_label",
        "cost_consultation_cad", "cost_testing_cad",
        "cost_cascade_cad", "cost_vus_followup_cad", "total_cost_per_proband_cad"
    )
    pheno_required <- c(
        "iteration_id", "phenotype_category", "strategy_id", "strategy_label", "n_probands",
        "total_cost_per_proband_cad", "diagnoses_per_proband", "cost_per_diagnosis_cad",
        "cost_consultation_cad", "cost_testing_cad", "cost_cascade_cad", "cost_vus_followup_cad"
    )
    trace_required <- c("iteration_id", "parameter_name", "parameter_value")

    assert_required_columns(iter_results_df, iter_required, "strategy_iteration_outcomes", exact = TRUE)
    assert_required_columns(cost_components_df, component_required, "iteration_cost_components", exact = TRUE)
    assert_required_columns(phenotype_results_df, pheno_required, "phenotype_iteration_outcomes", exact = TRUE)
    assert_required_columns(param_trace_df, trace_required, "iteration_parameter_trace", exact = TRUE)

    assert_no_na(iter_results_df, c("iteration_id", "strategy_id", "strategy_label"), "strategy_iteration_outcomes")

    saveRDS(list(
        iterations = iter_results_df,
        components = cost_components_df,
        phenotype = phenotype_results_df,
        parameters = param_trace_df
    ), OUTPUT_RESULTS)

    write_csv_validated(iter_results_df, OUTPUT_ITERATION_CSV, "strategy_iteration_outcomes")
    write_csv_validated(cost_components_df, OUTPUT_COMPONENT_CSV, "iteration_cost_components")
    write_csv_validated(phenotype_results_df, OUTPUT_PHENOTYPE_CSV, "phenotype_iteration_outcomes")
    write_csv_validated(param_trace_df, OUTPUT_PARAM_TRACE_CSV, "iteration_parameter_trace")

    cat("\nSimulation Complete.\n")
    cat("Results saved to:\n")
    cat("- ", OUTPUT_RESULTS, "\n")
    cat("- ", OUTPUT_ITERATION_CSV, "\n")
    cat("- ", OUTPUT_COMPONENT_CSV, "\n")
    cat("- ", OUTPUT_PHENOTYPE_CSV, "\n")
    cat("- ", OUTPUT_PARAM_TRACE_CSV, "\n")
} # end if (sys.nframe() == 0)
