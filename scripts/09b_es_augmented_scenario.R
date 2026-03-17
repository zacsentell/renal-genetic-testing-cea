# 09b_es_augmented_scenario.R
# Purpose: Augmented ES scenario analysis (§7.4)
#   Tests whether the efficiency frontier changes when ES is co-ordered with
#   phenotype-directed difficult-locus assays: PKD1 long-range PCR applied
#   universally (all phenotypes) and MUC1-VNTR assay applied to tubulointerstitial
#   probands only. Augmented Reflex applies the same upgraded ES detection at the
#   reflex step. Bundled assay costs are absorbed within the ES unit price.
#   Full PSA re-run required (detection changes affect per-proband diagnosis
#   probabilities; a post-hoc uplift is not valid).
# Author: Renal Genetics CEA Team

if (!require("dplyr")) install.packages("dplyr")
if (!require("tidyr")) install.packages("tidyr")
if (!require("readr")) install.packages("readr")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("ggrepel")) install.packages("ggrepel")
if (!requireNamespace("config", quietly = TRUE)) stop("Package 'config' is required.")

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(ggrepel)

# Source engine functions (uses sys.nframe() guard to avoid running base PSA)
source("scripts/04_cohort_generator.R")
source("scripts/05_simulation_engine.R")
source("scripts/utils_cea.R")

# ==============================================================================
# 1. Configuration
# ==============================================================================
cfg <- config::get()
PARAMS_DIR <- "data/params"
OUTPUT_DIR <- "outputs/results/supplement/es_augmented_scenario"
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

OUTPUT_ITER_CSV <- file.path(OUTPUT_DIR, "es_augmented_iteration_summary.csv")
OUTPUT_FRONTIER_PNG <- file.path(OUTPUT_DIR, "efficiency_frontier_es_augmented.png")
OUTPUT_PSA_CSV <- file.path(OUTPUT_DIR, "es_augmented_psa_robustness.csv")

N_ITER <- cfg$simulation$n_iter
N_PROBANDS <- cfg$simulation$n_probands
SEED_START <- cfg$simulation$seed_start
COST_CV <- cfg$simulation$cost_cv

STRATEGIES <- c("Panel", "ES_Augmented", "Reflex_Augmented")

STRATEGY_MAP_AUG <- c(
    "Panel"            = 1L,
    "ES_Augmented"     = 4L,
    "Reflex_Augmented" = 5L
)

STRATEGY_LABELS <- c(
    "Panel"            = "Panel",
    "ES_Augmented"     = "ES (Augmented)",
    "Reflex_Augmented" = "Reflex (Augmented)"
)

# ==============================================================================
# 2. Load Parameters
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
    panels = read_param("panels")
)

# Verify ES_Augmented is present in detection matrix
if (is.null(params$detection_betas$strategy_detection_matrix$ES_Augmented)) {
    stop("ES_Augmented not found in detection_betas. Re-run scripts/01d_load_detection.R first.")
}
cat("Parameters loaded. ES_Augmented detection matrix present.\n")

# ==============================================================================
# 3. PSA Loop
# ==============================================================================
cat(sprintf(
    "Starting Augmented ES scenario PSA: %d iterations, %d probands, %d strategies\n",
    N_ITER, N_PROBANDS, length(STRATEGIES)
))

PHENOTYPE_CATEGORIES <- names(params$phenotype_alphas)

iter_results_list <- vector("list", N_ITER * length(STRATEGIES))
list_idx <- 0

# PSA iteration seeds (same starting seed as base case for comparability)
set.seed(SEED_START)
iter_seeds <- sample.int(1e9, N_ITER)

for (i in 1:N_ITER) {
    if (i %% 100 == 0) cat(sprintf("  Iteration %d / %d\n", i, N_ITER))

    curr_seed <- iter_seeds[i]
    set.seed(curr_seed)

    # Sample parameters
    s_det <- sample_detection_matrix(params$detection_betas$strategy_detection_matrix)
    s_vus <- sample_vus_probs(params$vus_betas)

    curr_costs <- rapply(params$costs, function(x) {
        if (is.numeric(x) && length(x) == 1) sample_gamma_cost(x, cv = COST_CV) else x
    }, how = "replace")

    sampled_yields <- sapply(PHENOTYPE_CATEGORIES, function(p) {
        bp <- params$yield_betas[[p]]$params
        if (is.null(bp)) {
            return(0)
        }
        rbeta(1, shape1 = bp["shape1"], shape2 = bp["shape2"])
    })

    cohort_res <- generate_cohort(
        n_probands = N_PROBANDS,
        params = params,
        sampled_yields = sampled_yields,
        seed = curr_seed
    )
    cohort <- cohort_res$cohort

    curr_pheno_costs <- build_phenotype_cost_map(params$panels, curr_costs$tests$panel)

    for (strat in STRATEGIES) {
        list_idx <- list_idx + 1

        res <- run_simulation_step(
            cohort, strat, s_det, s_vus, curr_costs,
            params$panels, curr_pheno_costs,
            cascade_n_relatives = params$cascade$eligible_relatives_base
        )

        cost_per_diag <- if (res$diagnoses_per_proband > 0) {
            res$total_cost_per_proband_cad / res$diagnoses_per_proband
        } else {
            NA_real_
        }

        iter_results_list[[list_idx]] <- data.frame(
            iteration_id               = i,
            strategy_id                = STRATEGY_MAP_AUG[[strat]],
            strategy_label             = STRATEGY_LABELS[[strat]],
            total_cost_per_proband_cad = res$total_cost_per_proband_cad,
            diagnoses_per_proband      = res$diagnoses_per_proband,
            cost_per_diagnosis_cad     = cost_per_diag,
            pr_at_least_one_vus        = res$pr_at_least_one_vus,
            stringsAsFactors           = FALSE
        )
    }
}

iter_df <- bind_rows(iter_results_list)
cat("PSA complete.\n")

# ==============================================================================
# 4. Summary Statistics
# ==============================================================================
summary_df <- iter_df %>%
    group_by(strategy_id, strategy_label) %>%
    summarise(
        total_cost_per_proband_cad_mean = mean(total_cost_per_proband_cad),
        total_cost_per_proband_cad_ui_low = quantile(total_cost_per_proband_cad, 0.025),
        total_cost_per_proband_cad_ui_high = quantile(total_cost_per_proband_cad, 0.975),
        diagnoses_per_proband_mean = mean(diagnoses_per_proband),
        diagnoses_per_proband_ui_low = quantile(diagnoses_per_proband, 0.025),
        diagnoses_per_proband_ui_high = quantile(diagnoses_per_proband, 0.975),
        pr_at_least_one_vus_mean = mean(pr_at_least_one_vus),
        .groups = "drop"
    )

write_csv(iter_df, OUTPUT_ITER_CSV)
cat("Iteration summary written to:", OUTPUT_ITER_CSV, "\n")

# ==============================================================================
# 5. Incremental Analysis and Efficiency Frontier
# ==============================================================================
analysis_data <- summary_df %>%
    select(
        strategy_id,
        strategy_label,
        cost = total_cost_per_proband_cad_mean,
        effect = diagnoses_per_proband_mean,
        cost_ui_low = total_cost_per_proband_cad_ui_low,
        cost_ui_high = total_cost_per_proband_cad_ui_high,
        effect_ui_low = diagnoses_per_proband_ui_low,
        effect_ui_high = diagnoses_per_proband_ui_high
    )

results <- perform_incremental_analysis(analysis_data)

cat("\n--- Augmented ES Scenario: Dominance Analysis ---\n")
print(results %>% select(
    strategy_label, cost, effect, dominance_status,
    incremental_cost, incremental_effect, icer
))

# PSA robustness: proportion of iterations each strategy is non-dominated
psa_robustness <- iter_df %>%
    group_by(iteration_id) %>%
    arrange(cost = total_cost_per_proband_cad, .by_group = TRUE) %>%
    mutate(iter_data = list(cur_data())) %>%
    ungroup()

# Per-iteration frontier: compute non-dominated set for each iteration
compute_non_dominated_iter <- function(df) {
    # Sort by effect ascending; strategy is non-dominated if no other has
    # higher effect at lower cost (strict dominance) or lies on the frontier
    df <- df %>% arrange(diagnoses_per_proband)
    non_dom <- logical(nrow(df))
    min_cost_so_far <- Inf
    for (k in nrow(df):1) {
        if (df$total_cost_per_proband_cad[k] <= min_cost_so_far) {
            non_dom[k] <- TRUE
            min_cost_so_far <- df$total_cost_per_proband_cad[k]
        }
    }
    # Also check extended dominance via frontier
    frontier_idx <- which(non_dom)
    if (length(frontier_idx) >= 2) {
        # Check if any non-dominated strategy lies above the line
        frontier <- df[frontier_idx, ]
        for (ki in seq_along(frontier_idx)) {
            if (ki == 1 || ki == length(frontier_idx)) next
            slope_full <- (frontier$total_cost_per_proband_cad[length(frontier_idx)] -
                frontier$total_cost_per_proband_cad[1]) /
                (frontier$diagnoses_per_proband[length(frontier_idx)] -
                    frontier$diagnoses_per_proband[1])
            expected_cost <- frontier$total_cost_per_proband_cad[1] +
                slope_full * (frontier$diagnoses_per_proband[ki] -
                    frontier$diagnoses_per_proband[1])
            if (frontier$total_cost_per_proband_cad[ki] > expected_cost + 1e-6) {
                non_dom[frontier_idx[ki]] <- FALSE
            }
        }
    }
    df$non_dominated <- non_dom
    df
}

psa_nd <- iter_df %>%
    group_by(iteration_id) %>%
    group_modify(~ compute_non_dominated_iter(.x)) %>%
    ungroup()

robustness_tbl <- psa_nd %>%
    group_by(strategy_id, strategy_label) %>%
    summarise(pr_non_dominated = mean(non_dominated), .groups = "drop")

cat("\n--- PSA Robustness (Pr non-dominated) ---\n")
print(robustness_tbl)
cat(
    "Note: values sum to >100% because multiple strategies can be simultaneously",
    "non-dominated within the same iteration.\n"
)

write_csv(robustness_tbl, OUTPUT_PSA_CSV)

# ==============================================================================
# 6. Efficiency Frontier Figure
# ==============================================================================
plot_df <- results %>%
    left_join(summary_df %>% select(strategy_id, pr_at_least_one_vus_mean), by = "strategy_id") %>%
    mutate(
        is_dominated = dominance_status != "non_dominated",
        alpha_val = ifelse(is_dominated, 0.5, 1.0),
        vus_prob = pr_at_least_one_vus_mean,
        strategy_display = strategy_label,
        strategy_label_full = paste0(
            strategy_label, "\n(VUS: ", round(vus_prob * 100, 0), "%)"
        ),
        bubble_size = 3 + (vus_prob * 8)
    )

frontier_data <- plot_df %>%
    filter(dominance_status == "non_dominated") %>%
    arrange(effect)

p_frontier <- ggplot(plot_df, aes(x = effect, y = cost)) +
    geom_errorbar(
        aes(ymin = cost_ui_low, ymax = cost_ui_high, alpha = alpha_val),
        width = 0.005, color = "gray60", linewidth = 0.5
    ) +
    geom_errorbarh(
        aes(xmin = effect_ui_low, xmax = effect_ui_high, alpha = alpha_val),
        height = 50, color = "gray60", linewidth = 0.5
    ) +
    geom_line(
        data = frontier_data,
        aes(x = effect, y = cost),
        color = "black", linewidth = 1.0, linetype = "solid"
    ) +
    geom_point(
        aes(
            color = strategy_display,
            alpha = alpha_val,
            size  = bubble_size
        )
    ) +
    geom_label_repel(
        aes(label = strategy_label_full, alpha = alpha_val),
        size = 4.5, family = "sans", fontface = "bold",
        box.padding = 1.0, point.padding = 0.6,
        min.segment.length = 0.1, segment.size = 0.4,
        segment.color = "gray40",
        max.overlaps = 20,
        force = 3,
        force_pull = 0.4,
        fill = "white",
        label.size = 0.2,
        label.padding = unit(0.25, "lines")
    ) +
    scale_color_viridis_d(option = "D", begin = 0.15, end = 0.85, direction = 1) +
    scale_alpha_identity() +
    scale_size_identity() +
    scale_x_continuous(
        labels = scales::percent_format(accuracy = 1),
        expand = expansion(mult = 0.18)
    ) +
    scale_y_continuous(
        labels = scales::dollar_format(prefix = "$", suffix = ""),
        expand = expansion(mult = c(0.08, 0.12))
    ) +
    labs(
        title    = "Efficiency Frontier: Augmented ES Scenario",
        subtitle = "Bubble size represents VUS burden (larger = higher probability of VUS). Error bars show 95% UI.",
        x        = "Diagnostic Yield (Proportion of Probands Diagnosed)",
        y        = "Total Cost per Proband (CAD)",
        color    = NULL
    ) +
    theme_minimal(base_size = 13, base_family = "sans") +
    theme(
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.6),
        axis.ticks = element_line(color = "black", linewidth = 0.5),
        axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 13, face = "bold"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0),
        plot.subtitle = element_text(size = 11, color = "gray30", hjust = 0),
        legend.position = "none",
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        plot.margin = margin(15, 15, 10, 10)
    )

ggsave(OUTPUT_FRONTIER_PNG, p_frontier, width = 7, height = 5.5, dpi = 300, bg = "white")
cat("Efficiency frontier saved to:", OUTPUT_FRONTIER_PNG, "\n")

# ==============================================================================
# 7. Print key results for report
# ==============================================================================
cat("\n=== KEY RESULTS FOR §7.4 ===\n")
cat("\nMean outcomes:\n")
print(summary_df %>%
    mutate(
        yield_pct = sprintf("%.1f%%", diagnoses_per_proband_mean * 100),
        cost_fmt = sprintf("$%s", formatC(round(total_cost_per_proband_cad_mean), big.mark = ",")),
        yield_ui = sprintf("[%.1f%%–%.1f%%]", diagnoses_per_proband_ui_low * 100, diagnoses_per_proband_ui_high * 100),
        cost_ui = sprintf(
            "[$%s–$%s]",
            formatC(round(total_cost_per_proband_cad_ui_low), big.mark = ","),
            formatC(round(total_cost_per_proband_cad_ui_high), big.mark = ",")
        )
    ) %>%
    select(strategy_label, cost_fmt, cost_ui, yield_pct, yield_ui))

cat("\nIncremental analysis:\n")
print(results %>%
    select(strategy_label, dominance_status, incremental_cost, incremental_effect, icer) %>%
    mutate(icer_fmt = ifelse(is.na(icer), "—", sprintf("$%s", formatC(round(icer), big.mark = ",")))))

cat("\n=== SUCCESS ===\n")
cat("Outputs written to:", OUTPUT_DIR, "\n")
