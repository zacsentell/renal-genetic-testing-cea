# 09d_pkd1_full_detection_scenario.R
# Purpose: PKD1 full detection scenario analysis (§7.3)
#   Standard short-read WES detects only 7.14% of PKD1 variants in exons 1-32
#   due to competitive capture by six pseudogenes sharing 97.6-97.8% sequence
#   identity. Chang et al. (2022) showed that a modified xGEN probe design with
#   pseudogene-aware forced-alignment achieves complete coverage of all PKD1
#   coding regions at 35.8-68.6x mean depth, with 100% concordance to clinical
#   genetic testing. This scenario assigns full PKD1 variant detection (1.0) to
#   all Panel phenotypes and to ES, representing a near-future state in which
#   optimized capture and bioinformatic disambiguation make the entire PKD1 locus
#   accessible within standard exome workflows.
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
OUTPUT_DIR <- "outputs/results/scenario_analysis/pkd1_full_detection"
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

OUTPUT_ITER_CSV <- file.path(OUTPUT_DIR, "pkd1_full_iteration_outcomes.csv")
OUTPUT_SUMMARY_CSV <- file.path(OUTPUT_DIR, "pkd1_full_summary.csv")
OUTPUT_INCREMENTAL_CSV <- file.path(OUTPUT_DIR, "pkd1_full_incremental.csv")
OUTPUT_FRONTIER_PNG <- file.path(OUTPUT_DIR, "efficiency_frontier_pkd1_full.png")
OUTPUT_PSA_CSV <- file.path(OUTPUT_DIR, "pkd1_full_psa_robustness.csv")
OUTPUT_COMPARISON_CSV <- file.path(OUTPUT_DIR, "pkd1_full_comparison_vs_base.csv")

N_ITER <- cfg$simulation$n_iter
N_PROBANDS <- cfg$simulation$n_probands
SEED_START <- cfg$simulation$seed_start
COST_CV <- cfg$simulation$cost_cv

STRATEGIES <- c("Panel", "ES", "Panel_Reflex_ES")

STRATEGY_MAP <- c(
    "Panel"           = 1L,
    "ES"              = 2L,
    "Panel_Reflex_ES" = 3L
)

STRATEGY_LABELS <- c(
    "Panel"           = "Panel",
    "ES"              = "ES",
    "Panel_Reflex_ES" = "Reflex (Panel\u2192ES)"
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
    panels = read_param("panels"),
    cascade = read_param("cascade_params"),
    uptake = read_param("uptake_params")
)

# Verify base strategies present in detection matrix
stopifnot(
    !is.null(params$detection_betas$strategy_detection_matrix$Panel),
    !is.null(params$detection_betas$strategy_detection_matrix$ES)
)
cat("Parameters loaded. Panel and ES detection matrices present.\n")

# ==============================================================================
# 3. PKD1 Full Detection Override
# ==============================================================================
# Override sampled PKD1 detection probabilities to 1.0 for all Panel phenotypes
# and for ES (Universal). All other detection parameters (SNV/Indel, CNV/SV,
# MUC1) remain at their sampled values. This is a deterministic scenario
# assumption, not a sampled parameter.

override_pkd1_full <- function(s_det) {
    # Panel: phenotype-keyed, override each phenotype's PKD1
    for (pheno in names(s_det$Panel)) {
        if (!is.null(s_det$Panel[[pheno]]$Difficult_Locus$PKD1)) {
            s_det$Panel[[pheno]]$Difficult_Locus$PKD1 <- 1.0
        }
    }
    # ES: universal, override PKD1
    if (!is.null(s_det$ES$Universal$Difficult_Locus$PKD1)) {
        s_det$ES$Universal$Difficult_Locus$PKD1 <- 1.0
    }
    s_det
}

# ==============================================================================
# 4. PSA Loop
# ==============================================================================
cat(sprintf(
    "Starting PKD1 full detection scenario PSA: %d iterations, %d probands, %d strategies\n",
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

    # Apply PKD1 full detection override
    s_det <- override_pkd1_full(s_det)

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

    for (strat in STRATEGIES) {
        list_idx <- list_idx + 1

        res <- run_simulation_step(
            cohort, strat, s_det, s_vus, curr_costs,
            params$panels, curr_pheno_costs,
            cascade_n_relatives = params$cascade$eligible_relatives_base,
            uptake_reflex = uptake_reflex_i,
            uptake_cascade = uptake_cascade_i
        )

        cost_per_diag <- if (res$diagnoses_per_proband > 0) {
            res$total_cost_per_proband_cad / res$diagnoses_per_proband
        } else {
            NA_real_
        }

        iter_results_list[[list_idx]] <- data.frame(
            iteration_id               = i,
            strategy_id                = STRATEGY_MAP[[strat]],
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

write_csv(iter_df, OUTPUT_ITER_CSV)
cat("Iteration outcomes written to:", OUTPUT_ITER_CSV, "\n")

# ==============================================================================
# 5. Summary Statistics
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
        cost_per_diagnosis_cad_mean = mean(cost_per_diagnosis_cad, na.rm = TRUE),
        pr_at_least_one_vus_mean = mean(pr_at_least_one_vus),
        .groups = "drop"
    )

write_csv(summary_df, OUTPUT_SUMMARY_CSV)
cat("Summary written to:", OUTPUT_SUMMARY_CSV, "\n")

# ==============================================================================
# 6. Incremental Analysis and Efficiency Frontier
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

write_csv(results, OUTPUT_INCREMENTAL_CSV)
cat("Incremental analysis written to:", OUTPUT_INCREMENTAL_CSV, "\n")

cat("\n--- PKD1 Full Detection Scenario: Dominance Analysis ---\n")
print(results %>% select(
    strategy_label, cost, effect, dominance_status,
    incremental_cost, incremental_effect, icer
))

# ==============================================================================
# 7. PSA Robustness
# ==============================================================================
# Per-iteration frontier: compute non-dominated set for each iteration
compute_non_dominated_iter <- function(df) {
    df <- df %>% arrange(diagnoses_per_proband)
    non_dom <- logical(nrow(df))
    min_cost_so_far <- Inf
    for (k in nrow(df):1) {
        if (df$total_cost_per_proband_cad[k] <= min_cost_so_far) {
            non_dom[k] <- TRUE
            min_cost_so_far <- df$total_cost_per_proband_cad[k]
        }
    }
    # Check extended dominance via frontier
    frontier_idx <- which(non_dom)
    if (length(frontier_idx) >= 2) {
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
cat("PSA robustness written to:", OUTPUT_PSA_CSV, "\n")

# ==============================================================================
# 8. Comparison with Base Case
# ==============================================================================
BASE_CASE_CSV <- "outputs/results/base_case/iteration_level/strategy_iteration_outcomes.csv"
if (file.exists(BASE_CASE_CSV)) {
    base_iter <- read_csv(BASE_CASE_CSV, show_col_types = FALSE)

    base_summary <- base_iter %>%
        group_by(strategy_id) %>%
        summarise(
            base_cost = mean(total_cost_per_proband_cad),
            base_yield = mean(diagnoses_per_proband),
            .groups = "drop"
        )

    scenario_summary <- summary_df %>%
        select(
            strategy_id, strategy_label,
            scenario_cost = total_cost_per_proband_cad_mean,
            scenario_yield = diagnoses_per_proband_mean
        )

    # Join on strategy_id only (base case uses raw names, scenario uses display labels)
    comparison_df <- scenario_summary %>%
        inner_join(base_summary, by = "strategy_id") %>%
        mutate(
            delta_cost = scenario_cost - base_cost,
            delta_yield = scenario_yield - base_yield
        )

    write_csv(comparison_df, OUTPUT_COMPARISON_CSV)
    cat("Comparison with base case written to:", OUTPUT_COMPARISON_CSV, "\n")

    cat("\n--- Base Case vs PKD1 Full Detection ---\n")
    print(comparison_df %>%
        mutate(
            base_yield_pct = sprintf("%.2f%%", base_yield * 100),
            scenario_yield_pct = sprintf("%.2f%%", scenario_yield * 100),
            delta_yield_pct = sprintf("%+.2f pp", delta_yield * 100),
            delta_cost_fmt = sprintf("%+.0f CAD", delta_cost)
        ) %>%
        select(strategy_label, base_yield_pct, scenario_yield_pct, delta_yield_pct, delta_cost_fmt))
} else {
    cat("WARNING: Base case iteration outcomes not found at:", BASE_CASE_CSV, "\n")
    cat("Skipping base case comparison.\n")
}

# ==============================================================================
# 9. ES PKD1 Detection Threshold Analysis
# ==============================================================================
# Sweep ES PKD1 detection from the base rate to 1.0 to find the minimum
# detection rate at which ES becomes non-dominated (joins the frontier).
#
# Assumptions (deterministic, same approximation as GS uplift threshold):
#   - Panel yield is held at the base case mean. Panel already uses its bundled
#     long-range PCR assay (~99% PKD1 detection); improving ES detection does not
#     affect the Panel step.
#   - Reflex yield is held at the base case mean. Most cystic PKD1 probands are
#     diagnosed at the Panel step and never reach ES, so the second-order effect
#     of improved ES PKD1 detection on Reflex yield is negligible.
#   - ES yield scales linearly between the base case mean (at BASE_PKD1_DETECT)
#     and the full-detection scenario mean (at 1.0).
#   - ES cost is held at the base case mean (assay cost is not separately
#     modelled in this threshold sweep).

cat("\n--- ES PKD1 Detection Threshold Analysis ---\n")

OUTPUT_THRESHOLD_PKD1 <- file.path(OUTPUT_DIR, "pkd1_es_detection_threshold.csv")

BASE_CASE_SUMMARY_CSV <- "outputs/results/base_case/summary_tables/base_case_outcomes_by_strategy.csv"
BASE_PKD1_DETECT <- 0.125

if (file.exists(BASE_CASE_SUMMARY_CSV)) {
    base_summary_tbl <- read_csv(BASE_CASE_SUMMARY_CSV, show_col_types = FALSE)

    panel_base  <- base_summary_tbl %>% filter(strategy_id == 1)
    es_base     <- base_summary_tbl %>% filter(strategy_id == 2)
    reflex_base <- base_summary_tbl %>% filter(strategy_id == 3)

    panel_mean_cost   <- panel_base$total_cost_per_proband_cad_mean
    panel_mean_yield  <- panel_base$diagnoses_per_proband_mean
    es_mean_cost      <- es_base$total_cost_per_proband_cad_mean
    base_es_yield     <- es_base$diagnoses_per_proband_mean
    reflex_mean_cost  <- reflex_base$total_cost_per_proband_cad_mean
    reflex_mean_yield <- reflex_base$diagnoses_per_proband_mean

    full_es_yield <- summary_df %>%
        filter(strategy_id == 2) %>%
        pull(diagnoses_per_proband_mean)

    detect_seq <- seq(BASE_PKD1_DETECT, 1.0, by = 0.005)

    threshold_results_pkd1 <- lapply(detect_seq, function(d) {
        es_yield_d <- base_es_yield +
            (d - BASE_PKD1_DETECT) * (full_es_yield - base_es_yield) / (1.0 - BASE_PKD1_DETECT)

        frontier_input <- data.frame(
            strategy_label = c("Panel", "ES", "Reflex"),
            cost   = c(panel_mean_cost, es_mean_cost, reflex_mean_cost),
            effect = c(panel_mean_yield, es_yield_d, reflex_mean_yield),
            stringsAsFactors = FALSE
        )

        dom <- perform_incremental_analysis(frontier_input)
        es_row <- dom[dom$strategy_label == "ES", , drop = FALSE]

        data.frame(
            pkd1_detect_rate = d,
            es_yield         = es_yield_d,
            panel_yield      = panel_mean_yield,
            reflex_yield     = reflex_mean_yield,
            es_on_frontier   = es_row$dominance_status == "non_dominated",
            icpd_es_vs_panel = ifelse(
                es_row$dominance_status == "non_dominated" & !is.na(es_row$icer),
                es_row$icer, NA_real_
            ),
            stringsAsFactors = FALSE
        )
    })

    threshold_pkd1_df <- do.call(rbind, threshold_results_pkd1)
    write_csv(threshold_pkd1_df, OUTPUT_THRESHOLD_PKD1)
    cat("PKD1 ES detection threshold table written to:", OUTPUT_THRESHOLD_PKD1, "\n")

    min_frontier_pkd1 <- threshold_pkd1_df[threshold_pkd1_df$es_on_frontier, , drop = FALSE]
    if (nrow(min_frontier_pkd1) > 0) {
        min_row_pkd1 <- min_frontier_pkd1[1, , drop = FALSE]
        cat(sprintf(
            "Minimum ES PKD1 detection rate for ES on frontier: %.1f%% (ES yield %.1f%% vs Reflex %.1f%%; ICPD vs Panel $%s)\n",
            min_row_pkd1$pkd1_detect_rate * 100,
            min_row_pkd1$es_yield * 100,
            min_row_pkd1$reflex_yield * 100,
            formatC(round(min_row_pkd1$icpd_es_vs_panel), format = "f", digits = 0, big.mark = ",")
        ))
    } else {
        cat("ES does not enter the frontier at any tested detection rate (12.5%-100%).\n")
    }
} else {
    cat("WARNING: Base case summary not found at:", BASE_CASE_SUMMARY_CSV, "\n")
    cat("Skipping ES PKD1 detection threshold analysis.\n")
}

# ==============================================================================
# 10. Efficiency Frontier Figure
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
        title    = "Efficiency Frontier: PKD1 Full Detection Scenario",
        subtitle = "Bubble size represents VUS burden (larger = higher probability of VUS).\nError bars show 95% UI.",
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
# 11. Print key results for report
# ==============================================================================
cat("\n=== KEY RESULTS FOR \u00a77.3 ===\n")
cat("\nMean outcomes:\n")
print(summary_df %>%
    mutate(
        yield_pct = sprintf("%.1f%%", diagnoses_per_proband_mean * 100),
        cost_fmt = sprintf("$%s", formatC(round(total_cost_per_proband_cad_mean), big.mark = ",")),
        yield_ui = sprintf("[%.1f%%\u2013%.1f%%]", diagnoses_per_proband_ui_low * 100, diagnoses_per_proband_ui_high * 100),
        cost_ui = sprintf(
            "[$%s\u2013$%s]",
            formatC(round(total_cost_per_proband_cad_ui_low), big.mark = ","),
            formatC(round(total_cost_per_proband_cad_ui_high), big.mark = ",")
        )
    ) %>%
    select(strategy_label, cost_fmt, cost_ui, yield_pct, yield_ui))

cat("\nIncremental analysis:\n")
print(results %>%
    select(strategy_label, dominance_status, incremental_cost, incremental_effect, icer) %>%
    mutate(icer_fmt = ifelse(is.na(icer), "\u2014", sprintf("$%s", formatC(round(icer), big.mark = ",")))))

cat("\n=== SUCCESS ===\n")
cat("Outputs written to:", OUTPUT_DIR, "\n")
