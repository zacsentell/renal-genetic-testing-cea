# 09e_reflex_uptake_dsa.R
# Purpose: One-way deterministic sensitivity analysis (DSA) for the reflex
#   testing uptake probability (base case = 50%, range 20%-100%).
#
# Method: Full PSA re-simulation at each DSA uptake value. Unlike the cascade
#   family-size DSA (09c), analytical rescaling is not valid here because
#   changing uptake alters which probands receive ES, which changes diagnoses,
#   VUS outcomes, and all downstream costs non-linearly.
#
# Outputs:
#   outputs/results/uncertainty_sensitivity/reflex_uptake_dsa/reflex_uptake_dsa_icpd.csv
#   outputs/results/uncertainty_sensitivity/reflex_uptake_dsa/reflex_uptake_dsa_costs.csv
#   outputs/results/uncertainty_sensitivity/reflex_uptake_dsa/reflex_uptake_dsa_icpd_plot.png/.pdf

if (!require("dplyr")) install.packages("dplyr")
if (!require("tidyr")) install.packages("tidyr")
if (!require("readr")) install.packages("readr")
if (!require("ggplot2")) install.packages("ggplot2")
if (!requireNamespace("config", quietly = TRUE)) {
    stop("Package 'config' is required. Install it with install.packages('config').")
}

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

# Source engine functions (uses sys.nframe() guard to avoid running base PSA)
source("scripts/04_cohort_generator.R")
source("scripts/05_simulation_engine.R")
source("scripts/utils_schema.R")

# ==============================================================================
# 1. Configuration
# ==============================================================================
cfg <- config::get()
PARAMS_DIR <- "data/params"

N_ITER     <- cfg$simulation$n_iter
N_PROBANDS <- cfg$simulation$n_probands
SEED_START <- cfg$simulation$seed_start
COST_CV    <- cfg$simulation$cost_cv

# Load uptake parameters for DSA range
uptake_params <- readRDS(file.path(PARAMS_DIR, "uptake_params.rds"))

if (is.na(uptake_params$reflex$dsa_min) || is.na(uptake_params$reflex$dsa_max)) {
    stop("Reflex uptake DSA range (dsa_min, dsa_max) must be specified in uptake_parameters.csv")
}

base_uptake <- uptake_params$reflex$base_case
dsa_min     <- uptake_params$reflex$dsa_min
dsa_max     <- uptake_params$reflex$dsa_max
uptake_values <- seq(dsa_min, dsa_max, by = 0.10)

# Ensure base case is included in the sweep
if (!base_uptake %in% uptake_values) {
    uptake_values <- sort(unique(c(uptake_values, base_uptake)))
}

# Strategies: Panel and Reflex are needed for the ICPD comparison
STRATEGIES <- c("Panel", "ES", "Panel_Reflex_ES")
STRATEGY_MAP <- c("Panel" = 1L, "ES" = 2L, "Panel_Reflex_ES" = 3L)

OUTPUT_DIR <- "outputs/results/uncertainty_sensitivity/reflex_uptake_dsa"
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

OUTPUT_ICPD_CSV  <- file.path(OUTPUT_DIR, "reflex_uptake_dsa_icpd.csv")
OUTPUT_COSTS_CSV <- file.path(OUTPUT_DIR, "reflex_uptake_dsa_costs.csv")
OUTPUT_FIG_BASE  <- file.path(OUTPUT_DIR, "reflex_uptake_dsa_icpd_plot")

cat("--- Reflex Uptake DSA ---\n")
cat(sprintf("Base case: %.0f%%. Sweep: %.0f%% to %.0f%% (%d values)\n",
    base_uptake * 100, dsa_min * 100, dsa_max * 100, length(uptake_values)))

# ==============================================================================
# 2. Load Parameters
# ==============================================================================
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

PHENOTYPE_CATEGORIES <- names(params$phenotype_alphas)

# Derive iter_seeds identically to 05_simulation_engine.R
set.seed(SEED_START)
iter_seeds <- sample.int(1e9, N_ITER)

# ==============================================================================
# 3. DSA Loop: Re-simulate at Each Uptake Value
# ==============================================================================
all_results <- vector("list", length(uptake_values))

for (u_idx in seq_along(uptake_values)) {
    dsa_uptake <- uptake_values[u_idx]
    cat(sprintf("\n  DSA uptake = %.0f%% (%d/%d)\n",
        dsa_uptake * 100, u_idx, length(uptake_values)))

    iter_results_list <- vector("list", N_ITER * length(STRATEGIES))
    list_idx <- 0

    for (i in 1:N_ITER) {
        if (i %% 200 == 0) cat(sprintf("    Iteration %d / %d\n", i, N_ITER))

        curr_seed <- iter_seeds[i]
        set.seed(curr_seed)

        # Sample parameters in same order as 05_simulation_engine.R
        s_det <- sample_detection_matrix(params$detection_betas$strategy_detection_matrix)
        s_vus <- sample_vus_probs(params$vus_betas)

        curr_costs <- rapply(params$costs, function(x) {
            if (is.numeric(x) && length(x) == 1) sample_gamma_cost(x, cv = COST_CV) else x
        }, how = "replace")

        sampled_yields <- sapply(PHENOTYPE_CATEGORIES, function(p) {
            bp <- params$yield_betas[[p]]$params
            if (is.null(bp)) return(0)
            rbeta(1, shape1 = bp["shape1"], shape2 = bp["shape2"])
        })

        # Consume RNG for uptake sampling to maintain seed alignment,
        # but override reflex uptake with the DSA value
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

        # Override reflex uptake with DSA value (cascade remains as sampled)
        uptake_reflex_i <- dsa_uptake

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

            iter_results_list[[list_idx]] <- data.frame(
                uptake_reflex = dsa_uptake,
                iteration_id = i,
                strategy_id = STRATEGY_MAP[[strat]],
                strategy_label = strat,
                total_cost_per_proband_cad = res$total_cost_per_proband_cad,
                diagnoses_per_proband = res$diagnoses_per_proband,
                stringsAsFactors = FALSE
            )
        }
    }

    all_results[[u_idx]] <- bind_rows(iter_results_list)
}

dsa_data <- bind_rows(all_results)

cat(sprintf("\nDSA complete: %d total rows (%d uptake values x %d iterations x %d strategies)\n",
    nrow(dsa_data), length(uptake_values), N_ITER, length(STRATEGIES)))

# ==============================================================================
# 4. Compute ICPD: Reflex vs Panel at Each Uptake Level
# ==============================================================================
focal <- dsa_data %>%
    filter(strategy_label == "Panel_Reflex_ES") %>%
    select(uptake_reflex, iteration_id,
           cost_focal = total_cost_per_proband_cad,
           yield_focal = diagnoses_per_proband)

comp <- dsa_data %>%
    filter(strategy_label == "Panel") %>%
    select(uptake_reflex, iteration_id,
           cost_comp = total_cost_per_proband_cad,
           yield_comp = diagnoses_per_proband)

icpd_all <- focal %>%
    inner_join(comp, by = c("uptake_reflex", "iteration_id")) %>%
    mutate(
        delta_cost  = cost_focal - cost_comp,
        delta_yield = yield_focal - yield_comp,
        icpd = ifelse(delta_yield > 1e-9, delta_cost / delta_yield, NA_real_),
        comparison = "Panel_Reflex_ES vs Panel"
    )

n_invalid <- sum(is.na(icpd_all$icpd))
if (n_invalid > 0) {
    warning(sprintf(
        "%d ICPD values are NA (zero or negative delta_yield). Excluded from summaries.",
        n_invalid
    ))
}

# Summary: mean + 95% UI per uptake value
icpd_summary <- icpd_all %>%
    group_by(comparison, uptake_reflex) %>%
    summarise(
        icpd_mean    = mean(icpd, na.rm = TRUE),
        icpd_q025    = quantile(icpd, 0.025, na.rm = TRUE),
        icpd_q975    = quantile(icpd, 0.975, na.rm = TRUE),
        n_iterations = sum(!is.na(icpd)),
        .groups = "drop"
    ) %>%
    mutate(is_base_case = (uptake_reflex == base_uptake))

# ==============================================================================
# 5. Total Cost Summary by Strategy x Uptake
# ==============================================================================
cost_summary <- dsa_data %>%
    group_by(strategy_label, uptake_reflex) %>%
    summarise(
        total_cost_mean = mean(total_cost_per_proband_cad, na.rm = TRUE),
        total_cost_q025 = quantile(total_cost_per_proband_cad, 0.025, na.rm = TRUE),
        total_cost_q975 = quantile(total_cost_per_proband_cad, 0.975, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    mutate(is_base_case = (uptake_reflex == base_uptake))

# ==============================================================================
# 6. Write Outputs
# ==============================================================================
write_csv_validated(icpd_summary, OUTPUT_ICPD_CSV, "reflex_uptake_dsa_icpd")
write_csv_validated(cost_summary, OUTPUT_COSTS_CSV, "reflex_uptake_dsa_costs")

# ==============================================================================
# 7. Figure: ICPD vs Reflex Uptake
# ==============================================================================
plot_data <- icpd_summary

p <- ggplot(plot_data,
            aes(x = uptake_reflex * 100, y = icpd_mean)) +
    geom_ribbon(aes(ymin = icpd_q025, ymax = icpd_q975),
                alpha = 0.15, fill = "#440154") +
    geom_line(linewidth = 0.9, colour = "#440154") +
    geom_point(size = 2.8, colour = "#440154") +
    geom_vline(xintercept = base_uptake * 100, linetype = "dashed",
               colour = "grey45", linewidth = 0.6) +
    annotate("text",
             x     = base_uptake * 100 + 2,
             y     = max(plot_data$icpd_q975, na.rm = TRUE),
             label = paste0("Base case\n(", base_uptake * 100, "%)"),
             hjust = 0, vjust = 1.1,
             size  = 3, colour = "grey40") +
    scale_x_continuous(
        breaks = uptake_values * 100,
        labels = paste0(uptake_values * 100, "%")
    ) +
    scale_y_continuous(
        labels = scales::dollar_format(prefix = "$", suffix = "")
    ) +
    labs(
        title    = "One-Way DSA: Reflex Testing Uptake Probability",
        subtitle = paste0(
            "ICPD for Reflex (Panel\u2192ES) vs Panel by probability that\n",
            "panel-negative probands proceed to ES (base case = ",
            base_uptake * 100, "%)"
        ),
        x        = "Reflex ES uptake probability (%)",
        y        = "ICPD (CAD per additional diagnosis)",
        caption  = "Shaded band: 95% uncertainty interval (PSA, 1000 iterations). Dashed line: base case."
    ) +
    theme_bw(base_size = 11) +
    theme(
        plot.caption     = element_text(size = 8, colour = "grey50"),
        plot.subtitle    = element_text(size = 9)
    )

ggsave(paste0(OUTPUT_FIG_BASE, ".png"), p, dpi = 300, width = 7, height = 5.5)
ggsave(paste0(OUTPUT_FIG_BASE, ".pdf"), p, width = 7, height = 5.5)

# ==============================================================================
# 8. Console Summary
# ==============================================================================
cat("\nICPD summary at each reflex uptake value:\n")
print(
    icpd_summary %>%
        mutate(
            uptake_pct = paste0(uptake_reflex * 100, "%"),
            icpd_fmt = sprintf("$%s ($%s\u2013$%s)",
                formatC(icpd_mean,  format = "f", digits = 0, big.mark = ","),
                formatC(icpd_q025,  format = "f", digits = 0, big.mark = ","),
                formatC(icpd_q975,  format = "f", digits = 0, big.mark = ","))
        ) %>%
        select(comparison, uptake_pct, icpd_fmt),
    n = Inf
)

cat("\nReflex uptake DSA complete. Outputs written to:", OUTPUT_DIR, "\n")
cat("  -", OUTPUT_ICPD_CSV, "\n")
cat("  -", OUTPUT_COSTS_CSV, "\n")
cat("  -", paste0(OUTPUT_FIG_BASE, ".png"), "\n")
