# 09c_cascade_dsa.R
# Purpose: One-way deterministic sensitivity analysis (DSA) for the cascade
#   testing family-size assumption (eligible first-degree relatives per
#   diagnosed proband, base case = 2, range = 0–4).
#
# Method: Analytical derivation from existing PSA iteration-level cost
#   components. Because cascade cost = n_relatives * unit_cost * yield,
#   we can scale cost_cascade_cad to any n_relatives without re-running
#   the simulation:
#
#     total_cost_at_k = (total_cost_base - cost_cascade_base)
#                       + (k / base_n) * cost_cascade_base
#
#   This preserves full PSA uncertainty in unit costs (Gamma-sampled per
#   iteration) while sweeping the structural count assumption separately.
#
# Outputs (referenced in report §7.4):
#   outputs/results/uncertainty_sensitivity/cascade_dsa/cascade_dsa_icpd.csv
#   outputs/results/uncertainty_sensitivity/cascade_dsa/cascade_dsa_costs.csv
#   outputs/results/uncertainty_sensitivity/cascade_dsa/cascade_dsa_icpd_plot.png/.pdf

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
if (!requireNamespace("config", quietly = TRUE)) {
    stop("Package 'config' is required. Install it with install.packages('config').")
}
source("scripts/utils_schema.R")

# ==============================================================================
# 1. Configuration
# ==============================================================================
cascade_params <- readRDS("data/params/cascade_params.rds")
base_n <- cascade_params$eligible_relatives_base
min_n  <- cascade_params$eligible_relatives_min
max_n  <- cascade_params$eligible_relatives_max

if (base_n <= 0) stop("cascade$eligible_relatives_base must be > 0 for DSA scaling.")

relatives_values <- seq(min_n, max_n)

INPUT_COMPONENTS <- "outputs/results/base_case/cost_composition/iteration_cost_components.csv"
INPUT_OUTCOMES   <- "outputs/results/base_case/iteration_level/strategy_iteration_outcomes.csv"

OUTPUT_DIR <- "outputs/results/uncertainty_sensitivity/cascade_dsa"
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

OUTPUT_ICPD_CSV  <- file.path(OUTPUT_DIR, "cascade_dsa_icpd.csv")
OUTPUT_COSTS_CSV <- file.path(OUTPUT_DIR, "cascade_dsa_costs.csv")
OUTPUT_FIG_BASE  <- file.path(OUTPUT_DIR, "cascade_dsa_icpd_plot")

cat("--- Cascade DSA ---\n")
cat(sprintf("Base case: %g relatives. Sweep: %g to %g\n", base_n, min_n, max_n))

# ==============================================================================
# 2. Load and Validate Inputs
# ==============================================================================
if (!file.exists(INPUT_COMPONENTS)) stop("Cost components file not found: ", INPUT_COMPONENTS)
if (!file.exists(INPUT_OUTCOMES))   stop("Iteration outcomes file not found: ", INPUT_OUTCOMES)

comp_data <- read_csv(INPUT_COMPONENTS, show_col_types = FALSE)
iter_data <- read_csv(INPUT_OUTCOMES,   show_col_types = FALSE)

required_comp <- c(
    "iteration_id", "strategy_id", "strategy_label",
    "cost_cascade_cad", "total_cost_per_proband_cad"
)
required_iter <- c(
    "iteration_id", "strategy_id", "strategy_label",
    "diagnoses_per_proband"
)
assert_required_columns(comp_data, required_comp, "iteration_cost_components")
assert_required_columns(iter_data, required_iter, "strategy_iteration_outcomes")

# Merge cost components with yield
merged <- comp_data %>%
    select(iteration_id, strategy_id, strategy_label,
           cost_cascade_cad, total_cost_per_proband_cad) %>%
    left_join(
        iter_data %>% select(iteration_id, strategy_id, diagnoses_per_proband),
        by = c("iteration_id", "strategy_id")
    )

n_na_yield <- sum(is.na(merged$diagnoses_per_proband))
if (n_na_yield > 0) {
    stop(sprintf(
        "%d rows have NA diagnoses_per_proband after join. Check strategy_id alignment.",
        n_na_yield
    ))
}

# ==============================================================================
# 3. Analytical DSA Sweep
# ==============================================================================
# For each n_relatives = k:
#   total_cost_k = (total_cost_base - cost_cascade_base)
#                 + (k / base_n) * cost_cascade_base
#
# Derivation: cost_cascade_base = base_n * unit_cost * yield_indicator (per proband)
#             cost_cascade_k    = k      * unit_cost * yield_indicator
#          => cost_cascade_k    = (k / base_n) * cost_cascade_base  [exact]

dsa_iter <- merged %>%
    crossing(n_relatives = relatives_values) %>%
    mutate(
        total_cost_dsa = (total_cost_per_proband_cad - cost_cascade_cad) +
                         (n_relatives / base_n) * cost_cascade_cad,
        is_base_case = (n_relatives == base_n)
    )

# ==============================================================================
# 4. ICPD by Comparison × n_relatives × Iteration
# ==============================================================================
compute_icpd_series <- function(data, focal_label, comp_label) {
    focal <- data %>%
        filter(strategy_label == focal_label) %>%
        select(iteration_id, n_relatives,
               cost_focal = total_cost_dsa,
               yield_focal = diagnoses_per_proband)

    comp <- data %>%
        filter(strategy_label == comp_label) %>%
        select(iteration_id, n_relatives,
               cost_comp = total_cost_dsa,
               yield_comp = diagnoses_per_proband)

    focal %>%
        inner_join(comp, by = c("iteration_id", "n_relatives")) %>%
        mutate(
            delta_cost  = cost_focal - cost_comp,
            delta_yield = yield_focal - yield_comp,
            icpd        = ifelse(delta_yield > 1e-9, delta_cost / delta_yield, NA_real_),
            comparison  = sprintf("%s vs %s", focal_label, comp_label)
        )
}

icpd_es_panel  <- compute_icpd_series(dsa_iter, "ES",             "Panel")
icpd_reflex_es <- compute_icpd_series(dsa_iter, "Panel_Reflex_ES", "ES")

icpd_all <- bind_rows(icpd_es_panel, icpd_reflex_es)

n_invalid <- sum(is.na(icpd_all$icpd))
if (n_invalid > 0) {
    warning(sprintf(
        "%d ICPD values are NA (zero or negative delta_yield). These are excluded from summaries.",
        n_invalid
    ))
}

# Summary: mean + 95% UI per comparison × n_relatives
icpd_summary <- icpd_all %>%
    group_by(comparison, n_relatives) %>%
    summarise(
        icpd_mean    = mean(icpd, na.rm = TRUE),
        icpd_q025    = quantile(icpd, 0.025, na.rm = TRUE),
        icpd_q975    = quantile(icpd, 0.975, na.rm = TRUE),
        n_iterations = sum(!is.na(icpd)),
        .groups = "drop"
    ) %>%
    mutate(is_base_case = (n_relatives == base_n))

# ==============================================================================
# 5. Total Cost Summary by Strategy × n_relatives
# ==============================================================================
cost_summary <- dsa_iter %>%
    group_by(strategy_label, n_relatives) %>%
    summarise(
        total_cost_mean = mean(total_cost_dsa, na.rm = TRUE),
        total_cost_q025 = quantile(total_cost_dsa, 0.025, na.rm = TRUE),
        total_cost_q975 = quantile(total_cost_dsa, 0.975, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    mutate(is_base_case = (n_relatives == base_n))

# ==============================================================================
# 6. Write Outputs
# ==============================================================================
write_csv_validated(icpd_summary,  OUTPUT_ICPD_CSV,  "cascade_dsa_icpd")
write_csv_validated(cost_summary,  OUTPUT_COSTS_CSV, "cascade_dsa_costs")

# ==============================================================================
# 7. Figure: ICPD vs n_relatives
# ==============================================================================
comparison_labels <- c(
    "ES vs Panel"            = "ES vs Panel",
    "Panel_Reflex_ES vs ES"  = "Reflex (Panel\u2192ES) vs ES"
)

plot_data <- icpd_summary %>%
    mutate(
        comparison_label = recode(comparison, !!!comparison_labels),
        comparison_label = factor(comparison_label, levels = unname(comparison_labels))
    )

p <- ggplot(plot_data,
            aes(x = n_relatives, y = icpd_mean,
                colour = comparison_label, fill = comparison_label)) +
    geom_ribbon(aes(ymin = icpd_q025, ymax = icpd_q975),
                alpha = 0.15, colour = NA) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.8) +
    geom_vline(xintercept = base_n, linetype = "dashed",
               colour = "grey45", linewidth = 0.6) +
    annotate("text",
             x     = base_n + 0.12,
             y     = max(plot_data$icpd_q975, na.rm = TRUE),
             label = paste0("Base case\n(n = ", base_n, ")"),
             hjust = 0, vjust = 1.1,
             size  = 3, colour = "grey40") +
    scale_x_continuous(
        breaks = relatives_values,
        labels = relatives_values
    ) +
    scale_y_continuous(
        labels = scales::dollar_format(prefix = "$", suffix = "")
    ) +
    scale_colour_viridis_d(option = "D", begin = 0.15, end = 0.75,
                           name = "Comparison") +
    scale_fill_viridis_d(option = "D",   begin = 0.15, end = 0.75,
                         name = "Comparison") +
    labs(
        title    = "One-Way DSA: Cascade Testing Family Size",
        subtitle = paste0(
            "Incremental cost per additional diagnosis (ICPD) by number of eligible\n",
            "first-degree relatives tested per diagnosed proband (base case = ", base_n, ")"
        ),
        x        = "Eligible first-degree relatives per diagnosed proband (n)",
        y        = "ICPD (CAD per additional diagnosis)",
        caption  = "Shaded bands: 95% uncertainty interval (PSA, 1000 iterations). Dashed line: base case."
    ) +
    theme_bw(base_size = 11) +
    theme(
        legend.position  = "bottom",
        plot.caption     = element_text(size = 8, colour = "grey50"),
        plot.subtitle    = element_text(size = 9)
    )

ggsave(paste0(OUTPUT_FIG_BASE, ".png"), p, dpi = 300, width = 7, height = 5.5)
ggsave(paste0(OUTPUT_FIG_BASE, ".pdf"), p, width = 7, height = 5.5)

# ==============================================================================
# 8. Console Summary
# ==============================================================================
cat("\nICPD summary at each n_relatives value:\n")
print(
    icpd_summary %>%
        mutate(
            icpd_fmt = sprintf("$%s ($%s-$%s)",
                formatC(icpd_mean,  format = "f", digits = 0, big.mark = ","),
                formatC(icpd_q025,  format = "f", digits = 0, big.mark = ","),
                formatC(icpd_q975,  format = "f", digits = 0, big.mark = ","))
        ) %>%
        select(comparison, n_relatives, icpd_fmt),
    n = Inf
)

cat("\nCascade DSA complete. Outputs written to:", OUTPUT_DIR, "\n")
cat("  -", OUTPUT_ICPD_CSV, "\n")
cat("  -", OUTPUT_COSTS_CSV, "\n")
cat("  -", paste0(OUTPUT_FIG_BASE, ".png"), "\n")
