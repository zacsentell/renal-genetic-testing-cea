# 06b_incremental_analysis.R
# Purpose: Incremental analysis and efficiency frontier for base case (§7.5)
# Author: Renal Genetics CEA Team
# Date: 2025-12-30

library(dplyr)
library(readr)
library(ggplot2)
library(ggrepel)
source("scripts/utils_schema.R")

# ==============================================================================
# 1. Configuration
# ==============================================================================
INPUT_FILE <- "outputs/results/base_case/summary_tables/base_case_outcomes_by_strategy.csv"
OUTPUT_DIR <- "outputs/results/incremental_analysis/base_case"
OUTPUT_TABLE <- file.path(OUTPUT_DIR, "incremental_table_base_case.csv")
FIG_BASE <- file.path(OUTPUT_DIR, "efficiency_frontier_base_case")

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# ==============================================================================
# 2. Load Data
# ==============================================================================
cat("Loading base case outcomes from:", INPUT_FILE, "\n")
outcomes <- read_csv(INPUT_FILE, show_col_types = FALSE)

# Validate required columns
required_cols <- c(
    "strategy_id", "strategy_label",
    "total_cost_per_proband_cad_mean", "diagnoses_per_proband_mean",
    "total_cost_per_proband_cad_ui_low", "total_cost_per_proband_cad_ui_high",
    "diagnoses_per_proband_ui_low", "diagnoses_per_proband_ui_high",
    "pr_at_least_one_vus_mean"
)
stopifnot("Missing required columns" = all(required_cols %in% names(outcomes)))
cat("Loaded", nrow(outcomes), "strategies\n")

# Source shared utility
source("scripts/utils_cea.R")

# ==============================================================================
# 4. Perform Analysis
# ==============================================================================
cat("\n--- Performing Incremental Analysis ---\n")

# Prepare data for analysis (keep VUS data for later joining)
analysis_data <- outcomes %>%
    select(
        strategy_id, strategy_label,
        cost = total_cost_per_proband_cad_mean,
        effect = diagnoses_per_proband_mean,
        cost_ui_low = total_cost_per_proband_cad_ui_low,
        cost_ui_high = total_cost_per_proband_cad_ui_high,
        effect_ui_low = diagnoses_per_proband_ui_low,
        effect_ui_high = diagnoses_per_proband_ui_high
    )

# Keep VUS data separately for joining after incremental analysis
vus_data <- outcomes %>%
    select(strategy_id, vus_prob = pr_at_least_one_vus_mean)

# Run incremental analysis
results <- perform_incremental_analysis(analysis_data)

# Join VUS data back to results
results <- results %>%
    left_join(vus_data, by = "strategy_id")

# Format output table per §7.5 spec
output_table <- results %>%
    arrange(cost) %>%
    mutate(
        nnt_probands_per_additional_diagnosis_mean = ifelse(
            !is.na(incremental_effect) & incremental_effect > 0,
            1 / incremental_effect,
            NA_real_
        )
    ) %>%
    select(
        strategy_id,
        strategy_label,
        total_cost_per_proband_cad_mean = cost,
        total_cost_per_proband_cad_ui_low = cost_ui_low,
        total_cost_per_proband_cad_ui_high = cost_ui_high,
        diagnoses_per_proband_mean = effect,
        diagnoses_per_proband_ui_low = effect_ui_low,
        diagnoses_per_proband_ui_high = effect_ui_high,
        incremental_cost_vs_prev_cad_mean = incremental_cost,
        incremental_diagnoses_vs_prev_mean = incremental_effect,
        icpd_cad_per_additional_diagnosis_mean = icer,
        nnt_probands_per_additional_diagnosis_mean,
        dominance_status
    )

incremental_required <- c(
    "strategy_id", "strategy_label",
    "total_cost_per_proband_cad_mean", "total_cost_per_proband_cad_ui_low", "total_cost_per_proband_cad_ui_high",
    "diagnoses_per_proband_mean", "diagnoses_per_proband_ui_low", "diagnoses_per_proband_ui_high",
    "incremental_cost_vs_prev_cad_mean", "incremental_diagnoses_vs_prev_mean",
    "icpd_cad_per_additional_diagnosis_mean", "nnt_probands_per_additional_diagnosis_mean",
    "dominance_status"
)
assert_required_columns(output_table, incremental_required, "incremental_table_base_case", exact = TRUE)
assert_no_na(output_table, c("strategy_id", "strategy_label", "dominance_status"), "incremental_table_base_case")
assert_frontier_first_incremental_na(output_table, table_name = "incremental_table_base_case")

# ==============================================================================
# 5. Export Table
# ==============================================================================
write_csv_validated(output_table, OUTPUT_TABLE, "incremental_table_base_case")
cat("\nIncremental analysis table exported to:", OUTPUT_TABLE, "\n")

# Print summary
cat("\n=== Incremental Analysis Results ===\n\n")
print(
    as.data.frame(output_table %>%
        select(
            strategy_label, total_cost_per_proband_cad_mean,
            diagnoses_per_proband_mean, icpd_cad_per_additional_diagnosis_mean,
            dominance_status
        )),
    row.names = FALSE
)

# ==============================================================================
# 6. Efficiency Frontier Plot with VUS Burden
# ==============================================================================
cat("\n--- Generating Efficiency Frontier Plot ---\n")

# Clean strategy labels for display and create multi-line labels with VUS%
results <- results %>%
    mutate(
        strategy_display = case_when(
            strategy_label == "Panel_Reflex_ES" ~ "Panel reflex to ES",
            TRUE ~ strategy_label
        ),
        # Multi-line label with VUS probability
        strategy_label_full = paste0(
            strategy_display, "\n(VUS: ", round(vus_prob * 100, 0), "%)"
        )
    )

# Prepare plot data with VUS-based bubble sizing
# Note: Invert VUS for size so LOWER VUS = LARGER bubble (more desirable)
# Or: Keep direct mapping where HIGHER VUS = LARGER bubble (shows burden)
plot_data <- results %>%
    mutate(
        is_dominated = dominance_status != "non_dominated",
        alpha_val = ifelse(is_dominated, 0.5, 1.0),
        # Scale VUS prob to reasonable point sizes (range ~3 to 8)
        bubble_size = 3 + (vus_prob * 8)
    )

# Get frontier line data (non-dominated strategies)
frontier_data <- results %>%
    filter(dominance_status == "non_dominated") %>%
    arrange(effect)

# Create plot
p <- ggplot(plot_data, aes(x = effect, y = cost)) +
    # Error bars (behind points)
    geom_errorbar(
        aes(ymin = cost_ui_low, ymax = cost_ui_high, alpha = alpha_val),
        width = 0.005, color = "gray60", linewidth = 0.5
    ) +
    geom_errorbarh(
        aes(xmin = effect_ui_low, xmax = effect_ui_high, alpha = alpha_val),
        height = 50, color = "gray60", linewidth = 0.5
    ) +
    # Frontier line
    geom_line(
        data = frontier_data,
        aes(x = effect, y = cost),
        color = "black", linewidth = 1.0, linetype = "solid"
    ) +
    # Strategy points - bubble size encodes VUS burden
    geom_point(
        aes(
            color = strategy_display,
            alpha = alpha_val,
            size = bubble_size
        )
    ) +
    # Strategy labels with VUS% (using ggrepel for non-overlap)
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
    # Scales
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
    # Labels
    labs(
        title = "Efficiency Frontier: Cost vs Diagnostic Yield",
        subtitle = "Bubble size represents VUS burden (larger = higher probability of VUS).\nError bars show 95% UI.",
        x = "Diagnostic Yield (Proportion of Probands Diagnosed)",
        y = "Total Cost per Proband (CAD)",
        color = NULL
    ) +
    # Theme - larger base size for publication quality
    theme_minimal(base_size = 13, base_family = "sans") +
    theme(
        # Grid - remove all gridlines
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        # Axes
        axis.line = element_line(color = "black", linewidth = 0.6),
        axis.ticks = element_line(color = "black", linewidth = 0.5),
        axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 13, face = "bold"),
        # Title
        plot.title = element_text(size = 15, face = "bold", hjust = 0),
        plot.subtitle = element_text(size = 11, color = "gray30", hjust = 0),
        # Legend
        legend.position = "none", # Direct labels via ggrepel
        # Background
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        # Margins
        plot.margin = margin(15, 15, 10, 10)
    )

# ==============================================================================
# 7. Export Figure
# ==============================================================================
ggsave(paste0(FIG_BASE, ".png"), p, width = 7, height = 5.5, dpi = 300, bg = "white")
ggsave(paste0(FIG_BASE, ".svg"), p, width = 7, height = 5.5, bg = "white")

cat("Figure exported to:\n")
cat(" -", paste0(FIG_BASE, ".png"), "\n")
cat(" -", paste0(FIG_BASE, ".svg"), "\n")

# ==============================================================================
# 8. Verification
# ==============================================================================
cat("\n--- Verification ---\n")

# Check ICER monotonicity for frontier
frontier_check <- output_table %>%
    filter(dominance_status == "non_dominated") %>%
    arrange(total_cost_per_proband_cad_mean)

if (nrow(frontier_check) > 1) {
    icers <- na.omit(frontier_check$icpd_cad_per_additional_diagnosis_mean)
    if (length(icers) > 1) {
        is_monotonic <- all(diff(icers) >= 0)
        if (is_monotonic) {
            cat("PASS: Frontier ICERs are monotonically increasing\n")
        } else {
            warning("Frontier ICERs are NOT monotonically increasing!")
            print(frontier_check %>% select(strategy_label, icpd_cad_per_additional_diagnosis_mean))
        }
    }
}

cat("\nIncremental analysis complete.\n")
