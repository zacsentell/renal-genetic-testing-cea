# scripts/06e_phenotype_figure.R
# Purpose: Phenotype-stratified cost-effectiveness figures
#   Figure 1: Reflex cost vs yield by phenotype (canonical results figure)
#   Figure 2: Trajectory of Value – Reflex strategy across phenotypes
#   Figure 3: All-modes efficiency landscape (small multiples)
#   Figure 4: Most cost-effective strategy by phenotype at WTP anchors

library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(scales)
if (!requireNamespace("config", quietly = TRUE)) {
    stop("Package 'config' is required. Install it with install.packages('config').")
}
source("scripts/utils_schema.R")

# ==============================================================================
# IO Paths
# ==============================================================================
INPUT_PHENO_CSV <- "outputs/results/supplement/phenotype_stratified/base_case/phenotype_iteration_outcomes.csv"
OUTPUT_DIR <- "outputs/results/supplement/phenotype_stratified"

OUTPUT_FIG_REFLEX_PNG <- file.path(OUTPUT_DIR, "reflex_cost_yield_by_phenotype.png")
OUTPUT_FIG_REFLEX_SVG <- file.path(OUTPUT_DIR, "reflex_cost_yield_by_phenotype.svg")
OUTPUT_FIG_TRAJ_PNG <- file.path(OUTPUT_DIR, "trajectory_of_value.png")
OUTPUT_FIG_TRAJ_SVG <- file.path(OUTPUT_DIR, "trajectory_of_value.svg")
OUTPUT_FIG_ALL_PNG <- file.path(OUTPUT_DIR, "phenotype_stratified_all_modes.png")
OUTPUT_FIG_ALL_SVG <- file.path(OUTPUT_DIR, "phenotype_stratified_all_modes.svg")
OUTPUT_FIG_WTP_PNG <- file.path(OUTPUT_DIR, "phenotype_wtp_winner_matrix.png")
OUTPUT_FIG_WTP_SVG <- file.path(OUTPUT_DIR, "phenotype_wtp_winner_matrix.svg")
OUTPUT_WTP_LONG_CSV <- file.path(OUTPUT_DIR, "phenotype_wtp_winner_long.csv")
OUTPUT_WTP_WIDE_CSV <- file.path(OUTPUT_DIR, "phenotype_wtp_winner_wide.csv")

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# ==============================================================================
# Configuration
# ==============================================================================
cfg <- config::get()
WTP_MIN <- cfg$wtp$min
WTP_MAX <- cfg$wtp$max
WTP_STEP <- cfg$wtp$step
WTP_ANCHORS <- cfg$wtp$anchors

if (is.null(WTP_MIN) || is.null(WTP_MAX) || is.null(WTP_STEP)) {
    stop("Missing WTP configuration in config.yml (wtp:min/max/step).")
}
if (!is.numeric(WTP_STEP) || WTP_STEP <= 0) {
    stop("WTP step must be a positive numeric value.")
}
if (is.null(WTP_ANCHORS) || length(WTP_ANCHORS) == 0) {
    stop("Missing WTP anchors in config.yml (wtp:anchors).")
}
if (!is.numeric(WTP_ANCHORS)) {
    stop("WTP anchors must be numeric values.")
}

WTP_ANCHORS <- sort(unique(as.numeric(WTP_ANCHORS)))
if (any(WTP_ANCHORS < WTP_MIN | WTP_ANCHORS > WTP_MAX)) {
    stop("All WTP anchors must lie within [wtp:min, wtp:max].")
}
step_units <- (WTP_ANCHORS - WTP_MIN) / WTP_STEP
if (any(abs(step_units - round(step_units)) > 1e-9)) {
    stop("All WTP anchors must align to the configured wtp:step interval.")
}
if (length(WTP_ANCHORS) < 2) {
    stop("Provide at least two WTP anchors in config.yml for matrix display.")
}

cat("Using WTP anchors (CAD):", paste(WTP_ANCHORS, collapse = ", "), "\n")

# ==============================================================================
# Load and aggregate
# ==============================================================================
cat("Loading phenotype outcomes from:", INPUT_PHENO_CSV, "\n")
df <- read_csv(INPUT_PHENO_CSV, show_col_types = FALSE)

required_cols <- c(
    "iteration_id", "phenotype_category", "strategy_id", "strategy_label",
    "n_probands", "total_cost_per_proband_cad", "diagnoses_per_proband",
    "cost_per_diagnosis_cad"
)
assert_required_columns(df, required_cols, "phenotype_iteration_outcomes")
assert_no_na(
    df, c("iteration_id", "phenotype_category", "strategy_id", "strategy_label"),
    "phenotype_iteration_outcomes"
)

pheno_levels <- c("CKDu", "glomerular", "tubulointerstitial", "tubulopathies", "cystic")
pheno_labels <- c("CKDu", "Glomerular", "Tubulointerstitial", "Tubulopathies", "Cystic")

strategy_display_label <- function(x) {
    dplyr::case_when(
        x == "Panel_Reflex_ES" ~ "Reflex (Panel->ES)",
        TRUE ~ x
    )
}

strategy_display_levels <- c("Panel", "ES", "Reflex (Panel->ES)", "GS")
strategy_cols_wtp_all <- setNames(
    scales::viridis_pal(option = "D", begin = 0.15, end = 0.85, direction = 1)(
        length(strategy_display_levels)
    ),
    strategy_display_levels
)

df_summary <- df %>%
    group_by(strategy_label, phenotype_category) %>%
    summarise(
        yield_mean = mean(diagnoses_per_proband),
        cost_per_diag_mean = mean(cost_per_diagnosis_cad, na.rm = TRUE),
        cost_mean = mean(total_cost_per_proband_cad),
        .groups = "drop"
    ) %>%
    mutate(
        phenotype_factor = factor(phenotype_category, levels = pheno_levels, labels = pheno_labels)
    )

# ==============================================================================
# Figure 1: Reflex cost vs diagnostic yield by phenotype
# ==============================================================================
df_reflex_scatter <- df_summary %>%
    filter(strategy_label == "Panel_Reflex_ES") %>%
    arrange(yield_mean)

p_reflex <- ggplot(df_reflex_scatter, aes(x = yield_mean, y = cost_mean)) +
    geom_point(aes(fill = phenotype_factor), shape = 21, size = 6, color = "black", stroke = 1.2) +
    geom_label_repel(
        aes(label = paste0(phenotype_factor, "\n", percent(yield_mean, accuracy = 0.1))),
        box.padding = 0.7,
        point.padding = 0.5,
        size = 4,
        fontface = "bold",
        fill = "white",
        label.size = 0.2,
        min.segment.length = 0.1
    ) +
    scale_x_continuous(
        labels = percent_format(accuracy = 1),
        name   = "Diagnostic Yield (Proportion Diagnosed)",
        expand = expansion(mult = 0.15)
    ) +
    scale_y_continuous(
        labels = dollar_format(prefix = "$"),
        name   = "Total Cost per Proband (CAD)",
        expand = expansion(mult = c(0.05, 0.12))
    ) +
    scale_fill_viridis_d(option = "D", begin = 0.15, end = 0.85, guide = "none") +
    labs(
        title    = "Reflex Testing: Cost and Yield by Phenotype",
        subtitle = "Each point represents a phenotype category under the Panel -> ES reflex strategy"
    ) +
    theme_minimal(base_size = 13, base_family = "sans") +
    theme(
        panel.grid.minor  = element_blank(),
        panel.grid.major  = element_blank(),
        axis.line         = element_line(color = "black", linewidth = 0.6),
        axis.ticks        = element_line(color = "black", linewidth = 0.5),
        axis.text         = element_text(color = "black", size = 11),
        axis.title        = element_text(color = "black", size = 12, face = "bold"),
        plot.title        = element_text(size = 14, face = "bold", hjust = 0),
        plot.subtitle     = element_text(size = 10, color = "gray30", hjust = 0),
        panel.background  = element_rect(fill = "white", color = NA),
        plot.background   = element_rect(fill = "white", color = NA),
        plot.margin       = margin(15, 15, 10, 10)
    )

cat("Saving Figure 1 (Reflex cost vs yield) to:", OUTPUT_FIG_REFLEX_PNG, "\n")
ggsave(OUTPUT_FIG_REFLEX_PNG, plot = p_reflex, width = 8, height = 6, dpi = 300, bg = "white")
ggsave(OUTPUT_FIG_REFLEX_SVG, plot = p_reflex, width = 8, height = 6, bg = "white")

# ==============================================================================
# Figure 2: Trajectory of Value – Reflex strategy
# ==============================================================================
df_reflex <- df_summary %>%
    filter(strategy_label == "Panel_Reflex_ES") %>%
    arrange(yield_mean)

p_traj <- ggplot(df_reflex, aes(x = yield_mean, y = cost_per_diag_mean)) +
    geom_path(
        arrow     = arrow(length = unit(0.3, "cm"), type = "closed"),
        color     = "grey60",
        linewidth = 1,
        linetype  = "solid"
    ) +
    geom_point(aes(fill = yield_mean, size = yield_mean), shape = 21, color = "black", stroke = 1.2) +
    geom_text_repel(
        aes(label = paste0(phenotype_factor, "\n", percent(yield_mean, accuracy = 0.1))),
        box.padding = 0.6,
        point.padding = 0.6,
        size = 4.5,
        fontface = "bold",
        color = "black",
        bg.color = "white",
        bg.r = 0.15
    ) +
    scale_x_continuous(labels = percent_format(accuracy = 1), name = "Diagnostic Yield") +
    scale_y_continuous(
        labels = dollar_format(prefix = "$", suffix = " CAD"),
        name   = "Cost per Diagnosis",
        limits = c(0, NA)
    ) +
    scale_fill_viridis_c(option = "C", name = "Yield", guide = "none") +
    scale_size_continuous(range = c(6, 12), guide = "none") +
    annotate("text",
        x = max(df_reflex$yield_mean), y = -Inf, vjust = -0.5, hjust = 1.0,
        label = "High Efficiency Zone", size = 5, fontface = "bold", color = "#E65100"
    ) +
    labs(
        title    = "The Trajectory of Value",
        subtitle = "Reflex Strategy: Phenotype enrichment drives higher yield and lower cost per diagnosis"
    ) +
    theme_classic(base_size = 14) +
    theme(
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 18),
        plot.subtitle = element_text(face = "italic", size = 12)
    )

cat("Saving Figure 2 (Trajectory of Value) to:", OUTPUT_FIG_TRAJ_PNG, "\n")
ggsave(OUTPUT_FIG_TRAJ_PNG, plot = p_traj, width = 9, height = 7, dpi = 300, bg = "white")
ggsave(OUTPUT_FIG_TRAJ_SVG, plot = p_traj, width = 9, height = 7, bg = "white")

# ==============================================================================
# Figure 3: All-modes efficiency landscape (small multiples)
# ==============================================================================
strategy_cols <- c(
    "Panel" = "grey50",
    "Panel_Reflex_ES" = "#E65100",
    "ES" = "#1976D2",
    "GS" = "#7B1FA2"
)
strategy_labels_clean <- c(
    "Panel" = "Panel",
    "Panel_Reflex_ES" = "Reflex",
    "ES" = "Exome",
    "GS" = "Genome"
)

arrow_data <- df_summary %>%
    select(phenotype_factor, strategy_label, yield_mean, cost_per_diag_mean) %>%
    filter(strategy_label %in% c("Panel", "Panel_Reflex_ES")) %>%
    tidyr::pivot_wider(names_from = strategy_label, values_from = c(yield_mean, cost_per_diag_mean))

p_all <- ggplot() +
    geom_segment(
        data = arrow_data,
        aes(
            x = yield_mean_Panel, y = cost_per_diag_mean_Panel,
            xend = yield_mean_Panel_Reflex_ES, yend = cost_per_diag_mean_Panel_Reflex_ES
        ),
        arrow = arrow(length = unit(0.2, "cm")),
        color = "grey70",
        linewidth = 0.8
    ) +
    geom_point(
        data = df_summary,
        aes(x = yield_mean, y = cost_per_diag_mean, color = strategy_label),
        size = 4, alpha = 0.9
    ) +
    facet_wrap(~phenotype_factor, nrow = 1, scales = "fixed") +
    scale_color_manual(values = strategy_cols, labels = strategy_labels_clean, name = "Strategy") +
    scale_x_continuous(labels = percent_format(accuracy = 1), name = "Diagnostic Yield") +
    scale_y_continuous(
        labels = dollar_format(prefix = "$", scale = 1e-3, suffix = "k"),
        name   = "Cost per Diagnosis (CAD)"
    ) +
    theme_bw(base_size = 14) +
    theme(
        strip.background = element_rect(fill = "grey95"),
        strip.text       = element_text(face = "bold"),
        legend.position  = "bottom",
        panel.grid.minor = element_blank()
    ) +
    labs(
        title    = "Efficiency Landscape by Phenotype",
        subtitle = "Comparing strategies across clinical subgroups. Arrow indicates upgrade from Panel to Reflex."
    )

cat("Saving Figure 3 (All Modes) to:", OUTPUT_FIG_ALL_PNG, "\n")
ggsave(OUTPUT_FIG_ALL_PNG, plot = p_all, width = 14, height = 5, dpi = 300, bg = "white")
ggsave(OUTPUT_FIG_ALL_SVG, plot = p_all, width = 14, height = 5, bg = "white")

# ==============================================================================
# Figure 4 + output tables: optimal strategy by phenotype and WTP anchors
# ==============================================================================
cat("Computing phenotype-specific WTP winners...\n")

pheno_available <- intersect(pheno_levels, sort(unique(df$phenotype_category)))
strategy_lookup <- df %>%
    distinct(strategy_id, strategy_label) %>%
    arrange(strategy_id)

df_ce <- df %>%
    filter(!is.na(total_cost_per_proband_cad), !is.na(diagnoses_per_proband))

prob_rows <- vector("list", length(WTP_ANCHORS))
winner_rows <- vector("list", length(WTP_ANCHORS))
tie_groups <- 0L

for (k in seq_along(WTP_ANCHORS)) {
    wtp <- WTP_ANCHORS[k]

    wtp_data <- df_ce %>%
        mutate(
            wtp_threshold_cad = wtp,
            nmb = (wtp * diagnoses_per_proband) - total_cost_per_proband_cad
        )

    ties_k <- wtp_data %>%
        group_by(iteration_id, phenotype_category) %>%
        mutate(max_nmb = max(nmb, na.rm = TRUE)) %>%
        filter(abs(nmb - max_nmb) < 1e-12) %>%
        summarise(n_tied = n(), .groups = "drop") %>%
        filter(n_tied > 1)
    tie_groups <- tie_groups + nrow(ties_k)

    winners <- wtp_data %>%
        arrange(iteration_id, phenotype_category, desc(nmb), strategy_id) %>%
        group_by(iteration_id, phenotype_category) %>%
        slice(1) %>%
        ungroup()

    n_iter_by_pheno <- wtp_data %>%
        distinct(iteration_id, phenotype_category) %>%
        count(phenotype_category, name = "n_iter")

    pr_raw <- winners %>%
        count(phenotype_category, strategy_id, strategy_label, name = "n_optimal") %>%
        left_join(n_iter_by_pheno, by = "phenotype_category") %>%
        mutate(
            wtp_threshold_cad = wtp,
            pr_cost_effective = n_optimal / n_iter
        ) %>%
        select(wtp_threshold_cad, phenotype_category, strategy_id, strategy_label, pr_cost_effective)

    pr_complete <- tidyr::expand_grid(
        wtp_threshold_cad = wtp,
        phenotype_category = pheno_available,
        strategy_id = strategy_lookup$strategy_id
    ) %>%
        left_join(strategy_lookup, by = "strategy_id") %>%
        left_join(pr_raw, by = c("wtp_threshold_cad", "phenotype_category", "strategy_id", "strategy_label")) %>%
        mutate(pr_cost_effective = dplyr::coalesce(pr_cost_effective, 0))

    mean_nmb <- wtp_data %>%
        group_by(wtp_threshold_cad, phenotype_category, strategy_id, strategy_label) %>%
        summarise(mean_nmb = mean(nmb), .groups = "drop")

    winner_prob <- pr_complete %>%
        arrange(wtp_threshold_cad, phenotype_category, desc(pr_cost_effective), strategy_id) %>%
        group_by(wtp_threshold_cad, phenotype_category) %>%
        slice(1) %>%
        ungroup() %>%
        rename(
            winner_strategy_id = strategy_id,
            winner_strategy_label = strategy_label,
            winner_probability = pr_cost_effective
        ) %>%
        select(
            wtp_threshold_cad, phenotype_category,
            winner_strategy_id, winner_strategy_label, winner_probability
        )

    winner_nmb <- mean_nmb %>%
        rename(
            winner_strategy_id = strategy_id,
            winner_strategy_label = strategy_label,
            winner_mean_nmb_cad = mean_nmb
        )

    runner_up_nmb <- mean_nmb %>%
        inner_join(
            winner_prob %>% select(wtp_threshold_cad, phenotype_category, winner_strategy_id),
            by = c("wtp_threshold_cad", "phenotype_category")
        ) %>%
        filter(strategy_id != winner_strategy_id) %>%
        arrange(wtp_threshold_cad, phenotype_category, desc(mean_nmb), strategy_id) %>%
        group_by(wtp_threshold_cad, phenotype_category) %>%
        slice(1) %>%
        ungroup() %>%
        rename(
            runner_up_strategy_id = strategy_id,
            runner_up_strategy_label = strategy_label,
            runner_up_mean_nmb_cad = mean_nmb
        ) %>%
        select(
            wtp_threshold_cad, phenotype_category,
            runner_up_strategy_id, runner_up_strategy_label, runner_up_mean_nmb_cad
        )

    winner_rows[[k]] <- winner_prob %>%
        left_join(
            winner_nmb,
            by = c("wtp_threshold_cad", "phenotype_category", "winner_strategy_id", "winner_strategy_label")
        ) %>%
        left_join(runner_up_nmb, by = c("wtp_threshold_cad", "phenotype_category")) %>%
        mutate(nmb_margin_vs_runner_up_cad = winner_mean_nmb_cad - runner_up_mean_nmb_cad)

    prob_rows[[k]] <- pr_complete
}

probability_all <- bind_rows(prob_rows)
winner_long <- bind_rows(winner_rows) %>%
    mutate(
        phenotype = factor(phenotype_category, levels = pheno_levels, labels = pheno_labels),
        winner_strategy_display = factor(
            strategy_display_label(winner_strategy_label),
            levels = strategy_display_levels
        )
    ) %>%
    arrange(phenotype, wtp_threshold_cad)

strategy_levels_present <- strategy_display_levels[
    strategy_display_levels %in% unique(as.character(winner_long$winner_strategy_display))
]
winner_long <- winner_long %>%
    mutate(
        winner_strategy_display = factor(
            as.character(winner_strategy_display),
            levels = strategy_levels_present
        )
    )

assert_probabilities_sum_to_one(
    probability_all,
    c("wtp_threshold_cad", "phenotype_category"),
    "pr_cost_effective",
    "phenotype_wtp_strategy_probabilities"
)

winner_long_out <- winner_long %>%
    transmute(
        phenotype = as.character(phenotype),
        phenotype_category,
        wtp_threshold_cad,
        winner_strategy_id,
        winner_strategy_label,
        winner_probability,
        winner_mean_nmb_cad,
        runner_up_strategy_id,
        runner_up_strategy_label,
        runner_up_mean_nmb_cad,
        nmb_margin_vs_runner_up_cad
    )

assert_required_columns(
    winner_long_out,
    c(
        "phenotype", "phenotype_category", "wtp_threshold_cad",
        "winner_strategy_id", "winner_strategy_label", "winner_probability",
        "winner_mean_nmb_cad", "runner_up_strategy_id", "runner_up_strategy_label",
        "runner_up_mean_nmb_cad", "nmb_margin_vs_runner_up_cad"
    ),
    "phenotype_wtp_winner_long",
    exact = TRUE
)

winner_wide_out <- winner_long %>%
    transmute(
        phenotype = as.character(phenotype),
        wtp_col = paste0("wtp_", wtp_threshold_cad),
        winner_with_probability = paste0(
            winner_strategy_display,
            " (", percent(winner_probability, accuracy = 1), ")"
        )
    ) %>%
    pivot_wider(
        names_from = wtp_col,
        values_from = winner_with_probability
    )

assert_required_columns(winner_wide_out, c("phenotype"), "phenotype_wtp_winner_wide")

write_csv_validated(winner_long_out, OUTPUT_WTP_LONG_CSV, "phenotype_wtp_winner_long")
write_csv_validated(winner_wide_out, OUTPUT_WTP_WIDE_CSV, "phenotype_wtp_winner_wide")

wtp_values_sorted <- sort(unique(winner_long$wtp_threshold_cad))
wtp_labels <- setNames(dollar_format(prefix = "$", suffix = "")(wtp_values_sorted), wtp_values_sorted)

strategy_cols_wtp <- strategy_cols_wtp_all[strategy_levels_present]

p_wtp <- ggplot(
    winner_long,
    aes(
        x = factor(wtp_threshold_cad, levels = wtp_values_sorted),
        y = phenotype,
        fill = winner_strategy_display
    )
) +
    geom_tile(color = "gray97", linewidth = 1.2, width = 0.96, height = 0.96) +
    geom_text(
        aes(label = percent(winner_probability, accuracy = 1)),
        color = "white",
        size = 4.0,
        fontface = "bold"
    ) +
    scale_x_discrete(labels = wtp_labels) +
    scale_fill_manual(values = strategy_cols_wtp, drop = FALSE) +
    labs(
        title = "Most Cost-Effective Strategy by Phenotype",
        subtitle = "Tile color shows the winning strategy; labels show probability of optimality",
        x = "Willingness-to-Pay per Additional Diagnosis (CAD)",
        y = "Phenotype",
        fill = "Winning Strategy"
    ) +
    coord_fixed(ratio = 1) +
    theme_minimal(base_size = 12, base_family = "sans") +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.6),
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.ticks = element_line(color = "black", linewidth = 0.4),
        axis.text = element_text(color = "black", size = 10.5),
        axis.title = element_text(color = "black", size = 11.5, face = "bold"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0),
        plot.subtitle = element_text(size = 10, color = "gray30", hjust = 0, lineheight = 1.1),
        legend.position = "bottom",
        legend.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 9.5),
        legend.key.width = unit(1.2, "lines"),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        plot.margin = margin(12, 14, 8, 10)
    )

cat("Saving Figure 4 (WTP winner matrix) to:", OUTPUT_FIG_WTP_PNG, "\n")
ggsave(OUTPUT_FIG_WTP_PNG, plot = p_wtp, width = 10.5, height = 6.2, dpi = 300, bg = "white")
ggsave(OUTPUT_FIG_WTP_SVG, plot = p_wtp, width = 10.5, height = 6.2, bg = "white")

cat("Deterministic tie groups resolved by strategy_id ordering:", tie_groups, "\n")

cat("Done. All figures written to:", OUTPUT_DIR, "\n")
