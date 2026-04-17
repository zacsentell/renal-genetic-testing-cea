# 10_clinical_impact_analysis.R
# Purpose: Join gene-level clinical utility annotations to simulation diagnosis
#          counts and compute clinical impact proportions by strategy.
#
# Inputs:
#   - data/raw/clinical_utility/gene_clinical_utility.csv (curated)
#   - outputs/results/clinical_impact/iteration_gene_counts.csv (from 05)
#   - outputs/results/clinical_impact/gs_iteration_gene_counts.csv (from 05b)
#   - outputs/results/base_case/iteration_level/strategy_iteration_outcomes.csv
#   - outputs/results/scenario_analysis/gs_uplift/gs_baseline_iteration_outcomes.csv
#
# Outputs:
#   - outputs/results/clinical_impact/high_impact_outcomes_by_strategy.csv
#   - outputs/results/clinical_impact/gene_coverage_audit.csv
#   - outputs/results/clinical_impact/impact_by_phenotype.csv
#   - outputs/results/clinical_impact/gene_composition_by_impact.csv
#   - outputs/results/clinical_impact/clinical_impact_composite.png/.svg

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(patchwork)

source("scripts/utils_schema.R")

# ==============================================================================
# 1. Configuration
# ==============================================================================
UTILITY_TABLE_PATH <- "data/raw/clinical_utility/gene_clinical_utility.csv"
GENE_COUNTS_PATH <- "outputs/results/clinical_impact/iteration_gene_counts.csv"
GS_GENE_COUNTS_PATH <- "outputs/results/clinical_impact/gs_iteration_gene_counts.csv"
ITER_OUTCOMES_PATH <- "outputs/results/base_case/iteration_level/strategy_iteration_outcomes.csv"
GS_ITER_OUTCOMES_PATH <- "outputs/results/scenario_analysis/gs_uplift/gs_baseline_iteration_outcomes.csv"
SIM_GENE_PATH <- "outputs/parameters/joint_proportions_full.csv"
PHENO_ITER_PATH <- "outputs/results/supplement/phenotype_stratified/base_case/phenotype_iteration_outcomes.csv"

OUTPUT_DIR <- "outputs/results/clinical_impact"
OUTPUT_TABLE1 <- file.path(OUTPUT_DIR, "high_impact_outcomes_by_strategy.csv")
OUTPUT_TABLE2 <- file.path(OUTPUT_DIR, "gene_coverage_audit.csv")
OUTPUT_IMPACT_PHENO_CSV <- file.path(OUTPUT_DIR, "impact_by_phenotype.csv")
OUTPUT_GENE_COMP_CSV <- file.path(OUTPUT_DIR, "gene_composition_by_impact.csv")
FIG_COMPOSITE_BASE <- file.path(OUTPUT_DIR, "clinical_impact_composite")
LEGACY_FIGURES <- file.path(
    OUTPUT_DIR,
    c(
        "clinical_impact_by_strategy.png",
        "clinical_impact_by_strategy.svg",
        "impact_by_phenotype.png",
        "impact_by_phenotype.svg",
        "gene_composition_by_impact.png",
        "gene_composition_by_impact.svg"
    )
)

# Reference strategy for Figure 2 gene composition view.
# Panel_Reflex_ES is the provincially-reimbursed baseline strategy.
FIG2_REF_STRATEGY <- "Panel_Reflex_ES"

# Phenotype ordering used in Figures 1 and 2. Cystic first, then glomerular,
# then the two tubular buckets, with unknown-cause CKDu last.
# Values match the raw cohort phenotype labels (lowercase except CKDu);
# PHENOTYPE_DISPLAY provides the Title Case labels shown on axes and legends.
PHENOTYPE_LEVELS <- c(
    "cystic", "glomerular", "tubulointerstitial",
    "tubulopathies", "CKDu"
)
PHENOTYPE_DISPLAY <- c(
    "cystic"             = "Cystic",
    "glomerular"         = "Glomerular",
    "tubulointerstitial" = "Tubulointerstitial",
    "tubulopathies"      = "Tubulopathies",
    "CKDu"               = "CKDu"
)

# Fixed phenotype palette (ColorBrewer Dark2) re-used across Figures 1 and 2
# so readers can visually link phenotype attribution between views. Names
# use the display labels so manual scale_fill lookup works on the plot data.
PHENOTYPE_PALETTE <- c(
    "Cystic"             = "#1b9e77",
    "Glomerular"         = "#d95f02",
    "Tubulointerstitial" = "#7570b3",
    "Tubulopathies"      = "#e7298a",
    "CKDu"               = "#66a61e"
)

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
legacy_existing <- LEGACY_FIGURES[file.exists(LEGACY_FIGURES)]
if (length(legacy_existing) > 0) {
    removed_legacy <- file.remove(legacy_existing)
}

# Category columns
# Framework: three non-exclusive primary categories capturing distinct clinical
# actions triggered by the molecular diagnosis. `therapeutic` is further
# partitioned into three non-exclusive sub-categories for granularity.
PRIMARY_CATS <- c("therapeutic", "extrarenal_surveillance", "transplant_donor")
THERAPEUTIC_SUBS <- c("targeted_drug", "avoid_ineffective_therapy", "supportive_dietary")
ALL_CATS <- c("high_impact", PRIMARY_CATS, THERAPEUTIC_SUBS)

# ==============================================================================
# 2. Load and validate utility table
# ==============================================================================
cat("Loading clinical utility table...\n")

if (!file.exists(UTILITY_TABLE_PATH)) {
    stop(
        "Curated utility table not found: ", UTILITY_TABLE_PATH,
        "\nSee data/raw/clinical_utility/curation_methodology.md for provenance."
    )
}

utility <- read_csv(UTILITY_TABLE_PATH, show_col_types = FALSE)

utility_required <- c("gene_symbol", "high_impact", PRIMARY_CATS, THERAPEUTIC_SUBS)
assert_required_columns(utility, utility_required, "gene_clinical_utility")

# Validate: binary 0/1 for category columns
for (col in ALL_CATS) {
    vals <- utility[[col]]
    if (any(!is.na(vals) & !(vals %in% c(0L, 1L)))) {
        stop("Column '", col, "' in utility table contains non-binary values")
    }
}

# Validate: if therapeutic=1, at least one sub-category must be 1
therapeutic_rows <- utility %>% filter(therapeutic == 1)
if (nrow(therapeutic_rows) > 0) {
    has_sub <- therapeutic_rows$targeted_drug == 1 |
        therapeutic_rows$avoid_ineffective_therapy == 1 |
        therapeutic_rows$supportive_dietary == 1
    if (any(!has_sub, na.rm = TRUE)) {
        bad_genes <- therapeutic_rows$gene_symbol[!has_sub]
        stop(
            "Genes with therapeutic=1 but no sub-category: ",
            paste(bad_genes, collapse = ", ")
        )
    }
}

# Validate: high_impact = OR of primaries (invariant 3 of curation methodology)
or_primaries <- as.integer(utility$therapeutic == 1 |
    utility$extrarenal_surveillance == 1 |
    utility$transplant_donor == 1)
mismatch_rows <- which(utility$high_impact != or_primaries)
if (length(mismatch_rows) > 0) {
    bad_genes <- utility$gene_symbol[mismatch_rows]
    stop(
        "Invariant violated: high_impact must equal the OR of therapeutic, ",
        "extrarenal_surveillance, transplant_donor. Offending gene(s): ",
        paste(bad_genes, collapse = ", ")
    )
}

cat("  Utility table loaded:", nrow(utility), "genes\n")

# ==============================================================================
# 3. Load gene counts and iteration outcomes
# ==============================================================================
cat("Loading gene counts and iteration outcomes...\n")

gene_counts <- read_csv(GENE_COUNTS_PATH, show_col_types = FALSE)
gs_gene_counts <- read_csv(GS_GENE_COUNTS_PATH, show_col_types = FALSE)

# Combine base case and GS gene counts
all_gene_counts <- bind_rows(gene_counts, gs_gene_counts)

# Load iteration outcomes for cost-per-high-impact calculation
iter_outcomes <- read_csv(ITER_OUTCOMES_PATH, show_col_types = FALSE)
gs_iter_outcomes <- read_csv(GS_ITER_OUTCOMES_PATH, show_col_types = FALSE) %>%
    select(
        iteration_id, strategy_id, strategy_label,
        total_cost_per_proband_cad, diagnoses_per_proband
    )

all_iter_outcomes <- bind_rows(
    iter_outcomes %>% select(
        iteration_id, strategy_id, strategy_label,
        total_cost_per_proband_cad, diagnoses_per_proband
    ),
    gs_iter_outcomes
)

# Load phenotype-level iteration data for all-proband denominators (Fig 8.2B)
pheno_iter_raw <- read_csv(PHENO_ITER_PATH, show_col_types = FALSE) %>%
    select(iteration_id, phenotype_category, strategy_id, strategy_label, n_probands)

cat("  Gene counts:", nrow(all_gene_counts), "rows\n")
cat("  Strategies:", paste(unique(all_gene_counts$strategy_label), collapse = ", "), "\n")

# ==============================================================================
# 4. Join gene counts to utility table
# ==============================================================================
cat("Joining gene counts to utility table...\n")

# Check for genes in simulation but not in utility table.
# Missing genes would silently propagate NA through every per-iteration
# proportion and then get stripped at the mean() step, producing output
# CSVs that look fine but undercount clinical-impact categories. Fail
# fast so curation errors are caught at the first run instead.
sim_genes_in_counts <- unique(all_gene_counts$gene_symbol)
missing_genes <- setdiff(sim_genes_in_counts, utility$gene_symbol)

if (length(missing_genes) > 0) {
    stop(
        "Genes in simulation output not found in utility table: ",
        paste(missing_genes, collapse = ", "),
        "\nAdd rows to ", UTILITY_TABLE_PATH,
        " before rerunning the clinical impact analysis."
    )
}

# Left join: gene counts -> utility categories
counts_with_utility <- all_gene_counts %>%
    left_join(
        utility %>% select(gene_symbol, all_of(ALL_CATS)),
        by = "gene_symbol"
    )

# ==============================================================================
# 5. Compute per-iteration category proportions
# ==============================================================================
cat("Computing category proportions per iteration...\n")

# Two denominators are computed per iteration x strategy:
#   (a) n_diagnosed_total: total diagnosed probands (primary metric).
#       "Among diagnosed probands, what fraction have a high-impact result?"
#   (b) n_probands: all probands (used only for cost-per-high-impact-dx).
#       Cost is borne by all tested probands, so cost/high-impact-dx
#       uses the all-proband proportion.
iter_category_counts <- counts_with_utility %>%
    group_by(iteration_id, strategy_id, strategy_label, n_probands) %>%
    summarise(
        n_diagnosed_total = sum(n_diagnosed),
        across(
            all_of(ALL_CATS),
            ~ sum(n_diagnosed[. == 1], na.rm = FALSE),
            .names = "n_{.col}"
        ),
        .groups = "drop"
    )

iter_category_props <- iter_category_counts %>%
    mutate(
        across(
            all_of(paste0("n_", ALL_CATS)),
            ~ ifelse(n_diagnosed_total > 0, . / n_diagnosed_total, NA_real_),
            .names = "{gsub('n_', 'prop_of_diagnosed_', .col)}"
        ),
        # All-proband proportion: only for high_impact (cost calculation)
        prop_of_all_probands_high_impact = n_high_impact / n_probands
    )

# Join iteration-level cost for cost-per-high-impact calculation
iter_category_props <- iter_category_props %>%
    left_join(
        all_iter_outcomes %>% select(
            iteration_id, strategy_id,
            total_cost_per_proband_cad
        ),
        by = c("iteration_id", "strategy_id")
    ) %>%
    mutate(
        cost_per_high_impact_dx = ifelse(
            prop_of_all_probands_high_impact > 0,
            total_cost_per_proband_cad / prop_of_all_probands_high_impact,
            NA_real_
        )
    )

# ==============================================================================
# 6. Aggregate across iterations: mean and 95% UI
# ==============================================================================
cat("Aggregating across iterations...\n")

# Build long-form results for Table 1.
# Column naming encodes the denominator:
#   prop_of_diagnosed_* : numerator / diagnosed probands (primary metric)
#   cost_per_high_impact_dx_* : total cost per proband / (high-impact / all probands)
results_list <- list()

for (strat_id in unique(iter_category_props$strategy_id)) {
    strat_data <- iter_category_props %>% filter(strategy_id == strat_id)
    strat_label <- strat_data$strategy_label[1]

    # High impact aggregate (includes cost per high impact dx)
    results_list[[length(results_list) + 1]] <- data.frame(
        strategy_id = strat_id,
        strategy_label = strat_label,
        level = "aggregate",
        category = "high_impact",
        prop_of_diagnosed_mean = mean(strat_data$prop_of_diagnosed_high_impact, na.rm = TRUE),
        prop_of_diagnosed_ui_low = quantile(strat_data$prop_of_diagnosed_high_impact, 0.025, na.rm = TRUE),
        prop_of_diagnosed_ui_high = quantile(strat_data$prop_of_diagnosed_high_impact, 0.975, na.rm = TRUE),
        cost_per_high_impact_dx_mean = mean(strat_data$cost_per_high_impact_dx, na.rm = TRUE),
        cost_per_high_impact_dx_ui_low = quantile(strat_data$cost_per_high_impact_dx, 0.025, na.rm = TRUE),
        cost_per_high_impact_dx_ui_high = quantile(strat_data$cost_per_high_impact_dx, 0.975, na.rm = TRUE),
        stringsAsFactors = FALSE
    )

    # Primary categories
    for (cat in PRIMARY_CATS) {
        prop_col <- paste0("prop_of_diagnosed_", cat)
        results_list[[length(results_list) + 1]] <- data.frame(
            strategy_id = strat_id,
            strategy_label = strat_label,
            level = "primary",
            category = cat,
            prop_of_diagnosed_mean = mean(strat_data[[prop_col]], na.rm = TRUE),
            prop_of_diagnosed_ui_low = quantile(strat_data[[prop_col]], 0.025, na.rm = TRUE),
            prop_of_diagnosed_ui_high = quantile(strat_data[[prop_col]], 0.975, na.rm = TRUE),
            cost_per_high_impact_dx_mean = NA_real_,
            cost_per_high_impact_dx_ui_low = NA_real_,
            cost_per_high_impact_dx_ui_high = NA_real_,
            stringsAsFactors = FALSE
        )
    }

    # Therapeutic sub-categories
    for (cat in THERAPEUTIC_SUBS) {
        prop_col <- paste0("prop_of_diagnosed_", cat)
        results_list[[length(results_list) + 1]] <- data.frame(
            strategy_id = strat_id,
            strategy_label = strat_label,
            level = "therapeutic_sub",
            category = cat,
            prop_of_diagnosed_mean = mean(strat_data[[prop_col]], na.rm = TRUE),
            prop_of_diagnosed_ui_low = quantile(strat_data[[prop_col]], 0.025, na.rm = TRUE),
            prop_of_diagnosed_ui_high = quantile(strat_data[[prop_col]], 0.975, na.rm = TRUE),
            cost_per_high_impact_dx_mean = NA_real_,
            cost_per_high_impact_dx_ui_low = NA_real_,
            cost_per_high_impact_dx_ui_high = NA_real_,
            stringsAsFactors = FALSE
        )
    }
}

table1 <- do.call(rbind, results_list)

# ==============================================================================
# 7. Gene coverage audit table
# ==============================================================================
cat("Building gene coverage audit table...\n")

sim_genes <- read_csv(SIM_GENE_PATH, show_col_types = FALSE)
sim_gene_list <- sort(unique(sim_genes$gene))

# Full outer join of simulation genes and utility table genes
all_gene_symbols <- sort(unique(c(sim_gene_list, utility$gene_symbol)))

table2 <- data.frame(gene_symbol = all_gene_symbols, stringsAsFactors = FALSE) %>%
    mutate(in_simulation = gene_symbol %in% sim_gene_list) %>%
    left_join(
        utility %>%
            mutate(
                source = ifelse(is.na(knoers_table) | knoers_table == "",
                    "extension", "knoers_2022"
                ),
                in_utility_table = TRUE
            ) %>%
            select(gene_symbol, all_of(ALL_CATS), source, rationale, in_utility_table),
        by = "gene_symbol"
    ) %>%
    mutate(in_utility_table = ifelse(is.na(in_utility_table), FALSE, in_utility_table))

# ==============================================================================
# 8. Validate and write tables
# ==============================================================================
cat("Writing output tables...\n")

table1_required <- c(
    "strategy_id", "strategy_label", "level", "category",
    "prop_of_diagnosed_mean", "prop_of_diagnosed_ui_low",
    "prop_of_diagnosed_ui_high"
)
assert_required_columns(table1, table1_required, "high_impact_outcomes_by_strategy")
assert_no_na(
    table1, c("strategy_id", "strategy_label", "level", "category"),
    "high_impact_outcomes_by_strategy"
)

write_csv_validated(table1, OUTPUT_TABLE1, "high_impact_outcomes_by_strategy")

table2_required <- c("gene_symbol", "in_simulation", "in_utility_table")
assert_required_columns(table2, table2_required, "gene_coverage_audit")
write_csv_validated(table2, OUTPUT_TABLE2, "gene_coverage_audit")

# ==============================================================================
# 9. Figure: Composite clinical impact (three panels)
# ==============================================================================
cat("Generating composite clinical impact figure...\n")

# Strategy display labels
strategy_display <- c(
    "Panel" = "Panel",
    "ES" = "ES",
    "Panel_Reflex_ES" = "Reflex",
    "GS" = "GS"
)

IMPACT_CATS_FIG <- c(
    "high_impact",
    "therapeutic",
    "extrarenal_surveillance",
    "transplant_donor"
)

category_display <- c(
    "high_impact" = "High Impact (any)",
    "therapeutic" = "Therapeutic",
    "extrarenal_surveillance" = "Extrarenal Surveillance",
    "transplant_donor" = "Transplant / Donor"
)

category_display_panel_c <- c(
    "high_impact" = "High Impact",
    "therapeutic" = "Therapeutic",
    "extrarenal_surveillance" = "Extrarenal",
    "transplant_donor" = "Transplant"
)

phenotype_axis_display <- c(
    "cystic" = "Cystic",
    "glomerular" = "Glomerular",
    "tubulointerstitial" = "Tubulo-\ninterstitial",
    "tubulopathies" = "Tubulopathies",
    "CKDu" = "CKDu"
)

plot_theme <- theme_minimal(base_size = 11) +
    theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(color = "gray92"),
        axis.text.x = element_text(size = 9.5, face = "bold", color = "gray20"),
        axis.text.y = element_text(size = 9),
        axis.title.y = element_text(size = 10, margin = margin(r = 8)),
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 8),
        plot.title = element_text(face = "bold", size = 12, margin = margin(b = 3)),
        plot.subtitle = element_text(size = 9, color = "gray30", margin = margin(b = 8)),
        plot.margin = margin(t = 10, r = 12, b = 6, l = 10)
    )

# ---- Panel A: clinical impact by strategy ------------------------------------
panel_a_data <- table1 %>%
    filter(level %in% c("aggregate", "primary")) %>%
    mutate(
        category = factor(
            category,
            levels = IMPACT_CATS_FIG,
            labels = category_display[IMPACT_CATS_FIG]
        ),
        strategy = strategy_display[strategy_label],
        strategy = factor(strategy, levels = c("Panel", "ES", "Reflex", "GS"))
    )

p_panel_a <- ggplot(panel_a_data, aes(
    x = strategy, y = prop_of_diagnosed_mean,
    fill = category
)) +
    geom_col(position = position_dodge(width = 0.76), width = 0.66) +
    geom_errorbar(
        aes(ymin = prop_of_diagnosed_ui_low, ymax = prop_of_diagnosed_ui_high),
        position = position_dodge(width = 0.76),
        width = 0.2, linewidth = 0.35, color = "gray20"
    ) +
    scale_fill_viridis_d(
        option = "D", begin = 0.12, end = 0.88,
        direction = 1, name = "Impact category"
    ) +
    scale_y_continuous(
        labels = scales::percent_format(accuracy = 1),
        expand = expansion(mult = c(0, 0.05))
    ) +
    labs(
        x = NULL,
        y = "Proportion of diagnosed probands",
        title = "A. Clinical impact by strategy"
    ) +
    plot_theme +
    theme(legend.position = "bottom")

# ---- Panel B data prep: impact by phenotype ----------------------------------
counts_pheno_long <- all_gene_counts %>%
    left_join(
        utility %>% select(gene_symbol, all_of(IMPACT_CATS_FIG)),
        by = "gene_symbol"
    ) %>%
    pivot_longer(
        cols = all_of(IMPACT_CATS_FIG),
        names_to = "category",
        values_to = "cat_flag"
    )

pheno_denom_diagnosed <- all_gene_counts %>%
    group_by(iteration_id, strategy_id, strategy_label, phenotype) %>%
    summarise(n_diag_pheno = sum(n_diagnosed), .groups = "drop")

pheno_cat_iter <- counts_pheno_long %>%
    group_by(iteration_id, strategy_id, strategy_label, phenotype, category) %>%
    summarise(
        n_diag_cat = sum(n_diagnosed[cat_flag == 1], na.rm = TRUE),
        .groups = "drop"
    )

iter_pheno_diagnosed <- pheno_cat_iter %>%
    left_join(
        pheno_denom_diagnosed,
        by = c("iteration_id", "strategy_id", "strategy_label", "phenotype")
    ) %>%
    filter(n_diag_pheno > 0) %>%
    mutate(prop = n_diag_cat / n_diag_pheno)

agg_pheno_diagnosed <- iter_pheno_diagnosed %>%
    group_by(strategy_id, strategy_label, phenotype, category) %>%
    summarise(
        prop_of_diagnosed_mean = mean(prop),
        prop_of_diagnosed_ui_low = quantile(prop, 0.025),
        prop_of_diagnosed_ui_high = quantile(prop, 0.975),
        n_iter_nonzero = dplyr::n(),
        .groups = "drop"
    )

pheno_denom_all <- pheno_iter_raw %>%
    filter(strategy_id == 1) %>%
    select(iteration_id,
        phenotype = phenotype_category,
        n_probands_pheno = n_probands
    )

iter_pheno_all <- pheno_cat_iter %>%
    left_join(pheno_denom_all, by = c("iteration_id", "phenotype")) %>%
    filter(!is.na(n_probands_pheno), n_probands_pheno > 0) %>%
    mutate(prop = n_diag_cat / n_probands_pheno)

agg_pheno_all <- iter_pheno_all %>%
    group_by(strategy_id, strategy_label, phenotype, category) %>%
    summarise(
        prop_of_all_probands_mean = mean(prop),
        prop_of_all_probands_ui_low = quantile(prop, 0.025),
        prop_of_all_probands_ui_high = quantile(prop, 0.975),
        n_iter_nonzero = dplyr::n(),
        .groups = "drop"
    )

impact_pheno_out <- agg_pheno_diagnosed %>%
    left_join(
        agg_pheno_all %>%
            select(
                strategy_id, strategy_label, phenotype, category,
                prop_of_all_probands_mean,
                prop_of_all_probands_ui_low,
                prop_of_all_probands_ui_high
            ),
        by = c("strategy_id", "strategy_label", "phenotype", "category")
    ) %>%
    arrange(strategy_id, phenotype, category)
write_csv_validated(
    impact_pheno_out, OUTPUT_IMPACT_PHENO_CSV,
    "impact_by_phenotype"
)

# ---- Panel B: impact by phenotype among all probands -------------------------
panel_b_data <- agg_pheno_all %>%
    filter(strategy_label == FIG2_REF_STRATEGY) %>%
    mutate(
        phenotype = factor(
            phenotype_axis_display[phenotype],
            levels = phenotype_axis_display[PHENOTYPE_LEVELS]
        ),
        category = factor(
            category,
            levels = IMPACT_CATS_FIG,
            labels = category_display[IMPACT_CATS_FIG]
        )
    )

p_panel_b <- ggplot(panel_b_data, aes(
    x = phenotype,
    y = prop_of_all_probands_mean,
    fill = category
)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.72) +
    geom_errorbar(
        aes(ymin = prop_of_all_probands_ui_low, ymax = prop_of_all_probands_ui_high),
        position = position_dodge(width = 0.8),
        width = 0.22, linewidth = 0.35, color = "gray25"
    ) +
    scale_fill_viridis_d(
        option = "D", begin = 0.12, end = 0.88,
        direction = 1, name = "Impact category"
    ) +
    scale_y_continuous(
        labels = scales::percent_format(accuracy = 0.1),
        expand = expansion(mult = c(0, 0.05))
    ) +
    labs(
        x = NULL,
        y = "Proportion of all probands in phenotype",
        title = "B. Impact by phenotype among all probands",
        subtitle = "Reference strategy: Reflex (Panel to ES)."
    ) +
    plot_theme +
    theme(legend.position = "none")

# ---- Panel C data prep: top genes by impact category -------------------------
TOP_N_GENES <- 10

modal_phenotype <- all_gene_counts %>%
    group_by(gene_symbol, phenotype) %>%
    summarise(total = sum(n_diagnosed), .groups = "drop_last") %>%
    slice_max(total, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    transmute(gene_symbol, modal_phenotype = phenotype)

ref_counts <- all_gene_counts %>% filter(strategy_label == FIG2_REF_STRATEGY)
n_iter_fig2 <- length(unique(ref_counts$iteration_id))
n_probands_fig2 <- ref_counts$n_probands[1]

gene_full <- ref_counts %>%
    group_by(gene_symbol) %>%
    summarise(total_n = sum(n_diagnosed), .groups = "drop") %>%
    mutate(mean_prop = total_n / (n_iter_fig2 * n_probands_fig2)) %>%
    left_join(
        utility %>% select(gene_symbol, all_of(IMPACT_CATS_FIG)),
        by = "gene_symbol"
    ) %>%
    left_join(modal_phenotype, by = "gene_symbol")

top_gene_table <- gene_full %>%
    filter(high_impact == 1) %>%
    arrange(desc(mean_prop)) %>%
    slice_head(n = TOP_N_GENES) %>%
    mutate(gene_label = paste0(
        gene_symbol, " (",
        scales::percent(mean_prop, accuracy = 0.1), ")"
    ))

top_genes <- top_gene_table$gene_symbol
gene_levels_plot <- rev(top_gene_table$gene_symbol)
gene_labels_plot <- setNames(
    top_gene_table$gene_label,
    top_gene_table$gene_symbol
)[gene_levels_plot]

fig2_dots <- gene_full %>%
    filter(gene_symbol %in% top_genes) %>%
    pivot_longer(
        cols = all_of(IMPACT_CATS_FIG),
        names_to = "impact_category",
        values_to = "flag"
    ) %>%
    filter(flag == 1) %>%
    mutate(
        gene_symbol = factor(gene_symbol, levels = gene_levels_plot),
        impact_category = factor(
            impact_category,
            levels = IMPACT_CATS_FIG,
            labels = category_display_panel_c[IMPACT_CATS_FIG]
        )
    )

# ---- Supplementary table -----------------------------------------------------
gene_comp_out <- gene_full %>%
    pivot_longer(
        cols = all_of(IMPACT_CATS_FIG),
        names_to = "impact_category",
        values_to = "flag"
    ) %>%
    filter(flag == 1) %>%
    arrange(impact_category, desc(mean_prop)) %>%
    group_by(impact_category) %>%
    mutate(rank = row_number()) %>%
    ungroup() %>%
    mutate(in_top_n = gene_symbol %in% top_genes) %>%
    select(
        impact_category, rank, gene_symbol, modal_phenotype,
        mean_prop, in_top_n
    )
write_csv_validated(
    gene_comp_out, OUTPUT_GENE_COMP_CSV,
    "gene_composition_by_impact"
)

size_breaks <- unique(signif(
    quantile(top_gene_table$mean_prop, probs = c(0.25, 0.5, 0.75, 1)),
    3
))
size_breaks <- size_breaks[size_breaks > 0]
if (length(size_breaks) == 0) {
    size_breaks <- 0.001
}

# ---- Panel C: gene composition (top 10) --------------------------------------
p_panel_c <- ggplot(fig2_dots, aes(x = impact_category, y = gene_symbol)) +
    geom_point(aes(size = mean_prop, color = impact_category), alpha = 0.9) +
    scale_x_discrete(position = "top", expand = expansion(add = c(0.4, 0.4))) +
    scale_y_discrete(
        labels = gene_labels_plot,
        expand = expansion(add = c(0.5, 0.5))
    ) +
    scale_size_area(
        max_size = 13,
        breaks = size_breaks,
        labels = scales::percent_format(accuracy = 0.1),
        name = "Mean cohort proportion diagnosed"
    ) +
    scale_color_viridis_d(
        option = "D", begin = 0.12, end = 0.88,
        direction = 1, guide = "none"
    ) +
    labs(
        x = NULL,
        y = NULL,
        title = "C. Gene composition by impact category (Top 10)",
        subtitle = paste0(
            "Reference strategy: Reflex (Panel to ES). ",
            "Gene labels show High Impact cohort contribution."
        )
    ) +
    theme_minimal(base_size = 11) +
    theme(
        panel.grid.major = element_line(color = "gray94"),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(size = 9, face = "bold", color = "gray15"),
        axis.text.x.top = element_text(size = 9.5, face = "bold", color = "gray15"),
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 8),
        legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 12, margin = margin(b = 3)),
        plot.subtitle = element_text(size = 9, color = "gray30", margin = margin(b = 8)),
        plot.margin = margin(t = 10, r = 12, b = 6, l = 10)
    )

# ---- Composite figure export -------------------------------------------------
p_composite <- p_panel_a / p_panel_b / p_panel_c +
    plot_layout(heights = c(1.05, 1.05, 1.25)) +
    plot_annotation(
        caption = paste0(
            "Error bars: 95% UI across 1,000 simulation iterations. ",
            "Impact categories are non-exclusive and do not sum to 100%. ",
            "Panel B denominator is all probands within phenotype."
        )
    )

ggsave(paste0(FIG_COMPOSITE_BASE, ".png"), p_composite,
    width = 9.5, height = 14, dpi = 300, bg = "white"
)
ggsave(paste0(FIG_COMPOSITE_BASE, ".svg"), p_composite,
    width = 9.5, height = 14, bg = "white"
)

# ==============================================================================
# 10. Summary
# ==============================================================================
cat("\n=== CLINICAL IMPACT ANALYSIS COMPLETE ===\n")

cat("\nTable 1: High-impact outcomes by strategy (among diagnosed probands)\n")
hi_summary <- table1 %>%
    filter(category == "high_impact") %>%
    select(
        strategy_label, prop_of_diagnosed_mean, prop_of_diagnosed_ui_low,
        prop_of_diagnosed_ui_high, cost_per_high_impact_dx_mean
    )
for (i in seq_len(nrow(hi_summary))) {
    r <- hi_summary[i, ]
    cat(sprintf(
        "  %s: %.1f%% [%.1f%%, %.1f%%] | Cost/HI dx: $%.0f\n",
        r$strategy_label,
        r$prop_of_diagnosed_mean * 100, r$prop_of_diagnosed_ui_low * 100,
        r$prop_of_diagnosed_ui_high * 100,
        r$cost_per_high_impact_dx_mean
    ))
}

cat("\nTable 2: Gene coverage audit -", nrow(table2), "genes\n")
cat("  In simulation:", sum(table2$in_simulation), "\n")
cat("  In utility table:", sum(table2$in_utility_table), "\n")
cat("  In both:", sum(table2$in_simulation & table2$in_utility_table), "\n")

cat("\nOutputs:\n")
cat("  ", OUTPUT_TABLE1, "\n")
cat("  ", OUTPUT_TABLE2, "\n")
cat("  ", OUTPUT_IMPACT_PHENO_CSV, "\n")
cat("  ", OUTPUT_GENE_COMP_CSV, "\n")
cat("  ", paste0(FIG_COMPOSITE_BASE, ".png"), "\n")
