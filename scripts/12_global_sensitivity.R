# scripts/12_global_sensitivity.R
# Purpose: PRCC global sensitivity analysis for three strategy comparisons, plus composite figure.
# Author: Zachary Sentell

library(dplyr)
library(tidyr)
library(sensitivity)
library(readr)
library(ggplot2)
library(patchwork)

# ==============================================================================
# IO Paths
# ==============================================================================
INPUT_BASE_OUTCOMES   <- "outputs/results/base_case/iteration_level/strategy_iteration_outcomes.csv"
INPUT_PARAM_TRACE     <- "outputs/results/base_case/iteration_level/iteration_parameter_trace.csv"
INPUT_GS_SCENARIO_ITER <- "outputs/results/scenario_analysis/gs_uplift/gs_scenario_iteration_outcomes.csv"
OUTPUT_DIR            <- "outputs/results/uncertainty_sensitivity/prcc"

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

if (!file.exists(INPUT_BASE_OUTCOMES))    stop("Iteration outcomes not found: ",  INPUT_BASE_OUTCOMES)
if (!file.exists(INPUT_PARAM_TRACE))      stop("Parameter trace not found: ",     INPUT_PARAM_TRACE)
if (!file.exists(INPUT_GS_SCENARIO_ITER)) stop("GS scenario iteration outcomes not found: ", INPUT_GS_SCENARIO_ITER)

# ==============================================================================
# 1. Load and prepare data
# ==============================================================================
outcomes_df <- read_csv(INPUT_BASE_OUTCOMES, show_col_types = FALSE)
params_df   <- read_csv(INPUT_PARAM_TRACE, show_col_types = FALSE)
gs_iter_df  <- read_csv(INPUT_GS_SCENARIO_ITER, show_col_types = FALSE)

params_wide <- params_df %>%
    pivot_wider(names_from = parameter_name, values_from = parameter_value)

base_wide <- outcomes_df %>%
    select(iteration_id, strategy_label, total_cost_per_proband_cad, diagnoses_per_proband) %>%
    pivot_wider(names_from = strategy_label, values_from = c(total_cost_per_proband_cad, diagnoses_per_proband))

gs_scenario_only <- gs_iter_df %>%
    filter(grepl("^uplift_", scenario_label)) %>%
    transmute(
        iteration_id,
        total_cost_per_proband_cad_GS_scenario = total_cost_per_proband_cad,
        diagnoses_per_proband_GS_scenario = diagnoses_per_proband
    )

comparison_frame <- base_wide %>%
    left_join(gs_scenario_only, by = "iteration_id")

comparisons <- list(
    list(
        focal_cost_col = "total_cost_per_proband_cad_Panel_Reflex_ES",
        comp_cost_col  = "total_cost_per_proband_cad_Panel",
        focal_diag_col = "diagnoses_per_proband_Panel_Reflex_ES",
        comp_diag_col  = "diagnoses_per_proband_Panel",
        label          = "reflex_vs_panel",
        title_focal    = "Panel-first with exome reflex",
        title_comp     = "Phenotype-directed Panel"
    ),
    list(
        focal_cost_col = "total_cost_per_proband_cad_ES",
        comp_cost_col  = "total_cost_per_proband_cad_Panel",
        focal_diag_col = "diagnoses_per_proband_ES",
        comp_diag_col  = "diagnoses_per_proband_Panel",
        label          = "es_vs_panel",
        title_focal    = "Exome-first",
        title_comp     = "Phenotype-directed Panel"
    ),
    list(
        focal_cost_col = "total_cost_per_proband_cad_GS_scenario",
        comp_cost_col  = "total_cost_per_proband_cad_Panel_Reflex_ES",
        focal_diag_col = "diagnoses_per_proband_GS_scenario",
        comp_diag_col  = "diagnoses_per_proband_Panel_Reflex_ES",
        label          = "genome_vs_reflex",
        title_focal    = "Genome sequencing (+10% yield)",
        title_comp     = "Panel-first with exome reflex"
    )
)

# ==============================================================================
# 2. Run PRCC for each comparison
# ==============================================================================
for (comp in comparisons) {
    req <- c(comp$focal_cost_col, comp$comp_cost_col, comp$focal_diag_col, comp$comp_diag_col)
    missing <- setdiff(req, names(comparison_frame))
    if (length(missing) > 0) {
        warning("Skipping ", comp$label, " due to missing columns: ", paste(missing, collapse = ", "))
        next
    }

    analysis_df <- comparison_frame %>%
        mutate(
            delta_cost_cad = .data[[comp$focal_cost_col]] - .data[[comp$comp_cost_col]],
            delta_diagnoses = .data[[comp$focal_diag_col]] - .data[[comp$comp_diag_col]]
        ) %>%
        select(iteration_id, delta_cost_cad, delta_diagnoses) %>%
        inner_join(params_wide, by = "iteration_id")

    X <- analysis_df %>%
        select(-iteration_id, -delta_cost_cad, -delta_diagnoses) %>%
        select(where(is.numeric))

    names(X) <- make.names(names(X), unique = TRUE)

    col_vars <- apply(X, 2, function(x) var(x, na.rm = TRUE))
    X <- X[, !(is.na(col_vars) | col_vars == 0), drop = FALSE]

    if (ncol(X) < 2) {
        warning("Skipping ", comp$label, ": not enough non-constant parameters")
        next
    }

    nboot <- 500
    cat("Running PRCC for", comp$label, "(nboot=", nboot, ")\n", sep = "")
    pcc_cost <- pcc(X, analysis_df$delta_cost_cad, rank = TRUE, nboot = nboot, conf = 0.95)
    pcc_diag <- pcc(X, analysis_df$delta_diagnoses, rank = TRUE, nboot = nboot, conf = 0.95)

    prcc_cost <- pcc_cost$PRCC
    prcc_diag <- pcc_diag$PRCC

    prcc_table <- data.frame(
        parameter_name = rownames(prcc_cost),
        prcc_delta_cost_cad = prcc_cost$original,
        prcc_delta_cost_cad_ci_low = prcc_cost$`min. c.i.`,
        prcc_delta_cost_cad_ci_high = prcc_cost$`max. c.i.`,
        prcc_delta_diagnoses = prcc_diag$original,
        prcc_delta_diagnoses_ci_low = prcc_diag$`min. c.i.`,
        prcc_delta_diagnoses_ci_high = prcc_diag$`max. c.i.`,
        stringsAsFactors = FALSE
    ) %>%
        arrange(desc(abs(prcc_delta_cost_cad)))

    out_table <- file.path(OUTPUT_DIR, paste0("prcc_results_", comp$label, ".csv"))
    write_csv(prcc_table, out_table)

    cat("Saved PRCC results for", comp$label, "\n")
}

# ==============================================================================
# 3. Composite figure: ranked forest plot (primary) + heatmap (all comparisons)
# ==============================================================================
param_labels <- c(
    "prev_CKDu"               = "Prevalence: CKDu",
    "prev_cystic"             = "Prevalence: Cystic",
    "prev_glomerular"         = "Prevalence: Glomerular",
    "prev_tubulointerstitial" = "Prevalence: Tubulointerstitial",
    "prev_tubulopathies"      = "Prevalence: Tubulopathies",
    "yield_CKDu"              = "Monogenic probability: CKDu",
    "yield_cystic"            = "Monogenic probability: Cystic",
    "yield_glomerular"        = "Monogenic probability: Glomerular",
    "yield_tubulointerstitial"= "Monogenic probability: Tubulointerstitial",
    "yield_tubulopathies"     = "Monogenic probability: Tubulopathies",
    "prob_vus_panel_mean"     = "VUS probability: Panel",
    "prob_vus_es"             = "VUS probability: ES",
    "prob_vus_gs"             = "VUS probability: GS",
    "det_snv_es"              = "SNV/indel detection: ES",
    "det_cnv_es"              = "CNV/SV detection: ES",
    "det_snv_gs"              = "SNV/indel detection: GS",
    "det_cnv_gs"              = "CNV/SV detection: GS",
    "det_snv_panel_weighted"  = "SNV/indel detection: Panel",
    "det_cnv_panel_weighted"  = "CNV/SV detection: Panel",
    "det_pkd1_es"             = "PKD1 detection: ES",
    "det_pkd1_gs"             = "PKD1 detection: GS",
    "det_muc1_es"             = "MUC1 detection: ES",
    "det_muc1_gs"             = "MUC1 detection: GS",
    "det_pkd1_panel"          = "PKD1 detection: Cystic panel",
    "det_muc1_panel"          = "MUC1 detection: ADTKD panel",
    "cost_panel"              = "Panel unit cost",
    "cost_es"                 = "Exome unit cost",
    "cost_gs"                 = "Genome unit cost",
    "cost_familial"           = "Familial variant test cost",
    "cost_consult_pre"        = "Pre-test consultation cost",
    "cost_consult_post"       = "Post-test consultation cost",
    "uptake_reflex"           = "Reflex ES uptake probability",
    "uptake_cascade"          = "Cascade uptake probability"
)

get_label <- function(param_name) {
    clean <- gsub("\\.", "_", param_name)
    if (clean %in% names(param_labels)) return(param_labels[[clean]])
    tools::toTitleCase(gsub("_", " ", clean))
}

assign_category <- function(param_name) {
    if (grepl("^prev_",            param_name)) return("Phenotype prevalence")
    if (grepl("^cost_",            param_name)) return("Unit cost")
    if (grepl("det_pkd1|det_muc1", param_name)) return("Difficult-locus detection")
    if (grepl("^det_",             param_name)) return("Standard detection")
    if (grepl("^yield_",           param_name)) return("Monogenic probability")
    if (grepl("^prob_vus",         param_name)) return("VUS probability")
    if (grepl("^uptake_",          param_name)) return("Uptake probability")
    "Other"
}

comparison_meta <- list(
    list(label = "reflex_vs_panel", display = "Reflex vs Panel",  file = "prcc_results_reflex_vs_panel.csv"),
    list(label = "es_vs_panel",     display = "ES vs Panel",      file = "prcc_results_es_vs_panel.csv"),
    list(label = "genome_vs_reflex",display = "GS vs Reflex",     file = "prcc_results_genome_vs_reflex.csv")
)

prcc_results <- list()
for (m in comparison_meta) {
    path <- file.path(OUTPUT_DIR, m$file)
    if (!file.exists(path)) {
        message("PRCC result not found — skipping: ", path)
        next
    }
    df <- read_csv(path, show_col_types = FALSE)
    prcc_results[[m$label]] <- list(data = df, display = m$display)
}

if (length(prcc_results) == 0) stop("No PRCC result files found in ", OUTPUT_DIR)

primary_key <- if ("reflex_vs_panel" %in% names(prcc_results)) "reflex_vs_panel" else names(prcc_results)[1]
primary_df  <- prcc_results[[primary_key]]$data

PRCC_SHOW_THRESH  <- 0.10   # min |PRCC| to colour a heatmap cell (must also be sig)
PRCC_UNION_THRESH <- 0.20   # min |PRCC| for a parameter to enter the union set
N_PRIMARY         <- 6      # top-N from primary comparison to always include

is_sig <- function(lo, hi) !is.na(lo) & !is.na(hi) & sign(lo) == sign(hi)

top_primary <- primary_df %>%
    mutate(max_abs = pmax(abs(prcc_delta_cost_cad), abs(prcc_delta_diagnoses), na.rm = TRUE)) %>%
    arrange(desc(max_abs)) %>%
    slice_head(n = N_PRIMARY) %>%
    pull(parameter_name)

top_scenario <- unique(unlist(lapply(names(prcc_results), function(lbl) {
    df <- prcc_results[[lbl]]$data
    df %>%
        filter(
            (abs(prcc_delta_cost_cad)   >= PRCC_UNION_THRESH &
                 is_sig(prcc_delta_cost_cad_ci_low,   prcc_delta_cost_cad_ci_high)) |
            (abs(prcc_delta_diagnoses)  >= PRCC_UNION_THRESH &
                 is_sig(prcc_delta_diagnoses_ci_low, prcc_delta_diagnoses_ci_high))
        ) %>%
        pull(parameter_name)
})))

param_union <- unique(c(top_primary, top_scenario))

scenario_only <- setdiff(top_scenario, top_primary)
if (length(scenario_only) > 0) {
    scenario_ranks <- vapply(scenario_only, function(p) {
        max(sapply(names(prcc_results)[names(prcc_results) != primary_key], function(lbl) {
            df <- prcc_results[[lbl]]$data
            row <- df[df$parameter_name == p, ]
            if (nrow(row) == 0) return(0)
            max(abs(row$prcc_delta_cost_cad[1]), abs(row$prcc_delta_diagnoses[1]), na.rm = TRUE)
        }), na.rm = TRUE)
    }, numeric(1))
    scenario_only <- scenario_only[order(scenario_ranks, decreasing = TRUE)]
}
param_union <- unique(c(top_primary, scenario_only))

param_display  <- setNames(vapply(param_union, get_label,      character(1)), param_union)
param_category <- setNames(vapply(param_union, assign_category, character(1)), param_union)

primary_df_reflex <- primary_df %>% filter(parameter_name %in% param_union)

rank_order <- primary_df_reflex %>%
    mutate(max_abs = pmax(abs(prcc_delta_cost_cad), abs(prcc_delta_diagnoses), na.rm = TRUE),
           display_label = param_display[parameter_name]) %>%
    arrange(desc(max_abs)) %>%
    pull(display_label)
rank_order_full <- unique(c(rank_order, unname(param_display[scenario_only])))
y_levels <- rev(rank_order_full)

forest_df <- bind_rows(
    primary_df_reflex %>%
        transmute(parameter_name,
                  outcome  = "Incremental cost",
                  prcc     = prcc_delta_cost_cad,
                  ci_low   = prcc_delta_cost_cad_ci_low,
                  ci_high  = prcc_delta_cost_cad_ci_high),
    primary_df_reflex %>%
        transmute(parameter_name,
                  outcome  = "Incremental yield",
                  prcc     = prcc_delta_diagnoses,
                  ci_low   = prcc_delta_diagnoses_ci_low,
                  ci_high  = prcc_delta_diagnoses_ci_high)
) %>%
    mutate(
        display_label = factor(param_display[parameter_name], levels = y_levels),
        category      = param_category[parameter_name],
        outcome       = factor(outcome, levels = c("Incremental cost", "Incremental yield")),
        significant   = is_sig(ci_low, ci_high)
    )

comparison_order <- vapply(comparison_meta, `[[`, character(1), "display")
comparison_order <- comparison_order[comparison_order %in%
                                         sapply(prcc_results, `[[`, "display")]

heatmap_parts <- list()
for (lbl in names(prcc_results)) {
    info <- prcc_results[[lbl]]
    df   <- info$data
    heatmap_parts[[paste0(lbl, "_cost")]] <- df %>%
        transmute(parameter_name,
                  comparison = info$display,
                  outcome    = "Cost",
                  prcc       = prcc_delta_cost_cad,
                  ci_low     = prcc_delta_cost_cad_ci_low,
                  ci_high    = prcc_delta_cost_cad_ci_high)
    heatmap_parts[[paste0(lbl, "_yield")]] <- df %>%
        transmute(parameter_name,
                  comparison = info$display,
                  outcome    = "Yield",
                  prcc       = prcc_delta_diagnoses,
                  ci_low     = prcc_delta_diagnoses_ci_low,
                  ci_high    = prcc_delta_diagnoses_ci_high)
}

heatmap_df <- bind_rows(heatmap_parts) %>%
    filter(parameter_name %in% param_union) %>%
    complete(parameter_name = param_union,
             comparison     = comparison_order,
             outcome        = c("Cost", "Yield")) %>%
    mutate(
        significant = is_sig(ci_low, ci_high),
        material    = !is.na(prcc) & abs(prcc) >= PRCC_SHOW_THRESH & significant,
        fill_value  = ifelse(material, prcc, NA_real_),
        text_label  = ifelse(material,
                             ifelse(prcc >= 0,
                                    sprintf("+%.2f", prcc),
                                    sprintf("%.2f",  prcc)),
                             ""),
        text_color  = ifelse(material & abs(prcc) >= 0.60, "white", "grey15"),
        display_label = factor(
            vapply(parameter_name,
                   function(p) param_display[[gsub("\\.", "_", p)]],
                   character(1)),
            levels = y_levels
        ),
        comparison = factor(comparison, levels = comparison_order),
        outcome    = factor(outcome, levels = c("Cost", "Yield"))
    )

forest_plot <- ggplot(
    forest_df,
    aes(x = prcc, y = display_label,
        color = category, shape = outcome, group = outcome)
) +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = "grey45", linewidth = 0.5) +
    geom_errorbar(
        aes(xmin = ci_low, xmax = ci_high, alpha = significant),
        width       = 0,
        linewidth   = 0.65,
        orientation = "y",
        position    = position_dodge(width = 0.55)
    ) +
    geom_point(
        aes(alpha = significant),
        size     = 3.0,
        position = position_dodge(width = 0.55)
    ) +
    scale_alpha_manual(values = c("TRUE" = 1.0, "FALSE" = 0.22), guide = "none") +
    scale_shape_manual(
        values = c("Incremental cost" = 16, "Incremental yield" = 18),
        name   = "Outcome"
    ) +
    scale_color_viridis_d(
        option = "D",
        end    = 0.88,
        name   = "Parameter category"
    ) +
    scale_x_continuous(
        limits = c(-1, 1),
        breaks = c(-1, -0.5, 0, 0.5, 1),
        labels = c("-1", "-0.5", "0", "+0.5", "+1"),
        expand = expansion(mult = 0.02)
    ) +
    labs(
        title    = "A  Reflex vs Panel (primary comparison)",
        subtitle = paste0(
            "Parameters significant in any comparison, ranked by importance to primary decision; ",
            "faded = not significant here"
        ),
        x = "Partial rank correlation coefficient (PRCC)",
        y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
        panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank(),
        panel.grid.major.x = element_line(color = "grey92", linewidth = 0.35),
        axis.text.y        = element_text(size = 9.5, color = "black"),
        axis.text.x        = element_text(size = 9,   color = "black"),
        axis.title.x       = element_text(size = 10,  face = "bold",
                                          margin = margin(t = 6)),
        plot.title         = element_text(size = 12, face = "bold"),
        plot.subtitle      = element_text(size = 8.5, color = "grey40",
                                          margin = margin(b = 4)),
        legend.position    = "bottom",
        legend.box         = "vertical",
        legend.text        = element_text(size = 9),
        legend.title       = element_text(size = 9, face = "bold"),
        legend.key.size    = unit(0.45, "cm"),
        plot.margin        = margin(t = 6, r = 6, b = 6, l = 6)
    ) +
    guides(
        color = guide_legend(order = 1, nrow = 3,
                             override.aes = list(shape = 16, size = 3,
                                                 alpha = 1)),
        shape = guide_legend(order = 2,
                             override.aes = list(color = "grey30",
                                                 size = 3, alpha = 1))
    )

heatmap_plot <- ggplot(
    heatmap_df,
    aes(x = outcome, y = display_label, fill = fill_value)
) +
    geom_tile(color = "white", linewidth = 0.9) +
    geom_text(
        aes(label = text_label, color = text_color),
        size     = 2.7,
        fontface = "plain",
        show.legend = FALSE
    ) +
    scale_fill_gradient2(
        low      = "#2166AC",
        mid      = "white",
        high     = "#B2182B",
        midpoint = 0,
        limits   = c(-1, 1),
        na.value = "grey94",
        name     = "PRCC",
        breaks   = c(-1, -0.5, 0, 0.5, 1),
        labels   = c("-1", "-0.5", "0", "+0.5", "+1"),
        guide    = guide_colorbar(
            barwidth       = 0.7,
            barheight      = 5.5,
            title.position = "top",
            title.hjust    = 0.5
        )
    ) +
    scale_color_identity(guide = "none") +
    facet_wrap(~ comparison, nrow = 1, strip.position = "top") +
    labs(
        title    = "B  All comparisons",
        subtitle = "Coloured cells: |PRCC| >= 0.10, 95% CI excludes zero. Grey: not material or not significant.",
        x        = NULL,
        y        = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
        axis.text.y      = element_blank(),
        axis.ticks.y     = element_blank(),
        axis.text.x      = element_text(size = 9, color = "black"),
        strip.text       = element_text(size = 10, face = "bold"),
        strip.background = element_rect(fill = "grey92", color = NA),
        panel.grid       = element_blank(),
        plot.title       = element_text(size = 12, face = "bold"),
        plot.subtitle    = element_text(size = 8.5, color = "grey40",
                                        margin = margin(b = 4)),
        legend.position  = "right",
        legend.title     = element_text(size = 9.5, face = "bold"),
        legend.text      = element_text(size = 9),
        plot.margin      = margin(t = 6, r = 10, b = 6, l = 2)
    )

composite <- forest_plot + heatmap_plot +
    plot_layout(widths = c(3, 2))

out_png <- file.path(OUTPUT_DIR, "prcc_tornado_composite.png")
out_pdf <- file.path(OUTPUT_DIR, "prcc_tornado_composite.pdf")

ggsave(out_png, composite, width = 16, height = 9.0, dpi = 300, bg = "white")
ggsave(out_pdf, composite, width = 16, height = 9.0, bg = "white")

cat("Saved composite PRCC figure:\n  ", out_png, "\n  ", out_pdf, "\n")
cat("Parameter set (", length(param_union), " parameters):\n",
    paste0("  ", param_union, " -> ", param_display, "\n"), sep = "")
cat("\nGlobal sensitivity analysis complete.\n")
