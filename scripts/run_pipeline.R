# run_pipeline.R
# Purpose: Canonical pipeline runner for the renal genetics CEA workflow.

if (!requireNamespace("config", quietly = TRUE)) {
    stop("Package 'config' is required. Install it with install.packages('config').")
}
library(readr)

cfg <- config::get()
scripts <- cfg$pipeline$scripts

run_script <- function(script_path) {
    if (!file.exists(script_path)) stop("Script not found: ", script_path)
    cat("\n=== Running", script_path, "===\n")
    status <- system2("Rscript", script_path)
    if (!identical(status, 0L)) {
        stop("Script failed: ", script_path, " (exit code ", status, ")")
    }
}

for (script_path in scripts) {
    run_script(script_path)
}

required_outputs <- c(
    "outputs/results/base_case/iteration_level/strategy_iteration_outcomes.csv",
    "outputs/results/base_case/summary_tables/base_case_outcomes_by_strategy.csv",
    "outputs/results/base_case/cost_composition/base_case_cost_components_by_strategy.csv",
    "outputs/results/base_case/cost_composition/mean_cost_composition_by_strategy_base_case.png",
    "outputs/results/incremental_analysis/base_case/incremental_table_base_case.csv",
    "outputs/results/incremental_analysis/base_case/efficiency_frontier_base_case.png",
    "outputs/results/supplement/cystic_subgroup/cystic_outcomes_by_strategy.csv",
    "outputs/results/supplement/cystic_subgroup/cystic_cost_components_by_strategy.csv",
    "outputs/results/supplement/cystic_subgroup/incremental_table_cystic.csv",
    "outputs/results/supplement/cystic_subgroup/mean_cost_composition_by_strategy_cystic.png",
    "outputs/results/supplement/cystic_subgroup/efficiency_frontier_cystic.png",
    "outputs/results/base_case/vus_burden/outcome_proportions_by_strategy.csv",
    "outputs/results/base_case/vus_burden/outcome_proportions_stacked.png",
    "outputs/results/uncertainty_sensitivity/probabilistic_sensitivity_analysis/psa_decision_robustness.csv",
    "outputs/results/uncertainty_sensitivity/probabilistic_sensitivity_analysis/psa_cost_diagnosis_scatter.png",
    "outputs/results/uncertainty_sensitivity/global_sensitivity_analysis_prcc/prcc_results_reflex_vs_panel.csv",
    "outputs/results/uncertainty_sensitivity/global_sensitivity_analysis_prcc/prcc_results_es_vs_panel.csv",
    "outputs/results/uncertainty_sensitivity/global_sensitivity_analysis_prcc/prcc_results_genome_vs_reflex.csv",
    "outputs/results/uncertainty_sensitivity/willingness_to_pay/ceac_table.csv",
    "outputs/results/uncertainty_sensitivity/willingness_to_pay/ceac_plot.png",
    "outputs/results/scenario_analysis/gs_uplift/gs_baseline_iteration_outcomes.csv",
    "outputs/results/scenario_analysis/gs_uplift/gs_scenario_iteration_outcomes.csv",
    "outputs/results/scenario_analysis/gs_uplift/gs_scenario_summary.csv",
    "outputs/results/scenario_analysis/gs_uplift/gs_scenario_incremental.csv",
    "outputs/results/scenario_analysis/gs_uplift/gs_scenario_comparison.png"
)

artifact_manifest <- data.frame(
    path = required_outputs,
    exists = file.exists(required_outputs),
    stringsAsFactors = FALSE
)

artifact_manifest$file_size_bytes <- ifelse(
    artifact_manifest$exists,
    file.info(artifact_manifest$path)$size,
    NA
)

artifact_manifest$md5 <- NA_character_
existing <- artifact_manifest$path[artifact_manifest$exists]
if (length(existing) > 0) {
    artifact_manifest$md5[artifact_manifest$exists] <- as.character(tools::md5sum(existing))
}

manifest_path <- "outputs/results/pipeline_artifact_manifest.csv"
if (!dir.exists(dirname(manifest_path))) dir.create(dirname(manifest_path), recursive = TRUE)
write_csv(artifact_manifest, manifest_path)
cat("\nPipeline manifest written to", manifest_path, "\n")

missing <- artifact_manifest$path[!artifact_manifest$exists]
if (length(missing) > 0) {
    stop("Pipeline completed with missing required outputs:\n", paste(missing, collapse = "\n"))
}

cat("\nPipeline run complete. All required outputs were produced.\n")
