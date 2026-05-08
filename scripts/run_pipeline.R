# scripts/run_pipeline.R
# Purpose: Canonical pipeline runner. Sources scripts in order and verifies required outputs.
# Author: Zachary Sentell

if (!requireNamespace("config", quietly = TRUE)) {
    stop("Package 'config' is required. Install it with install.packages('config').")
}
library(readr)

cfg <- config::get()
scripts <- cfg$pipeline$scripts
required_outputs <- cfg$pipeline$required_outputs

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
