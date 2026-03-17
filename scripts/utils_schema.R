# utils_schema.R
# Purpose: Shared schema and output validation utilities.

library(readr)

assert_required_columns <- function(df, required_cols, table_name, exact = FALSE) {
    missing <- setdiff(required_cols, names(df))
    if (length(missing) > 0) {
        stop(sprintf(
            "%s missing required columns: %s",
            table_name,
            paste(missing, collapse = ", ")
        ))
    }

    if (exact) {
        extras <- setdiff(names(df), required_cols)
        if (length(extras) > 0) {
            stop(sprintf(
                "%s has unexpected columns: %s",
                table_name,
                paste(extras, collapse = ", ")
            ))
        }
    }

    invisible(TRUE)
}

assert_no_na <- function(df, cols, table_name) {
    offenders <- cols[sapply(cols, function(col) any(is.na(df[[col]])))]
    if (length(offenders) > 0) {
        stop(sprintf(
            "%s has NA values in required non-missing columns: %s",
            table_name,
            paste(offenders, collapse = ", ")
        ))
    }

    invisible(TRUE)
}

assert_probabilities_sum_to_one <- function(df, group_cols, prob_col, table_name, tol = 1e-6) {
    grouped <- aggregate(df[[prob_col]], by = df[group_cols], FUN = sum)
    bad <- abs(grouped$x - 1) > tol
    if (any(bad)) {
        stop(sprintf(
            "%s probability sums are not 1 for %d group(s)",
            table_name,
            sum(bad)
        ))
    }

    invisible(TRUE)
}

assert_frontier_first_incremental_na <- function(df,
                                                 status_col = "dominance_status",
                                                 cost_col = "total_cost_per_proband_cad_mean",
                                                 incr_cols = c(
                                                     "incremental_cost_vs_prev_cad_mean",
                                                     "incremental_diagnoses_vs_prev_mean",
                                                     "icpd_cad_per_additional_diagnosis_mean"
                                                 ),
                                                 table_name = "incremental_table") {
    frontier <- df[df[[status_col]] == "non_dominated", , drop = FALSE]
    if (nrow(frontier) == 0) {
        stop(sprintf("%s has no non_dominated strategies", table_name))
    }

    first_idx <- which.min(frontier[[cost_col]])
    first_frontier <- frontier[first_idx, , drop = FALSE]

    for (col in incr_cols) {
        if (!is.na(first_frontier[[col]][1])) {
            stop(sprintf(
                "%s expected NA for first frontier strategy column '%s'",
                table_name,
                col
            ))
        }
    }

    invisible(TRUE)
}

write_csv_validated <- function(df, path, table_name = basename(path)) {
    if (!dir.exists(dirname(path))) {
        dir.create(dirname(path), recursive = TRUE)
    }
    write_csv(df, path)
    message(sprintf("Wrote %s: %s", table_name, path))
    invisible(path)
}
