# scripts/00_prisma_figure.R
# Generate PRISMA 2020 flow diagram for synthetic cohort creation
# Input: data/raw/pubmed_query_results.xlsx (PRISMA_log sheet)
# Output: outputs/figures/prisma/prisma_2020_flowdiagram.{png,svg}
#         outputs/figures/prisma/caption.txt

suppressPackageStartupMessages({
    library(readxl)
    library(dplyr)
    library(stringr)
    library(PRISMA2020)
    # Check for high-res save capability
    has_rsvg <- requireNamespace("rsvg", quietly = TRUE)
})

# Ensure output directory exists
output_dir <- "outputs/figures/prisma"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

message("Reading screening log...")
log_path <- "data/raw/pubmed_query_results.xlsx"
if (!file.exists(log_path)) stop("Input file not found: ", log_path)

log <- read_excel(log_path, sheet = "PRISMA_log")

# --- 1. Compute Box Counts from Log ---

# Database results (identification)
database_results <- nrow(log)
duplicates <- sum(log$duplicate == 1, na.rm = TRUE)

# Initial screening (Title/Abstract)
records_screened <- database_results - duplicates

# TA exclusion count and breakdown by reason.
# Record 32312792 (CAKUT-only paper) was tagged "Not adult predominant cohort" at the TA
# stage but is correctly classified as a single-phenotype cohort; it is remapped here.
ta_excludes <- log %>%
    filter(ta_decision == "exclude") %>%
    mutate(reason_display = case_when(
        fulltext_exclusion_reason == "Not renal genetic testing study" ~
            "Out of scope",
        fulltext_exclusion_reason %in% c("Single phenotype", "Not adult predominant cohort") ~
            "Single-phenotype cohort",
        TRUE ~ fulltext_exclusion_reason
    ))

ta_reason_counts <- ta_excludes %>%
    count(reason_display, name = "n") %>%
    arrange(desc(n))

records_excluded <- sum(log$ta_decision == "exclude", na.rm = TRUE)
records_excluded_str <- paste0(
    ta_reason_counts$reason_display, ", ", ta_reason_counts$n, collapse = "; "
)

# Full-text retrieval
dbr_sought_reports <- sum(log$fulltext_sought == 1, na.rm = TRUE)
dbr_notretrieved_reports <- sum(log$fulltext_sought == 1 & log$fulltext_retrieved == 0, na.rm = TRUE)
dbr_assessed <- dbr_sought_reports - dbr_notretrieved_reports

# Full-text exclusions and reasons
ft_excludes <- log %>%
    filter(fulltext_sought == 1, fulltext_retrieved == 1, fulltext_decision == "exclude")

# Harmonize exclusion reasons to match PRISMA display labels
ft_excludes_clean <- ft_excludes %>%
    mutate(reason_clean = case_when(
        fulltext_exclusion_reason == "Not adult predominant cohort" ~ "Not adult-predominant cohort",
        fulltext_exclusion_reason == "Review commentary protocol"   ~ "Review, commentary, or protocol",
        fulltext_exclusion_reason == "Does not meet N proband criteria" ~ "Sample size <150 probands",
        fulltext_exclusion_reason == "Single phenotype" ~ "Single-phenotype cohort",
        TRUE ~ fulltext_exclusion_reason
    ))

# Construct reason summary string
reason_counts <- ft_excludes_clean %>%
    count(reason_clean, name = "n") %>%
    arrange(desc(n))

dbr_excluded_str <- if (nrow(reason_counts) == 0) {
    "0"
} else {
    # Format: "Reason, n; Reason, n"
    paste0(reason_counts$reason_clean, ", ", reason_counts$n, collapse = "; ")
}

# Included studies/reports
included <- log %>%
    filter(fulltext_decision == "include")

new_reports <- nrow(included)
new_studies <- n_distinct(included$study_id, na.rm = TRUE)

# --- 2. Populate PRISMA Template Correctly ---

message("Populating PRISMA template...")

# Load empty template from package
template_path <- system.file("extdata", "PRISMA.csv", package = "PRISMA2020")
prisma_data <- read.csv(template_path, stringsAsFactors = FALSE)

# Helper function to update 'n' by matching 'boxtext' pattern
update_n <- function(df, pattern, value) {
    idx <- which(grepl(pattern, df$boxtext, ignore.case = TRUE))
    if (length(idx) == 0) {
        warning("Pattern not found: ", pattern)
        return(df)
    }
    if (length(idx) > 1) idx <- idx[1]

    df$n[idx] <- value
    return(df)
}

# Helper function to update 'boxtext' (label) by matching 'boxtext' pattern
update_label <- function(df, pattern, new_label) {
    idx <- which(grepl(pattern, df$boxtext, ignore.case = TRUE))
    if (length(idx) == 0) {
        return(df)
    }
    if (length(idx) > 1) idx <- idx[1]
    df$boxtext[idx] <- new_label
    return(df)
}

# Update counts based on specific text mapping
# To hide boxes (Registers, various exclusions), we set n = NA
prisma_filled <- prisma_data %>%
    update_n("^Databases$", database_results) %>%
    update_n("^Registers$", NA) %>% # HIDE Registers
    update_n("Duplicate records", duplicates) %>%
    update_n("ineligible by automation", NA) %>% # HIDE Automation
    update_n("removed for other reasons", NA) %>% # HIDE Other
    update_n("Records screened", records_screened) %>%
    update_n("Records excluded", records_excluded) %>%
    update_n("Reports sought for retrieval", dbr_sought_reports) %>%
    update_n("Reports not retrieved", dbr_notretrieved_reports) %>%
    update_n("Reports assessed for eligibility", dbr_assessed) %>%
    update_n("Reports excluded", dbr_excluded_str) %>%
    update_n("New studies included in review", new_studies) %>%
    update_n("Reports of new included studies", new_reports)

# Update Labels for clarity
# Embed TA exclusion breakdown in the "Records excluded" box label so the figure
# shows reason counts rather than just the total (the n field carries the total).
ta_label_str <- paste0(
    "Records excluded\n   ",
    paste0(ta_reason_counts$reason_display, ": ", ta_reason_counts$n, collapse = "\n   ")
)
prisma_filled <- prisma_filled %>%
    update_label("Records excluded", ta_label_str) %>%
    update_label("New studies included in review", "Studies (cohorts) included for parameterization") # %>%
# update_label("Reports of new included studies", paste0("reports (n = ", new_reports, ")")) # Not strictly used as label but usually n is appended.
# Actually, 'Reports of new included studies' is the text for the second line in the final box.
# If we want a combined box text like "Studies ... (n=6); reports (n=6)", PRISMA2020 might split them.
# Let's inspect the template behavior. box10 has 2 lines.
# We can change the top line to: "Studies (cohorts) included for parameterization"

# Process data for plotting
flow_data <- PRISMA_data(prisma_filled)

# --- 3. Generate High-Quality Plots ---

message("Generating plots with enhanced quality...")

# Generate Plot Object with larger font
plot_obj <- PRISMA_flowdiagram(
    flow_data,
    interactive = FALSE,
    previous = FALSE,
    other = FALSE,
    fontsize = 14
)

# Function to save using high-res if possible
save_high_res <- function(plot_obj, filename_base) {
    svg_path <- file.path(output_dir, paste0(filename_base, ".svg"))
    png_path <- file.path(output_dir, paste0(filename_base, ".png"))

    if (requireNamespace("DiagrammeRsvg", quietly = TRUE) && requireNamespace("rsvg", quietly = TRUE)) {
        svg_code <- DiagrammeRsvg::export_svg(plot_obj)
        writeLines(svg_code, svg_path)
        message("Saved SVG: ", svg_path)
        rsvg::rsvg_png(charToRaw(svg_code), file = png_path, width = 2400)
        message("Saved High-Res PNG (300dpi class): ", png_path)
    } else {
        message("Falling back to standard PRISMA_save")
        PRISMA_save(plot_obj, filename = png_path, filetype = "PNG", overwrite = TRUE)
        PRISMA_save(plot_obj, filename = svg_path, filetype = "SVG", overwrite = TRUE)
    }
}

tryCatch(
    {
        save_high_res(plot_obj, "prisma_2020_flowdiagram")
    },
    error = function(e) {
        message("Error in custom save: ", e$message)
        PRISMA_save(plot_obj, filename = file.path(output_dir, "prisma_2020_flowdiagram.png"), filetype = "PNG", overwrite = TRUE)
    }
)

# --- 4. Write Caption ---

caption_text <- "Title: Study selection for synthetic cohort parameterization

Caption:
PRISMA 2020 flow diagram depicting the study selection process for synthetic cohort parameterization.
Database search (PubMed/MEDLINE; title-only; 2020/01/01 to 2025/09/12) identified 179 records.
Title/abstract screening excluded 166 records: 160 out of scope (not a renal genetic testing cohort study) and 6 single-phenotype cohorts.
Full-text review excluded 7 reports: 3 not adult-predominant, 3 reviews or protocols, 1 sample size <150 probands.
Six cohorts were retained and provided phenotype prevalence, diagnostic yields, and gene/variant architecture parameters for the cost-effectiveness model.
'Registers' and 'Automation tools' were not used and are omitted for clarity."

writeLines(caption_text, file.path(output_dir, "caption.txt"))
message("Saved caption to: ", file.path(output_dir, "caption.txt"))
message("Done.")
