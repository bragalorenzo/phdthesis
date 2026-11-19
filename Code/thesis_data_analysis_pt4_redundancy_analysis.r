######################################################################
#--- PhD Thesis Data Analysis ---#--- Part 3: Redundancy Analysis ---#
######################################################################

# Global and pollen-zone models + exploratory biplot
#
# Reproducibility rules:
#   - Original spreadsheets/CSVs are never overwritten.
#   - Samples are joined by ID, not assumed to share an implicit row order.
#   - Pollen-zone assignments come only from the depth intervals in "zonation".
#   - Response exclusions, climate predictors, and plotting choices are explicit.
#   - No permutation tests are run: temporal dependence has not been assessed.
#   - For the thesis: check temporal structure, response residuals and predictor
#     distributions before drawing inferential conclusions.
###############################################################################

library(tidyverse)
library(vegan)

set.seed(777)
# Setup #####
# Select the record: "inf1" or "cpz2"

record <- "inf1"

workdir <- paste0(
    "G:/My Drive/Work/PhD Thesis/Data Analysis/",
    toupper(record),
    " Data"
)

input_file <- paste0(record, "_compendium.xlsx")

setwd(workdir)

# Palette and pollen-zone order
if (tolower(record) == "inf1") {
    retro_palette <- c(
        "MOD4b" = "#ce4a4a",
        "MOD4a" = "#eaaf41",
        "MOD3"  = "#3cb68e",
        "MOD2"  = "#6688c3",
        "MOD1b" = "#b25da6",
        "MOD1a" = "#48a56a"
    )

    zone_levels <- c(
        "MOD4b", "MOD4a", "MOD3",
        "MOD2", "MOD1b", "MOD1a"
    )
} else if (tolower(record) == "cpz2") {
    retro_palette <- c(
        "MAL4"  = "#ce4a4a",
        "MAL3"  = "#eaaf41",
        "MAL2"  = "#6688c3",
        "MAL1b" = "#b25da6",
        "MAL1a" = "#48a56a"
    )

    zone_levels <- c(
        "MAL4", "MAL3", "MAL2",
        "MAL1b", "MAL1a"
    )
} else {
    stop("Record not found.")
}

# RDA settings
settings <- list(
    input_file = input_file,
    counts_sheet = "counts_aggregated",
    output_dir = "RDA_results_ecogroups_11-11",
    climate_dir = file.path(
        workdir,
        paste0(toupper(record), " PaleoView")
    ),
    climate_files = c(
        precipitation = "annual_mean_precipitation.csv",
        temperature = "annual_mean_temperature.csv",
        precipitation_seasonality = "annual_precipitation_seasonality.csv",
        relative_humidity = "annual_relative_humidity.csv",
        summer_precipitation = "precipitation_summer_jja.csv"
    ),

    # Climate variables used in all models
    climate_vars = c(
        "temperature",
        "summer_precipitation"
    ),

    # Response matrix exclusions
    exclude_response = c("nd"),

    # Climate data processing
    climate_step_years = 10,
    relative_humidity_zero_as_na = TRUE,
    standardize_predictors = TRUE,

    # Independent within-zone exploratory RDAs
    run_zone_models = TRUE,

    # Plot settings (do not affect fitted models)
    make_plot = TRUE,
    scaling = 2,
    n_group_labels = 5,
    groups_to_label = NULL,
    environment_arrow_mult = 1.2,
    group_arrow_mult = 10,
    plot_title = "RDA – Main ecological groups",
    environment_labels = c(
        temperature = "Mean annual temperature",
        precipitation = "Annual precipitation",
        precipitation_seasonality = "Precipitation seasonality",
        relative_humidity = "Relative humidity",
        summer_precipitation = "Summer precipitation"
    ),
    group_labels = c(
        limno = "Limnophytes",
        hydro_root = "Rooting hydrophytes",
        hydro_float = "Floating hydrophytes",
        wet_meadow = "Wet meadow",
        meadow = "Grassland",
        mediterranean_tree = "Mediterranean trees/shrubs",
        api = "API",
        conifer = "Conifers",
        oakwood = "Mixed oakwood",
        potentially_anthropic = "Potential anthropogenic indicators"
    ),
    retro_colors = retro_palette
)

# Data #####

counts <- rio::import(
    settings$input_file,
    sheet = settings$counts_sheet
)

metadata <- rio::import(
    settings$input_file,
    sheet = "sample_metadata"
)

zones <- rio::import(
    settings$input_file,
    sheet = "zonation"
)

taxon_lookup <- rio::import(
    settings$input_file,
    sheet = "taxon_lookup"
)

# 2. SMALL REUSABLE FUNCTIONS -------------------------------------------------

# Fail with a readable message rather than silently producing a wrong model.
assert <- function(condition, message) {
    if (!isTRUE(condition)) stop(message, call. = FALSE)
}

# Find the pollen zone containing one sample depth.
# Intervals are [topDepth, bottomDepth), except the deepest endpoint is included.
zone_from_depth <- function(depth, zone_def) {
    idx <- which(
        depth >= zone_def$topDepth &
            (depth < zone_def$bottomDepth |
                (depth == max(zone_def$bottomDepth) &
                    depth == zone_def$bottomDepth))
    )
    assert(
        length(idx) == 1L,
        paste0("Depth ", depth, " cm matches ", length(idx), " zones")
    )
    zone_def$pollenZone[idx]
}

# Fit one Hellinger-RDA. The response Y must ALREADY be Hellinger-transformed.
# X contains ONLY the selected numerical climatic predictors.
fit_climate_rda <- function(Y, X, standardize = TRUE) {
    assert(nrow(Y) == nrow(X), "Y/X row-count mismatch")
    assert(nrow(Y) > ncol(X) + 1L, "Too few samples for this many predictors")
    assert(all(vapply(X, is.numeric, logical(1))), "Non-numeric predictor")
    assert(all(is.finite(as.matrix(X))), "Missing/non-finite predictor")

    x_sd <- vapply(X, sd, numeric(1))
    assert(all(is.finite(x_sd) & x_sd > 0), "Constant predictor found")

    X_model <- if (standardize) {
        as.data.frame(scale(X))
    } else {
        as.data.frame(X)
    }

    assert(
        qr(cbind(Intercept = 1, as.matrix(X_model)))$rank == ncol(X_model) + 1L,
        "Predictor matrix is rank deficient"
    )

    # Formula RHS (.) is safe here: X_model contains ONLY selected predictors.
    model <- vegan::rda(Y ~ ., data = X_model)
    list(model = model, predictors = X_model)
}

# Summary for one fitted model; VIF is exported separately.
model_metrics <- function(model, n_samples, n_predictors) {
    r2 <- vegan::RsquareAdj(model)
    tibble(
        n_samples = n_samples,
        n_predictors = n_predictors,
        residual_df = n_samples - n_predictors - 1L,
        R2 = r2$r.squared,
        R2_adj = r2$adj.r.squared,
        max_VIF = max(vegan::vif.cca(model))
    )
}

# Plot scores from the exact model supplied: no refitting in the plot function.
plot_rda <- function(model, sample_zones, zone_order, settings) {
    assert(
        length(model$CCA$eig) >= 2L,
        "The biplot needs at least two constrained RDA axes"
    )

    get_scores <- function(display, id) {
        vegan::scores(
            model,
            display = display, choices = 1:2,
            scaling = settings$scaling
        ) |>
            as.data.frame() |>
            tibble::rownames_to_column(id)
    }

    sites <- get_scores("sites", "sample") |>
        left_join(sample_zones |> select(sample, pollenZone), by = "sample") |>
        mutate(pollenZone = factor(pollenZone, levels = zone_order))
    assert(!anyNA(sites$pollenZone), "Plot: samples without pollen zones")

    environment <- get_scores("bp", "variable") |>
        mutate(
            x_plot = RDA1 * settings$environment_arrow_mult,
            y_plot = RDA2 * settings$environment_arrow_mult,
            label = coalesce(
                unname(settings$environment_labels[variable]),
                str_replace_all(variable, "_", " ")
            )
        )

    groups <- get_scores("species", "group") |>
        mutate(
            length_2d = sqrt(RDA1^2 + RDA2^2),
            x_plot = RDA1 * settings$group_arrow_mult,
            y_plot = RDA2 * settings$group_arrow_mult,
            label = dplyr::coalesce(
                unname(settings$group_labels[group]),
                str_replace_all(group, "_", " ")
            )
        )

    # Select response labels by plotted-arrow length OR explicit group names.
    if (is.null(settings$groups_to_label)) {
        groups <- groups |>
            slice_max(length_2d, n = settings$n_group_labels, with_ties = FALSE)
    } else {
        missing_groups <- setdiff(settings$groups_to_label, groups$group)
        assert(
            length(missing_groups) == 0L,
            paste("Unknown plotted groups:", paste(missing_groups, collapse = ", "))
        )
        groups <- groups |> filter(group %in% settings$groups_to_label)
    }

    eig <- model$CCA$eig[1:2]
    axis_labels <- sprintf(
        "RDA%d (%.1f%% of total variation)",
        1:2, 100 * eig / model$tot.chi
    )
    r2 <- vegan::RsquareAdj(model)
    subtitle <- sprintf(
        "R2 = %.3f; adjusted R2 = %.3f",
        r2$r.squared, r2$adj.r.squared
    )
    labels_data <- bind_rows(
        groups |>
            transmute(
                x = x_plot,
                y = y_plot,
                label,
                face = "plain"
            ),
        environment |>
            transmute(
                x = x_plot,
                y = y_plot,
                label,
                face = "bold"
            ),
        sites |>
            transmute(
                x = RDA1,
                y = RDA2,
                label = "",
                face = "plain"
            )
    )

    ggplot() +
        geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) +
        geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.3) +
        # Dashed: ecological response-group scores (visually enlarged).
        geom_segment(
            data = groups,
            aes(x = 0, y = 0, xend = x_plot, yend = y_plot),
            arrow = grid::arrow(length = grid::unit(0.13, "cm")),
            linewidth = 0.45, linetype = "dashed", colour = "grey45"
        ) +
        # Solid: environmental biplot scores (separate graphical multiplier).
        geom_segment(
            data = environment,
            aes(x = 0, y = 0, xend = x_plot, yend = y_plot),
            arrow = grid::arrow(length = grid::unit(0.22, "cm")),
            linewidth = 0.75, colour = "grey20"
        ) +
        geom_point(
            data = sites, aes(RDA1, RDA2, fill = pollenZone),
            shape = 21, size = 2.8, stroke = 0.35,
            colour = "black", alpha = 0.9
        ) +
        ggrepel::geom_label_repel(
            data = labels_data,
            aes(
                x = x,
                y = y,
                label = label,
                fontface = face
            ),
            seed = 777,
            size = 3,
            colour = "grey20",
            fill = "white",
            label.size = 0.15,
            label.padding = grid::unit(0.12, "lines"),
            label.r = grid::unit(0.1, "lines"),
            box.padding = 0.4,
            point.padding = 0.15,
            force = 2,
            force_pull = 0.5,
            max.iter = 20000,
            max.time = 3,
            max.overlaps = Inf,
            min.segment.length = 0,
            segment.colour = "grey60",
            segment.size = 0.25
        ) +
        scale_fill_manual(values = retro_palette, drop = FALSE) +
        coord_equal() +
        labs(
            title = settings$plot_title, subtitle = subtitle,
            x = axis_labels[1], y = axis_labels[2], fill = "Pollen Zone",
            caption = paste0(
                "Hellinger distance RDA; graphical arrow multipliers: ",
                "climate x", settings$environment_arrow_mult,
                ", ecological groups x", settings$group_arrow_mult,
                "."
            )
        ) +
        theme_classic() +
        theme(
            legend.position = "right",
            plot.caption = element_text(size = 8, colour = "grey40"),
            plot.caption.position = "plot"
        )
}

# 3. IMPORT AND VALIDATE INPUTS ------------------------------------------------

assert(file.exists(settings$input_file), "Missing input spreadsheet")
climate_paths <- file.path(settings$climate_dir, settings$climate_files)
assert(
    all(file.exists(climate_paths)),
    paste("Missing PaleoView files:", paste(climate_paths[!file.exists(climate_paths)],
        collapse = ", "
    ))
)
assert(
    all(settings$climate_vars %in% names(settings$climate_files)),
    "Selected climate variable is absent from climate_files"
)
assert(!anyDuplicated(settings$climate_vars), "Repeated climate variable")

counts_table <- rio::import(settings$input_file, sheet = settings$counts_sheet)
metadata <- rio::import(settings$input_file, sheet = "sample_metadata") |>
    filter(excludeSample == "no")
zones <- rio::import(settings$input_file, sheet = "zonation")

assert(
    all(c("sample", "ageCalBP", "depthCm") %in% names(metadata)),
    "Metadata lacks sample, ageCalBP, or depthCm"
)
assert("sample" %in% names(counts_table), "Counts sheet lacks sample ID column")
assert(
    all(c("topDepth", "bottomDepth", "pollenZone") %in% names(zones)),
    "Zonation sheet has unexpected column names"
)

metadata <- metadata |>
    mutate(sample = as.character(sample)) |>
    arrange(desc(ageCalBP))
counts_table <- counts_table |>
    mutate(sample = as.character(sample))

assert(
    is.numeric(metadata$ageCalBP) && is.numeric(metadata$depthCm),
    "ageCalBP and depthCm must be numeric"
)
assert(
    !anyNA(metadata$ageCalBP) && !anyNA(metadata$depthCm),
    "Missing age/depth in selected samples"
)
assert(
    !anyDuplicated(metadata$sample) && !anyDuplicated(counts_table$sample),
    "Duplicated sample ID"
)
assert(
    setequal(metadata$sample, counts_table$sample),
    "Counts/metadata have different included sample IDs"
)

# Zonation is defined by the depth ranges, NOT a legacy metadata$pollenZone.
assert(
    is.numeric(zones$topDepth) && is.numeric(zones$bottomDepth),
    "Zonation depth limits must be numeric"
)
zone_def <- zones |>
    mutate(pollenZone = as.character(pollenZone)) |>
    arrange(topDepth)
assert(all(zone_def$topDepth < zone_def$bottomDepth), "Invalid depth intervals")
assert(!anyDuplicated(zone_def$pollenZone), "Duplicate zone names")
assert(
    isTRUE(all.equal(
        head(zone_def$bottomDepth, -1),
        tail(zone_def$topDepth, -1)
    )),
    "Zonation intervals are not contiguous"
)
zone_order <- rev(zone_def$pollenZone) # Oldest/deepest -> youngest/shallowest

sample_zones <- metadata |>
    transmute(
        sample, depthCm, ageCalBP,
        pollenZone = vapply(depthCm, zone_from_depth,
            character(1),
            zone_def = zone_def
        )
    )
zone_sizes <- sample_zones |>
    mutate(pollenZone = factor(pollenZone, levels = zone_order)) |>
    count(pollenZone, .drop = FALSE)
assert(all(zone_sizes$n > 0L), "At least one pollen zone has no samples")
cat("\nSAMPLES BY ZONE (oldest -> youngest)\n")
print(zone_sizes, n = "inf")

# 4. PREPARE RESPONSE MATRIX --------------------------------------------------

# Sample ID is a KEY, never a taxon.
counts_all <- counts_table |>
    tibble::column_to_rownames("sample")
assert(
    all(vapply(counts_all, is.numeric, logical(1))),
    "Some count columns are not numeric"
)
counts_all <- as.matrix(counts_all)[sample_zones$sample, , drop = FALSE]
assert(
    all(is.finite(counts_all)) && all(counts_all >= 0),
    "Counts must be finite and nonnegative"
)
assert(all(rowSums(counts_all) > 0), "Some samples have zero pollen sum")

missing_exclusions <- setdiff(settings$exclude_response, colnames(counts_all))
assert(
    length(missing_exclusions) == 0L,
    paste("Response columns not found:", paste(missing_exclusions, collapse = ", "))
)

# Exclusion changes the proportional denominator: record how much was retained.
counts_rda <- counts_all[
    , !colnames(counts_all) %in% settings$exclude_response,
    drop = FALSE
]
assert(
    ncol(counts_rda) > 0L && all(rowSums(counts_rda) > 0),
    "No positive counts remain in at least one sample"
)

response_coverage <- sample_zones |>
    mutate(
        total_counts = rowSums(counts_all),
        retained_counts = rowSums(counts_rda),
        retained_fraction = retained_counts / total_counts
    )
cat("\nFRACTION OF COUNTS RETAINED IN THE RESPONSE\n")
print(summary(response_coverage$retained_fraction))

# Recalculate Hellinger AFTER choosing/aggregating the response taxa/groups.
counts_hellinger <- vegan::decostand(counts_rda, method = "hellinger")
assert(
    identical(rownames(counts_hellinger), sample_zones$sample),
    "Counts are not in metadata order"
)

# 5. PALEOVIEW: IMPORT, QC, AGE MATCHING --------------------------------------

climate_raw <- setNames(
    lapply(climate_paths, rio::import),
    names(settings$climate_files)
)

climate_series <- purrr::imap(climate_raw, function(df, variable) {
    assert(
        all(c("Year (BP)", "Area Mean") %in% names(df)),
        paste("Unexpected columns:", variable)
    )
    assert(
        is.numeric(df[["Year (BP)"]]) && is.numeric(df[["Area Mean"]]),
        paste("Non-numeric PaleoView data:", variable)
    )
    assert(
        !anyNA(df[["Year (BP)"]]) && !anyDuplicated(df[["Year (BP)"]]),
        paste("Missing or duplicate PaleoView years:", variable)
    )
    tibble(AgeBP = df[["Year (BP)"]], !!variable := df[["Area Mean"]]) |>
        arrange(desc(AgeBP))
})

# Require identical time grids. A join must not conceal a shifted series.
reference_ages <- climate_series[[1]]$AgeBP
assert(
    all(diff(reference_ages) == -settings$climate_step_years),
    "PaleoView time step differs from the configured value"
)
for (series in climate_series) {
    assert(
        identical(series$AgeBP, reference_ages),
        "PaleoView CSVs do not have exactly matching age grids"
    )
}
climate_data <- purrr::reduce(climate_series, dplyr::left_join, by = "AgeBP")

# Diagnostic only: report zero and missing values, including consecutive runs.
longest_zeros <- function(x) {
    runs <- rle(!is.na(x) & x == 0)
    if (any(runs$values)) max(runs$lengths[runs$values]) else 0L
}
climate_qc_summary <- climate_data |>
    pivot_longer(-AgeBP, names_to = "variable", values_to = "value") |>
    group_by(variable) |>
    summarise(
        n_zeros = sum(value == 0, na.rm = TRUE),
        n_missing = sum(is.na(value)),
        longest_zero_run = longest_zeros(value),
        .groups = "drop"
    )
cat("\nPALEOVIEW ZERO/MISSING-VALUE CHECK\n")
print(climate_qc_summary, n = Inf)

# Preserve climate_data untouched. Only RH zeros are flagged as missing in QC.
climate_data_qc <- climate_data
if (settings$relative_humidity_zero_as_na &&
    "relative_humidity" %in% names(climate_data_qc)) {
    climate_data_qc <- climate_data_qc |>
        mutate(relative_humidity = na_if(relative_humidity, 0))
}

assert(
    all(sample_zones$ageCalBP >= min(climate_data_qc$AgeBP) &
        sample_zones$ageCalBP <= max(climate_data_qc$AgeBP)),
    "One or more pollen ages fall outside the PaleoView age range"
)

# Nearest-neighbour in time; ties select the older age because rows descend.
nearest_idx <- vapply(sample_zones$ageCalBP, function(age) {
    which.min(abs(climate_data_qc$AgeBP - age))
}, integer(1))

climate_aligned <- sample_zones |>
    mutate(
        AgeBP_climate = climate_data_qc$AgeBP[nearest_idx],
        age_difference = ageCalBP - AgeBP_climate
    ) |>
    bind_cols(climate_data_qc[nearest_idx, -1, drop = FALSE])

assert(
    all(abs(climate_aligned$age_difference) <=
        settings$climate_step_years / 2 + 1e-8),
    "Pollen-climate age matching exceeds half a climate step"
)

# X contains only the chosen climate variables, in pollen-sample row order.
env_rda <- climate_aligned |>
    select(all_of(settings$climate_vars))
assert(
    all(is.finite(as.matrix(env_rda))),
    "Selected predictors contain NA/Inf (check RH gap or file QC)"
)
assert(nrow(counts_hellinger) == nrow(env_rda), "Y/X row-count mismatch")
cat("\nAGE-MATCHING ERROR (years; signed)\n")
print(summary(climate_aligned$age_difference))
cat("\nSELECTED PREDICTOR CORRELATIONS\n")
print(round(cor(env_rda), 3))


# 6. GLOBAL RDA ----------------------------------------------------------------

fit_global <- fit_climate_rda(
    Y = counts_hellinger, X = env_rda,
    standardize = settings$standardize_predictors
)
rda_climate <- fit_global$model

global_summary <- model_metrics(
    rda_climate, nrow(counts_hellinger), ncol(env_rda)
)
global_vif <- tibble(
    variable = names(vegan::vif.cca(rda_climate)),
    VIF = unname(vegan::vif.cca(rda_climate))
)
cat("\nGLOBAL RDA\n")
print(global_summary)
print(global_vif)
print(summary(rda_climate))

# Overall significance of the RDA
anova(rda_climate, permutations = 9999)

# Significance of individual constrained axes
anova(rda_climate, by = "axis", permutations = 9999)

# Significance of individual climatic predictors
anova(rda_climate, by = "margin", permutations = 9999)

# 7. SEPARATE WITHIN-ZONE RDAs (OPTIONAL) -------------------------------------

rda_by_zone <- list()
rda_zone_summary <- tibble()
vif_by_zone <- tibble()

if (settings$run_zone_models) {
    # Same selected climate predictors in every zone. No pooled/centroid bins.
    rda_by_zone <- setNames(vector("list", length(zone_order)), zone_order)
    for (z in zone_order) {
        idx <- which(sample_zones$pollenZone == z)
        zone_fit <- fit_climate_rda(
            Y = counts_hellinger[idx, , drop = FALSE],
            X = env_rda[idx, , drop = FALSE],
            standardize = settings$standardize_predictors
        )
        rda_by_zone[[z]] <- zone_fit$model
    }

    rda_zone_summary <- purrr::imap_dfr(rda_by_zone, function(model, z) {
        n <- sum(sample_zones$pollenZone == z)
        model_metrics(model, n, length(settings$climate_vars)) |>
            mutate(pollenZone = z, .before = 1)
    })
    vif_by_zone <- purrr::imap_dfr(rda_by_zone, function(model, z) {
        vif <- vegan::vif.cca(model)
        tibble(pollenZone = z, variable = names(vif), VIF = unname(vif))
    })
    cat("\nZONE-WISE RDA SUMMARY\n")
    print(rda_zone_summary, n = Inf)
    cat("\nZONE-WISE VIF\n")
    print(vif_by_zone, n = Inf)
}

# 8. EXPLORATORY BIPLOT -------------------------------------------------------

if (settings$make_plot) {
    p_rda <- plot_rda(rda_climate, sample_zones, zone_order, settings)
    print(p_rda)
}


# 9. EXPORT + REPRODUCIBILITY -------------------------------------------------

dir.create(settings$output_dir, recursive = TRUE, showWarnings = FALSE)

rio::export(global_summary, file.path(settings$output_dir, "global_rda_summary.csv"))
rio::export(global_vif, file.path(settings$output_dir, "global_rda_vif.csv"))
rio::export(zone_sizes, file.path(settings$output_dir, "zone_sample_sizes.csv"))
rio::export(
    response_coverage,
    file.path(settings$output_dir, "response_retained_fraction.csv")
)
rio::export(
    climate_qc_summary,
    file.path(settings$output_dir, "paleoview_quality_check.csv")
)

if (settings$run_zone_models) {
    rio::export(
        rda_zone_summary,
        file.path(settings$output_dir, "zone_rda_summary.csv")
    )
    rio::export(
        vif_by_zone,
        file.path(settings$output_dir, "zone_rda_vif.csv")
    )
}
if (settings$make_plot) {
    ggsave(file.path(settings$output_dir, "rda_biplot.pdf"),
        p_rda,
        width = 8.2, height = 6.4
    )
    ggsave(file.path(settings$output_dir, "rda_biplot.png"),
        p_rda,
        width = 8.2, height = 6.4, dpi = 300
    )
}

# Preserve precise model inputs, settings and fitted model objects.
saveRDS(
    list(
        settings = settings, sample_zones = sample_zones,
        counts_hellinger = counts_hellinger, env_rda = env_rda,
        global_model = rda_climate, zone_models = rda_by_zone
    ),
    file.path(settings$output_dir, "rda_analysis.rds")
)

input_paths <- c(settings$input_file, unname(climate_paths))
input_manifest <- tibble(
    path = input_paths,
    md5 = unname(tools::md5sum(input_paths))
)
rio::export(input_manifest, file.path(settings$output_dir, "input_manifest.csv"))
writeLines(
    capture.output(sessionInfo()),
    file.path(settings$output_dir, "R_session_info.txt")
)

cat("\nFinished. Outputs in: ", settings$output_dir, "\n", sep = "")

# DCA #####

counts_dca <- rio::import(
    "cpz2_compendium.xlsx",
    sheet = "counts_aggregated"
) |>
    select(-nd) |>
    column_to_rownames("sample") |>
    as.matrix()

counts_dca <- counts_dca[
    , colSums(counts_dca) > 0,
    drop = FALSE
]

dca <- vegan::decorana(counts_dca)

dca
