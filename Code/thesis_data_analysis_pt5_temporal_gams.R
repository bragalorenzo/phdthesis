################################################################
#--- PhD Thesis Data Analysis ---#--- Part 5: Temporal GAMs ---#
################################################################

# Savestates #####
save.image("inf1_roc_workspace.RData")
load("cpz2_cta_workspace.RData")

# Libraries #####

library(tidyverse)
library(vegan)
library(mgcv)

# Setup #####

record <- "inf1"

workdir <- paste0(
    "G:/My Drive/Work/PhD Thesis/Data Analysis/",
    toupper(record),
    " Data"
)

input_file <- paste0(
    record, "_compendium.xlsx"
)

setwd(workdir)

## Palette
if (tolower(record) == "inf1") {
    # INF1 palette
    retro_palette <- c(
        "MOD4b" = "#ce4a4a",
        "MOD4a" = "#eaaf41",
        "MOD3"  = "#3cb68e",
        "MOD2"  = "#6688c3",
        "MOD1b" = "#b25da6",
        "MOD1a" = "#48a56a"
    )
    # Chronological order of pollen zones: recent -> old
    zone_levels <- c("MOD4b", "MOD4a", "MOD3", "MOD2", "MOD1b", "MOD1a")
} else if (tolower(record) == "cpz2") {
    # CPZ2 palette
    retro_palette <- c(
        "MAL4" = "#ce4a4a",
        "MAL3" = "#eaaf41",
        "MAL2" = "#6688c3",
        "MAL1b" = "#b25da6",
        "MAL1a" = "#48a56a"
    )
    # Chronological order of pollen zones: recent -> old
    zone_levels <- c(
        "MAL4",
        "MAL3",
        "MAL2",
        "MAL1b",
        "MAL1a"
    )
} else {
    message("Record not found.")
}

# Data #####

counts_harmonised <- rio::import(input_file, sheet = "counts_harmonised")

metadata <- rio::import(input_file, sheet = "sample_metadata")

zones <- rio::import(input_file, sheet = "zonation")

taxon_lookup <- rio::import(input_file, sheet = "taxon_lookup")

# GAM data preparation #####

set.seed(777)

# Filter and order metadata
metadata_gam <- metadata |>
    filter(excludeSample == "no") |>
    arrange(desc(ageCalBP))

sample_ids <- as.character(metadata_gam$sample)

# Align counts with metadata
gam_counts <- counts_harmonised |>
    select(-"Indeterminable") |> # Remove indeterminable
    mutate(sample = as.character(sample)) |>
    tibble::column_to_rownames("sample")

stopifnot(
    !anyDuplicated(sample_ids),
    !anyDuplicated(rownames(gam_counts)),
    !anyNA(metadata_gam$ageCalBP),
    all(sample_ids %in% rownames(gam_counts))
)

pollen_counts <- as.matrix(
    gam_counts[sample_ids, , drop = FALSE]
)

storage.mode(pollen_counts) <- "double"

stopifnot(
    !anyNA(pollen_counts),
    all(is.finite(pollen_counts)),
    all(pollen_counts >= 0),
    all(rowSums(pollen_counts) > 0)
)

# Hellinger transformation
pollen_hell <- vegan::decostand(
    pollen_counts,
    method = "hellinger"
)

# Global GAM screening #####

# Select nonconstant taxa
taxa_to_fit <- colnames(pollen_hell)[
    apply(pollen_hell, 2, sd) > 0
]

# Fit one GAM per taxon
gam_fits <- setNames(
    lapply(taxa_to_fit, function(tx) {
        dat <- data.frame(
            y = as.numeric(pollen_hell[, tx]),
            ageCalBP = metadata_gam$ageCalBP
        )

        tryCatch(
            mgcv::gam(
                y ~ s(ageCalBP, k = 6),
                data = dat,
                family = gaussian(),
                method = "REML"
            ),
            error = function(e) {
                message("GAM failed [", tx, "]: ", conditionMessage(e))
                NULL
            }
        )
    }),
    taxa_to_fit
)

# Extract model statistics
gam_results <- purrr::imap_dfr(gam_fits, function(model, tx) {
    if (is.null(model)) {
        return(tibble(
            taxon = tx,
            edf = NA_real_,
            p_value = NA_real_,
            R2_adj = NA_real_,
            dev_expl = NA_real_
        ))
    }

    sm <- summary(model)

    tibble(
        taxon = tx,
        edf = unname(sm$s.table[1, "edf"]),
        p_value = unname(sm$s.table[1, "p-value"]),
        R2_adj = sm$r.sq,
        dev_expl = sm$dev.expl
    )
}) |>
    mutate(p_adj = p.adjust(p_value, method = "BH")) |>
    arrange(p_adj)

# Identify significant taxa
significant_taxa <- gam_results |>
    filter(p_adj < 0.05) |>
    pull(taxon)

# Screening summary
gam_summary <- gam_results |>
    summarise(
        n_taxa = n(),
        n_fitted = sum(!is.na(p_value)),
        n_BH_selected = sum(p_adj < 0.05, na.rm = TRUE),
        n_EDF_above_1_5 = sum(
            edf > 1.5 & p_adj < 0.05,
            na.rm = TRUE
        )
    )

print(gam_summary)
print(gam_results |> filter(p_adj < 0.05))

# Taxon frequency and abundance #####

taxon_frequency <- tibble(
    taxon = colnames(pollen_counts),
    n_present = colSums(pollen_counts > 0),
    total_count = colSums(pollen_counts),
    max_count = apply(pollen_counts, 2, max)
) |>
    mutate(
        occurrence_pct = 100 * n_present / nrow(pollen_counts),
        max_count_fraction = max_count / pmax(total_count, 1)
    )

# Frequency statistics for significant taxa
significant_taxa_frequency <- taxon_frequency |>
    filter(taxon %in% significant_taxa) |>
    left_join(
        gam_results |>
            select(taxon, edf, R2_adj, p_adj),
        by = "taxon"
    ) |>
    arrange(p_adj)

print(significant_taxa_frequency, n = Inf)

# GAM predictions #####

# Select taxa for detailed analysis
selected_taxa <- significant_taxa

stopifnot(
    length(selected_taxa) > 0,
    all(selected_taxa %in% names(gam_fits)),
    all(!vapply(gam_fits[selected_taxa], is.null, logical(1)))
)

# Prediction grid
age_grid <- tibble(
    ageCalBP = seq(
        min(metadata_gam$ageCalBP),
        max(metadata_gam$ageCalBP),
        length.out = 300
    )
)

# Observed values
gam_observed <- bind_cols(
    tibble(ageCalBP = metadata_gam$ageCalBP),
    as_tibble(
        pollen_hell[, selected_taxa, drop = FALSE],
        .name_repair = "minimal"
    )
) |>
    pivot_longer(
        -ageCalBP,
        names_to = "taxon",
        values_to = "observed"
    )

# Model predictions
gam_predictions <- purrr::map_dfr(selected_taxa, function(tx) {
    pred <- predict(
        gam_fits[[tx]],
        newdata = age_grid,
        type = "response",
        se.fit = TRUE
    )
    tibble(
        taxon = tx,
        ageCalBP = age_grid$ageCalBP,
        fit = as.numeric(pred$fit),
        lower = as.numeric(pred$fit - 1.96 * pred$se.fit),
        upper = as.numeric(pred$fit + 1.96 * pred$se.fit)
    )
})

# GAM curves #####

plot_ncol <- 3

p_gam <- ggplot() +
    geom_rect(
        data = zones,
        aes(
            xmin = topAge,
            xmax = bottomAge,
            ymin = -Inf,
            ymax = Inf,
            fill = pollenZone
        ),
        alpha = 0.12,
        inherit.aes = FALSE
    ) +
    geom_ribbon(
        data = gam_predictions,
        aes(x = ageCalBP, ymin = lower, ymax = upper),
        fill = "grey80", alpha = 0.6
    ) +
    geom_line(
        data = gam_predictions,
        aes(ageCalBP, fit),
        colour = "#6688c3", linewidth = 0.9
    ) +
    geom_point(
        data = gam_observed,
        aes(ageCalBP, observed),
        size = 1.3, alpha = 0.6
    ) +
    facet_wrap(
        ~ factor(taxon, levels = selected_taxa),
        scales = "free_y",
        ncol = plot_ncol
    ) +
    scale_fill_manual(
        name = "Pollen Zone",
        values = retro_palette
    ) +
    scale_x_reverse() +
    labs(
        x = "Age (cal yr BP)",
        y = "Hellinger-transformed pollen abundance"
    ) +
    theme_classic()

p_gam

# GAM first derivatives #####

set.seed(777)

gam_derivatives <- purrr::map_dfr(selected_taxa, function(tx) {
    d <- gratia::derivatives(
        gam_fits[[tx]],
        select = "s(ageCalBP)",
        data = age_grid,
        type = "central",
        interval = "simultaneous",
        level = 0.95,
        n_sim = 10000
    ) |>
        as_tibble()

    # Compatibility with older gratia versions
    if (!".lower_ci" %in% names(d) && "lower" %in% names(d)) {
        d <- rename(d, .lower_ci = lower, .upper_ci = upper)
    }

    if (!".derivative" %in% names(d) && "derivative" %in% names(d)) {
        d <- rename(d, .derivative = derivative)
    }

    stopifnot(all(
        c("ageCalBP", ".derivative", ".lower_ci", ".upper_ci") %in% names(d)
    ))

    d |>
        mutate(
            taxon = tx,
            change = case_when(
                .lower_ci > 0 ~ "Positive derivative",
                .upper_ci < 0 ~ "Negative derivative",
                TRUE ~ "No supported change"
            )
        )
})

# GAM derivatives plot #####

p_derivatives <- ggplot(
    gam_derivatives,
    aes(ageCalBP, .derivative)
) +
    geom_rect(
        data = zones,
        aes(
            xmin = topAge,
            xmax = bottomAge,
            ymin = -Inf,
            ymax = Inf,
            fill = pollenZone
        ),
        alpha = 0.12,
        inherit.aes = FALSE
    ) +
    geom_ribbon(
        aes(ymin = .lower_ci, ymax = .upper_ci),
        fill = "grey80", alpha = 0.6
    ) +
    geom_line(
        colour = "#ce4a4a",
        linewidth = 0.9
    ) +
    geom_hline(
        yintercept = 0,
        linetype = "dashed"
    ) +
    facet_wrap(
        ~ factor(taxon, levels = selected_taxa),
        scales = "free_y",
        ncol = plot_ncol
    ) +
    scale_fill_manual(
        name = "Pollen Zone",
        values = retro_palette
    ) +
    scale_x_reverse() +
    labs(
        x = "Age (cal yr BP)",
        y = "First derivative (per BP year)"
    ) +
    theme_classic()

p_derivatives

# Significant change intervals #####

change_intervals <- gam_derivatives |>
    arrange(taxon, desc(ageCalBP)) |>
    group_by(taxon) |>
    mutate(
        interval_id = cumsum(
            change != lag(change, default = first(change))
        )
    ) |>
    group_by(taxon, interval_id, change) |>
    summarise(
        oldest_BP = max(ageCalBP),
        youngest_BP = min(ageCalBP),
        duration_years = oldest_BP - youngest_BP,
        .groups = "drop"
    ) |>
    filter(change != "No supported change") |>
    mutate(
        direction_to_present = if_else(
            change == "Positive derivative",
            "Decrease",
            "Increase"
        )
    ) |>
    select(-interval_id) |>
    arrange(taxon, desc(oldest_BP))

print(change_intervals, n = Inf)

# GAM diagnostics #####

set.seed(777)

k_diagnostics <- purrr::map_dfr(selected_taxa, function(tx) {
    kc <- as.data.frame(mgcv::k.check(gam_fits[[tx]]))

    tibble(
        taxon = tx,
        k_prime = kc[["k'"]][1],
        edf = kc[["edf"]][1],
        k_index = kc[["k-index"]][1],
        p_value_k = kc[["p-value"]][1]
    )
})

print(k_diagnostics, n = Inf)

# Diagnostic plots #####

pdf(
    file = paste0(
        "./Figures/",
        record, "_GAM_diagnostics_",
        Sys.Date(), ".pdf"
    ),
    width = 10,
    height = 8,
    onefile = TRUE
)

par(oma = c(0, 0, 3, 0))

for (tx in selected_taxa) {
    mgcv::gam.check(gam_fits[[tx]])

    mtext(
        paste("GAM diagnostics -", tx),
        side = 3,
        outer = TRUE,
        line = 1,
        font = 2,
        cex = 1.2
    )
}

dev.off()

# Export GAM tables #####

rio::export(
    list(
        "summary" = gam_summary,
        "screening" = gam_results,
        "frequency" = taxon_frequency,
        "significant_taxa" = significant_taxa_frequency,
        "predictions" = gam_predictions,
        "derivatives" = gam_derivatives,
        "change_intervals" = change_intervals,
        "k_diagnostics" = k_diagnostics
    ),
    file = paste0(
        "./outputs/",
        record, "_GAM_",
        Sys.Date(),
        ".xlsx"
    )
)

# Export GAM figures #####

plot_height <- 2.8 * ceiling(
    length(selected_taxa) / plot_ncol
) + 1

gam_figures <- list(
    "curves" = p_gam,
    "derivatives" = p_derivatives,
    "curves_derivatives" = p_gam / p_derivatives
)

figure_heights <- c(
    curves = plot_height,
    derivatives = plot_height,
    curves_derivatives = 2 * plot_height
)

for (fig in names(gam_figures)) {
    svg(
        filename = paste0(
            "./Figures/",
            record, "_GAM_", fig, "_",
            Sys.Date(), ".svg"
        ),
        width = 12,
        height = figure_heights[[fig]]
    )

    print(gam_figures[[fig]])
    dev.off()

    ggsave(
        filename = paste0(
            "./Figures/",
            record, "_GAM_", fig, "_",
            Sys.Date(), ".png"
        ),
        plot = gam_figures[[fig]],
        width = 12,
        height = figure_heights[[fig]],
        dpi = 300
    )
}

# Residual autocorrelation (ACF) #####

chronological_order <- order(
    metadata_gam$ageCalBP,
    decreasing = TRUE
)

pdf(
    file = paste0(
        "./Figures/",
        record, "_GAM_residual_ACF_",
        Sys.Date(), ".pdf"
    ),
    width = 12,
    height = 9
)

par(mfrow = c(3, 2))

for (tx in selected_taxa) {
    res <- residuals(
        gam_fits[[tx]],
        type = "pearson"
    )

    acf(
        res[chronological_order],
        main = tx,
        cex.main = 0.9,
        lag.max = 10
    )
}

dev.off()

# Save GAM workspace #####

saveRDS(
    list(
        record = record,
        models = gam_fits,
        screening = gam_results,
        frequency = taxon_frequency,
        predictions = gam_predictions,
        derivatives = gam_derivatives,
        change_intervals = change_intervals,
        k_diagnostics = k_diagnostics
    ),
    file = paste0(
        "./outputs/",
        record, "_GAM_analysis_",
        Sys.Date(), ".rds"
    )
)

# Avena-Triticum residuals check: oscillatory trend? #####

tx <- grep("Avena", names(gam_fits), value = TRUE)

stopifnot(length(tx) == 1)

gam_residuals <- tibble(
    ageCalBP = metadata_gam$ageCalBP,
    residual = as.numeric(
        residuals(gam_fits[[tx]], type = "pearson")
    )
) |>
    arrange(desc(ageCalBP))

p_residuals <- ggplot(
    gam_residuals,
    aes(ageCalBP, residual)
) +
    geom_hline(
        yintercept = 0,
        linetype = "dashed",
        colour = "grey50"
    ) +
    geom_line(
        colour = "#6688c3",
        linewidth = 0.7
    ) +
    geom_point(
        colour = "#6688c3",
        size = 2
    ) +
    scale_x_reverse() +
    labs(
        title = tx,
        x = "Age (cal yr BP)",
        y = "Pearson residuals"
    ) +
    theme_classic()

p_residuals

svg(
    filename = paste0(
        "./Figures/",
        record, "_avena-triticum_residuals_",
        Sys.Date(), ".svg"
    ),
    width = 16,
    height = 9
)
p_residuals
dev.off()

ggsave(
    filename = paste0(
        "./Figures/",
        record, "_avena-triticum_residuals_",
        Sys.Date(), ".png"
    ),
    plot = p_residuals,
    width = 16,
    height = 9,
    dpi = 300
)

rio::export(
    gam_residuals,
    file = paste0(
        "./outputs/",
        record, "_avena-triticum_residuals_",
        Sys.Date(),
        ".xlsx"
    )
)
