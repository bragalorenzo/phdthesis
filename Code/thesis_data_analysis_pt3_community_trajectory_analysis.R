################################################################################
#--- PhD Thesis Data Analysis ---#--- Part 3: Community Trajectory Analysis ---#
################################################################################

# Savestates #####
save.image("inf1_cta_workspace.RData")
load("cpz2_cta_workspace.RData")

# Libraries #####

# remotes::install_github("emf-creaf/ecotraj")

library(tidyverse)
library(vegan)
library(ecotraj)

# Setup #####

record <- "cpz2"

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

## Clean input data

# Sort metadata from oldest to youngest
metadata_traj <- metadata |>
    filter(excludeSample == "no") |>
    arrange(desc(ageCalBP))

# Put sample IDs in row names ONLY for the numeric community matrix

counts_traj <- counts_harmonised |>
    select(-"Indeterminable") |> # Remove Indeterminable pollen
    tibble::column_to_rownames("sample")

counts_traj <- counts_traj[
    metadata_traj$sample, ,
    drop = FALSE
]

## Check for enough community states

n <- nrow(metadata_traj)
stopifnot(n >= 3) # Enough observations for CTA

## Attribute samples to pollen zone

zones <- zones |>
    arrange(topDepth) |>
    mutate(pollenZone = as.character(pollenZone))

zone_breaks <- c(zones$topDepth[1], zones$bottomDepth)
stopifnot(
    all(diff(zone_breaks) > 0),
    all(metadata_traj$depthCm >= min(zone_breaks)),
    all(metadata_traj$depthCm <= max(zone_breaks)),
    !any(metadata_traj$depthCm %in%
        zone_breaks[-c(1, length(zone_breaks))])
)

metadata_traj <- metadata_traj |>
    mutate(
        pollenZone = as.character(cut(
            depthCm,
            breaks = zone_breaks,
            labels = zones$pollenZone,
            include.lowest = TRUE
        )),
    )

stopifnot(!anyNA(metadata_traj$pollenZone))

# Character levels for chronological plotting, oldest -> youngest.
zone_levels <- rev(zones$pollenZone)

zone_summary <- metadata_traj |>
    group_by(pollenZone) |>
    summarise(
        n = n(),
        ageOld = max(ageCalBP),
        ageYoung = min(ageCalBP),
        duration = ageOld - ageYoung,
        depthTop = min(depthCm),
        depthBottom = max(depthCm),
        maxAgeGap = if (n() >= 2) max(-diff(ageCalBP)) else NA_real_,
        .groups = "drop"
    )

# Global trajectory #####

## Hellinger distance matrix

counts_hellinger <- vegan::decostand(counts_traj, method = "hellinger")
dist_hellinger <- vegan::vegdist(counts_hellinger, method = "euclidean")

## Times & surveys for global record

times <- -metadata_traj$ageCalBP # Increasing from past to present
surveys <- seq_len(n)

## Compute trajectory

trajectory <- ecotraj::defineTrajectories(
    d = dist_hellinger,
    sites = rep("record", n),
    surveys = surveys,
    times = times
)

stopifnot(inherits(trajectory, "trajectories"))

## Global trajectory metrics

trajectory_metrics <- ecotraj::trajectoryMetrics(trajectory)
trajectory_lengths <- ecotraj::trajectoryLengths(trajectory)
trajectory_speeds <- ecotraj::trajectorySpeeds(trajectory)
trajectory_angles <- ecotraj::trajectoryAngles(trajectory, stats = FALSE)
trajectory_variation <- ecotraj::trajectoryInternalVariation(
    trajectory,
    relativeContributions = TRUE
)

## Segment and relative speed data frame

segment_cols <- paste0("S", seq_len(n - 1))

segment_df <- data.frame(
    sampleFrom = metadata_traj$sample[-n], # Remove youngest sample
    sampleTo = metadata_traj$sample[-1], # Remove oldest sample
    ageFrom = metadata_traj$ageCalBP[-n],
    ageTo = metadata_traj$ageCalBP[-1],
    zoneFrom = metadata_traj$pollenZone[-n],
    zoneTo = metadata_traj$pollenZone[-1],
    length = as.numeric(unlist(
        trajectory_lengths[1, segment_cols, drop = FALSE],
        use.names = FALSE
    )),
    speed = as.numeric(unlist(
        trajectory_speeds[1, segment_cols, drop = FALSE],
        use.names = FALSE
    ))
) |>
    mutate(
        ageMidpoint = (ageFrom + ageTo) / 2,
        speed100 = speed * 100,
        crossesZone = zoneFrom != zoneTo,
        cumulativeLength = cumsum(length)
    )

## Angles

angle_df <- data.frame(
    sample = metadata_traj$sample[2:(n - 1)],
    ageCalBP = metadata_traj$ageCalBP[2:(n - 1)],
    angle = as.numeric(unlist(
        trajectory_angles[1, , drop = FALSE],
        use.names = FALSE
    ))
)

## Per-state relative contribution to global internal variation (LCBD-like)

rc_cols <- paste0("rc_", seq_len(n))

internal_df <- data.frame(
    sample = metadata_traj$sample,
    ageCalBP = metadata_traj$ageCalBP,
    pollenZone = metadata_traj$pollenZone,
    contribution = as.numeric(unlist(
        trajectory_variation[1, rc_cols, drop = FALSE],
        use.names = FALSE
    ))
) |>
    mutate(
        contributionPct = round(contribution * 100, 2)
    )

# Pollen zone subtrajectories #####

# Each subtrajectory requires at least three samples for all metrics
stopifnot(all(zone_summary$n >= 3))

## Compute subtrajectories

subtrajectories <- ecotraj::defineTrajectories(
    d = dist_hellinger,
    sites = metadata_traj$pollenZone,
    surveys = surveys,
    times = times
)

## Calculate subtrajectory metrics

subtraj_metrics <- ecotraj::trajectoryMetrics(subtrajectories) |>
    rename(pollenZone = trajectory) |>
    select(
        pollenZone, length, mean_speed, mean_angle,
        directionality, internal_ss, internal_variance
    )

## Combine metrics with age, depth and pollen zone

subtraj_summary <- zone_summary |>
    left_join(subtraj_metrics, by = "pollenZone")

# Taxon-wise decomposition of internal variance #####

# Multivariate sample variance is the sum of individual taxon variances

taxon.variance <- function(group_ids) {
    group_ids <- as.character(group_ids)

    purrr::map(unique(group_ids), function(group_name) {
        x <- as.matrix(counts_hellinger[
            group_ids == group_name, ,
            drop = FALSE
        ])

        v <- apply(x, 2, stats::var)
        total <- sum(v)

        data.frame(
            group = group_name,
            pollenTaxon = names(v),
            taxonVariance = as.numeric(v),
            relativeContribution = if (total > 0) v / total else 0
        )
    }) |>
        purrr::list_rbind()
}

## Decompose taxon variance for the global record

taxon_variance_global <- data.frame(
    pollenTaxon = colnames(counts_hellinger),
    taxonVariance = apply(counts_hellinger, 2, stats::var)
) |>
    mutate(
        relativeContribution = taxonVariance / sum(taxonVariance),
        contributionPct = relativeContribution * 100
    ) |>
    arrange(desc(relativeContribution))

# Check consistency with global trajectory internal variance

stopifnot(
    abs(
        sum(taxon_variance_global$taxonVariance) -
            trajectory_variation$internal_variance[1]
    ) < 1e-8
)

## Decompose taxon variance for each pollen zone

taxon_variance_zones <- taxon.variance(
    metadata_traj$pollenZone
) |>
    rename(pollenZone = group)


## Select the 10 taxa contributing most to global variance

top_taxa_global <- taxon_variance_global |>
    slice_max(
        order_by = relativeContribution,
        n = 10,
        with_ties = FALSE
    ) |>
    arrange(relativeContribution) |>
    mutate(
        pollenTaxon = factor(pollenTaxon, levels = pollenTaxon)
    )

# Independent check: sum of individual taxon variances must reproduce the
# internal variance returned by ecotraj for every zone and main zone.
check_zone <- taxon_variance_zones |>
    group_by(pollenZone) |>
    summarise(varSum = sum(taxonVariance), .groups = "drop") |>
    left_join(
        subtraj_summary |> select(pollenZone, internal_variance),
        by = "pollenZone"
    )

check_main <- taxon_variance_global |>
    summarise(varSum = sum(taxonVariance), .groups = "drop")

stopifnot(
    all(abs(check_zone$varSum - check_zone$internal_variance) < 1e-8),
    all(abs(check_main$varSum - check_main$internal_variance) < 1e-8)
)

# Export CTA results #####

rio::export(
    list(
        "global_segments" = segment_df,
        "global_angles" = angle_df,
        "global_contributions" = internal_df,
        "subtraj_metrics" = subtraj_summary,
        "sub_taxon_variance" = taxon_variance_zones,
        "global_taxon_variance" = taxon_variance_global
    ),
    file = paste0(
        "./Outputs/",
        record,
        "_cta_",
        Sys.Date(),
        ".xlsx"
    )
)

# Visualisation #####
# 3 figures total
# Figure 1/3 = PCoA trajectory
# Figure 2A/3 = LCBD-like through time
# Figure 2B/3 = Top 10 taxa contributing to global variance
# Figure 3A/3 = Internal variance by pollen zone
# Figure 3B/3 = Top 10 taxa contributing to within-zone variance

## PCoA

pcoa <- stats::cmdscale(dist_hellinger, k = 2, eig = TRUE)
pcoa_df <- as.data.frame(pcoa$points, .name_repair = "minimal") |>
    stats::setNames(c("PCoA1", "PCoA2")) |>
    mutate(
        sample = metadata_traj$sample,
        ageCalBP = metadata_traj$ageCalBP,
        pollenZone = factor(metadata_traj$pollenZone, levels = zone_levels)
    )

axis_percent <- round(
    100 * pcoa$eig[1:2] / sum(pmax(pcoa$eig, 0)),
    1
)

## Order factors before plotting

pcoa_df <- pcoa_df |>
    mutate(pollenZone = factor(as.character(pollenZone), levels = zone_levels))

internal_df <- internal_df |>
    mutate(pollenZone = factor(as.character(pollenZone), levels = zone_levels))

subtraj_summary <- subtraj_summary |>
    mutate(pollenZone = factor(as.character(pollenZone), levels = zone_levels))

# Global taxon variance ------------------------------------------------------

taxon_variance_global <- data.frame(
    pollenTaxon = colnames(counts_hellinger),
    taxonVariance = apply(counts_hellinger, 2, stats::var)
) |>
    mutate(
        relativeContribution = taxonVariance / sum(taxonVariance),
        contributionPct = relativeContribution * 100
    ) |>
    arrange(desc(relativeContribution))

stopifnot(
    abs(
        sum(taxon_variance_global$taxonVariance) -
            trajectory_variation$internal_variance[1]
    ) < 1e-8
)

top_taxa_global <- taxon_variance_global |>
    slice_max(
        order_by = relativeContribution,
        n = 10,
        with_ties = FALSE
    ) |>
    arrange(relativeContribution) |>
    mutate(
        pollenTaxon = factor(pollenTaxon, levels = pollenTaxon)
    )

# Figure 1/3 : PCoA

p_pcoa <- ggplot(pcoa_df, aes(PCoA1, PCoA2)) +
    geom_path(colour = "grey65", linewidth = 0.5) +
    geom_point(aes(colour = pollenZone), size = 2.2) +
    geom_point(
        data = pcoa_df[1, , drop = FALSE],
        aes(shape = "Start"),
        size = 3.8, stroke = 1, colour = "black"
    ) +
    geom_point(
        data = pcoa_df[nrow(pcoa_df), , drop = FALSE],
        aes(shape = "End"),
        size = 3.8, stroke = 1, colour = "black"
    ) +
    scale_shape_manual(
        values = c("Start" = 0, "End" = 2),
        breaks = c("Start", "End")
    ) +
    scale_colour_manual(
        values = retro_palette,
        breaks = zone_levels,
        drop = FALSE
    ) +
    coord_equal() +
    labs(
        x = paste0("PCoA 1 (", axis_percent[1], "%)"),
        y = paste0("PCoA 2 (", axis_percent[2], "%)"),
        colour = "Pollen zone",
        shape = "Trajectory"
    ) +
    theme_classic() +
    theme(legend.position = "right")

svg(
    filename = paste0(
        "./Figures/",
        record,
        "_cta_pcoa_",
        Sys.Date(),
        ".svg"
    ),
    width = 16, height = 9
)
p_pcoa
dev.off()

ggsave(
    filename = paste0(
        "./Figures/",
        record,
        "_cta_pcoa_",
        Sys.Date(),
        ".png"
    ),
    plot = p_pcoa,
    width = 16, height = 9
)

# Figure 2A/3 : sample contribution to trajectory variance

p_lcbd <- ggplot(internal_df, aes(x = contribution, y = ageCalBP)) +
    geom_rect(
        data = zones,
        aes(
            xmin = -Inf, xmax = Inf,
            ymin = topAge, ymax = bottomAge,
            fill = pollenZone
        ),
        alpha = 0.12,
        inherit.aes = FALSE
    ) +
    geom_vline(
        xintercept = 1 / n,
        linetype = "dashed",
        colour = "#f46379"
    ) +
    geom_path(colour = "grey25", linewidth = 0.7) +
    geom_point(
        shape = 21,
        fill = "white",
        size = 1.5,
        na.rm = TRUE
    ) +
    scale_fill_manual(
        name = "Pollen zone",
        values = retro_palette,
        breaks = names(retro_palette)
    ) +
    scale_y_reverse() +
    expand_limits(y = internal_df$ageCalBP) +
    labs(
        title = "A",
        x = "Relative contribution to variance",
        y = "Age (cal yr BP)"
    ) +
    theme_classic()

p_lcbd

# Figure 2B/3 : Pollen zone variance

p_zone_variance <- ggplot(
    subtraj_summary,
    aes(x = pollenZone, y = internal_variance, fill = pollenZone)
) +
    geom_col(width = 0.72, show.legend = FALSE) +
    geom_text(
        aes(label = paste0("n=", n)),
        vjust = -0.45,
        size = 3.1
    ) +
    scale_fill_manual(
        values = retro_palette,
        breaks = zone_levels,
        drop = FALSE
    ) +
    labs(
        title = "B",
        x = "",
        y = "Internal variance"
    ) +
    theme_classic()

p_zone_variance

# Figure 2/3: panel assembly

library(patchwork)

p_lcbd_zone_var <- p_lcbd + p_zone_variance +
    plot_layout(
        widths = c(0.5, 0.5)
    )

p_lcbd_zone_var

svg(
    filename = paste0(
        "./Figures/",
        record,
        "_cta_lcbd_zone_variance_",
        Sys.Date(),
        ".svg"
    ),
    width = 16, height = 9
)
p_lcbd_zone_var
dev.off()

ggsave(
    filename = paste0(
        "./Figures/",
        record,
        "_cta_lcbd_zone_variance_",
        Sys.Date(),
        ".png"
    ),
    plot = p_lcbd_zone_var,
    width = 16, height = 9
)

# Figure 3A/3 : top 10 taxa contributing to global variance

# Taxon labels for plotting
## INF1
taxon_labels <- c(
    "Ornithogalum umbellatum-Typ B" = "Ornithogalum umbellatum type",
    "Scirpus sl" = "Scirpus s.l.",
    "Sparganium emersum type P" = "Sparganium emersum type",
    "Sparganium erectum type P" = "Sparganium erectum type"
)

## CPZ2
taxon_labels <- c(
    "Poaceae wild grass" = "Poaceae wild\ngrass group",
    "Aster type M" = "Aster type",
    "Scilla type M" = "Scilla type",
    "Avena-Triticum group A" = "Avena-Triticum group",
    "Ranunculus aquatilis group P" = "Ranunculus aquatilis group",
    "Cirsium palustre type P" = "Cirsium palustre type",
    "Urtica dioica type P" = "Urtica dioica type",
    "Pulsatilla vulgaris group P" = "Pulsatilla vulgaris group",
    "Ephedra dystachya-Typ B" = "Ephedra distachya type"
)

p_taxon_variance_global <- ggplot(
    top_taxa_global,
    aes(x = pollenTaxon, y = contributionPct)
) +
    geom_col(fill = "#3cb68e", width = 0.75) +
    scale_y_continuous() +
    labs(
        title = "A",
        y = "Contribution to global internal variance (%)",
        x = ""
    ) +
    scale_x_discrete(
        labels = function(x) {
            labels <- sub("___.*$", "", x)
            dplyr::recode(
                labels, !!!taxon_labels,
                .default = labels
            )
        }
    ) +
    theme_classic() +
    theme(axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5,
    ))

p_taxon_variance_global

# Figure 3B/3 : top 10 taxa per zone

top_taxa <- taxon_variance_zones |>
    group_by(pollenZone) |>
    slice_max(
        order_by = relativeContribution,
        n = 10,
        with_ties = FALSE
    ) |>
    ungroup() |>
    mutate(
        pollenZone = factor(pollenZone, levels = zone_levels),
        taxon_zone = paste(pollenTaxon, pollenZone, sep = "___")
    ) |>
    arrange(pollenZone, relativeContribution) |>
    mutate(
        taxon_zone = factor(taxon_zone, levels = unique(taxon_zone)),
        contributionPct = relativeContribution * 100
    )

p_taxon_variance <- ggplot(
    top_taxa,
    aes(x = taxon_zone, y = contributionPct)
) +
    geom_col(
        aes(fill = pollenZone),
        width = 0.75
    ) +
    scale_fill_manual(
        values = retro_palette,
        breaks = zone_levels,
        guide = "none"
    ) +
    facet_wrap(
        ~pollenZone,
        scales = "free_x",
        ncol = 2
    ) +
    scale_x_discrete(
        labels = function(x) {
            labels <- sub("___.*$", "", x)
            dplyr::recode(
                labels, !!!taxon_labels,
                .default = labels
            )
        }
    ) +
    scale_y_continuous(
        expand = expansion(mult = c(0, 0.08))
    ) +
    labs(
        title = "B",
        x = NULL,
        y = "Contribution to within-zone variance (%)"
    ) +
    theme_classic() +
    theme(
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(
            angle = 90,
            hjust = 1,
            vjust = 0.5,
            size = 9
        )
    )

p_taxon_variance

taxa_variance <- p_taxon_variance_global + p_taxon_variance +
    plot_layout(widths = c(0.4, 0.6))

taxa_variance

# Save

svg(
    filename = paste0(
        "./Figures/",
        record,
        "_cta_taxa_variance",
        Sys.Date(),
        ".svg"
    ),
    width = 16, height = 9
)
taxa_variance
dev.off()

ggsave(
    filename = paste0(
        "./Figures/",
        record,
        "_cta_taxa_variance",
        Sys.Date(),
        ".png"
    ),
    plot = taxa_variance,
    width = 16, height = 9
)

# Plotting other metrics #####

p_length <- plot_age_metric(
    segment_df, "length", "Segment length", "Hellinger distance", "#507c78"
) +
    geom_point(
        data = filter(segment_df, crossesZone),
        colour = "firebrick", shape = 4, size = 2.1
    ) +
    theme_classic()

p_speed <- plot_age_metric(
    segment_df, "speed100", "Segment speed",
    "Hellinger distance / 100 yr", "#b36f50"
)

p_angle <- plot_age_metric(
    angle_df, "angle", "Turning angle", "Degrees", "#8b7194"
) +
    scale_x_continuous(limits = c(0, 180), breaks = seq(0, 180, 45))
