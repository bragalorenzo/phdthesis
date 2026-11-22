#######################################################################
#--- L.B. PhD Thesis ---#--- Part 3: Community Trajectory Analysis ---#
#######################################################################

# Savestates #####
save.image("inf1_cta_workspace.RData")
load("cpz2_cta_workspace.RData")

# Libraries #####

library(tidyverse)
library(vegan)
library(ecotraj)
library(patchwork)

## If ecotraj is not installed

remotes::install_github("emf-creaf/ecotraj")

# Setup and Data #####

set.seed(777)

record <- "inf1" # Manually change to swap records

load(
    paste0(
        "./Github/phdthesis/data/", # Path
        record, "_harmonised_data.Rdata" # File name
    )
)

# Helper functions #####

## Extract taxon variance

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

## Palette
if (tolower(record) == "inf1") {
    # INF1 palette
    palette <- c(
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
    palette <- c(
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
    arrange(desc(ageCalBP))

# Put sample IDs in row names ONLY for the numeric community matrix

counts_traj <- counts_harmonised |>
    select(-"Indeterminable") |> # Remove Indeterminable pollen
    tibble::column_to_rownames("sample")

counts_traj <- counts_traj[
    metadata_traj$sample, ,
    drop = FALSE
]

## Check for enough community states (3)

n <- nrow(metadata_traj)
stopifnot(n >= 3) # Enough observations for CTA

## Attribute samples to pollen zone

zones <- zones |>
    arrange(topDepth) |>
    mutate(pollenZone = as.character(pollenZone))

zone_breaks <- c(zones$topDepth[1], zones$bottomDepth)

metadata_traj <- metadata_traj |>
    mutate(
        pollenZone = as.character(cut(
            depthCm,
            breaks = zone_breaks,
            labels = zones$pollenZone,
            include.lowest = TRUE
        )),
    )

# Character levels for chronological plotting, oldest -> youngest.
zone_levels <- rev(zones$pollenZone)

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

## Global trajectory metrics #####

### General metrics

trajectory_metrics <- ecotraj::trajectoryMetrics(trajectory)

### Explore segment lengths and speeds

trajectory_lengths <- ecotraj::trajectoryLengths(trajectory)

trajectory_speeds <- ecotraj::trajectorySpeeds(trajectory)

segment_df <- data.frame(
    sampleFrom = metadata_traj$sample[-n], # Remove youngest sample
    sampleTo = metadata_traj$sample[-1], # Remove oldest sample
    ageFrom = metadata_traj$ageCalBP[-n],
    ageTo = metadata_traj$ageCalBP[-1],
    zoneFrom = metadata_traj$pollenZone[-n],
    zoneTo = metadata_traj$pollenZone[-1],
    length = as.numeric(unlist(
        trajectory_lengths[1, paste0("S", seq_len(n - 1)), drop = FALSE],
        use.names = FALSE
    )),
    speed = as.numeric(unlist(
        trajectory_speeds[1, paste0("S", seq_len(n - 1)), drop = FALSE],
        use.names = FALSE
    ))
) |>
    mutate(
        ageMidpoint = (ageFrom + ageTo) / 2,
        speed100 = speed * 100,
        crossesZone = zoneFrom != zoneTo,
        cumulativeLength = cumsum(length)
    )

### Explore angles

trajectory_angles <- ecotraj::trajectoryAngles(trajectory)

angle_df <- data.frame(
    sample = metadata_traj$sample[2:(n - 1)],
    ageCalBP = metadata_traj$ageCalBP[2:(n - 1)],
    angle = as.numeric(unlist(
        trajectory_angles[1, , drop = FALSE],
        use.names = FALSE
    ))
)

### State-wise (sample-wise) relative contribution to variance

trajectory_variation <- ecotraj::trajectoryInternalVariation(
    trajectory,
    relativeContributions = TRUE
)

internal_df <- data.frame(
    sample = metadata_traj$sample,
    ageCalBP = metadata_traj$ageCalBP,
    pollenZone = metadata_traj$pollenZone,
    contribution = as.numeric(unlist(
        trajectory_variation[1, paste0("rc_", seq_len(n)), drop = FALSE],
        use.names = FALSE
    ))
) |>
    mutate(
        contributionPct = round(contribution * 100, 2)
    )

# Pollen zone subtrajectories #####

# Each subtrajectory requires at least three samples for all metrics

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
# Species Contribution to Beta Diversity (SCBD)

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

## Decompose taxon variance for each pollen zone

taxon_variance_zones <- taxon.variance(
    metadata_traj$pollenZone
) |>
    rename(pollenZone = group)

## Run PCoA to plot whole trajectory in reduced space

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

# Visualisation #####

## Order pollen zone factors before plotting

pcoa_df <- pcoa_df |>
    mutate(pollenZone = factor(as.character(pollenZone), levels = zone_levels))

internal_df <- internal_df |>
    mutate(pollenZone = factor(as.character(pollenZone), levels = zone_levels))

subtraj_summary <- subtraj_summary |>
    mutate(pollenZone = factor(as.character(pollenZone), levels = zone_levels))

## Figure 1: trajectory PCoA

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
        values = palette,
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

p_pcoa

## Figure 2A: community state (sample) contribution to trajectory variance

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
        values = palette,
        breaks = names(palette)
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

## Figure 2B: internal variance of each pollen zone

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
        values = palette,
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

## Figure 2: complete

p_lcbd_zone_var <- p_lcbd + p_zone_variance +
    plot_layout(
        widths = c(0.5, 0.5)
    )

p_lcbd_zone_var

## Figure 3A: top 10 taxa contributing to global variance

### Clean taxon labels for plotting

if (record == tolower("inf1")) {
    # INF1
    taxon_labels <- c(
        "Ornithogalum umbellatum-Typ B" = "Ornithogalum umbellatum type",
        "Scirpus sl" = "Scirpus s.l.",
        "Sparganium emersum type P" = "Sparganium emersum type",
        "Sparganium erectum type P" = "Sparganium erectum type"
    )
} else if (record == tolower("cpz2")) {
    # CPZ2
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
}

### Plot

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

## Figure 3B: top 10 taxa per zone

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
        values = palette,
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

## Figure 3: complete

taxa_variance <- p_taxon_variance_global + p_taxon_variance +
    plot_layout(widths = c(0.4, 0.6))

taxa_variance
