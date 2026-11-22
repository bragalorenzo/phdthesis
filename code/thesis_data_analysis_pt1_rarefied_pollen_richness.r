#########################################################################
#--- L.B. PhD Thesis ---#--- Part 1: Rarefied Palynological Richness ---#
#########################################################################

# Libraries ######

library(tidyverse)
library(vegan)
library(patchwork)

# Settings and data #####

set.seed(777)

record <- "inf1" # Manually change to swap records

load(
    paste0(
        "./phdthesis/data/", # Path
        record, "_rarefaction.Rdata" # File name
    )
)

# Helper functions #####

assign_zone <- function(depth, zones_table) {
    zone_idx <- which(depth >= zones_table$topDepth & depth <= zones_table$bottomDepth)

    if (length(zone_idx) == 1) {
        return(as.character(zones_table$pollenZone[zone_idx]))
    } else {
        return(NA_character_)
    }
}

# Rarefied palynological richness #####

## Run analysis

rarefied_richness <- vegan::rarefy(
    counts,
    sample = min(rowSums(counts))
)

## Organise results

diversity <- metadata |>
    select(-excludeSample) |>
    mutate(
        sample = rownames(diversity),
        depthCm = as.numeric(depthCm),
        ageCalBP = as.numeric(ageCalBP),
        rarefiedRichness = round(rarefied_richness, digits = 0),
        pollenZone = map_chr(depthCm, assign_zone, zones_table = zones)
    )

# Visualisation ####

## Clean data for plotting

diversity_plot <- diversity |>
    mutate(
        pollenZone = factor(
            pollenZone,
            levels = zone_levels
        )
    ) |>
    arrange(ageCalBP)

## Shared scale for richness

richness_limits <- range(
    diversity_plot$rarefiedRichness,
    na.rm = TRUE
) + c(-2, 2)

## Panel A: global boxplot for richness

p_global <- ggplot(
    diversity_plot,
    aes(y = rarefiedRichness, x = factor(1))
) +
    geom_boxplot(
        width = 0.35,
        outlier.shape = NA,
        fill = "white"
    ) +
    geom_jitter(
        aes(colour = pollenZone),
        width = 0.08,
        alpha = 0.7,
        size = 1.2
    ) +
    scale_colour_manual(values = palette) +
    scale_y_continuous(limits = richness_limits) +
    labs(
        x = toupper(record),
        y = "Rarefied palynological richness",
        title = "A"
    ) +
    theme_classic() +
    theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none"
    )

## Panel B: boxplot for each pollen zone

p_zones <- ggplot(
    diversity_plot,
    aes(y = rarefiedRichness, x = pollenZone)
) +
    geom_boxplot(
        aes(fill = pollenZone),
        width = 0.55,
        outlier.shape = NA,
        alpha = 0.75
    ) +
    geom_jitter(
        aes(colour = pollenZone),
        width = 0.10,
        alpha = 0.75,
        size = 1.2
    ) +
    scale_fill_manual(values = palette) +
    scale_colour_manual(values = palette) +
    scale_y_continuous(limits = richness_limits) +
    scale_x_discrete(limits = rev(zone_levels)) +
    labs(
        x = "",
        y = "",
        title = "B"
    ) +
    theme_classic() +
    theme(
        legend.position = "none"
    )

## Panel C: temporal behaviour of richness
p_time <- ggplot(
    diversity_plot,
    aes(x = ageCalBP, y = rarefiedRichness)
) +
    geom_path(
        linewidth = 0.4,
        colour = "grey25"
    ) +
    geom_point(
        aes(colour = pollenZone),
        size = 1.7
    ) +
    scale_colour_manual(values = palette) +
    scale_y_continuous(limits = richness_limits) +
    scale_x_reverse() +
    labs(
        x = "Age (cal yr BP)",
        y = "",
        title = "C"
    ) +
    theme_classic() +
    theme(
        legend.position = "none"
    )

## Full plot

right_panels <- (p_zones / p_time) +
    plot_layout(
        heights = c(7, 3)
    )

richness_figure <- (p_global + right_panels) +
    plot_layout(
        widths = c(2, 6)
    )

richness_figure
