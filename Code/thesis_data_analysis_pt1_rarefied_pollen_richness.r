#########################################################################
#--- PhD Thesis Data Analysis ---#--- Part 1: Palynological Richness ---#
#########################################################################

library(tidyverse)
library(vegan)
library(rioja)

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

counts <- rio::import(input_file, sheet = "counts")

metadata <- rio::import(input_file, sheet = "sample_metadata")

zones <- rio::import(input_file, sheet = "zonation")

taxon_lookup <- rio::import(input_file, sheet = "taxon_lookup")

## Transpose

counts_t <- counts |>
    t() |>
    as.data.frame()

colnames(counts_t) <- counts_t[1, ]

counts_t <- counts_t[-1, ]

## Clean

metadata <- counts_t |>
    filter(excludeSample == "no") |>
    select(
        excludeSample,
        depthCm,
        ageCalBP
    )

counts_clean <- counts_t |>
    filter(excludeSample == "no") |>
    select(-c(excludeSample, depthCm, ageCalBP)) |>
    mutate(across(everything(), as.numeric))

## Structure

head(counts_clean)
str(counts_clean)

# Rarefied palynological richness ####
# Rarefaction
rarefied_richness <- vegan::rarefy(
    counts_clean,
    sample = min(rowSums(counts_clean))
)

summary(rarefied_richness)

# Function to assign pollen zone from depth
assign_zone <- function(depth, zones_table) {
    zone_idx <- which(depth >= zones_table$topDepth & depth <= zones_table$bottomDepth)

    if (length(zone_idx) == 1) {
        return(as.character(zones_table$pollenZone[zone_idx]))
    } else {
        return(NA_character_)
    }
}

# Diversity table
diversity <- metadata |>
    select(-excludeSample) |>
    mutate(
        sample = rownames(diversity),
        depthCm = as.numeric(depthCm),
        ageCalBP = as.numeric(ageCalBP),
        rarefiedRichness = round(rarefied_richness, digits = 0),
        pollenZone = map_chr(depthCm, assign_zone, zones_table = zones)
    )

head(diversity)
table(diversity$pollenZone, useNA = "ifany")

# Save data

rio::export(
    diversity,
    file = paste0(
        "./Outputs/",
        record,
        "_palynological_richness_",
        Sys.Date(),
        ".xlsx"
    )
)

# Visualisation ####

library(patchwork)

diversity_plot <- diversity |>
    mutate(
        pollenZone = factor(
            pollenZone,
            levels = zone_levels
        )
    ) |>
    arrange(ageCalBP)

## Shared scale for richness

richness_limits <- range(diversity_plot$rarefiedRichness, na.rm = TRUE) + c(-2, 2)

# Panel A: global boxplot for richness

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
    scale_colour_manual(values = retro_palette) +
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

# Panel B: boxplot for each pollen zone

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
    scale_fill_manual(values = retro_palette) +
    scale_colour_manual(values = retro_palette) +
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

# Panel C: temporal behaviour of richness
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
    scale_colour_manual(values = retro_palette) +
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

richness_panels <- (p_global + right_panels) +
    plot_layout(
        widths = c(2, 6)
    )

richness_panels

## Save

svg(
    filename = paste0(
        "./Figures/",
        record,
        "_palynological_richness_",
        Sys.Date(),
        ".svg"
    ),
    width = 16, height = 9
)
richness_panels
dev.off()

ggsave(
    filename = paste0(
        "./Figures/",
        record,
        "_palynological_richness_",
        Sys.Date(),
        ".png"
    ),
    plot = richness_panels,
    width = 16, height = 9
)

## Save (poster)

svg("palynological_richness_roc.svg", width = 30, height = 13.23)
full_figure
dev.off()

ggsave("palynological_richness_roc.png")

# In case of palynological richness + RoC plot #####

richness_group <- richness_panels +
    plot_annotation(
        title = "Rarefied pollen richness",
        theme = theme(
            plot.title = element_text(
                hjust = 0.5,
                size = 14
            )
        )
    )

full_figure <- wrap_elements(full = richness_group) + roc_plot +
    plot_layout(
        widths = c(8, 2)
    )

full_figure

## Save

svg("./Figures/palynological_richness_roc.svg", width = 16, height = 9)
full_figure
dev.off()

ggsave("palynological_richness_roc.png")

## Save (poster)

svg("palynological_richness_roc.svg", width = 30, height = 13.23)
full_figure
dev.off()

ggsave("palynological_richness_roc.png")
