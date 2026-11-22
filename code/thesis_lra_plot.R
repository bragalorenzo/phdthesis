#########################################################################
#--- PhD Thesis Data Analysis ---#--- Part 1: Palynological Richness ---#
#########################################################################

library(tidyverse)
library(vegan)

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

counts <- rio::import(input_file, sheet = "counts")

metadata <- rio::import(input_file, sheet = "sample_metadata")

zones <- rio::import(input_file, sheet = "zonation")

taxon_lookup <- rio::import(input_file, sheet = "taxon_lookup")

lra_taxon_lookup <- rio::import(paste0(record, "_lra_taxon_lookup.xlsx"))

lra_output <- rio::import(paste0(
    record,
    "_lra_full_biggest_radius_",
    "2026-11-16", # Date
    ".xlsx"
)) # INF1

# Unexplainable wave of intense emotions  #####

# Select outermost ring and aggregate by land cover
lra_landcover <- lra_output |>
    filter(outerRadius == max(outerRadius)) |>
    pivot_longer(
        cols = starts_with("mean."),
        names_to = "pollenTaxon",
        values_to = "cover"
    ) |>
    left_join(
        lra_taxon_lookup |> select(pollenTaxon, landType),
        by = "pollenTaxon"
    ) |>
    mutate(
        landType = replace_na(landType, "Unclassified")
    ) |>
    summarise(
        cover = sum(cover, na.rm = TRUE),
        .by = c(pollenZone, landType)
    ) |>
    mutate(
        pollenZone = factor(pollenZone, levels = zone_levels)
    )

unique(lra_taxon_lookup$landType)

# Plot ####

## Palette

landcover_palette <- c(
    "Conifers" = "#82A55B",
    "Deciduous forest" = "#466D53",
    "Riparian forest" = "#418C83",
    "Wetland meadow" = "#79B9CE",
    "Meadow/Pasture" = "#D9CD75",
    "Pasture" = "#E9A858",
    "Mediterranean" = "#9273A7",
    "Anthropic" = "#C77365",
    "Mesic woody shrubs" = "#9A9D88",
    "Unclassified" = "#BDBDBD"
)

## The actual plot

zone_totals <- lra_landcover |>
    summarise(
        total = sum(cover),
        .by = pollenZone
    )

p_lra <- ggplot(lra_landcover, aes(x = pollenZone, y = cover)) +
    geom_col(
        aes(fill = landType),
        width = 0.8
    ) +
    geom_col(
        data = zone_totals,
        aes(y = total),
        colour = "black",
        fill = NA,
        width = 0.8,
        show.legend = FALSE
    ) +
    scale_fill_manual(values = landcover_palette) +
    scale_y_continuous() +
    scale_x_discrete(
        limits = rev(zone_levels),
        labels = rev(zone_levels)
    ) +
    labs(
        x = "Pollen zone",
        y = "Reconstructed land cover (%)",
        fill = "Land cover"
    ) +
    theme_classic() +
    theme(
        legend.position = "right"
    )

p_lra

## Save
svg(
    filename = paste0(
        "./Figures/",
        record,
        "_lra_stacked_barplot_",
        Sys.Date(),
        ".svg"
    ),
    width = 16, height = 9
)
p_lra
dev.off()

ggsave(
    filename = paste0(
        "./Figures/",
        record,
        "_lra_stacked_barplot_",
        Sys.Date(),
        ".png"
    ),
    plot = p_lra,
    width = 16, height = 9
)

rio::export(
    lra_landcover,
    file = paste0(
        workdir,
        "/Outputs/",
        record,
        "_landcover_class_",
        Sys.Date(),
        ".xlsx"
    )
)

# Export a taxon-category land cover by zone table #####

lra_taxon_cover <- lra_output |>
    filter(outerRadius == max(outerRadius)) |>
    pivot_longer(
        cols = starts_with("mean."),
        names_to = "pollenTaxon",
        values_to = "cover"
    ) |>
    left_join(
        lra_taxon_lookup |>
            select(pollenTaxon, landType, order),
        by = "pollenTaxon"
    ) |>
    mutate(
        pollenZone = factor(pollenZone, levels = zone_levels),
        landType = replace_na(landType, "Unclassified"),
        pollenTaxon = sub("^mean\\.", "", pollenTaxon)
    ) |>
    arrange(pollenZone, order, pollenTaxon) |>
    select(pollenZone, pollenTaxon, cover, landType)

rio::export(
    lra_taxon_cover,
    file = paste0(
        workdir,
        "/Outputs/",
        record,
        "_landcover_taxon_and_class_",
        Sys.Date(),
        ".xlsx"
    )
)
