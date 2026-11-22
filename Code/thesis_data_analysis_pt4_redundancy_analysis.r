#############################################################
#--- L.B. PhD Thesis ---#--- Part 3: Redundancy Analysis ---#
#############################################################

# Libraries #####

library(tidyverse)
library(vegan)

# Setup and Data #####

set.seed(777)

record <- "inf1" # Manually change to swap records

load(
    paste0(
        "./phdthesis/data/", # Path
        record, "_rda_data.Rdata" # File name
    )
)

# Response matrix #####

## Prepare pollen counts
counts_rda <- counts_aggregated |>
    column_to_rownames("sample") |>
    select(-nd) |> # Remove Indeterminable pollen grains
    as.matrix()

## Hellinger transformation (samples ordered as in metadata)
counts_hellinger <- vegan::decostand(
    counts_rda,
    method = "hellinger"
)

# Predictor matrix #####

## Align climatic data to pollen sample ages (nearest neighbour)
climate_aligned <- metadata |>
    mutate(
        AgeBP_climate = climate_data$AgeBP[
            vapply(
                ageCalBP,
                \(age) which.min(abs(climate_data$AgeBP - age)),
                integer(1)
            )
        ]
    ) |>
    left_join(
        climate_data,
        by = c("AgeBP_climate" = "AgeBP")
    )

## Select and scale climatic predictors

climate_selected <- c( # Only two chosen because of multicollinearity
    "temperature",
    "summer_precipitation"
)

env_rda <- climate_aligned |>
    select(climate_selected) |>
    scale() |>
    as.data.frame()

# RDA #####

rda_climate <- vegan::rda(counts_hellinger ~ ., data = env_rda)

## Model summary
summary(rda_climate)
vegan::RsquareAdj(rda_climate)
vegan::vif.cca(rda_climate)

## Permutation tests
anova(rda_climate, permutations = 9999)
anova(rda_climate, by = "axis", permutations = 9999)
anova(rda_climate, by = "margin", permutations = 9999)

# Visualisation #####

## Explained variance
constrained_pct <- 100 * sum(rda_climate$CCA$eig) / rda_climate$tot.chi
axis_pct <- 100 * rda_climate$CCA$eig[1:2] / sum(rda_climate$CCA$eig)

## Samples' scores

samples <- vegan::scores(
    rda_climate,
    display = "sites", scaling = 2
) |>
    as.data.frame() |>
    rownames_to_column("sample") |>
    left_join(
        metadata |> select(sample, pollenZone),
        by = "sample"
    )

## Climatic variables' scores

environment <- vegan::scores(
    rda_climate,
    display = "bp", scaling = 2
) |>
    as.data.frame() |>
    rownames_to_column("variable")

### Correct labels

environment_labels <- c(
    temperature = "Mean annual temperature",
    precipitation = "Annual precipitation",
    precipitation_seasonality = "Precipitation seasonality",
    relative_humidity = "Relative humidity",
    summer_precipitation = "Summer precipitation"
)

environment <- environment |>
    mutate(label = coalesce(
        unname(environment_labels[variable]), variable
    ))


## Ecological groups' scores

# Pollen taxa contained in each group are visible in taxon_lookup
# ecoMacrogroup column

### Labels

group_labels <- c(
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
)

groups <- vegan::scores(
    rda_climate,
    display = "species", scaling = 2
) |>
    as.data.frame() |>
    rownames_to_column("group") |>
    mutate(
        length_2d = sqrt(RDA1^2 + RDA2^2),
        label = coalesce(unname(group_labels[group]), group)
    ) |>
    slice_max(length_2d, n = 5, with_ties = FALSE)

## Plot

p_rda <- ggplot() +
    geom_hline(yintercept = 0, colour = "grey80") +
    geom_vline(xintercept = 0, colour = "grey80") +

    # Ecological groups (arrow length multiplier = 10)
    geom_segment(
        data = groups,
        aes(
            x = 0, y = 0,
            xend = RDA1 * 10,
            yend = RDA2 * 10
        ),
        arrow = grid::arrow(length = grid::unit(0.15, "cm")),
        colour = "grey50",
        linetype = "dashed"
    ) +
    geom_text(
        data = groups,
        aes(
            x = RDA1 * 10,
            y = RDA2 * 10,
            label = label
        ),
        colour = "grey40",
        size = 3,
        vjust = 1
    ) +

    # Climatic variables (arrow multiplier = 1.2)
    geom_segment(
        data = environment,
        aes(
            x = 0, y = 0,
            xend = RDA1 * 1.2,
            yend = RDA2 * 1.2
        ),
        arrow = grid::arrow(length = grid::unit(0.2, "cm")),
        colour = "black",
        linewidth = 0.7
    ) +
    geom_text(
        data = environment,
        aes(
            x = RDA1 * 1.2,
            y = RDA2 * 1.2,
            label = label
        ),
        size = 3,
        vjust = -1,
        hjust = 1
    ) +

    # Samples
    geom_point(
        data = samples,
        aes(RDA1, RDA2, fill = pollenZone),
        shape = 21, size = 3, colour = "black"
    ) +
    scale_fill_manual(values = palette) +
    coord_equal() +
    labs(
        title = sprintf(
            "Constrained variance: %.1f%% of total",
            constrained_pct
        ),
        x = sprintf("RDA1 (%.1f%%)", axis_pct[1]),
        y = sprintf("RDA2 (%.1f%%)", axis_pct[2]),
        fill = "Pollen zone"
    ) +
    theme_classic() +
    theme(legend.position = "right")

p_rda
