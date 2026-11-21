##################################################################
#--- L.B. PhD Thesis  ---#--- Part 2: Rate of Change Analysis ---#
##################################################################

# Libraries #####

library(tidyverse)
library(RRatepol)
library(patchwork)

## If RRatepol is not installed

remotes::install_github("HOPE-UIB-BIO/R-Ratepol-package")

# Setup and Data #####

set.seed(777)

record <- "inf1" # Manually change to swap records

load(
    paste0(
        "./phdthesis/data/", # Path
        record, "_harmonised_data.Rdata" # File name
    )
)

## RoC bin sizes
if (tolower(record) == "inf1") {
    # Breaks for y axis of RoC plots
    y_breaks <- seq(500, 4000, 500)
    # Bin sizes for base run and sensitivity analyses
    roc_bin_size <- 100
    bin_s1 <- 150
    bin_s2 <- 200
} else if (tolower(record) == "cpz2") {
    # Breaks for y axis of RoC plots
    y_breaks <- seq(500, 7000, 500)
    # Bin sizes for base run and sensitivity analyses
    roc_bin_size <- 300
    bin_s1 <- 200
    bin_s2 <- 400
} else {
    message("Record not found.")
}

# Rate of Change ####

## Prepare data for RRatepol input

community <- counts_harmonised |>
    select(-"Indeterminable") |> # Remove Indeterminable
    rename(sample_id = sample)

ages <- metadata |>
    select(sample, ageCalBP) |>
    rename(sample_id = sample, age = ageCalBP) |>
    mutate(age = as.numeric(age))

## Check age resolution to set RoC bin size (beginning of script)

age_resolution <- ages |>
    arrange(age) |>
    summarise(
        min    = min(diff(age)),
        q25    = quantile(diff(age), 0.25),
        median = median(diff(age)),
        mean   = mean(diff(age)),
        q75    = quantile(diff(age), 0.75),
        max    = max(diff(age))
    )

age_resolution

ages |>
    arrange(age) |>
    mutate(dt = c(NA, diff(age))) |>
    ggplot(aes(x = age, y = dt)) +
    geom_line() +
    geom_point() +
    scale_x_reverse() +
    labs(
        x = "Age (cal yr BP)",
        y = "Temporal resolution (yr)"
    ) +
    theme_classic()

## Run

### Decide bin size pre-run (based on age resolution)

roc <- RRatepol::estimate_roc(
    data_source_community = community,
    data_source_age = ages,
    working_units = "MW",
    bin_size = roc_bin_size,
    number_of_shifts = 5,
    dissimilarity_coefficient = "chord",
    time_standardisation = 100,
    smooth_method = "age.w",
    rand = 1000
)

## Detect peaks

roc_peaks <- RRatepol::detect_peak_points(
    data_source = roc,
    sel_method = "trend_non_linear"
)

roc_trend <- roc_peaks |>
    mutate(
        ROC_smooth = predict(
            mgcv::gam(ROC ~ s(Age, k = 3), data = roc_peaks),
            newdata = roc_peaks
        )
    )

RRatepol::plot_roc(
    data_source = roc_peaks,
    peaks = TRUE,
    trend = "trend_non_linear"
)

## Sensitivity analysis with different bin sizes #####

### First sensitivity analysis

roc_s1 <- RRatepol::estimate_roc(
    data_source_community = community,
    data_source_age = ages,
    working_units = "MW",
    bin_size = bin_s1,
    number_of_shifts = 5,
    dissimilarity_coefficient = "chord",
    time_standardisation = 100,
    smooth_method = "age.w",
    rand = 1000
)

roc_s1_peaks <- RRatepol::detect_peak_points(
    data_source = roc_s1,
    sel_method = "trend_non_linear"
)

### Second sensitivity analysis

roc_s2 <- RRatepol::estimate_roc(
    data_source_community = community,
    data_source_age = ages,
    working_units = "MW",
    bin_size = bin_s2,
    number_of_shifts = 5,
    dissimilarity_coefficient = "chord",
    time_standardisation = 100,
    smooth_method = "age.w",
    rand = 1000
)

roc_s2_peaks <- RRatepol::detect_peak_points(
    data_source = roc_s2,
    sel_method = "trend_non_linear"
)

### Fast visualisation

roc_plot <- RRatepol::plot_roc(
    data_source = roc_peaks,
    peaks = TRUE,
    trend = "trend_non_linear"
) +
    labs(title = paste("RoC", roc_bin_size))

roc_s1_plot <- RRatepol::plot_roc(
    data_source = roc_s1_peaks,
    peaks = TRUE,
    trend = "trend_non_linear"
) +
    labs(title = paste("RoC", bin_s1))

roc_s2_plot <- RRatepol::plot_roc(
    data_source = roc_s2_peaks,
    peaks = TRUE,
    trend = "trend_non_linear"
) +
    labs(title = paste("RoC", bin_s2))

roc_sensitivity_check <- roc_plot + roc_s1_plot + roc_s2_plot

roc_sensitivity_check

# Visualisation with zone-based palette #####

## GAM trends (smooth line on plot) #####

roc_trend <- roc_peaks |>
    mutate(
        ROC_smooth = predict(
            mgcv::gam(ROC ~ s(Age, k = 3), data = roc_peaks),
            newdata = roc_peaks
        )
    )

roc_s1_trend <- roc_s1_peaks |>
    mutate(
        ROC_smooth = predict(
            mgcv::gam(ROC ~ s(Age, k = 3), data = roc_s1_peaks),
            newdata = roc_s1_peaks
        )
    )

roc_s2_trend <- roc_s2_peaks |>
    mutate(
        ROC_smooth = predict(
            mgcv::gam(ROC ~ s(Age, k = 3), data = roc_s2_peaks),
            newdata = roc_s2_peaks
        )
    )


## Panel A: selected bin size RoC plot

roc_plot <- ggplot(roc_peaks, aes(x = ROC, y = Age)) +
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
    geom_ribbon(
        aes(
            xmin = ROC_dw,
            xmax = ROC_up,
            y = Age
        ),
        fill = "grey80",
        alpha = 0.5,
        orientation = "y",
        inherit.aes = FALSE
    ) +
    geom_path(colour = "grey25", linewidth = 0.7) +
    geom_path(
        data = roc_trend |> arrange(Age),
        aes(x = ROC_smooth, y = Age),
        colour = "#6688c3",
        linewidth = 0.9,
        inherit.aes = FALSE
    ) +
    geom_point(
        data = subset(roc_peaks, Peak),
        colour = "#39d353",
        size = 2.5
    ) +
    scale_fill_manual(
        name = "Pollen zones",
        values = palette,
        breaks = names(palette)
    ) +
    guides(
        fill = guide_legend(
            override.aes = list(alpha = 0.5)
        )
    ) +
    theme(
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_text(face = "bold"),
        legend.key.size = unit(0.5, "cm")
    ) +
    scale_y_reverse(breaks = y_breaks) +
    coord_cartesian(ylim = c(max(metadata$ageCalBP) + 100, 300)) +
    labs(
        title = paste0("A – MW = ", roc_bin_size, " years"),
        x = NULL,
        y = "Age (cal yr BP)"
    ) +
    theme_classic()

## Panel B: RoC sensitivity 1

roc_s1_plot <- ggplot(roc_s1_peaks, aes(x = ROC, y = Age)) +
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
    geom_ribbon(
        aes(
            xmin = ROC_dw,
            xmax = ROC_up,
            y = Age
        ),
        fill = "grey80",
        alpha = 0.5,
        orientation = "y",
        inherit.aes = FALSE
    ) +
    geom_path(colour = "grey25", linewidth = 0.7) +
    geom_path(
        data = roc_s1_trend |> arrange(Age),
        aes(x = ROC_smooth, y = Age),
        colour = "#6688c3",
        linewidth = 0.9,
        inherit.aes = FALSE
    ) +
    geom_point(
        data = subset(roc_s1_peaks, Peak),
        colour = "#39d353",
        size = 2.5
    ) +
    scale_fill_manual(
        name = "Pollen zones",
        values = palette,
        breaks = names(palette)
    ) +
    guides(
        fill = guide_legend(
            override.aes = list(alpha = 0.5)
        )
    ) +
    theme(
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_text(face = "bold"),
        legend.key.size = unit(0.5, "cm")
    ) +
    scale_y_reverse(breaks = y_breaks) +
    coord_cartesian(ylim = c(max(metadata$ageCalBP) + 100, 300)) +
    labs(
        title = paste0("B – MW = ", bin_s1, " years"),
        x = "RoC score",
        y = "Age (cal yr BP)"
    ) +
    theme_classic() +
    theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
    )

## Panel C: RoC sensitivity 2

roc_s2_plot <- ggplot(roc_s2_peaks, aes(x = ROC, y = Age)) +
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
    geom_ribbon(
        aes(
            xmin = ROC_dw,
            xmax = ROC_up,
            y = Age
        ),
        fill = "grey80",
        alpha = 0.5,
        orientation = "y",
        inherit.aes = FALSE
    ) +
    geom_path(colour = "grey25", linewidth = 0.7) +
    geom_path(
        data = roc_s2_trend |> arrange(Age),
        aes(x = ROC_smooth, y = Age),
        colour = "#6688c3",
        linewidth = 0.9,
        inherit.aes = FALSE
    ) +
    geom_point(
        data = subset(roc_s2_peaks, Peak),
        colour = "#39d353",
        size = 2.5
    ) +
    scale_fill_manual(
        name = "Pollen zones",
        values = palette,
        breaks = names(palette)
    ) +
    guides(
        fill = guide_legend(
            override.aes = list(alpha = 0.5)
        )
    ) +
    theme(
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_text(face = "bold"),
        legend.key.size = unit(0.5, "cm")
    ) +
    scale_y_reverse(breaks = y_breaks) +
    coord_cartesian(ylim = c(max(metadata$ageCalBP) + 100, 300)) +
    labs(
        title = paste0("C – MW = ", bin_s2, " years"),
        x = NULL,
        y = "Age (cal yr BP)"
    ) +
    theme_classic() +
    theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
    )

## Combined figure

roc_sensitivity <- (roc_plot + roc_s1_plot + roc_s2_plot) +
    plot_layout(ncol = 3, guides = "collect") &
    theme(legend.position = "right")

roc_sensitivity
