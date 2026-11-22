########################
#--- REVEALS + LOVE ---#
########################

library(tidyverse)

# Preparation ####

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

# Regional data preparation ######

regional_taxon_lookup <- rio::import("fimon_regional.xlsx",
    sheet = "lra_taxon_lookup"
)

regional_counts <- rio::import("fimon_regional.xlsx", sheet = "counts")

rpp_fsp_lookup <- rio::import("rpp_fsp_lookup.xlsx")

# Normalise non-breaking spaces in RPP taxon names

rpp_fsp_lookup <- rpp_fsp_lookup |>
    mutate(
        taxon = gsub("\u00A0", " ", taxon, fixed = TRUE)
    )

## Build harmonisation table based on RPP IDs

harmon_regional <- regional_taxon_lookup |>
    filter(
        rppID != "abs"
    ) |>
    mutate(rppID = as.integer(rppID)) |>
    group_by(rppID) |>
    summarise(
        includedTaxa = paste(unique(pollenTaxon), collapse = " + "),
        .groups = "drop"
    ) |>
    left_join(
        rpp_fsp_lookup |>
            select(
                rppTaxon = taxon,
                rppID,
                fsp,
                rpp,
                rppSd
            ),
        by = "rppID"
    ) |>
    select(rppTaxon, includedTaxa, rppID, fsp, rpp, rppSd)

harmon_local <- taxon_lookup |>
    filter(
        rppID != "abs"
    ) |>
    mutate(rppID = as.integer(rppID)) |>
    group_by(rppID) |>
    summarise(
        includedTaxa = paste(unique(pollenTaxon), collapse = " + "),
        .groups = "drop"
    ) |>
    left_join(
        rpp_fsp_lookup |>
            select(
                rppTaxon = taxon,
                rppID,
                fsp,
                rpp,
                rppSd
            ),
        by = "rppID"
    ) |>
    select(rppTaxon, includedTaxa, rppID, fsp, rpp, rppSd)

## Aggregate counts by rppID

# Local record

local_rpp <- taxon_lookup |>
    filter(!is.na(rppID), rppID != "abs") |>
    mutate(rppID = as.integer(rppID)) |>
    filter(
        pollenTaxon %in% names(counts_clean),
        rppID %in% rpp_fsp_lookup$rppID
    )

local_groups <- split(local_rpp$pollenTaxon, local_rpp$rppID)

local_counts_rpp <- do.call(cbind, lapply(local_groups, function(taxa) {
    rowSums(counts_clean[, taxa, drop = FALSE])
}))

colnames(local_counts_rpp) <- rpp_fsp_lookup$taxon[
    match(as.integer(names(local_groups)), rpp_fsp_lookup$rppID)
]

local_counts_rpp <- bind_cols(
    metadata |>
        transmute(
            depthCm = as.numeric(depthCm),
            ageCalBP = as.numeric(ageCalBP)
        ),
    local_counts_rpp
)

glimpse(local_counts_rpp)

# Regional record

regional_rpp <- regional_taxon_lookup |>
    filter(!is.na(rppID), rppID != "abs") |>
    mutate(rppID = as.integer(rppID)) |>
    filter(
        pollenTaxon %in% names(regional_counts),
        rppID %in% rpp_fsp_lookup$rppID
    )

regional_groups <- split(regional_rpp$pollenTaxon, regional_rpp$rppID)

regional_counts_rpp <- do.call(cbind, lapply(regional_groups, function(taxa) {
    rowSums(regional_counts[, taxa, drop = FALSE], na.rm = TRUE)
}))

colnames(regional_counts_rpp) <- rpp_fsp_lookup$taxon[
    match(as.integer(names(regional_groups)), rpp_fsp_lookup$rppID)
]

regional_counts_rpp <- bind_cols(
    regional_counts |> select(ageCalBP),
    as_tibble(regional_counts_rpp)
)

glimpse(regional_counts_rpp)

# Record-exclusive taxa #####

# Common taxa
intersect(
    names(local_counts_rpp)[-c(1, 2)],
    names(regional_counts_rpp)[-1]
)

# Local-exclusive taxa
setdiff(
    names(local_counts_rpp)[-c(1, 2)],
    names(regional_counts_rpp)[-1]
)

# Regional-exclusive taxa
setdiff(
    names(regional_counts_rpp)[-1],
    names(local_counts_rpp)[-c(1, 2)]
)

# Define taxa to exclude from reconstruction #####

exclude_taxa <- c(
    "Abies alba", # Upland
    "Betula (mainly B. pubescens, B. pendula)",
    "Picea abies", # Upland
    "Fagus sylvatica", # Upland
    "Pinus (mainly P. sylvestris)", # Upland
    "Fagus", # Upland
    setdiff(
        names(regional_counts_rpp)[-1],
        names(local_counts_rpp)[-c(1, 2)]
    ) # Regional record-exclusive taxa
)

local_counts_rpp <- local_counts_rpp |>
    select(-any_of(exclude_taxa))

regional_counts_rpp <- regional_counts_rpp |>
    select(-any_of(exclude_taxa))

# Chronological zonation #####

zones <- zones |>
    arrange(topAge)

age_breaks <- c(zones$topAge, last(zones$bottomAge))

# Local
local_zoned <- local_counts_rpp |>
    mutate(
        pollenZone = cut(
            ageCalBP,
            breaks = age_breaks,
            labels = zones$pollenZone,
            include.lowest = TRUE,
            right = FALSE
        )
    )

# Regional
regional_zoned <- regional_counts_rpp |>
    mutate(
        pollenZone = cut(
            ageCalBP,
            breaks = age_breaks,
            labels = zones$pollenZone,
            include.lowest = TRUE,
            right = FALSE
        )
    ) |>
    filter(!is.na(pollenZone))

# Checks
table(local_zoned$pollenZone)
table(regional_zoned$pollenZone)

# Chronological aggregation of counts #####

local_pool <- local_zoned |>
    select(-depthCm, -ageCalBP, -any_of(exclude_taxa)) |>
    group_by(pollenZone) |>
    summarise(across(everything(), sum), .groups = "drop")

# Regional record
regional_pool <- regional_zoned |>
    select(-ageCalBP, -any_of(exclude_taxa)) |>
    group_by(pollenZone) |>
    summarise(across(everything(), sum), .groups = "drop")

# Check exclusions
stopifnot(
    !any(exclude_taxa %in% names(local_pool)),
    !any(exclude_taxa %in% names(regional_pool))
)
# Checks
stopifnot(
    nrow(local_pool) == nrow(zones),
    nrow(regional_pool) == nrow(zones),
    identical(
        as.character(local_pool$pollenZone),
        as.character(regional_pool$pollenZone)
    )
)

# Checks #####

table(local_zoned$pollenZone)
table(regional_zoned$pollenZone)

dim(local_pool)
dim(regional_pool)

stopifnot(
    identical(
        as.character(local_pool$pollenZone),
        as.character(regional_pool$pollenZone)
    )
)

# Pre-LRA data cleaning #####

## Remove special chars & spaces from names (required by disqover)
names(local_pool) <- gsub("[.,()/+[:space:]-]+", "_", names(local_pool))
names(regional_pool) <- gsub("[.,()/+[:space:]-]+", "_", names(regional_pool))

rio::export(
    list(
        "local_agg_counts" = local_pool,
        "regional_agg_counts" = regional_pool
    ),
    file = paste0(
        workdir,
        "/Outputs/",
        record,
        "_lra_input_counts_",
        Sys.Date(),
        ".xlsx"
    )
)

## Create output directory

output_dir <- paste(
    toupper(record),
    "LRA",
    Sys.Date()
)

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

## Data

local_data <- local_pool
regional_data <- regional_pool
rings <- rio::import("rings.xlsx")

## Create run parameter list

params_data <- harmon_local |>
    select(
        rppID,
        modelledLocalTaxa = includedTaxa
    ) |>
    left_join(
        harmon_regional |>
            select(
                rppID,
                modelledRegionalTaxa = includedTaxa
            ),
        by = "rppID"
    ) |>
    left_join(
        rpp_fsp_lookup |>
            select(
                rppID,
                rppTaxon = taxon,
                fsp,
                rpp,
                rppSd,
                source = study
            ),
        by = "rppID"
    ) |>
    mutate(
        cleanRppTaxon = gsub(
            "[.,()/+[:space:]-]+", "_", rppTaxon
        )
    ) |>
    filter(cleanRppTaxon %in% names(local_pool)[-1]) |>
    select(
        rppTaxon,
        cleanRppTaxon,
        modelledLocalTaxa,
        modelledRegionalTaxa,
        fsp,
        rpp,
        rppSd,
        source
    )

rio::export(params_data, "params.xlsx")

# Local basin diameter in m (100 INF1, 50 CPZ2)

local_basindiameter <- 100

# Area is 67 hectares (67/100 km2) Fimon, 6 km2 Albano

regional_basinarea <- 0.67

params_run <- list(
    ppes = data.frame(
        Taxon = params_data$cleanRppTaxon,
        fallspeed = params_data$fsp,
        ppe = params_data$rpp,
        ppe.error = params_data$rppSd
    ),
    fallspeed = data.frame(
        Taxon = params_data$cleanRppTaxon,
        fallspeed = params_data$fsp
    ),
    dBasinRegional = round(
        sqrt(
            regional_basinarea / pi
        ) * 2 * 1000,
        digits = 0
    ), # diameter in m
    dBasinLocal = local_basindiameter,
    tBasinRegional = "lake",
    tBasinLocal = "peatland",
    dwm = "lsm unstable",
    repeats = 10,
    n = 1,
    error = TRUE
)

# Match regional taxa to local taxa ####

local_taxa <- setdiff(names(local_data), "pollenZone")

missing_taxa <- setdiff(local_taxa, names(regional_data))

regional_data[missing_taxa] <- 0

regional_data <- regional_data[c("pollenZone", local_taxa)]

stopifnot(identical(names(local_data), names(regional_data)))

# Check that original regional counts are unchanged

shared_taxa <- intersect(
    names(regional_pool)[-1],
    names(regional_data)[-1]
)

stopifnot(
    identical(
        regional_pool[shared_taxa],
        regional_data[shared_taxa]
    ),
    identical(
        names(local_data),
        names(regional_data)
    ),
    identical(
        as.character(local_data$pollenZone),
        as.character(regional_data$pollenZone)
    ),
    all(regional_data$Phragmites_australis == 0)
)

stopifnot(
    identical(
        names(local_data)[-1],
        params_run$ppes$Taxon
    ),
    identical(
        params_run$ppes$Taxon,
        params_run$fallspeed$Taxon
    )
)

# LOVE run ####

library(disqover)

## Set seed for reproducibility

set.seed(777)

## Pass each time interval in local-regional couples

for (i in 1:nrow(local_data)) {
    loc_input <- cbind(age = as.numeric(i), local_data[i, -1])
    reg_input <- cbind(age = as.numeric(i), regional_data[i, -1])

    bin <- LOVEoptim(
        pollenLoc = loc_input,
        pollenReg = reg_input,
        lParams = params_run,
        dfRadii = rings,
        maxRingLove = 40
    )

    # Save results of the current time interval

    binc <- as.data.frame(bin$cover) * 100

    binc$pollenZone <- zones$pollenZone[i]
    binc$innerRadius <- rings$inner[-41]
    binc$outerRadius <- rings$outer[-41]

    rio::export(
        binc,
        paste0(
            output_dir,
            "/",
            record,
            "_lra_",
            zones$pollenZone[i],
            "_",
            Sys.Date(),
            ".xlsx"
        )
    )
}

## Optional: fuse all the age interval outputs in a single data frame

setwd(output_dir)

names <- paste0(
    record,
    "_lra_",
    zones$pollenZone,
    "_",
    Sys.Date(),
    ".xlsx"
)

for (j in seq_along(names)) {
    lraresult <- rio::import(names[j])
    lraresult$ageInterval <- local_data$pollenZone[j]
    lraresult$innerRadius <- rings$inner[-41] # Remove last ring [not modelled]
    lraresult$outerRadius <- rings$outer[-41] # Remove last ring [not modelled]

    if (j == 1) {
        output <- bind_rows(lraresult)
    } else {
        output <- bind_rows(
            output,
            lraresult
        )
    }
}

# Data frame with modelled land cover means per age interval for each ring #####

output_clean <- output |>
    select(c(
        ageInterval,
        innerRadius,
        outerRadius,
        starts_with("mean.") # Remove land cover quantile data, keep only means
    )) |>
    filter(outerRadius > 100) # Remove unmodelled radii (0-100 m) rows (empty)

rio::export(
    output_clean,
    paste0(record, "_modelled_landcover_ring_means_", Sys.Date(), ".xlsx")
)

## Keep only the biggest radius

output_biggest_radius <- output_clean |>
    filter(outerRadius == max(outerRadius))

rio::export(
    output_biggest_radius,
    paste0(
        record,
        "_lra_full_biggest_radius_",
        Sys.Date(),
        ".xlsx"
    )
)

# Mean percentages of modelled 7210 indicator taxa for each age interval ####

target_taxa <- c(
    "mean.Cladium_mariscus",
    "mean.Potentilla",
    "mean.Filipendula",
    "mean.Alnus",
    "mean.Fraxinus_ornus",
    "mean.Salix",
    "mean.Phragmites_australis",
    "mean.Carex",
    "mean.Galium_type",
    "mean.Scirpus_sl",
    "mean.Peucedanum_palustre_group"
)

output_long <- output_clean |>
    select(ageInterval, all_of(target_taxa)) |>
    pivot_longer(
        cols = starts_with("mean."),
        names_to = "Taxon",
        values_to = "LandCover"
    ) |>
    mutate(
        indicatorTaxon = gsub("mean\\.", "", Taxon),
        indicatorTaxon = gsub("_", " ", indicatorTaxon),
        indicatorTaxon = ifelse(indicatorTaxon == "Scirpus sl", "Scirpus s.l.", indicatorTaxon),
        isCladium = ifelse(Taxon == "mean.Cladium_mariscus", "Cladium", "Other")
    )

output_long$ageInterval <- factor(
    output_long$ageInterval,
    c(
        "1000-400",
        "1500-1000",
        "2000-1500",
        "2500-2000",
        "3000-2500",
        "3500-3000",
        "4000-3500",
        "4200-4000"
    )
)

output_indicators <- output_long |>
    group_by(ageInterval, indicatorTaxon) |>
    summarise(
        meanLandCover = round(mean(LandCover, na.rm = TRUE), digits = 1), # Round percentages to 1d decimal place
        .groups = "drop"
    )

rio::export(output_indicators, paste0("modelled_landcover_fullrange_indicators_means_", Sys.Date(), ".xlsx"))

# Data frame with modelled land cover means per age interval (largest ring) #####

output_means <- output_clean |>
    filter(outerRadius == 2000) |>
    group_by(ageInterval)

rio::export(output_means, paste0("modelled_landcover_2km_ring_", Sys.Date(), ".xlsx"))

# Run the chunk before to get the output_indicators data frame, then run this chunk to export it as an Excel file
rio::export(output_indicators, paste0("modelled_landcover_2km_ring_indicators_", Sys.Date(), ".xlsx"))

## Checks

# Taxa nei tre oggetti
taxa_loc <- colnames(local_data)[-1]
taxa_reg <- colnames(regional_data)[-1]
taxa_params <- params_run$ppes$Taxon

# Confronto
setdiff(taxa_loc, taxa_params)
setdiff(taxa_params, taxa_loc)

setdiff(taxa_loc, taxa_reg)
setdiff(taxa_reg, taxa_loc)

# Numero di taxa
c(
    local = length(taxa_loc),
    regional = length(taxa_reg),
    parameters = length(taxa_params)
)
