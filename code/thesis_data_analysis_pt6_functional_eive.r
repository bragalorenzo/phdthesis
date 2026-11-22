#######################
#--- EIVE matching ---#
#######################

library(tidyverse)

# Load data

eive_data <- rio::import("G:/My Drive/Work/PhD Thesis/Data Analysis/EIVE_Paper_1.0_SM_08.xlsx", sheet = "mainTable")

record <- "cpz2"

record_dir <- paste0(
    "G:/My Drive/Work/PhD Thesis/Data Analysis/", toupper(record), " Data"
)

proto_matches <- rio::import(paste0(
    record_dir, "/", record, "_proto_match_results.csv"
))

# Match #####

## Extract individual flora species from PROTO
proto_species <- proto_matches |>
    filter(Table == "Flora species score") |>
    select(
        speciesName = Name,
        gbifID = `GBIF ID`,
        pollenMatches = `Matched with`,
        Matches
    ) |>
    distinct()

## Check
glimpse(proto_species)
n_distinct(proto_species$speciesName)

## Normalise scientific names to genus + specific epithet
normalise_binomial <- function(x) {
    x |>
        str_replace_all("×", "") |>
        str_squish() |>
        str_extract("^[A-Z][A-Za-z-]+\\s+[a-z][a-z-]+") |>
        str_to_lower()
}

## Prepare EIVE species
eive_species <- eive_data |>
    filter(TaxonRank == "Species") |>
    mutate(
        binomial = normalise_binomial(TaxonConcept)
    ) |>
    filter(!is.na(binomial)) |>
    add_count(binomial, name = "nCandidates")

## Prepare PROTO species
proto_species <- proto_species |>
    mutate(
        binomial = normalise_binomial(speciesName)
    )

## Match PROTO species with EIVE
proto_eive <- proto_species |>
    left_join(
        eive_species,
        by = "binomial",
        relationship = "many-to-many"
    ) |>
    mutate(
        matchStatus = case_when(
            is.na(binomial) ~ "invalid_name",
            is.na(TaxonConcept) ~ "unmatched",
            nCandidates > 1 ~ "ambiguous",
            str_detect(
                speciesName,
                regex("\\b(subsp|ssp|var)\\.", ignore_case = TRUE)
            ) ~ "infraspecific_review",
            TRUE ~ "matched"
        )
    )

## Matching summary

proto_eive |>
    distinct(gbifID, speciesName, matchStatus) |>
    count(matchStatus) |>
    mutate(
        percentage = round(100 * n / n_distinct(proto_species$gbifID), 1)
    )

## Unmatched species
proto_eive |>
    filter(matchStatus == "unmatched") |>
    select(speciesName, gbifID, binomial) |>
    distinct() |>
    arrange(speciesName)

## Ambiguous matches
proto_eive |>
    filter(matchStatus %in% c(
        "ambiguous",
        "infraspecific_review"
    )) |>
    select(speciesName, TaxonConcept, matchStatus) |>
    distinct()

# Check known synonyms in EIVE
synonyms_to_check <- c(
    "anagallis arvensis",
    "myosoton aquaticum",
    "sedum compressum"
)

eive_data |>
    mutate(
        binomial = normalise_binomial(TaxonConcept)
    ) |>
    filter(binomial %in% synonyms_to_check) |>
    select(
        TaxonConcept,
        TaxonRank,
        starts_with("EIVEres-")
    )

# Correct matches #####

# Prepare EIVE
eive_species <- eive_data |>
    filter(
        TaxonRank == "Species",
        !str_detect(TaxonConcept, fixed("×"))
    ) |>
    mutate(
        binomial_eive = normalise_binomial(TaxonConcept)
    ) |>
    filter(!is.na(binomial_eive)) |>
    add_count(binomial_eive, name = "nCandidates")

# Prepare PROTO species and resolve synonyms (INF1)
proto_eive <- proto_species |>
    select(speciesName, gbifID, pollenMatches, Matches) |>
    mutate(
        binomial_proto = normalise_binomial(speciesName),
        binomial_eive = case_when(
            binomial_proto == "lysimachia arvensis" ~
                "anagallis arvensis",
            binomial_proto == "stellaria aquatica" ~
                "myosoton aquaticum",
            TRUE ~ binomial_proto
        )
    ) |>
    left_join(
        eive_species,
        by = "binomial_eive",
        relationship = "many-to-many"
    ) |>
    mutate(
        matchStatus = case_when(
            is.na(TaxonConcept) ~ "unmatched",
            nCandidates > 1 ~ "ambiguous",
            binomial_proto != binomial_eive ~ "synonym",
            TRUE ~ "exact"
        )
    )

# cpz2 synonyms

candidates <- tibble::tribble(
    ~speciesName, ~candidate,
    "Aetheorhiza bulbosa", "Sonchus bulbosus",
    "Centaurea cyanus", "Cyanus segetum",
    "Daucus broteroi", "Daucus broteri",
    "Dianthus virgineus", "Dianthus sylvestris subsp. longicaulis",
    "Lotus hirsutus", "Dorycnium hirsutum",
    "Salsola squarrosa", "Salsola controversa",
    "Senecio vernalis", "Senecio leucanthemifolius subsp. vernalis",
    "Soda inermis", "Salsola soda",
    "Visnaga daucoides", "Ammi visnaga"
)

eive_names <- tolower(as.character(eive_data$TaxonConcept))

candidates$in_eive <- vapply(
    candidates$candidate,
    function(x) {
        any(startsWith(eive_names, tolower(x)))
    },
    logical(1)
)

print(candidates, n = Inf)

proto_eive <- proto_species |>
    select(speciesName, gbifID, pollenMatches, Matches) |>
    mutate(
        binomial_proto = normalise_binomial(speciesName),
        binomial_eive = case_when(
            binomial_proto == "aetheorhiza bulbosa" ~
                "sonchus bulbosus",
            binomial_proto == "centaurea cyanus" ~
                "cyanus segetum",
            binomial_proto == "daucus broteroi" ~
                "daucus broteri",
            binomial_proto == "dianthus virgineus" ~
                "dianthus sylvestris",
            binomial_proto == "lotus hirsutus" ~
                "dorycnium hirsutum",
            binomial_proto == "senecio vernalis" ~
                "senecio leucanthemifolius",
            binomial_proto == "soda inermis" ~
                "salsola soda",
            binomial_proto == "visnaga daucoides" ~
                "ammi visnaga",
            TRUE ~ binomial_proto
        )
    ) |>
    left_join(
        eive_species,
        by = c("binomial_eive" = "binomial"),
        relationship = "many-to-many"
    ) |>
    mutate(
        matchStatus = case_when(
            is.na(TaxonConcept) ~ "unmatched",
            nCandidates > 1 ~ "ambiguous",
            binomial_proto != binomial_eive ~ "synonym",
            TRUE ~ "exact"
        )
    )

# Check results
proto_eive |>
    distinct(gbifID, speciesName, matchStatus) |>
    count(matchStatus)

# Check for unexpected duplicates
proto_eive |>
    count(gbifID) |>
    filter(n > 1)

proto_eive |>
    filter(binomial_proto %in% c(
        "dianthus virgineus",
        "senecio vernalis"
    )) |>
    distinct(speciesName, TaxonConcept, nCandidates, matchStatus)

# Value coverage by EIVE #####

# EIVE indicators
indicator_cols <- paste0(
    "EIVEres-", c("M", "N", "R", "L", "T")
)

# Coverage by indicator
indicator_coverage <- proto_eive |>
    select(gbifID, all_of(indicator_cols)) |>
    pivot_longer(
        cols = all_of(indicator_cols),
        names_to = "indicator",
        values_to = "value"
    ) |>
    group_by(indicator) |>
    summarise(
        nAvailable = sum(!is.na(value)),
        nMissing = sum(is.na(value)),
        coveragePct = round(
            100 * nAvailable / n(), 2
        ),
        .groups = "drop"
    )

print(indicator_coverage)

proto_eive |>
    filter(if_any(
        all_of(indicator_cols),
        is.na
    )) |>
    select(
        speciesName,
        TaxonConcept,
        matchStatus,
        all_of(indicator_cols)
    )

# Create pollen type x EIVE indicator table #####

# EIVE indicators
indicator_cols <- paste0(
    "EIVEres-", c("M", "N", "R", "L", "T")
)

# Original pollen taxa from PROTO
pollen_index <- proto_matches |>
    filter(Table == "Pollen taxon score") |>
    transmute(
        pollenType = str_squish(Type),
        pollenTaxon = Name,
        pollenID = `GBIF ID`
    ) |>
    distinct()

# Expand species-pollen associations
species_links <- proto_eive |>
    select(
        speciesName,
        gbifID,
        TaxonConcept,
        matchStatus,
        pollenMatches,
        all_of(indicator_cols)
    ) |>
    separate_rows(
        pollenMatches,
        sep = "\\s*\\|\\s*"
    ) |>
    rename(linkLabel = pollenMatches) |>
    mutate(
        pollenType = str_squish(
            str_match(
                linkLabel,
                "\\(([^,()]+),\\s*GBIF\\s*\\d+\\)\\s*$"
            )[, 2]
        ),
        pollenID = as.integer(
            str_match(
                linkLabel,
                "GBIF\\s*(\\d+)\\)\\s*$"
            )[, 2]
        )
    )

# Unparsed links
unparsed <- species_links |>
    filter(is.na(pollenType) | is.na(pollenID))

# Non-unique identifiers within pollen types
duplicate_keys <- pollen_index |>
    count(pollenType, pollenID) |>
    filter(n > 1)

# Links without a corresponding pollen taxon
unlinked <- species_links |>
    anti_join(
        pollen_index,
        by = c("pollenType", "pollenID")
    )

# Check parsed associations against PROTO counts
pair_check <- species_links |>
    count(gbifID, name = "nParsed") |>
    full_join(
        proto_species |>
            select(gbifID, Matches),
        by = "gbifID"
    ) |>
    filter(
        coalesce(nParsed, 0L) != coalesce(Matches, 0L)
    )

# Diagnostics
tibble(
    check = c(
        "Unparsed links",
        "Duplicate pollen keys",
        "Unlinked associations",
        "Species with inconsistent pair counts"
    ),
    n = c(
        nrow(unparsed),
        nrow(duplicate_keys),
        nrow(unlinked),
        nrow(pair_check)
    )
)

stopifnot(
    nrow(unparsed) == 0,
    nrow(duplicate_keys) == 0,
    nrow(unlinked) == 0,
    nrow(pair_check) == 0
)

# Create pollen-species table

pollen_species <- species_links |>
    left_join(
        pollen_index,
        by = c("pollenType", "pollenID"),
        relationship = "many-to-one"
    ) |>
    select(
        pollenTaxon,
        pollenType,
        pollenID,
        speciesName,
        gbifID,
        TaxonConcept,
        matchStatus,
        all_of(indicator_cols)
    )

# Ensure each pollen-species pair is unique
duplicate_pairs <- pollen_species |>
    count(pollenType, pollenID, gbifID) |>
    filter(n > 1)

stopifnot(nrow(duplicate_pairs) == 0)

# Check dataset
glimpse(pollen_species)

# Long format + indicator dispersions

pollen_eive_long <- pollen_species |>
    pivot_longer(
        cols = all_of(indicator_cols),
        names_to = "indicator",
        values_to = "value"
    ) |>
    mutate(
        indicator = str_remove(
            indicator, "^EIVEres-"
        )
    ) |>
    group_by(
        pollenTaxon,
        pollenType,
        pollenID,
        indicator
    ) |>
    summarise(
        nSpecies = n_distinct(gbifID),
        nEIVE = sum(!is.na(value)),
        median = if (all(is.na(value))) {
            NA_real_
        } else {
            median(value, na.rm = TRUE)
        },
        IQR = if (all(is.na(value))) {
            NA_real_
        } else {
            IQR(value, na.rm = TRUE)
        },
        min = if (all(is.na(value))) {
            NA_real_
        } else {
            min(value, na.rm = TRUE)
        },
        max = if (all(is.na(value))) {
            NA_real_
        } else {
            max(value, na.rm = TRUE)
        },
        .groups = "drop"
    )

# Create final matrix of pollen taxa x EIVE indicators

pollen_eive_matrix <- pollen_eive_long |>
    pivot_wider(
        id_cols = c(
            pollenTaxon,
            pollenType,
            pollenID,
            nSpecies
        ),
        names_from = indicator,
        values_from = c(
            nEIVE,
            median,
            IQR,
            min,
            max
        ),
        names_glue = "{indicator}_{.value}"
    ) |>
    arrange(pollenType, pollenTaxon)

glimpse(pollen_eive_matrix)

# Export
rio::export(
    pollen_species,
    paste0(
        record_dir, "/",
        record, "_eive_species_links.csv"
    )
)

rio::export(
    pollen_eive_matrix,
    paste0(
        record_dir, "/",
        record, "_eive_pollen_taxa.csv"
    )
)

# Results #####

# Number of species per pollen taxon
pollen_eive_matrix |>
    count(nSpecies) |>
    arrange(nSpecies)

# Pollen taxa with the greatest ecological heterogeneity
pollen_eive_matrix |>
    arrange(desc(M_IQR)) |>
    select(
        pollenTaxon,
        pollenType,
        nSpecies,
        M_median,
        M_IQR,
        M_min,
        M_max
    ) |>
    print(n = 30)

# Weighting of ecological values by species and pollen taxon #####

## Function to weight EIVE values by species and pollen taxon #####
eive.weight <- function(eive, hellinger) {
    h <- hellinger[!is.na(hellinger) & hellinger > 0]

    if (length(h) == 0L) {
        return(eive * 0)
    }

    eive * mean(h)
}

# Import count data & prepare for downstream zone analysis

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

## Data #####

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

## Start processing data for EIVE weighting #####

library(vegan)

counts_hellinger <- counts_clean |>
    vegan::decostand(method = "hellinger")

# Check: Unique pollen taxa associated with species
proto_taxa <- pollen_species |>
    distinct(
        pollenTaxon,
        pollenType,
        pollenID
    )

# Pollen taxa in PROTO but not in counts
unmatched_taxa <- proto_taxa |>
    filter(
        !pollenTaxon %in% colnames(counts_hellinger)
    )

print(unmatched_taxa, n = Inf)

# Prepare pollen taxon correspondence

taxon_map <- pollen_species |>
    distinct(pollenTaxon, pollenType, pollenID) |>
    mutate(
        terminalName = str_squish(
            str_split_i(pollenTaxon, "\\s*>\\s*", -1)
        ),
        countTaxon = case_when(
            pollenTaxon %in% colnames(counts_hellinger) ~
                pollenTaxon,
            pollenType == "Genus" ~
                str_extract(terminalName, "^[A-Z][A-Za-z-]+"),
            pollenType == "Species" ~
                str_extract(
                    terminalName,
                    "^[A-Z][A-Za-z-]+\\s+[a-z][a-z-]+"
                ),
            pollenType == "Family" ~
                paste0(terminalName, " undiff."),
            pollenType == "Punt" ~
                paste0(terminalName, " P"),
            pollenType == "Moore" ~
                paste0(terminalName, " M"),
            pollenType == "Beug" ~
                str_replace(terminalName, "-Typ$", "-Typ B"),
            pollenType == "Florenzano" ~
                str_remove(terminalName, "\\s*\\(undiff\\.\\)$"),
            TRUE ~ terminalName
        )
    ) |>
    mutate(
        # Handle undifferentiated genera
        countTaxon = if_else(
            pollenType == "Genus" &
                !countTaxon %in% colnames(counts_hellinger) &
                paste0(countTaxon, " undiff.") %in%
                    colnames(counts_hellinger),
            paste0(countTaxon, " undiff."),
            countTaxon
        )
    ) |>
    select(
        pollenTaxon,
        pollenType,
        pollenID,
        countTaxon
    )

unresolved_taxa <- taxon_map |>
    filter(
        !countTaxon %in% colnames(counts_hellinger)
    )

print(unresolved_taxa, n = Inf)

## Correct

# Manual nomenclatural reconciliation

taxon_map <- taxon_map |>
    mutate(
        countTaxon = case_when(
            pollenTaxon == "Sinapis type" ~
                "Brassica type (Sinapis type M)",
            pollenTaxon == "Sorbus-Gruppe" ~
                "Rosaceae woody undiff.",
            pollenTaxon == "Amaranthaceae" ~
                "Chenopodiaceae undiff.",
            str_detect(pollenTaxon, fixed("Quercus robur-pubescens-Typ")) ~
                "Deciduous Quercus",
            pollenTaxon == "Quercus  >  Quercus cerris-Typ" ~
                "Deciduous Quercus",
            pollenTaxon == "Ostrya-Typ" ~
                "Ostrya carpinifolia/Carpinus orientalis",
            pollenTaxon == "Lysimachia Tourn. ex L." ~
                "Anagallis",
            TRUE ~ countTaxon
        )
    )

unresolved_taxa <- taxon_map |>
    filter(!countTaxon %in% colnames(counts_hellinger))

print(unresolved_taxa, n = Inf)

# Match #####

# Indicator columns
indicator_cols <- paste0(
    "EIVEres-", c("M", "N", "R", "L", "T")
)

# Hellinger matrix
h <- as.matrix(counts_hellinger)

# Sample identifiers, in the same order as counts_hellinger
samples <- rownames(counts_t)[counts_t$excludeSample == "no"]

stopifnot(nrow(h) == length(samples))

# Species-level EIVE values
species_eive <- proto_eive |>
    select(
        gbifID,
        speciesName,
        all_of(indicator_cols)
    ) |>
    distinct()

# Connect PROTO associations to pollen count names
species_links <- pollen_species |>
    distinct(
        gbifID,
        pollenTaxon,
        pollenType,
        pollenID
    ) |>
    left_join(
        taxon_map,
        by = c("pollenTaxon", "pollenType", "pollenID"),
        relationship = "many-to-one"
    ) |>
    distinct(gbifID, countTaxon)

stopifnot(
    !anyNA(species_links$countTaxon),
    all(species_links$countTaxon %in% colnames(h))
)

# Create weighted matrix

# Initialise output
weighted_cols <- paste0(indicator_cols, "_weighted")

weighted_matrix <- matrix(
    NA_real_,
    nrow = nrow(h) * nrow(species_eive),
    ncol = length(indicator_cols),
    dimnames = list(NULL, weighted_cols)
)

# Species x sample weighting
is_present <- integer(nrow(h) * nrow(species_eive))

for (s in seq_len(nrow(species_eive))) {
    # Retrieve pollen taxa associated with the species
    taxa <- species_links$countTaxon[
        species_links$gbifID == species_eive$gbifID[s]
    ]

    # Retrieve original EIVE values
    eive <- unlist(
        species_eive[s, indicator_cols],
        use.names = FALSE
    )

    for (i in seq_len(nrow(h))) {
        idx <- (s - 1) * nrow(h) + i

        h_values <- h[i, taxa]

        # Presence of at least one matched pollen taxon
        is_present[idx] <- as.integer(any(h_values > 0))

        # EIVE weighting
        weighted_matrix[idx, ] <- eive.weight(
            eive = eive,
            hellinger = h_values
        )
    }

    species_eive_weighted <- expand.grid(
        sample = samples,
        gbifID = species_eive$gbifID,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    ) |>
        left_join(
            species_eive,
            by = "gbifID",
            relationship = "many-to-one"
        ) |>
        mutate(isPresent = is_present) |>
        bind_cols(
            as_tibble(weighted_matrix)
        )
}

glimpse(species_eive_weighted)

# Export

rio::export(
    species_eive_weighted,
    paste0(
        record_dir, "/",
        record, "_eive_species_weighted_",
        Sys.Date(),
        ".xlsx"
    )
)

# Savestate

species_eive_weighted <- rio::import(
    paste0(
        record_dir, "/Outputs/", record, "_eive_species_weighted_",
        Sys.Date(),
        ".xlsx"
    )
)

# Viz #####

# Summarise weighted EIVE by sample

eive_temporal <- species_eive_weighted |>
    filter(isPresent == 1) |>
    pivot_longer(
        cols = ends_with("_weighted"),
        names_to = "indicator",
        values_to = "value"
    ) |>
    group_by(sample, indicator) |>
    summarise(
        median = median(value, na.rm = TRUE),
        q25 = quantile(value, 0.25, na.rm = TRUE),
        q75 = quantile(value, 0.75, na.rm = TRUE),
        .groups = "drop"
    ) |>
    left_join(
        metadata |>
            tibble::rownames_to_column("sample") |>
            select(sample, ageCalBP),
        by = "sample"
    ) |>
    mutate(
        indicator = stringr::str_remove(
            indicator, "^EIVEres-|_weighted$"
        ),
        ageCalBP = as.numeric(ageCalBP)
    )

# Whole-record weighted EIVE

ggplot(
    eive_temporal,
    aes(x = ageCalBP, y = median)
) +
    geom_ribbon(
        aes(ymin = q25, ymax = q75),
        fill = "grey70",
        alpha = 0.35
    ) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 1) +
    facet_wrap(
        ~indicator,
        ncol = 1,
        scales = "free_y"
    ) +
    scale_x_reverse() +
    labs(
        x = "Age (cal yr BP)",
        y = "Weighted EIVE contribution"
    ) +
    theme_classic()

# Unweighted EIVE by sample #####

eive_unweighted <- species_eive_weighted |>
    filter(isPresent == 1) |>
    pivot_longer(
        cols = all_of(indicator_cols),
        names_to = "indicator",
        values_to = "value"
    ) |>
    group_by(sample, indicator) |>
    summarise(
        median = median(value, na.rm = TRUE),
        q25 = quantile(value, 0.25, na.rm = TRUE),
        q75 = quantile(value, 0.75, na.rm = TRUE),
        .groups = "drop"
    ) |>
    left_join(
        eive_temporal |>
            distinct(sample, ageCalBP),
        by = "sample"
    )

eive_unweighted |>
    select(sample, indicator, median) |>
    pivot_wider(
        names_from = indicator,
        values_from = median
    ) |>
    select(-sample) |>
    cor(use = "pairwise.complete.obs") |>
    round(3)

# Indicator labels
eive_unweighted <- eive_unweighted |>
    mutate(
        indicator = factor(
            indicator,
            levels = paste0("EIVEres-", c("M", "N", "R", "L", "T")),
            labels = c(
                "Moisture (M)",
                "Nutrients (N)",
                "Soil reaction (R)",
                "Light (L)",
                "Temperature (T)"
            )
        )
    )

# Temporal EIVE plot
# Add pollen zones to the temporal dataset

eive_temporal_plot <- eive_unweighted |>
    left_join(
        sample_zones |>
            select(sample, pollenZone) |>
            distinct(),
        by = "sample"
    )

# Temporal EIVE plot

p_eive_temporal <- ggplot(
    eive_temporal_plot,
    aes(x = ageCalBP, y = median)
) +
    geom_ribbon(
        aes(ymin = q25, ymax = q75),
        fill = "grey75",
        alpha = 0.45
    ) +
    geom_line(
        colour = "grey35",
        linewidth = 0.5
    ) +
    geom_point(
        aes(colour = pollenZone),
        size = 1.8
    ) +
    scale_colour_manual(
        values = retro_palette,
        breaks = levels(eive_violin$pollenZone),
        drop = FALSE
    ) +
    facet_wrap(
        ~indicator,
        ncol = 1,
        scales = "free_y"
    ) +
    scale_x_reverse() +
    labs(
        x = "Age (cal yr BP)",
        y = "Median EIVE",
        colour = "Pollen zone"
    ) +
    theme_classic() +
    theme(
        strip.background = element_blank(),
        strip.text = element_text(
            face = "bold",
            size = 11
        ),
        panel.spacing = unit(0.9, "lines"),
        legend.position = "right"
    )

p_eive_temporal

svg(
    filename = paste0(
        "./Figures/",
        record,
        "_eive_temporal_",
        Sys.Date(),
        ".svg"
    ),
    width = 16, height = 9
)
p_eive_temporal
dev.off()

ggsave(
    filename = paste0(
        "./Figures/",
        record,
        "_eive_temporal_",
        Sys.Date(),
        ".png"
    ),
    plot = p_eive_temporal,
    width = 16, height = 9
)

# By zone #####

# Assign pollen zones to samples
sample_zones <- metadata |>
    tibble::rownames_to_column("sample") |>
    mutate(
        depthCm = as.numeric(depthCm),
        ageCalBP = as.numeric(ageCalBP)
    ) |>
    left_join(
        zones,
        by = join_by(
            depthCm >= topDepth,
            depthCm <= bottomDepth
        )
    )

# Prepare species-level EIVE by zone
eive_violin <- species_eive_weighted |>
    filter(isPresent == 1) |>
    inner_join(
        sample_zones |> select(sample, pollenZone),
        by = "sample"
    ) |>
    distinct(pollenZone, gbifID, .keep_all = TRUE) |>
    pivot_longer(
        cols = all_of(indicator_cols),
        names_to = "indicator",
        values_to = "EIVE"
    ) |>
    mutate(
        pollenZone = factor(pollenZone, levels = zone_levels),
        indicator = factor(
            indicator,
            levels = paste0("EIVEres-", c("M", "N", "R", "L", "T")),
            labels = c(
                "Moisture (M)",
                "Nutrients (N)",
                "Soil reaction (R)",
                "Light (L)",
                "Temperature (T)"
            )
        )
    )

# EIVE distributions by pollen zone
p_eive_violin <- ggplot(
    eive_violin_clean,
    aes(x = pollenZone, y = EIVE, fill = pollenZone)
) +
    geom_violin(
        trim = FALSE,
        alpha = 0.65,
        colour = "grey30"
    ) +
    geom_boxplot(
        width = 0.12,
        outlier.shape = NA,
        fill = "white",
        alpha = 0.8
    ) +
    stat_summary(
        fun = median,
        geom = "line",
        aes(group = 1),
        colour = "black",
        linewidth = 0.6
    ) +
    stat_summary(
        fun = median,
        geom = "point",
        colour = "black",
        size = 1.8
    ) +
    scale_fill_manual(values = retro_palette) +
    facet_wrap(
        ~indicator,
        ncol = 1
    ) +
    labs(
        x = "Pollen zone",
        y = "EIVE"
    ) +
    theme_classic() +
    theme(
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.spacing = unit(0.8, "lines")
    )

p_eive_violin

svg(
    filename = paste0(
        "./Figures/",
        record,
        "_eive_by_zone_violins_",
        Sys.Date(),
        ".svg"
    ),
    width = 12, height = 16
)
p_eive_violin
dev.off()

ggsave(
    filename = paste0(
        "./Figures/",
        record,
        "_eive_by_zone_violins_",
        Sys.Date(),
        ".png"
    ),
    plot = p_eive_violin,
    width = 12, height = 16
)

eive_zone_summary <- eive_violin |>
    filter(is.finite(EIVE), !is.na(pollenZone)) |>
    group_by(indicator, pollenZone) |>
    summarise(
        n = n(),
        median = median(EIVE),
        Q1 = quantile(EIVE, 0.25),
        Q3 = quantile(EIVE, 0.75),
        min = min(EIVE),
        max = max(EIVE),
        .groups = "drop"
    ) |>
    mutate(
        across(
            c(median, Q1, Q3, min, max),
            ~ round(.x, 2)
        )
    )

print(eive_zone_summary, n = Inf)

rio::export(
    eive_zone_summary,
    "CPZ2_EIVE_zone_summary.xlsx"
)
