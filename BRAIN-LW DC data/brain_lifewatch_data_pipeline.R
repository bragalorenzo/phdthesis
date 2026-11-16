###############################################
### --- BRAIN to LifeWatch data pipeline ---###
###############################################

# Preparation ####

library(tidyverse)
library(taxize)
library(rgbif)
library(parzer)

# Batch processing of datasets ####

# Load list of datasets

brainsites <- rio::import(paste0("./brainsites_", Sys.Date(), ".xlsx"))
datasetlist <- rio::import(paste0("./Tank/datasetlist_", Sys.Date(), ".xlsx"))
brainbackbone <- rio::import("brain_taxonomic_backbone.xlsx")
# Write backbone backup
rio::export(brainbackbone, "brain_taxonomic_backbone_backup.xlsx")

## Custom functions

### Extract list of taxa from dataset

get_brain_taxa <- function(x) {
  stopifnot(is.data.frame(x))
  # Signpost variables
  temp_beginning <- which(x[, 1] == "Trees, shrubs and lianas") + 1
  temp_end <- which(x[, 1] == "Indeterminables") - 1
  # Extract taxa from dataset
  taxalist <- x[temp_beginning:temp_end, 1] |>
    as.data.frame()
  # Rename column
  colnames(taxalist) <- "verbatimIdentification"
  # Remove "Herbs" and NAs
  taxalist_clean <- taxalist |>
    filter(
      verbatimIdentification != "Herbs",
      !is.na(verbatimIdentification)
    )
  # Return
  return(taxalist_clean)
}

### Scientific name normalisation

sciname_norm <- function(x) {
  x <- as.character(x)
  x <- tolower(x) # case-insensitive
  x <- gsub("\\s+", " ", x) # spazi multipli -> uno solo
  x <- trimws(x) # toglie spazi iniziali/finali
  x
}

### Get taxonomic rank name from classification list

get_rank_name <- function(class_list, target_rank) {
  if (is.null(class_list)) {
    return(NA)
  }
  for (item in class_list) {
    if (!is.null(item$rank) && toupper(item$rank) == toupper(target_rank)) {
      if (!is.null(item$name)) {
        return(item$name)
      }
    }
  }
  return(NA)
}

### Get taxonomic field from usage object

get_usage_field <- function(usage_obj, field) {
  if (!is.null(usage_obj) && !is.null(usage_obj[[field]])) {
    return(usage_obj[[field]])
  } else {
    return(NA)
  }
}

### Type, group, etc. keyword removal

keyword_pattern <- "(?i)\\b(type|group|deciduous)\\b"

## Batch process

for (i in seq_len(nrow(datasetlist))) {
  # Read the updated backbone for each iteration
  brainbackbone <- rio::import("brain_taxonomic_backbone.xlsx")

  # Check if backbone has duplicated taxa
  # If yes, stop the process and print the duplicated taxa
  duplicate_taxa <- brainbackbone |>
    mutate(scient_norm = sciname_norm(verbatimIdentification)) |>
    count(scient_norm) |>
    filter(n > 1)

  if (nrow(duplicate_taxa) > 0) {
    stop(
      "Duplicated taxa found in BRAIN backbone: ",
      paste(duplicate_taxa$scient_norm, collapse = ", ")
    )
  }

  # Import the relevant data sheet
  if (datasetlist$sheetNames[i] == "default") {
    ## TRUE: import the default sheet (first sheet) of the Excel file
    curr_dataset_raw <- rio::import(paste0(
      "./Tank/", datasetlist$filename[i]
    ))
  } else {
    ## FALSE: more than one data sheet. Import the specified sheet of the Excel file
    curr_dataset_raw <- rio::import(
      paste0(
        "./Tank/", datasetlist$filename[i]
      ),
      sheet = datasetlist$sheetNames[i]
    )
  }

  message(paste("Processing dataset ID:", datasetlist$brainId[i]))
  message(paste("Name:", datasetlist$filename[i]))
  message(paste(i, "/", max(datasetlist$processNo)))

  # Extract the list of taxa
  curr_dataset <- get_brain_taxa(curr_dataset_raw)

  # Extract site metadata

  curr_metadata <- brainsites[
    tolower(brainsites$ID) %in% datasetlist$brainId[i],
  ]

  # Initialise the processed dataset's dataframe
  curr_dataset_processed <- data.frame(
    catalogNumber = seq_len(nrow(curr_dataset)),
    decimalLatitude = parzer::parse_lat(
      curr_metadata$Latitude
    ),
    decimalLongitude = parzer::parse_lon(
      curr_metadata$Longitude
    ),
    verbatimLocality = curr_metadata$Name,
    locality = datasetlist$locality[i],
    stateProvince = curr_metadata$Region,
    county = curr_metadata$Province,
    eventID = curr_metadata$ID,
    eventDate = Sys.Date(),
    eventRemarks = "null",
    verbatimIdentification = curr_dataset$verbatimIdentification,
    scientificName = NA,
    scientificNameAuthorship = NA,
    acceptedNameUsage = NA,
    acceptedNameUsageID = NA,
    kingdom = "Plantae",
    phylum = NA,
    class = NA,
    order = NA,
    family = NA,
    genus = NA,
    identificationRemarks = NA,
    FossilSpecimen = "Fossil pollen",
    taxonRemarks = NA,
    isBRAIN = 0
  )

  # Populate taxonomic data columns
  for (j in seq_len(nrow(curr_dataset_processed))) {
    ## Step 1: check if the taxon contains "type", "group", etc.
    if (grepl(
      keyword_pattern,
      curr_dataset_processed$verbatimIdentification[j],
      perl = TRUE
    )) {
      ## Extract corresponding keyword(s)
      matched_keywords <- regmatches(
        curr_dataset_processed$verbatimIdentification[j],
        gregexpr(keyword_pattern,
          curr_dataset_processed$verbatimIdentification[j],
          perl = TRUE
        )
      )[[1]]
      ## Add the keyword(s) in identificationRemarks
      curr_dataset_processed[j, "identificationRemarks"] <- paste(
        unique(matched_keywords),
        collapse = ", "
      )
      ## Remove the keyword(s) from verbatimIdentification and clean whitespaces
      curr_dataset_processed[j, "verbatimIdentification"] <- trimws(
        gsub(
          "\\s+", " ",
          sub(keyword_pattern, "",
            curr_dataset_processed$verbatimIdentification[j],
            perl = TRUE
          )
        )
      )
    } else {
      ## No keywords found, leave everything untouched
      curr_dataset_processed[j, "identificationRemarks"] <- NA
    }
    ## Check if verbatim taxon is already in backbone
    message(paste(
      "Processing taxon:",
      curr_dataset_processed$verbatimIdentification[j]
    ))
    if (
      sciname_norm(curr_dataset_processed$verbatimIdentification[j]) %in%
        sciname_norm(brainbackbone$verbatimIdentification)) {
      ## TRUE: populate dataset with taxon data from backbone
      message(paste(
        "Taxon",
        curr_dataset_processed$verbatimIdentification[j],
        "is in BRAIN backbone. Populating..."
      ))
      curr_dataset_processed$isBRAIN[j] <- 1 # It's in BRAIN backbone
      ### Extract data from backbone
      taxondata_brainbackbone <- brainbackbone[
        which(sciname_norm(brainbackbone$verbatimIdentification) ==
          sciname_norm(curr_dataset_processed$verbatimIdentification)[j]),
        c(
          "scientificName",
          "scientificNameAuthorship",
          "acceptedNameUsage",
          "acceptedNameUsageID",
          "kingdom",
          "phylum",
          "class",
          "order",
          "family",
          "genus",
          "taxonRemarks"
        )
      ]

      ## Check if the backbone lookup returned duplicates

      if (nrow(taxondata_brainbackbone) != 1) {
        stop(
          "Backbone lookup for '",
          curr_dataset_processed$verbatimIdentification[j],
          "' returned ",
          nrow(taxondata_brainbackbone),
          " rows."
        )
      }

      ### Populate taxonomic dataset fields with the backbone data
      curr_dataset_processed[j, c(
        "scientificName",
        "scientificNameAuthorship",
        "acceptedNameUsage",
        "acceptedNameUsageID",
        "kingdom",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "taxonRemarks"
      )] <- taxondata_brainbackbone
    } else {
      ## FALSE: retrieve the data from GBIF
      message(paste(
        "Taxon",
        curr_dataset_processed$verbatimIdentification[j],
        "is not in BRAIN backbone. Retrieving data..."
      ))
      curr_dataset_processed$isBRAIN[j] <- 0
      taxondata_gbif <- taxize::get_gbifid_(
        curr_dataset_processed$verbatimIdentification[j]
      ) |>
        as.data.frame()
      ### Convert it to CoL
      taxondata_col <- taxondata_gbif[1, ] |>
        select(ends_with(".usagekey")) |>
        rgbif::gbif_to_col()
      ### Check if a taxon was successfully retrieved
      if (taxondata_col$diagnostics$matchType == "NONE") {
        ## TRUE: no taxon was found in CoL. Leave the fields as NA and print a message
        message("Taxon was not found. Check manually.")
      } else {
        ## FALSE: a taxon was found in CoL. Populate the fields with the retrieved data
        curr_dataset_processed[j, "scientificName"] <- get_usage_field(taxondata_col$usage, "canonicalName")
        curr_dataset_processed[j, "scientificNameAuthorship"] <- get_usage_field(taxondata_col$usage, "authorship")
        curr_dataset_processed[j, "acceptedNameUsage"] <- get_usage_field(taxondata_col$usage, "canonicalName")
        key_val <- get_usage_field(taxondata_col$usage, "key")
        ### Check if the key value is not NA before constructing the URL
        curr_dataset_processed[j, "acceptedNameUsageID"] <- if (!is.na(key_val)) {
          ### TRUE: key value is not NA, construct the URL
          paste0("https://www.catalogueoflife.org/data/taxon/", key_val)
        } else {
          "invalid_key"
        }
        curr_dataset_processed[j, "kingdom"] <- get_rank_name(taxondata_col$classification, "kingdom")
        curr_dataset_processed[j, "phylum"] <- get_rank_name(taxondata_col$classification, "phylum")
        curr_dataset_processed[j, "class"] <- get_rank_name(taxondata_col$classification, "class")
        curr_dataset_processed[j, "order"] <- get_rank_name(taxondata_col$classification, "order")
        curr_dataset_processed[j, "family"] <- get_rank_name(taxondata_col$classification, "family")
        curr_dataset_processed[j, "genus"] <- get_rank_name(taxondata_col$classification, "genus")
      }

      ### Then, update the backbone with the new taxon data

      newtaxon_backbone <- data.frame(
        verbatimIdentification = curr_dataset_processed$verbatimIdentification[j],
        scientificName = curr_dataset_processed$scientificName[j],
        scientificNameAuthorship = curr_dataset_processed$scientificNameAuthorship[j],
        acceptedNameUsage = curr_dataset_processed$acceptedNameUsage[j],
        acceptedNameUsageID = curr_dataset_processed$acceptedNameUsageID[j],
        kingdom = get_rank_name(taxondata_col$classification, "kingdom"),
        phylum = get_rank_name(taxondata_col$classification, "phylum"),
        class = get_rank_name(taxondata_col$classification, "class"),
        order = get_rank_name(taxondata_col$classification, "order"),
        family = get_rank_name(taxondata_col$classification, "family"),
        genus = get_rank_name(taxondata_col$classification, "genus"),
        taxonRemarks = curr_dataset_processed$taxonRemarks[j]
      )

      #### Append the new taxon to the backbone & overwrite the backbone file
      brainbackbone <- bind_rows(brainbackbone, newtaxon_backbone)

      rio::export(
        brainbackbone,
        "brain_taxonomic_backbone.xlsx"
      )
    }

    ## Get new species & export info file on them
    newtaxa <- curr_dataset_processed |>
      filter(isBRAIN == 0)
    rio::export(
      newtaxa,
      paste0("./Backbone Updates/", curr_metadata$ID, "_newtaxa", Sys.Date(), ".xlsx")
    )
    ## Export processed dataset
    rio::export(
      curr_dataset_processed,
      paste0("./Processed/", curr_metadata$ID, "_lwprocessed", Sys.Date(), ".xlsx")
    )
  }
}

# New taxa added to backbone #####

backup <- rio::import("brain_taxonomic_backbone_backup.xlsx")

newtaxa <- anti_join(
  x = brainbackbone,
  y = backup,
  by = "verbatimIdentification"
)

rio::export(newtaxa, paste0("./Backbone Updates/newtaxa_", Sys.Date(), ".xlsx"))
