# Numerical ecology of Holocene pollen records

**Reproducible R workflows for analysing long-term vegetation and biodiversity dynamics in two Italian pollen records: INF1 and CPZ2.**

This repository documents the data-processing and numerical-ecology workflows developed for the doctoral research project *Furthering the relationships between palaeoecology and biodiversity: numerical ecology methods and development of informatics tools in palaeopalynological case studies* (University of Modena and Reggio Emilia, Italy).

The analyses combine pollen-assemblage data, age–depth information, environmental reconstructions and taxonomic reference tables to investigate compositional change, temporal turnover, local vegetation dynamics and the relationship between palaeoecological records and contemporary biodiversity information.

> **Repository status:** Documentation scaffold. Analysis scripts, exact package versions and executable run instructions will be added and checked against the final research code. This repository should not yet be considered a fully reproducible release.

## Case studies

| Record | Study site | Approximate temporal coverage | Samples used in the numerical analyses | Pollen zones (oldest → youngest) |
| --- | --- | --- | ---: | --- |
| **INF1** | Teor–Modeano palaeochannel, Friuli Venezia Giulia, northeastern Italy | ca. 4,200–400 cal yr BP | 73 | MOD1a → MOD1b → MOD2 → MOD3 → MOD4a → MOD4b |
| **CPZ2** | Malafede Pool, Castelporziano Presidential Estate, central Italy | ca. 7,000–400 cal yr BP | 29 | MAL1a → MAL1b → MAL2 → MAL3 → MAL4 |

The records differ in depositional setting, temporal resolution and ecological context. Consequently, the analytical approaches are comparable at the methodological level, but their parameters and interpretation are **not assumed to be interchangeable**.

## Analytical scope

The repository is organised around the following components. Not every component is necessarily applied in the same form to both records; record-specific implementations and settings belong with their respective scripts.

| Component | Purpose | Key methods or considerations |
| --- | --- | --- |
| Data preparation and harmonisation | Validate counts, metadata, chronologies and taxonomic identifiers | Sample filtering, chronological alignment and taxon lookup tables |
| Pollen stratigraphy and zonation | Describe major compositional changes | Pollen diagrams and stratigraphically constrained clustering (CONISS) |
| Palynological richness | Compare sample-level richness while accounting for count size | Rarefaction |
| Rate of Change (RoC) | Quantify temporal compositional turnover | Chord distance; moving-window settings specific to each record |
| Community Trajectory Analysis (CTA) | Characterise assemblage trajectories and compositional heterogeneity | Principal Coordinates Analysis (PCoA), sample-wise contributions to variance and within-zone variation |
| Climate–composition relationships | Explore the association between pollen composition and palaeoclimate | Redundancy Analysis (RDA), including temperature and summer precipitation |
| Temporal structure and nonlinear trends | Assess temporal patterns in assemblages or individual taxa | Distance-based Moran's Eigenvector Maps (dbMEM) and Generalised Additive Models (GAMs), where applied |
| Ecological indicator analyses | Explore ecological information associated with floristically matched pollen taxa | Ellenberg-like indicator values and explicit taxonomic matching assumptions, where applied |
| Landscape Reconstruction Algorithm (LRA) | Estimate local vegetation/land-cover composition | REVEALS/LOVE-related workflows, pollen productivity estimates and spatial buffers |

### Record-specific settings worth preserving

- **INF1:** the primary RoC analysis uses a 100-year moving window, with sensitivity analyses at alternative window widths. Its final sample set excludes eight subsamples that did not meet the minimum pollen-count threshold.
- **CPZ2:** the primary RoC analysis uses a 200-year moving window, reflecting its coarser temporal resolution.
- **Both records:** the LRA workflow uses **unweighted rings** and retains estimates for the **largest modelled ring** for zone-level visualisations. The exact inputs, taxonomic exclusions and model parameters must be verified in the corresponding scripts.
- **Chronology:** pollen-zone labels must be interpreted in the *oldest-to-youngest* order given in the table above. Plot order must not be inferred from alphabetic sorting.

## Repository layout

```text
.
├── README.md                  # Project overview and reproducibility information
├── .gitignore                 # Local R files, private data and generated output
├── scripts/
│   ├── INF1/                  # INF1 analysis scripts and run-order documentation
│   └── CPZ2/                  # CPZ2 analysis scripts and run-order documentation
├── R/                         # Shared helper functions, if needed
├── data/
│   ├── README.md              # Input inventory and data-access policy
│   ├── raw/                   # Source data (not committed by default)
│   └── processed/             # Derived analysis inputs (not committed by default)
├── outputs/                   # Tables and figures (not committed by default)
└── docs/
    ├── CODE_STYLE.md          # Script documentation conventions
    └── script-template.R      # Header template for new analysis scripts
```

Scripts will be grouped by **record and analytical step**, with a documented execution order. Shared functions should be introduced only where they genuinely reduce duplication or make methodological assumptions clearer.

## Input data and provenance

Depending on the analysis, the workflows use:

1. **Pollen counts:** sample-by-taxon count matrices, with sample identifiers consistent across tables.
2. **Chronological metadata:** sample depths, calibrated ages, age–depth model outputs and pollen-zone assignments.
3. **Taxonomic lookup tables:** harmonised pollen taxon names, ecological groupings and, when required, mappings to species or pollen productivity categories.
4. **Palaeoclimate reconstructions:** externally derived climate time series (including PaleoView outputs for relevant analyses).
5. **Landscape reconstruction inputs:** regional pollen records, pollen productivity estimates and LRA parameter tables.
6. **Ecological indicators:** species-level indicator information and the matching tables used to relate it to the pollen data, where relevant.

**Data availability is not identical to code availability.** Raw counts, third-party datasets, site metadata and derived datasets will be distributed only where permissions and source licenses allow it. See [`data/README.md`](data/README.md) for the input inventory to complete before release.

## Software and reproducibility

The workflow is developed in **R**, using packages such as `tidyverse`, `vegan`, `ggplot2`, `rioja`, `ecoTraj`, `RRatepol`, `mgcv`, `adespatial` and `rio`, as appropriate to the analysis. Additional packages or external tools used for LRA and taxonomic matching will be documented alongside the corresponding scripts.

The intended release procedure is:

1. Obtain the repository and the datasets permitted for distribution.
2. Install the documented R version and restore package versions from `renv.lock` **once that file is included**.
3. Place non-distributed input datasets in the locations specified in `data/README.md`.
4. Run the scripts in the order specified in `scripts/INF1/README.md` or `scripts/CPZ2/README.md`.
5. Compare generated tables and figures with the expected-output inventory supplied in each case-study directory.

Random seeds, taxonomic exclusions, temporal windows, transformations and model parameters should be written explicitly in the code. **An executable, end-to-end reproduction is not claimed until the scripts, input inventory and environment lockfile have been validated together.**

## Coding and documentation conventions

Every analytical script should state its **purpose, required inputs, generated outputs, package dependencies, main parameters and methodological assumptions** before performing the analysis. Methods-related comments should explain *why* a choice was made, rather than merely restating the R command.

The project uses readable, linear scripts for record-specific analyses and simple `ggplot2` code for figures. See [`docs/CODE_STYLE.md`](docs/CODE_STYLE.md) for the naming, commenting and validation conventions.

## Interpretation and limitations

The two records are separate palaeoecological case studies. Differences in sedimentary setting, sampling density, chronological uncertainty, taxonomic resolution and source-data coverage can influence downstream analyses. Outputs should be interpreted with these constraints in mind.

Species-level matches derived from pollen taxa are **taxonomic associations**, not direct observations of fossil species. Where a descriptive pollen–flora matching score is used, it must not be interpreted as a probability of past presence, ecological suitability or identification certainty.

## Citation and reuse

If you use the code or ideas implemented here, please cite the relevant peer-reviewed publications and, when available, the final archived release of this repository. A formal citation file and permanent identifier will be added when release metadata have been confirmed.

**License:** to be selected before public release. Code and third-party datasets may require different licenses; no reuse license should be inferred from their presence in the repository.

## Contact

For methodological questions, dataset provenance or reproducibility queries, please use the repository's GitHub Issues once the repository is published.
