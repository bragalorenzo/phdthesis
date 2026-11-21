# Numerical ecology of Holocene pollen records

**Reproducible R workflows for my PhD Thesis**

 *Furthering the relationships between palaeoecology and biodiversity: numerical ecology methods and development of informatics tools in palaeopalynological case studies* (University of Modena and Reggio Emilia, Italy).

The analyses combine pollen-assemblage data, age–depth information, environmental reconstructions and taxonomic reference tables to investigate compositional change, temporal turnover, local vegetation dynamics and the relationship between palaeoecological records and contemporary biodiversity information.

## Case studies

| Record | Study site | Approximate temporal coverage | Samples used in the numerical analyses | Pollen zones (oldest → youngest) |
| --- | --- | --- | ---: | --- |
| **INF1** | Teor–Modeano palaeochannel, Friuli Venezia Giulia, northeastern Italy | ca. 4,200–400 cal yr BP | 73 | MOD1a → MOD1b → MOD2 → MOD3 → MOD4a → MOD4b |
| **CPZ2** | Malafede Pool, Castelporziano Presidential Estate, central Italy | ca. 7,000–400 cal yr BP | 29 | MAL1a → MAL1b → MAL2 → MAL3 → MAL4 |

The records differ in depositional setting, temporal resolution and ecological context. The analytical approaches thus share methodology only.

## Analytical scope

The repository is organised around the following components. Not every component is necessarily applied in the same form to both records; record-specific implementations and settings belong with their respective scripts.

| Component | Purpose |
| --- | --- |
| Palynological richness (Rarefaction analysis) | Compare sample-level richness standardising pollen count |
| Rate of Change (RoC) analysis| Quantify compositional turnover in relation to time|
| Community Trajectory Analysis (CTA) | Characterise assemblage trajectories and compositional heterogeneity with pollen zone-wise internal variance, sample-wise and taxon-wise contributions to record- and pollen zone-wise variance |
| Redundancy Analysis (RDA)| Mean annual temperature and summer precipitation as drivers of change of INF1 and CPZ2 communities |
| Generalised Additive Models (GAMs) | Assess staitistically significant temporal patterns in individual taxa |
| Ecological indicator analyses with EIVE| Explore ecological information associated with floristically matched pollen taxa |
| Landscape Reconstruction Algorithm (LRA) | Estimate local vegetation/land-cover composition |

Analyses share a common starting record-specific .Rds file, containing data frames with data used for each one.

## Citation and reuse

This repository is explicitly meant for consultation only by the thesis's referees and evaluating committee.
