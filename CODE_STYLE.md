# Code style and documentation conventions

This project prioritises **traceability of scientific decisions** and readability over elaborate abstractions.

## One script, one analytical goal

Group scripts by record (`scripts/INF1/` or `scripts/CPZ2/`) and by step. Use a numeric prefix to document execution order once the actual workflow is confirmed (for example, `03_rate_of_change.R`). Avoid introducing general-purpose wrappers for code that is clearer as a short linear sequence.

## Script headers

Every script begins with a concise header specifying:

- **Purpose** — the question or output addressed.
- **Inputs** — exact relative paths, table names and essential columns.
- **Outputs** — exact relative paths and what they contain.
- **Dependencies** — required packages, external tools and upstream scripts.
- **Parameters** — non-default settings, random seeds and rationale.
- **Caveats** — exclusions, data availability and known limitations.

Use [`script-template.R`](script-template.R) as a starting point.

## Naming and directories

- Use consistent `snake_case` for new objects and script names. Preserve existing scientifically meaningful identifiers (e.g. `INF1`, `CPZ2`, `MOD1a`).
- Keep raw files in `data/raw/`, transformed analysis inputs in `data/processed/` and generated figures/tables in `outputs/`.
- Use **relative** file paths. Do not hardcode machine-specific paths or secrets.
- Treat sample order, age order and factor-level order explicitly. For INF1, the correct order is `MOD1a`, `MOD1b`, `MOD2`, `MOD3`, `MOD4a`, `MOD4b`; for CPZ2, `MAL1a`, `MAL1b`, `MAL2`, `MAL3`, `MAL4`.

## Data checks

Before model fitting, check that sample IDs align across data frames, expected columns are present, counts are non-negative and age/depth variables use stated units. Where appropriate, document missing values, zero-valued predictors, outliers and filtering criteria. Do not silently drop samples or taxa.

## Methodological transparency

- Explain why a transformation or distance measure was chosen.
- Keep all run parameters visible near the code that uses them.
- Call `set.seed()` when stochastic operations affect results and record the value.
- Distinguish diagnostics, sensitivity analyses and final results.
- Never label a taxonomic matching score as identification certainty or past-presence probability.
- Do not assume that parameters are transferable between INF1 and CPZ2.

## Figure scripts

Use direct `ggplot2` code where possible. Keep zone palettes, axis direction, factor order, plot dimensions and file formats explicit. Export publication-ready figures with reproducible settings.

## Reproducibility checklist

Before publishing an analysis script, verify that it (1) runs from the repository root, (2) has documented inputs, (3) produces the expected output, (4) works with the declared package versions, (5) includes complete parameter settings, and (6) does not commit private or restricted data.
