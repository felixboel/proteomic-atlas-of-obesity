# Proteomic Atlas of Obesity (Shiny webtool)

Interactive Shiny web application accompanying the study **“Multi-organ proteomic atlas of obesity regression in male mice”** (Boel et al.).

The tool enables exploration of tissue- and time-specific proteomic shifts across mouse organs during obesity and its regression.

---

## Hosted application

- **Public instance:** *TBD*

---

## Repository contents

- `app/` — Shiny application  
  - `app/app.R` — app entry point  
  - `app/R/` — helper scripts (`config.R`, `data_load.R`, `helpers.R`)  
  - `app/www/` — static assets (CSS, fonts)
- `Dockerfile` — container build for deployment
- `renv.lock`, `renv/` — pinned R package versions for reproducibility
- `data/` — placeholder for processed data (not included during embargo)

---

## Webtool features

### Pathway viewer
- Navigate Reactome pathway hierarchy using an interactive pathway tree.
- Bubble plot summarizing regulated proteins across organs and timepoints.
- Adjustable filtering by BH-adjusted p-value cutoff and log fold-change cutoff.

### Gene viewer
- Query one or more gene symbols (one per line).
- Organ selection interface with color-coded organ toggles.
- Multi-panel interactive plots showing gene-level changes across selected organs.

### General
- Interactive Plotly visualizations (zoom, hover tooltips, export).
- Consistent organ color palette across all views.

---

## Data availability

Processed data files are **not included in this repository during peer review**.

Upon publication, processed data will be made publicly available via PRIDE (accession: PXD066875), and this repository will be updated accordingly.

---

## Running locally

### 1) Restore packages using renv
From R:

```r
renv::restore()
```

### 2) Run the Shiny app
```r
shiny::runApp("app")
```
Requires the processed data files in data/processed/ (see Data availability above).

---

## Running locally
Build image
```bash
docker build -t proteomic-atlas-of-obesity .
```
Run container:
```bash
docker run --rm -p 3838:3838 proteomic-atlas-of-obesity
```
Requires the processed data files in data/processed/ (see Data availability above).

---

## Citation
If you use this webtool or the underlying dataset, please cite:

Boel F., et al. Multi-organ proteomic atlas of obesity regression in male mice. (Journal, Year) — DOI: TBD






