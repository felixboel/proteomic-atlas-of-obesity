# Proteomic Atlas of Obesity

Interactive R Shiny web application for exploring tissue- and time-specific proteomic changes during obesity and weight-loss regression in male mice.

The application provides interactive exploration of the multi-organ proteomic atlas described in:

**Boel, F. et al.** *Multi-organ proteomic atlas of obesity regression in male mice.* **Nature Metabolism** (2026).  
https://doi.org/10.1038/s42255-026-01599-5

---

## Hosted application

https://computproteomics.bmb.sdu.dk/app_direct/obesity-atlas/

---

## Data availability

Raw mass spectrometry proteomics data are available through the PRIDE repository under accession:

**PXD066875**

https://www.ebi.ac.uk/pride/archive/projects/PXD066875

Processed data used by the application are included in this repository under:

```text
data/processed/
```

---

## Webtool features

### Pathway Viewer

- Browse the Reactome pathway hierarchy using an interactive pathway tree.
- Visualize significantly regulated proteins across tissues and timepoints using interactive bubble plots.
- Adjust BH-adjusted P-value and absolute log2 fold-change thresholds.

### Gene Viewer

- Query one or more gene symbols (one per line).
- Select tissues interactively.
- Visualize estimated log2 fold changes across obesity and regression timepoints.
- Profiles are shown relative to age-matched lean controls with ±1 SE.

---

## Study design

Male mice were fed a high-fat diet (HFD) for 18 weeks to induce obesity and were subsequently switched to a low-fat diet (LFD) to induce weight loss.

- **OBE** — after 18 weeks HFD, immediately before switching to LFD
- **STR** — 2 weeks after switching to LFD
- **MTR** — 6 weeks after switching to LFD
- **LTR** — 12 weeks after switching to LFD

Differential protein expression was estimated using robust ridge regression (MSqRob2), where each timepoint was compared against age-matched lean control mice maintained on LFD throughout.

---

## Run locally with R

**Requirements**

- R 4.5.2
- renv

From the project root:

```bash
R -e "install.packages('renv', repos='https://cloud.r-project.org')"
R -e "renv::restore(prompt = FALSE)"
R -e "shiny::runApp('app', host = '0.0.0.0', port = 3838)"
```

Then open:

http://localhost:3838

---

## Run with Docker

Pull the published Docker image:

```bash
docker pull felixboel/proteomic-atlas-of-obesity:latest
```

Run the container:

```bash
docker run --rm -p 3838:3838 felixboel/proteomic-atlas-of-obesity:latest
```

Then open:

http://localhost:3838

Versioned Docker images are available.

---

## Citation

If you use this resource or the underlying dataset in academic work, please cite:

**Boel, F. et al.** *Multi-organ proteomic atlas of obesity regression in male mice.* **Nature Metabolism** (2026).  
https://doi.org/10.1038/s42255-026-01599-5

---

## License

Released under the MIT License. See [LICENSE](LICENSE) for details.
