# tcga_survivome

Interactive cohort building and survival analysis across the TCGA pan-cancer cohort
(33 cancer types, 11,160 patients), stratified by clinical, therapy and genomic-summary
features. It is built to run as a static site on GitHub Pages.

## Status

- [x] Data pipeline: one harmonized row per patient (`data/processed/`)
- [x] Cohort builder UI (DuckDB-WASM): filters on any variable, live counts
- [x] Kaplan–Meier, log-rank, numbers at risk, medians and landmark survival (JS)
- [x] Cox models with Efron ties and optional stratification by cancer type (JS)
- [ ] Proportional-hazards checks (possibly webR + `survival`)
- [ ] GitHub repository, Actions build + Pages deploy

## Running the site locally

```sh
python3 -m http.server -d site 8765    # then open http://localhost:8765/
```

The site is plain HTML, CSS and ES modules in `site/` with no build step. `site/data` is a
symlink to `data/processed/`. DuckDB-WASM loads from the jsDelivr CDN, so the first load
needs a network connection.

| File | Role |
|---|---|
| `site/js/stats.js` | Kaplan–Meier, log-rank, Cox; no dependencies |
| `site/js/kmplot.js` | SVG survival plot with at-risk table and hover readout |
| `site/js/db.js` | DuckDB-WASM setup and queries |
| `site/js/app.js` | Filters, grouping, tables, Cox panel, link and download actions |

## Tests

The statistics are checked against R's `survival` package (`survfit`, `survdiff`, `coxph`);
how each one is computed is described in [docs/methods.md](docs/methods.md).

```sh
Rscript tests/make_fixtures.R      # only when the data or the reference cases change
node --test tests/stats.test.mjs
```

## Rebuilding the data

Requires R ≥ 4.3 with `data.table`, `readxl`, `arrow`, `jsonlite`, `httr2`, `TCGAbiolinks`.

```sh
Rscript pipeline/01_download.R          # GDC PanCanAtlas files, Xena, cBioPortal (~40 MB)
Rscript pipeline/02_download_therapy.R  # GDC BCR Biotab drug/radiation tables (~15 min)
Rscript pipeline/03_build.R             # -> data/processed/
```

Outputs:

| File | Contents |
|---|---|
| `survivome.parquet` | Analysis table, 1 row per patient (loaded by the site) |
| `survivome.csv.gz` | Same table as CSV |
| `data_dictionary.json` | Label, group, type, unit, source, description, coverage per column |
| `qc_report.md` | Events, follow-up and feature coverage by cancer type |

## What's in the table

| Group | Features |
|---|---|
| Survival | OS, DSS, DFI, PFI (event + time in days) from TCGA-CDR |
| Clinical | cancer type, age, sex, race, ethnicity, genetic ancestry, stage (I–IV + detail), histology, grade, dx year |
| Therapy | chemo / hormone / targeted / immunotherapy / radiation (yes vs. "no record"), neoadjuvant history |
| Mutation | TMB (+ class), MSI (MANTIS, MSIsensor, status) |
| DNA repair | HRD score and components (LOH, TAI, LST), HRD class |
| Chromosomal instability | aneuploidy score, 39 arm-level calls, fraction genome altered, segment counts, tumor break load, ploidy, whole-genome doubling, subclonal fraction |
| Tumor | purity, hypoxia scores (Buffa, Ragnum, Winter) |
| Subtypes | PanCanAtlas selected / mRNA / methylation / CNA / integrative subtypes |
| Immune | immune subtype C1–C6, leukocyte fraction, TCR/BCR diversity |
| Viral | HPV, HBV, HCV, EBV, MCPyV read counts |
| QC | AWG do-not-use and pathology exclusion flags, redaction |

## Known limitations

- **Therapy data is sparse.** A "no record" value means nothing was recorded, not that
  the patient was confirmed untreated. Drug records exist for ~4,300 patients; none for LAML.
- **Survival endpoints vary in quality by cancer type.** TCGA-CDR recommends PFI over
  OS given short follow-up. Some types have very few OS events (PRAD 10, TGCT 4, PCPG 6,
  THCA 16, THYM 9).
- **Stage** is unavailable for GBM, LGG, LAML, PCPG, PRAD, SARC and THYM (no AJCC staging).
- **MSI** is missing for metastatic samples, so most SKCM cases lack it.
- **Immune subtypes** cover primary tumors only; DLBC and THYM are not assigned.
- **TMB** uses cBioPortal's MC3-derived values, about 1.3× the PanCanAtlas mutation-load values
  (both included). The TMB-H ≥ 10 cutoff depends on which you use.
- One tumor sample per patient: primary (-01) > primary blood (-03) > metastatic (-06) > recurrent.
- MRD isn't available: TCGA has a single pre-treatment sample per patient.

## Data sources and citations

All data are open access. Please cite the original publications:

- Liu J et al. An integrated TCGA pan-cancer clinical data resource (TCGA-CDR). *Cell* 2018.
- Ellrott K et al. Scalable open science approach for mutation calling (MC3). *Cell Syst* 2018.
- Thorsson V et al. The immune landscape of cancer. *Immunity* 2018.
- Knijnenburg TA et al. Genomic and molecular landscape of DNA damage repair deficiency. *Cell Rep* 2018.
- Taylor AM et al. Genomic and functional approaches to understanding cancer aneuploidy. *Cancer Cell* 2018.
- Bonneville R et al. Landscape of microsatellite instability across 39 cancer types. *JCO PO* 2017.
- Hoadley KA et al. Cell-of-origin patterns dominate the molecular classification of 10,000 tumors. *Cell* 2018.
- Cerami E et al. / Gao J et al. cBioPortal. *Cancer Discov* 2012 / *Sci Signal* 2013.
- Colaprico A et al. TCGAbiolinks. *Nucleic Acids Res* 2016.

Raw file URLs are listed in `pipeline/sources.csv`.
