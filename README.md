# CAW4AMR — Causal Analytic Workflow for AMR Surveillance

Public demonstration repository for the manuscript:
*The Causal Analytic Workflow: A Practical Guide to Structured Data
Processing for Antimicrobial Resistance Surveillance* (JAC-2026-0340).

## What this does

Processes raw neonatal infection surveillance data through an Extract-Transform-Load
pipeline to classify organisms as likely pathogens or contaminants, then aggregates
results across four levels: organism → specimen → episode → patient.

## Organism classification: two layers

1. **Taxonomy layer (AMR package)** — the R [`AMR`](https://amr-for-r.org) package resolves
   free-text organism names (including typos, synonyms, and reclassified names) into
   stable `mo` codes backed by LPSN/GBIF taxonomy. This runs at pipeline startup
   via `as.mo()`.

2. **Clinical curation layer (CSV dictionary)** — a version-controlled reference table
   (`algorithms/patho_grouping_v1.0_*.csv`) classifies each organism's pathogenicity
   as high or low, informed by the ACORN2/ANISA frameworks and expert clinical opinion.
   The pipeline joins taxonomy and curation on the `mo` code — never on free-text strings.

## Quick start

```r
# Clone and run the demo (R ≥ 4.0 required)
source("scripts/run_acorndemo_pipeline.R")
```

Missing R packages (`dplyr`, `readxl`, `AMR`) are auto-installed.

## Repository structure

```
├── README.md
├── algorithms/                          ← pathogenicity dictionary
│   └── patho_grouping_v1.0_*.csv
├── data/raw/                            ← synthetic demo data
│   └── acorndemo_data.xlsx
├── data/output/                         ← generated outputs (git-ignored)
├── scripts/
│   ├── run_acorndemo_pipeline.R         ← master driver
│   ├── initial_acorndemo_data_extract.R
│   ├── initial_acorndemo_data_transform.R
│   ├── initial_acorndemo_data_load.R
│   └── fcns_patho_grouping_v1.R         ← pathogenicity classification functions
└── .gitignore
```

## Dependencies

`dplyr`, `readxl`, `AMR` (auto-installed)

## Research question

What are the predominant pathogens responsible for early- and late-onset
neonatal infections in a tertiary NICU setting, and which isolated
organisms are more likely to be contaminants?

## Output column definitions

### acorndemo_organisms.csv (one row per isolated organism)

| Column | Type | Description |
|--------|------|-------------|
| `patient_id` | chr | De-identified patient code |
| `episode_id` | chr | De-identified episode code |
| `orgname` | chr | Organism name as reported by the laboratory |
| `specimen_id` | chr | De-identified specimen code |
| `likely_cause` | chr | Classification: `Yes` (likely pathogen), `No` (likely contaminant), `No growth`, `No organism name / culture result`, `Ungrouped` |
| `likely_contaminant` | chr | Inverse of `likely_cause`: `Yes` (likely contaminant), `No` (likely pathogen), etc. |

### acorndemo_episodes.csv (one row per clinical episode)

| Column | Type | Description |
|--------|------|-------------|
| `patient_id` | chr | De-identified patient code |
| `episode_id` | chr | De-identified episode code |
| `med_device_vent` | chr | Ventilator / intubation status |
| `med_device_pcv` | chr | Peripheral IV catheter status |
| `med_device_cvc` | chr | Central IV catheter status |
| `likely_causes_e` | chr | Episode-level summary of likely pathogens |
| `likely_contaminants_e` | chr | Episode-level summary of likely contaminants |

### acorndemo_specimens.csv (one row per specimen)

| Column | Type | Description |
|--------|------|-------------|
| `specimen_id` | chr | De-identified specimen code |
| `likely_causes_s` | chr | Specimen-level summary of likely pathogens |
| `likely_contaminants_s` | chr | Specimen-level summary of likely contaminants |

### acorndemo_persons.csv (one row per patient)

| Column | Type | Description |
|--------|------|-------------|
| `patient_id` | chr | De-identified patient code |
| `likely_causes_p` | chr | Patient-level summary across all episodes |
| `likely_contaminants_p` | chr | Patient-level summary across all episodes |

## Technical overview

### How the classification rules work

1. **High pathogenicity** (e.g. *Acinetobacter baumannii*, *Klebsiella pneumoniae*,
   *Escherichia coli*) → always classified as **likely cause** (`Yes`).

2. **Low pathogenicity** (e.g. *Staphylococcus epidermidis*, *Staphylococcus sp*)
   → normally classified as **likely contaminant** (`No`), unless the patient
   had an indwelling medical device (ventilator, peripheral IV, or central IV
   catheter) at the time of the episode, in which case the organism is upgraded
   to **likely cause** (`Yes`). This reflects the device-associated infection
   pathway common in NICUs.

3. **"No growth"** and **"No organism name / culture result"** are passed through
   as special categories.

4. Organisms **not found** in the dictionary are returned as `Ungrouped`.

### Hierarchy

Classifications propagate through four levels:

```
organism → specimen → episode → patient
```

- At the organism level: `assign_likely_cause()` classifies each isolate.
- At the specimen level: `assign_likely_causes_s()` collects all causes.
- At the episode level: `assign_likely_causes_e()` aggregates across specimens, with
  a special rule that removes *Staphylococcus sp* from the cause list if multiple
  organisms are present (it is assumed to be a contaminant in that context).
- At the patient level: `assign_likely_causes_p()` aggregates across all episodes.

The full set of 16 functions is in `scripts/fcns_patho_grouping_v1.R`.

## License

MIT. Data files in this repository are synthetic/dictionary only — no patient data.