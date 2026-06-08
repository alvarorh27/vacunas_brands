# Vaccine Effectiveness by Brand and Vaccination Regimen (COVID-19)

This repository contains the analysis pipeline for a study of COVID-19 vaccine
effectiveness by brand, vaccination regimen (one, two, or three doses), and
time since last immunization, using population-level health records from
Andalusia, Spain. The pipeline builds an analytic dataset from raw
administrative and clinical extracts, characterizes the study population, and
fits regression and generalized additive models (GAMs) to estimate the
association between vaccination status and a composite clinical outcome
(hospitalization, ICU admission, or 30-day mortality following infection).

## Repository structure

```
analysis/   R Markdown notebooks implementing the analysis pipeline (see below)
Rfiles/     Intermediate and final analytic datasets (.rds; not version-controlled)
data/       Exported odds-ratio (OR) tables in CSV format
figures/    Generated plots and exported tables (PNG/JPG), organized by type:
              Forest_plot/  forest plots of GAM-based odds ratios
              GAM_plot/     GAM smooth-term plots (probability vs. age / days
                            since last immunization), including breakpoint analyses
              OR/           odds-ratio summary tables exported as images
              Table1/       descriptive "Table 1" exported as an image
tables/     Descriptive tables exported as Word documents (.docx)
```

Raw data extracts are read from a shared, access-restricted location
(`/opt/datos_compartidos/vacunas_data/`) and are not distributed with this
repository. The `Rfiles/`, `data/`, and rendered HTML outputs are excluded from
version control via `.gitignore`.

## Analysis workflow

The notebooks in `analysis/` are designed to be run in sequence, each
consuming the dataset produced by the previous stage:

```
01_dataset.Rmd  →  04_eda.Rmd  →  02_tables.Rmd
                                  03_modeling.Rmd
```

1. **`01_dataset.Rmd` — Cohort construction and dataset assembly**
   Loads and harmonizes the raw extracts (PCR test results, target population
   registry, comorbidity records [PUPA], hospital discharge records [CMBD],
   and vaccination records), then:
   - Identifies COVID-19 infection episodes and reinfections (using a 90-day
     window to distinguish reinfections from persistent positivity).
   - Derives the study population, applying five sequential exclusion
     criteria (out-of-period PCR dates; invalid sex; implausible age;
     more than three vaccine doses; children under 18 with more than two doses).
   - Classifies vaccination history by brand, by EU/non-EU origin, by mRNA vs.
     non-mRNA platform, and by 1-/2-/3-dose regimen and brand combination.
   - Derives outcome variables: hospitalization, ICU admission, 30-day
     mortality, and a composite `clinical_event` outcome.
   - Assigns a predominant SARS-CoV-2 variant period to each infection.
   - **Output:** `Rfiles/final_dataset_brands.rds` (one row per infection
     episode, with demographic, clinical, comorbidity, and vaccination
     variables).

2. **`04_eda.Rmd` — Exploratory data analysis and additional exclusions**
   Examines the distribution of vaccine-shot combinations, age, and time since
   last immunization/infection in the assembled dataset. Defines grouped bands
   for "days since last immunization" (used as a stratification variable in
   later models) based on the empirical 95th percentile.
   - **Input:** `Rfiles/final_dataset_brands.rds`
   - **Output:** `Rfiles/final_dataset_brands_excluded.rds` (the analytic
     dataset used by all downstream descriptive and modeling notebooks).

3. **`02_tables.Rmd` — Descriptive statistics ("Table 1")**
   Produces the descriptive characteristics table stratified by the composite
   clinical outcome, including:
   - Normality testing (Kolmogorov–Smirnov) for continuous variables.
   - Construction of a labeled "Table 1" via `gtsummary`/`tableone`, with
     translated factor levels and cleaned variable labels for publication.
   - Export of the table as a formatted Word document (`tables/*.docx`) and as
     an image (`figures/Table1/`).
   - **Input/Output:** reads `final_dataset_brands_excluded.rds`; writes to
     `tables/` and `figures/Table1/`.

4. **`03_modeling.Rmd` — Regression and GAM-based effectiveness analyses**
   Implements the core modeling workflows used to estimate vaccine
   effectiveness against the composite clinical outcome, stratified by
   vaccination regimen (`vac_1shot`, `vac_2shot`, `vac_3shot`):
   - **Logistic regression with odds-ratio extraction**: helper functions fit
     logistic models and format coefficients as odds ratios with confidence
     intervals, exported to `data/*.csv`.
   - **GAM models with one spline (age)**: smooth-term models adjusted for
     age, sex, comorbidities, and SARS-CoV-2 variant period; results are
     summarized as forest plots (`figures/Forest_plot/fp_gam_1spline_*`) and
     odds-ratio tables (`figures/OR/`).
   - **GAM models with two splines (age and days since last
     immunization/infection)**: probability curves plotted against each
     smooth term (`figures/GAM_plot/GAM_*_2spline_*`).
   - **Piecewise/breakpoint segmentation**: identifies change points in the
     relationship between time since last immunization and outcome
     probability (`figures/GAM_plot/*_breakpoints.png`).
   - **Forest plots**: both "reduced" (key covariates only) and "all
     variables" model specifications are produced for each regimen.
   - **Interactive GAM visualizations**: `plotly`-based HTML widgets with
     tooltips highlighting minimum/maximum predicted probabilities, saved as
     self-contained HTML files in `figures/GAM_plot/`.
   - **Input:** `final_dataset_brands_excluded.rds`
   - **Output:** forest plots, GAM plots, OR tables, and interactive HTML
     widgets under `figures/`, and OR tables in `data/`.

## Key derived variables

- `vac_1shot`, `vac_2shot`, `vac_3shot` — vaccination regimen by number of
  doses and brand/category combination.
- `group_vaccines`, `european_vaccines` — brand groupings (single-brand,
  mixed mRNA, mixed non-mRNA, EU vs. non-EU origin).
- `ultima_inmun_infecc_dias` / `ultima_inmun_infecc_grouped` — days (and
  grouped bands) between the most recent immunizing event (vaccination or
  prior infection) and the index infection.
- `clinical_event` — composite outcome combining hospitalization, ICU
  admission, and 30-day mortality.
- `mortality30d` — death within 30 days of the index infection, used as the
  basis for the composite outcome.
- `variant_period` — predominant circulating SARS-CoV-2 variant at the time
  of infection.

## Modeling methods

- **Logistic regression** for odds-ratio estimation of vaccination effects on
  the composite clinical outcome.
- **Generalized additive models (GAMs)**, fit with `mgcv::gam()` using
  penalized regression splines (`s()`), in two specifications:
  - *one spline* on age, and
  - *two splines* on age and on days since last immunization/infection,
  both adjusted for sex, comorbidities, and circulating variant period.
- **Piecewise linear segmentation / breakpoint analysis** of GAM smooth terms
  to characterize how protection changes over time since last immunization.
- **Forest plots** summarizing odds ratios (and 95% confidence intervals)
  across covariates, for both reduced and fully adjusted models, by
  vaccination regimen.

## Outputs

| Location | Content |
|---|---|
| `Rfiles/*.rds` | Intermediate and final analytic datasets |
| `tables/*.docx` | Descriptive "Table 1", exported as Word documents |
| `data/*.csv` | Odds-ratio tables from logistic and GAM models |
| `figures/Table1/` | Descriptive table exported as an image |
| `figures/Forest_plot/` | Forest plots of adjusted odds ratios by regimen |
| `figures/GAM_plot/` | GAM smooth-term plots, breakpoint plots, and interactive HTML widgets |
| `figures/OR/` | Odds-ratio summary tables exported as images |

## Dependencies

The analysis is implemented in R (R Markdown notebooks) and relies primarily
on:

- **Data wrangling:** `data.table`, `dtplyr`, `dplyr`, `tidyr`, `purrr`,
  `magrittr`, `lubridate`, `stringr`, `stringi`, `janitor`, `zoo`, `vroom`,
  `rio`, `feather`
- **Modeling:** `mgcv` (GAMs), `broom`, `caret`, `parallel`
- **Descriptive statistics and tables:** `tableone`, `gtsummary`,
  `summarytools`, `skimr`, `flextable`, `officer`, `kableExtra`, `htmlTable`,
  `knitr`
- **Visualization:** `ggplot2`, `htmlwidgets`, `webshot`

## Notes on data access

Raw data extracts (PCR results, vaccination records, hospital discharge
records, population registries) originate from restricted-access health
administration sources and are not included in this repository. Paths to
these sources are referenced directly in `01_dataset.Rmd` for provenance and
reproducibility within the originating institution; they cannot be resolved
outside that environment.
