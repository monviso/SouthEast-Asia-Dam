# Southeast Asian dam ECOSTRESS WST modelling workflow

This repository contains the scripts and input data used to analyse ECOSTRESS-derived water surface temperature (WST) observations for Southeast Asian reservoirs and downstream outlet domains.

The workflow supports the manuscript:

**Predicting Water Temperature in Southeast Asian Reservoirs and its Impact on Downstream River Thermal Regimes**

The analysis uses sparse, high-quality ECOSTRESS WST observations, dam metadata, and hourly meteorological predictors to test whether reservoir and downstream outlet temperatures show different lagged responses to meteorological forcing.

## Repository structure

```text
.
├── dat/
│   ├── sites_complete_HQ_obs_ECOSTRESS.csv
│   ├── site_full.shp
│   ├── site_full.shx
│   ├── site_full.dbf
│   ├── site_full.prj
│   
│
├── script/
│   ├── Step2_get_climate_predictors_v2_fix_missing_coords.R
│   ├── Step3A_main_GAM_lag_LODO_noHour.R
│   ├── Step3B_sampling_coverage_hour_sensitivity.R
│   ├── Step3C_LODO_error_diagnostics.R
│   ├── Step3D_plot_significant_GAM_effects_ggplot_v5.R
│   └── plot_ecostress_dam_figures_from_table.R
|   └── Plots_paper_rev.R
│
└── README.md
```

## Input data

### `dat/sites_complete_HQ_obs_ECOSTRESS.csv`

Main ECOSTRESS high-quality observation table. The table contains retained water-surface-temperature observations for reservoir and outlet domains after ECOSTRESS quality filtering and water masking.

Important fields include:

| Field | Description |
|---|---|
| `dam` / `dam_key` | Dam/site identifier. |
| `date.time` / `datetime_utc` | ECOSTRESS acquisition time. |
| `obs_domain` | Observation domain: paired, lake-only, or outlet-only. |
| `Tm_l` | Reservoir/lake WST in degrees Celsius. |
| `Tm_o` | Downstream outlet WST in degrees Celsius. |
| `T_diff` | Outlet minus reservoir WST, calculated as `Tm_o - Tm_l`. |
| `Latitude_l`, `Longitude_l` | Lake-domain coordinates. |
| `Latitude_o`, `Longitude_o` | Outlet-domain coordinates. |
| `Latitude_ref`, `Longitude_ref` | Reference dam/site coordinates, where available. |
| `year`, `month`, `doy` | Derived acquisition-date variables. |
| `H_UTC`, `H_LOC` | UTC and local acquisition hour. |
| `LSTerr_l`, `LSTerr_o` | ECOSTRESS uncertainty/error fields for lake and outlet domains, where available. |
| `Dam_hgt` | Dam height, in metres. |
| `Res_area_kmsq` | Reservoir area, in square kilometres. |
| `Res_capaci_kmcube` | Reservoir capacity, in cubic kilometres. |
| `Country` | Country of the dam. |

The exact column set may vary slightly depending on the archived table version. If file names differ, update the input paths at the top of each script.

### `dat/site_full.shp`

Dam-site spatial layer used for mapping and site metadata. All shapefile sidecar files must remain in the same directory, including `.shx`, `.dbf`, `.prj`, and any additional files produced with the shapefile.

## Software requirements

The workflow was developed in R. The main R packages used are:

```r
mgcv
dplyr
tidyr
readr
purrr
lubridate
stringr
sf
terra
ggplot2
patchwork
viridis
scales
jsonlite
```

Install missing packages with:

```r
install.packages(c(
  "mgcv", "dplyr", "tidyr", "readr", "purrr", "lubridate",
  "stringr", "sf", "terra", "ggplot2", "patchwork", "viridis",
  "scales", "jsonlite"
))
```

## Running the workflow

Run the scripts from the repository root, or edit the path settings at the top of each script so that they point to the local `dat/`, `script/`, and output directories.

The recommended run order is:

```text
1. Step2_get_climate_predictors_v2_fix_missing_coords.R
2. Step3A_main_GAM_lag_LODO_noHour.R
3. Step3B_sampling_coverage_hour_sensitivity.R
4. Step3C_LODO_error_diagnostics.R
5. Step3D_plot_significant_GAM_effects_ggplot_v5.R
6. plot_ecostress_dam_figures_from_table.R
7. Plots_paper_rev
```

### Step 2: meteorological predictors

```text
script/Step2_get_climate_predictors_v2_fix_missing_coords.R
```

This script reads the ECOSTRESS HQ observation table and retrieves or prepares hourly meteorological predictors for each observation using the NASA POWER hourly point API.

Main outputs include:

| Output | Description |
|---|---|
| `dt0_nasa_power_zerolag.csv` | Observation table with zero-lag meteorological predictors. |
| `lag_long_nasa_power_0_48h.csv` | Long-format hourly meteorological lag table from 0 to 48 h. |
| `gam_lag_data_nasa_power.rds` | Model-ready lagged predictor object used by Step 3A. |
| `climate_hourly_nasa_power.rds` | Cached hourly meteorological data. |

The script also checks and reconciles missing coordinates where possible. Observations without sufficient coordinate information are reported and excluded from meteorological extraction.

### Step 3A: main GAM lag models and leave-one-dam-out validation

```text
script/Step3A_main_GAM_lag_LODO_noHour.R
```

This is the main modelling script. It fits generalized additive models for:

1. Reservoir/lake WST;
2. Downstream outlet WST;
3. Outlet–reservoir WST difference, defined as `Tm_o - Tm_l`.

The script compares models with no meteorological lag and with 4 h, 12 h, 24 h, and 48 h lag windows. The primary models include meteorological forcing, dam characteristics, and seasonal timing, but exclude local solar hour from the main structure.

The script also performs leave-one-dam-out cross-validation (LODO-CV), where each dam is held out in turn and predicted using a model fitted to the remaining dams.

Main outputs include:

| Output | Description |
|---|---|
| `model_comparison_zerolag.csv` | Comparison of zero-lag models. |
| `model_comparison_lag.csv` | Comparison of lagged models. |
| `model_comparison_all.csv` | Combined model-comparison table. |
| `best_model_by_response.csv` | Selected lag/model for each response. |
| `model_sample_sizes.csv` | Sample sizes used by each model. |
| `fits/` | Saved fitted GAM objects. |
| `lodo_predictions/` | Held-out dam predictions. |
| `lodo_by_dam/` | LODO-CV performance summarised by dam. |
| `model_summaries/` | GAM summaries and smooth statistics. |

### Step 3B: sampling coverage and sensitivity analyses

```text
script/Step3B_sampling_coverage_hour_sensitivity.R
```

This script evaluates whether model results are sensitive to uneven ECOSTRESS sampling across dams, months, and local acquisition hours.

It produces summaries of observation coverage and compares model structures with and without local solar hour. It also evaluates sampling-balance sensitivity using dam-month-hour groupings.

Main outputs include:

| Output | Description |
|---|---|
| Coverage tables | Observation counts by dam, month, day of year, local hour, and dam-month-hour combinations. |
| Evenness metrics | Metrics describing sampling imbalance. |
| Hour-sensitivity models | Comparison of no-hour, with-hour, no-cyclic, and cyclic-only model structures. |
| Sampling-balance summaries | Model performance under full and balanced sampling structures. |

### Step 3C: leave-one-dam-out error diagnostics

```text
script/Step3C_LODO_error_diagnostics.R
```

This optional diagnostic script uses the LODO-CV predictions from Step 3A to explore where prediction errors are largest.

It summarises errors by dam, month, local hour, observation source, and other available metadata. It is intended for diagnostic interpretation and supporting material, rather than primary model selection.

### Step 3D: significant GAM effect plots

```text
script/Step3D_plot_significant_GAM_effects_ggplot_v5.R
```

This script produces publication-style plots of significant smooth effects from the selected GAMs.

It uses the fitted models from Step 3A and generates figures for significant one-dimensional and two-dimensional smooth terms, including meteorological-lag response surfaces.

Main outputs include:

| Output | Description |
|---|---|
| Significant-effect figures | PNG/PDF figures of significant GAM smooths. |
| `significant_terms_selected_models_ggplot_v5.csv` | Table of selected significant terms. |
| `panel_status_ggplot_v5.csv` | Plotting status for each smooth term. |
| `shared_2d_colour_limits_ggplot_v5.csv` | Shared colour-scale diagnostics for 2D lag surfaces. |

### Figure script: site maps and observation-summary plots

```text
script/plot_ecostress_dam_figures_from_table.R
```

This script produces figures summarising the ECOSTRESS observation dataset, including dam-location maps and WST distribution plots.

It uses the main ECOSTRESS HQ table and the `site_full.shp` spatial layer.

## Main model responses

The main response variables are:

| Response | Definition |
|---|---|
| Reservoir WST | `Tm_l`, mean lake/reservoir-domain water surface temperature. |
| Outlet WST | `Tm_o`, mean downstream outlet-domain water surface temperature. |
| WST difference | `T_diff = Tm_o - Tm_l`. Negative values indicate cooler outlet WST than reservoir WST. |

## Notes on interpretation

ECOSTRESS-derived values are interpreted as **water surface temperature** after water masking. They should not be interpreted as direct measurements of bulk water temperature or release temperature from dam intakes.

The models are intended primarily for inference and feasibility testing. They assess whether sparse ECOSTRESS observations contain information on lagged reservoir and outlet thermal dynamics. They should not be treated as a fully operational forecasting system without additional local calibration and validation.

## Reproducibility notes

- Run scripts in the order listed above.
- Keep all shapefile sidecar files together in the `dat/` folder.
- Check path variables at the top of each script before running.
- NASA POWER retrieval requires internet access unless cached meteorological files are already available.
- Output directories are created automatically by the scripts where possible.
- If using a different filename for the ECOSTRESS HQ observation table, update the corresponding input path in each script.

## Citation

If using this code or dataset, please cite the associated manuscript and dataset record:

Redana, M. et al. *Predicting Water Temperature in Southeast Asian Reservoirs and its Impact on Downstream River Thermal Regimes.*
