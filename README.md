# Future Droughts R Scripts

This repository contains R scripts to compute evapotranspiration and SPEI from gridded climate data.

## Scripts

### functions_future_droughts.r
Main workflow script for raster-based drought analysis.

What it covers:
- Evapotranspiration functions
- `eto_pixel`: FAO-56 Penman-Monteith reference ET (ETo) per pixel
- `pet_pixel`: Priestley-Taylor-like PET variant per pixel
- `compute_et`: parallel wrapper for SpatRaster ET/PET computation

- Reference-period parameter estimation for SPEI
- `pwm_to_abc`: log-logistic PWM parameter estimation (a, b, c)
- `get_all_abc`: computes monthly parameter sets (12 months x 3 params)
- `compute_params`: parallel wrapper returning 36 parameter layers (`m01_a` ... `m12_c`)

- SPEI computation
- `compute_spei`: computes monthly SPEI from `tp`, `pev`, and 36-layer parameter raster
- Returns a SpatRaster aligned with input precipitation/evapotranspiration rasters

Notes:
- Inputs are expected to have matching geometry (extent, resolution, CRS).
- SPEI inputs should be monthly and have layer count divisible by 12.
- Parameter raster is expected to have 36 layers ordered by month and parameter.

### bkp.r
Legacy/backup script with early versions of ET and SPEI utilities.

What it includes:
- Scalar ET helper functions (`atm_press`, `psychometric_constant`, `sat_vapour_pressure`, `slope_vapour_pressure_curve`, `get_eto`)
- Early SPEI parameter and computation helpers (`param_loglogist`, `get_spei`)
- Exploratory raster/data loading snippets

Use this mainly as historical reference; prefer `functions_future_droughts.r` for the current workflow.

## Repository Tracking Rules

The `.gitignore` is configured to:
- Track all `.r` files
- Track everything under `files/`
- Exclude everything under `data/`
- Ignore common system/editor artifacts
