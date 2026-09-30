# R Source Code

This folder contains cleaned, reusable pieces migrated from the original thesis workspace.

## Files

- `gw_qxgb.R` — core geographically weighted quantile XGBoost estimator and spatial kernel helpers.
- `gw_qxgb_cv.R` — leave-one-location-out bandwidth selection with golden-section search.
- `metrics.R` — Pinball Loss and quantile Pseudo R² helpers.

## About GW-QRF

The original thesis workspace references a separate local file named `gqrf.R`, but that source file is not currently present in the connected thesis repository.

For that reason, a GW-QRF implementation has **not** been recreated from memory. The original file should be added when available.
