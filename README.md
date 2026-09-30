# Locally Weighted Tree-Based Quantile Regression for Spatial Heterogeneity

Portfolio repository for my master's thesis:

**Locally Weighted Tree-Based Quantile Regression for Spatial Heterogeneity Analysis**

This research develops geographically weighted tree-based quantile models to capture spatial heterogeneity, nonlinear relationships, and distributional differences across conditional quantiles.

## Proposed Framework

The thesis proposes two local tree-based quantile modeling approaches:

- **GW-QRF** — Geographically Weighted Quantile Regression Forest
- **GW-QXGB** — Geographically Weighted Quantile XGBoost

The core idea is to assign spatial kernel weights around each target location and fit a local tree-based quantile model.

## Research Motivation

Geographically Weighted Quantile Regression (GWQR) can model both spatial and quantile heterogeneity, but its linear structure can be restrictive when relationships are nonlinear or involve complex interactions.

Global tree-based quantile models can capture nonlinear effects, but they do not directly estimate a location-specific model for every observation.

This work combines local spatial weighting with tree-based quantile learning.

## Evaluation Design

The thesis evaluates the proposed framework through:

- simulation experiments,
- four public empirical spatial datasets,
- multiple conditional quantiles,
- comparisons with QR, GWQR, QRF, and QXGBoost,
- Pinball Loss and Pseudo R²,
- residual spatial autocorrelation using Moran's I,
- local SHAP interpretation,
- bootstrap-based stability analysis.

## Main Findings

Across many of the evaluated settings, the geographically weighted tree-based approaches improved predictive performance and reduced residual spatial autocorrelation relative to global alternatives.

Local SHAP analysis also showed that feature contributions can vary across both geographic locations and quantile levels.

## Repository Structure

```text
spatial-quantile-machine-learning/
├── README.md
├── R/
│   ├── gw_qxgb.R
│   ├── gw_qxgb_cv.R
│   ├── metrics.R
│   └── README.md
├── simulation/
│   └── README.md
├── empirical/
│   └── README.md
├── figures/
│   └── README.md
├── data/
│   └── README.md
└── .gitignore
```

## Included Code

### `R/gw_qxgb.R`

A cleaned version of the thesis GW-QXGB implementation. It:

- builds a model matrix from an R formula,
- computes spatial distances,
- constructs global / Gaussian / adaptive bisquare weights,
- trains a local XGBoost model for each location,
- supports quantile prediction,
- returns fitted values and residuals.

### `R/gw_qxgb_cv.R`

A cleaned bandwidth-selection routine using leave-one-location-out prediction and golden-section search.

### `R/metrics.R`

Utility functions for Pinball Loss and Pseudo R².

## GW-QRF Source Status

The original thesis workspace references a separate local `gqrf.R` implementation, but that file is not currently present in the connected thesis repository. I have therefore **not fabricated a replacement implementation**.

The GW-QRF source can be added here later when the original file is available.

## Reproducibility Notes

The original research workspace contains exploratory scripts, local file paths, intermediate objects, and experiment-specific code. This public repository is intentionally organized as a cleaner portfolio version.

Empirical and simulation scripts will be migrated incrementally after removing machine-specific paths and confirming the provenance of each dataset.

## Thesis

The full thesis is available through my Notion portfolio:

[Data Analytics & Science Portfolio](https://app.notion.com/p/Data-Analytics-Science-Portfolio-160bcf9053b6800198faddd8f6e6a8ab)

## Tools

- R
- XGBoost
- Spatial statistics
- Quantile regression
- Spatial kernel weighting
- SHAP
- Parallel computing
- Moran's I
- Bootstrap
- Simulation

---
Portfolio repository maintained by **Yi-Ling Dai**.
