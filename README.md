# Analysis Code — *Impact of dynamic atmospheric correction choice on SWOT sea surface height over continental shelves: a case study in the Southwestern Atlantic*

**Authors:** Cornejo-Guzmán S.E., Saraceno M., Ruiz-Etcheverry L.A., Birol F., Lyard F.

---

## Overview

This repository contains the scripts used to perform the tidal and atmospheric correction assessment, satellite/in-situ validation, spectral analysis, and spatiotemporal characterization of coastally trapped waves presented in the manuscript above.

The analysis combines **MATLAB** scripts (validation statistics, spectral analysis, CEOF decomposition) and **Python** scripts (tidal constituent extraction and tidal elevation prediction).

---

## Repository Contents

### Python scripts

| File | Description |
|---|---|
| `extract_insitu_harmonics_pyfes.py` | Performs harmonic analysis on in-situ SSH time series using **pyFES** to extract tidal constituents at each station. |
| `tidesFES2022bconstituents.py` | Extracts tidal harmonic constituents (amplitude and phase) from the **FES2022b** global model at each station using **pyTMD**. |
| `tidesTPXO10v2constituents.py` | Extracts constituents from the **TPXO10-atlas-v2** model using **pyTMD**. |
| `tidesEOT20constituents.py` | Extracts constituents from the **EOT20** global ocean tide model using **pyTMD**. |
| `tidesGOT55constituents.py` | Extracts constituents from **GOT5.5**, combining ocean and load tides via complex addition. |
| `tidesfromFES2022b.py` | Generates high-resolution (1-minute) tidal elevation time series for each station over the SWOT CalVal period using **FES2022b** via **pyTMD**. |

### MATLAB scripts

| File | Description |
|---|---|
| `evaluate_tidal_rss_misfit.m` | Computes Root Sum Square (RSS) misfits between modelled and in-situ tidal constituents for all five tidal models and five stations. Produces the bar chart figure and an Excel summary table. |
| `evaluate_swot_insitu_stats.m` | Computes Pearson correlation (*r*), RMSE, and R² between SWOT and in-situ SSHA for each DAC configuration. Applies iterative variance-based quality control and exports time-series figures and a statistics summary. |
| `compute_swot_insitu_spectra.m` | Computes variance-preserving Power Spectral Density (PSD) estimates for each station and a regional composite using Welch's method with a Tukey window. Includes red-noise significance testing via a first-order autoregressive model. |
| `analyze_swot_ceof_hovmoller.m` | Applies multi-step high-pass filtering (adapted from Dinápoli et al., 2025), performs Complex EOF (CEOF) decomposition, and constructs Hovmöller diagrams for the default and SIROCCO-corrected SSHA fields. |
| `EOF_.m` | Helper function implementing SVD-based EOF decomposition, used by `analyze_swot_ceof_hovmoller.m`. |

---

## Dependencies

### Python
- [`pyFES`](https://github.com/CNES/aviso-fes) — harmonic analysis
- [`pyTMD`](https://github.com/tsutterley/pyTMD) — tidal model extraction and prediction
- `numpy`, `pandas`

### MATLAB
- MATLAB R2021b or later (Signal Processing Toolbox for `filtfilt`, `butter`)

### Tidal model files (not included)
The following model directories must be provided locally:

| Model | Directory |
|---|---|
| FES2022b | `tide_models/fes2022b/` |
| TPXO10-atlas-v2 | `tide_models/` |
| EOT20 | `tide_models/` |
| GOT5.5 | `tide_models/GOT5.5/` |

---

## Stations

All scripts operate on the five in-situ stations used in the study:

| Station | Longitude | Latitude |
|---|---|---|
| Carina | 67.22° W | 52.76° S |
| Vega-Pleyade | 67.75° W | 53.30° S |
| Lander | 66.67° W | 52.14° S |
| Puerto Deseado | 65.92° W | 47.75° S |
| San Matías (SBE26) | 63.78° W | 41.18° S |

---

## License

Scripts are provided for reproducibility purposes. Please cite the associated manuscript if you use this code.
