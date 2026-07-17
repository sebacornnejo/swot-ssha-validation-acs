# Analysis Code — *On the Accuracy of SWOT Sea Surface Height Over Large Continental Shelves*

**Authors:** Cornejo-Guzmán S.E., Saraceno M., Ruiz-Etcheverry L.A., Birol F., Lyard F.

---

## Overview

This repository contains the scripts used in the manuscript above to validate SWOT L3 sea surface height over the Argentine Continental Shelf against five in situ records located under pass 007 during the 1-day CalVal phase. The analysis covers:

- the assessment of five tide models (FES2022b, TPXO10-atlas-v2, EOT20, GOT5.5, and the regional MSAS) against in situ harmonic constituents;
- the intercomparison of inverse-barometer (IB) references built from different atmospheric pressure products;
- the cross-spectral assessment (coherence, phase lag, variance ratio) of five dynamic atmospheric correction (DAC) products against a common ERA5-MSLP IB reference, on both the SWOT and the in situ sides;
- the spatiotemporal characterization of coastal trapped waves (Hovmöller diagrams, Lomb-Scargle spectral analysis, phase velocities) in the 2–10-day band;
- the validation of SWOT-derived geostrophic velocities against a bottom-mounted Aquadopp current profiler.

The repository combines **MATLAB** scripts (validation statistics, cross-spectral analysis, Hovmöller/spectral analysis, figure generation) and **Python** scripts (tidal constituent extraction and prediction, and the SWOT–ADCP computation engine).

---

## Repository Contents

### `figures/` — manuscript figures

Each script saves its figure under the same name as the manuscript figure it produces (`Fig01.png` … `Fig08.png`, `FigS1.png` … `FigS3.png`; `Fig05.m` saves its three panels as `Fig05a/b/c`). Scripts named `FigS*_SI.m` generate the Supporting Information (SI) figures.

| File | Description |
|---|---|
| `Fig01.m` | Study-area map: SWOT passes 007/022 swaths and conventional altimetry SLA over the Southwestern Atlantic, with the in situ stations and the 200-m isobath. |
| `Fig02.m` | Total Root Sum Square (RSS) misfit of the five tide models against in situ constituents (pyFES harmonic analysis), grouped bar chart per station plus summary table. |
| `Fig03.m` | Skill of the inverse-barometer series built from different pressure products (ERA5 MSLP, CFSv2, ECMWF IFS HRES, ERA5 sp) against the SWOT L2 IB reference: correlation, RMSE, centred RMSD, and variability ratio per station. |
| `Fig04.m` | Cross-spectral synthesis of SWOT vs in situ SSHA per correction: mean squared coherence (2–10 d), coherence-weighted phase lag, and variance ratio. |
| `Fig05.m` | R² bubble maps of SWOT vs in situ SSHA: (a) no correction, (b) tide + IB, (c) the five DAC products per station. |
| `Fig06.m` | Latitude–time Hovmöller of high-pass-filtered SSHA (IB-only) with coastal-trapped-wave phase-velocity estimates; manual axes/date-tick layout. |
| `Fig07.m` | Same pipeline as `Fig06.m` in a two-panel subplot layout, with extended per-latitude Lomb-Scargle diagnostics (peak-rank medians and spectral energy partition) printed to the console. |
| `Fig08.m` | SWOT-derived geostrophic velocities vs Aquadopp currents under the tide + IB correction: time series and spatial/temporal-filter sensitivity (plots the output of `compute_swot_adcp_inputs.py`). |
| `FigS1_SI.m` | Supporting Information Figure S1. Cross-spectral assessment of the five DACs against the ERA5-MSLP IB reference per station: coherence of the correction difference, weighted phase lag, and fractional variance reduction (in situ and SWOT sides). |
| `FigS2_SI.m` | Supporting Information Figure S2. Hovmöller diagrams of DAC-corrected SSHA (six panels: IB-only plus five DACs) showing the degradation of the coastal-trapped-wave signal. |
| `FigS3_SI.m` | Supporting Information Figure S3. Counterpart of `Fig08.m` under the tide + full DAC correction, same axes for direct contrast. |
| `compute_swot_insitu_stats.m` | SWOT vs in situ SSHA validation statistics: variance-based quality control, common-period centring, and per-DAC r / RMSE / R² (the R² values plotted in `Fig05.m`), plus per-station SSHA comparison time series and a statistics summary table. |
| `compute_swot_adcp_inputs.py` | Computation engine for `Fig08.m`/`FigS3_SI.m`: ADCP tide removal (pyFES), Ekman-depth estimate from ERA5 winds, barotropic reference velocity, and the spatial × temporal × barotropic-fraction sweep against SWOT geostrophic velocities. Exports a `.mat` file consumed by the MATLAB plotting scripts. |

### `tides/` — tidal constituent extraction and prediction

| File | Description |
|---|---|
| `extract_insitu_harmonics_pyfes.py` | Harmonic analysis of the in situ SSH time series using **pyFES** to extract tidal constituents at each station. |
| `tidesFES2022bconstituents.py` | Extracts tidal constituents (amplitude and phase) from the **FES2022b** global model at each station using **pyTMD**. |
| `tidesTPXO10v2constituents.py` | Extracts constituents from the **TPXO10-atlas-v2** model using **pyTMD**. |
| `tidesEOT20constituents.py` | Extracts constituents from the **EOT20** global ocean tide model using **pyTMD**. |
| `tidesGOT55constituents.py` | Extracts constituents from **GOT5.5**, combining ocean and load tides via complex addition. |
| `tidesfromFES2022b.py` | Generates high-resolution (1-minute) tidal elevation time series for each station over the SWOT CalVal period using **FES2022b** via **pyTMD**. |

---

## Dependencies

### Python
- [`pyFES`](https://github.com/CNES/aviso-fes) — harmonic analysis and tide removal
- [`pyTMD`](https://github.com/tsutterley/pyTMD) — tidal model extraction and prediction
- `numpy`, `pandas`, `scipy`, `xarray`
- `tqdm` (optional, progress bars)

### MATLAB
- MATLAB R2021b or later (Signal Processing Toolbox for `filtfilt`, `butter`, `mscohere`, `cpsd`)
- [`m_map`](https://www.eoas.ubc.ca/~rich/map.html) mapping toolbox (map figures)

### Tidal model files (not included)
The following model directories must be provided locally:

| Model | Directory |
|---|---|
| FES2022b | `tide_models/fes2022b/` |
| TPXO10-atlas-v2 | `tide_models/` |
| EOT20 | `tide_models/` |
| GOT5.5 | `tide_models/GOT5.5/` |

Input data (SWOT L3 swaths, in situ series, ADCP, ERA5 fields) are not distributed with this repository; see the Open Research section of the manuscript for the data sources.

---

## Stations

All scripts operate on the five in situ stations used in the study:

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
