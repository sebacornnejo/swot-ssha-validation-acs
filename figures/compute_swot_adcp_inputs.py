# -*- coding: utf-8 -*-
"""
===============================================================================
 SWOT vs ADCP COMPUTATION ENGINE (PYTHON) -> INPUT FOR THE MATLAB FIGURE
===============================================================================
 What it does:
   1. Loads the ADCP, removes outliers and REMOVES THE TIDE with pyfes
      (u_tide, v_tide).
   2. Computes the Ekman depth from the wind.
   3. Builds the S1 (barotropic) scheme for each percentage 'pct'
      (average over layers with depth > Ekman and <= max_depth*pct).
   4. For the two SWOT cases (NEAR and N9 = 9-pixel average):
        SPATIAL (0..150 km) x TEMPORAL (2..20 days) x pct sweep
        - Gaussian filter applied to SWOT and to the ADCP
        - both series brought to a COMMON FINE GRID (6 h)
        - the COMMON-PERIOD MEAN is removed (true anomalies)
        - r (corr) and RMSE are computed
   5. Saves 'swot_figure_inputs.mat' with everything MATLAB needs to plot
      (optimal series + sensitivity matrices + metadata).

 Note on figure/metric consistency: the sweep, the optimum selection, the
 sensitivity matrices and the r/RMSE values are all computed ON THE SAME
 6-h fine grid that is plotted, so the r and RMSE shown in the figure
 correspond exactly to the curves drawn. (The smoothing reduces the
 effective degrees of freedom; this is a sensitivity analysis, not a
 formal test.)
===============================================================================
"""

import os
import re
import warnings
import numpy as np
import pandas as pd
import matplotlib.dates as mdates
from scipy import signal
from scipy.interpolate import interp1d
from scipy.io import savemat

try:
    from tqdm import tqdm
except ImportError:

    def tqdm(iterable, **kwargs):
        return iterable


try:
    import pyfes
except ImportError as e:
    raise ImportError(
        "pyfes is required to remove the tide from the ADCP series. "
        "Install it in the Python environment."
    ) from e

warnings.filterwarnings("ignore")

# ==============================================================================
# 0. CONFIGURATION
# ==============================================================================
print("=== SWOT vs ADCP ENGINE: NEAR + N9 | SPATIAL x TEMPORAL x S1 SWEEP ===")

# --- Station / physics ---
LAT_STATION = -52.14
LON_STATION = -66.67
DEPTHS_TO_DISCARD = [5.0, 9.0]

# --- Sweep ranges ---
TEMPORAL_RANGE = np.arange(2, 21, 1)  # days (Gaussian filter)
PCT_RANGE = np.arange(0.50, 1.00, 0.01)  # barotropic percentage (S1 scheme)

# *** Spatial values shown in the figure (1 colour per value) ***
#     If you change this, update the 'cmap_dac' palette in MATLAB as well.
SPATIAL_LIST_KM = [0, 30, 60, 90, 120, 150]

# --- Files ---
ADCP_FILE = "./ADCP/SWOT02.csv"
WIND_FILES = ["./DATA/WINDs/uv_wind_Y2023.nc", "./DATA/WINDs/uv_wind_Y2024.nc"]
BASE_DIR = "./outputs/"
OUTPUT_MAT = "swot_figure_inputs.mat"
OUTPUT_CSV = "summary_optima_near_n9.csv"

MIN_POINTS = 15  # minimum number of valid points in the common period

# m/s -> cm/s conversion (internal series are in m/s; the figure uses cm/s)
TO_CMS = 100.0


# SWOT filename resolution by mode and km
def resolve_swot_file(mode, km):
    """Returns the SWOT CSV filename for the given mode (near / n9) and km."""
    if mode == "near":
        if km == 0:
            return "series_sla_velgeo_lander_original.csv"
        return f"series_sla_velgeo_lander_filtrado_{km}km.csv"
    elif mode == "n9":
        if km == 0:
            return "series_sla_velgeo_lander_filtrado_0km_n09_v3.csv"
        return f"series_sla_velgeo_lander_filtrado_{km}km_n09_v3.csv"
    raise ValueError(mode)


# ==============================================================================
# 1. FILTERING FUNCTIONS
# ==============================================================================
def remove_outliers_percentile(data, p_low=0.1, p_high=99.9):
    if np.all(np.isnan(data)):
        return data
    low, high = np.nanpercentile(data, [p_low, p_high])
    clean = data.copy()
    clean[(data < low) | (data > high)] = np.nan
    return clean


def apply_butterworth_filter(data, cutoff_days, dt_sec, order=2):
    fs = 1.0 / dt_sec
    nyquist = 0.5 * fs
    cutoff_freq = 1.0 / (cutoff_days * 86400.0)
    if cutoff_freq / nyquist >= 0.99:
        return data
    b, a = signal.butter(order, cutoff_freq / nyquist, btype="low")
    series = pd.Series(data)
    if series.isna().any():
        series = (
            series.interpolate(limit_direction="both")
            .fillna(method="bfill")
            .fillna(method="ffill")
        )
    try:
        return signal.filtfilt(b, a, series.values)
    except Exception:
        return data


def apply_gaussian_filter(data, window_days, dt_sec):
    """Gaussian filter. Window = window_days; sigma = window/6 (covers +-3 sigma).
    Interpolates interior NaNs before filtering; leaves NaNs at the edges."""
    window_points = int(window_days * 86400 / dt_sec)
    if window_points < 3:
        return np.asarray(data, dtype=float)
    std_points = window_points / 6.0
    series = pd.Series(np.asarray(data, dtype=float))
    if series.isna().any():
        series = (
            series.interpolate(limit_direction="both")
            .fillna(method="bfill")
            .fillna(method="ffill")
        )
    return (
        series.rolling(window=window_points, win_type="gaussian", center=True)
        .mean(std=std_points)
        .values
    )


# ==============================================================================
# 2. HELPER FUNCTIONS (ADCP / TIDE / EKMAN / S1)
# ==============================================================================
def tidal_analysis_removal(dates, u, v):
    """Removes the tide (u_tide, v_tide) by harmonic analysis with pyfes >= 2026.x
    (wave_table_factory + WaveTableInterface API, no leap seconds)."""
    u_int = np.asarray(
        pd.Series(u).interpolate(limit=5).bfill().ffill().values, dtype="float64"
    )
    v_int = np.asarray(
        pd.Series(v).interpolate(limit=5).bfill().ffill().values, dtype="float64"
    )

    # The official FES 2026.x example uses datetime64[us]
    timestamps = dates.values.astype("datetime64[us]")

    constituents = [
        "Mm",
        "Mf",
        "Mtm",
        "Msqm",
        "2Q1",
        "Sigma1",
        "Q1",
        "Rho1",
        "O1",
        "MP1",
        "M11",
        "M12",
        "M13",
        "Chi1",
        "Pi1",
        "P1",
        "S1",
        "K1",
        "Psi1",
        "Phi1",
        "Theta1",
        "J1",
        "OO1",
        "MNS2",
        "Eps2",
        "2N2",
        "Mu2",
        "2MS2",
        "N2",
        "Nu2",
        "M2",
        "MKS2",
        "Lambda2",
        "L2",
        "2MN2",
        "T2",
        "S2",
        "R2",
        "K2",
        "MSN2",
        "Eta2",
        "2SM2",
        "MO3",
        "2MK3",
        "M3",
        "MK3",
        "N4",
        "MN4",
        "M4",
        "SN4",
        "MS4",
        "MK4",
        "S4",
        "SK4",
        "R4",
        "2MN6",
        "M6",
        "MSN6",
        "2MS6",
        "2MK6",
        "2SM6",
        "MSK6",
        "S6",
        "M8",
        "MSf",
        "Ssa",
        "Sa",
    ]

    # FES 2026.x API: factory + DARWIN engine
    wt = pyfes.wave_table_factory(pyfes.DARWIN, constituents)
    f_nodal, vu_nodal = wt.compute_nodal_modulations(timestamps)  # no leap seconds
    w_u = wt.harmonic_analysis(u_int, f_nodal, vu_nodal)
    w_v = wt.harmonic_analysis(v_int, f_nodal, vu_nodal)
    u_tide = wt.tide_from_tide_series(timestamps, w_u)  # no leap seconds
    v_tide = wt.tide_from_tide_series(timestamps, w_v)
    return u - u_tide, v - v_tide


def calculate_ekman_depth_series(dates, lat, wind_files):
    try:
        import xarray as xr

        ds = xr.open_mfdataset(wind_files, combine="by_coords")
        w_point = ds.sel(longitude=LON_STATION, latitude=LAT_STATION, method="nearest")
        t_wind = pd.to_datetime(w_point["valid_time"].values)
        u10, v10 = w_point["u10"].values, w_point["v10"].values
        dt_wind = (t_wind[1] - t_wind[0]).total_seconds()
        u10_filt = apply_butterworth_filter(u10, 2, dt_wind)
        v10_filt = apply_butterworth_filter(v10, 2, dt_wind)
        w_speed = np.sqrt(u10_filt**2 + v10_filt**2)
        t_adcp_nums = mdates.date2num(dates)
        t_wind_nums = mdates.date2num(t_wind)
        w_speed_on_adcp = interp1d(
            t_wind_nums, w_speed, bounds_error=False, fill_value="extrapolate"
        )(t_adcp_nums)
        lat_rad = np.abs(np.deg2rad(lat))
        ekman_depth = (7.6 * w_speed_on_adcp) / np.sqrt(np.sin(lat_rad))
        return np.clip(ekman_depth, 10, 150)
    except Exception as e:
        print(f"  Warning: Ekman depth computation failed ({e}). Using constant 35 m.")
        return np.full(len(dates), 35.0)


def compute_reference_velocity(
    u_res_dict, v_res_dict, ekman_series, params, avail_depths
):
    """S1 (barotropic) scheme: for each time step, averages the layers with
    depth > Ekman(t) and <= max_depth * pct."""
    n_times = len(ekman_series)
    max_depth = max(avail_depths)
    sorted_depths = sorted(avail_depths)
    u_ref = np.full(n_times, np.nan)
    v_ref = np.full(n_times, np.nan)
    limit = max_depth * params["pct"]
    for i in range(n_times):
        d_ek = ekman_series[i]
        v_bins = [d for d in sorted_depths if d > d_ek and d <= limit]
        if v_bins:
            u_vals = [u_res_dict[d][i] for d in v_bins]
            v_vals = [v_res_dict[d][i] for d in v_bins]
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                u_ref[i] = np.nanmean(u_vals)
                v_ref[i] = np.nanmean(v_vals)
    return u_ref, v_ref


def to_matlab_datenum(t_seconds, base_time):
    """Converts seconds-since-base_time to MATLAB datenum."""
    t_dt = base_time + pd.to_timedelta(np.asarray(t_seconds, dtype=float), unit="s")
    days = (t_dt - pd.Timestamp("1970-01-01")) / pd.Timedelta(days=1)
    return 719529.0 + np.asarray(days, dtype=float)


# ==============================================================================
# 3. ADCP LOADING AND PROCESSING (shared by NEAR and N9)
# ==============================================================================
print("\n--- 1. LOADING ADCP DATA ---")
if not os.path.exists(ADCP_FILE):
    raise FileNotFoundError(f"Missing {ADCP_FILE}")

df = pd.read_csv(ADCP_FILE, delimiter=";")
df.columns = df.columns.str.strip()
df["DateTime"] = pd.to_datetime(
    df["DateTime"], format="%d-%m-%Y %H:%M:%S", errors="coerce"
)
df.dropna(subset=["DateTime"], inplace=True)
df = df.iloc[24:-1080].reset_index(drop=True)

# --- Diagnostics: show what is actually loaded ---
print(
    f"Columns detected ({len(df.columns)}): "
    f"{list(df.columns)[:12]}{' ...' if len(df.columns) > 12 else ''}"
)
print(f"Rows after DateTime cleaning: {len(df)}")
speed_cols = [c for c in df.columns if c.startswith("Speed#")]
print(f"'Speed#' columns found: {len(speed_cols)} -> {speed_cols[:6]}")
if len(speed_cols) == 0:
    raise RuntimeError(
        "No 'Speed#' columns. The delimiter is probably wrong "
        "(try delimiter=',' in the ADCP pd.read_csv) or the names differ."
    )

u_res_dict, v_res_dict, valid_depths_list = {}, {}, []
print(f"Discarding levels: {DEPTHS_TO_DISCARD}")
n_fail = 0
for col in speed_cols:
    try:
        d = float(col.split("(")[1].split("m")[0])
    except Exception as e:
        print(f"  [skip] could not parse depth from '{col}': {e}")
        continue
    if d < 3.0:
        continue
    if any(np.isclose(d, disc, atol=0.1) for disc in DEPTHS_TO_DISCARD):
        continue
    dir_col = col.replace("Speed", "Dir")
    if dir_col not in df.columns:
        print(f"  [skip] missing direction column for {col} (expected '{dir_col}')")
        continue
    try:
        s = df[col].values
        rad = np.deg2rad(df[dir_col].values)
        u_raw, v_raw = s * np.sin(rad), s * np.cos(rad)
        u_cl = remove_outliers_percentile(u_raw)
        v_cl = remove_outliers_percentile(v_raw)
        u_r, v_r = tidal_analysis_removal(df["DateTime"], u_cl, v_cl)  # pyfes 2026.x
        u_res_dict[d], v_res_dict[d] = u_r, v_r
        valid_depths_list.append(d)
    except Exception as e:
        n_fail += 1
        print(f"  [FAIL] depth {d} m: {type(e).__name__}: {e}")

valid_depths_list.sort()
print(
    f"Valid depths ({len(valid_depths_list)}): {valid_depths_list}  | failures: {n_fail}"
)
if not valid_depths_list:
    raise RuntimeError(
        "valid_depths_list is empty. Check the [FAIL]/[skip] messages above."
    )

ekman_series = calculate_ekman_depth_series(df["DateTime"], LAT_STATION, WIND_FILES)
dt_adcp = (df["DateTime"].iloc[1] - df["DateTime"].iloc[0]).total_seconds()
base_time = df["DateTime"].min()
t_adcp_sec = (df["DateTime"] - base_time).dt.total_seconds().values

# ==============================================================================
# 4. S1 CACHE (raw)  +  FILTERED ADCP CACHE per (pct, days)
# ==============================================================================
print("\n--- 2. S1 (barotropic) PRECOMPUTATION per pct ---")
s1_raw_store = {}  # pct -> (u_ref_raw, v_ref_raw)
for pct in tqdm(PCT_RANGE, desc="S1 cache"):
    u_r, v_r = compute_reference_velocity(
        u_res_dict, v_res_dict, ekman_series, {"pct": pct}, valid_depths_list
    )
    s1_raw_store[round(float(pct), 4)] = (u_r, v_r)

print("--- 3. FILTERED ADCP (Gaussian) PRECOMPUTATION per (pct, days) ---")
# ADCP filtering only depends on (pct, days) -> reused by NEAR and N9.
adcp_filt_cache = {}
for pct in tqdm(PCT_RANGE, desc="Filtered ADCP"):
    key_pct = round(float(pct), 4)
    u_r, v_r = s1_raw_store[key_pct]
    for days in TEMPORAL_RANGE:
        u_f = apply_gaussian_filter(u_r, days, dt_adcp)
        v_f = apply_gaussian_filter(v_r, days, dt_adcp)
        adcp_filt_cache[(key_pct, int(days))] = (u_f, v_f)


# ==============================================================================
# 5. SWEEP PER MODE (NEAR / N9) -> SPATIAL x TEMPORAL x pct
# ==============================================================================
def process_mode(mode):
    print(f"\n=== MODE: {mode.upper()} ===")
    num_days = len(TEMPORAL_RANGE)
    num_spat = len(SPATIAL_LIST_KM)
    sens_matrix = np.full((num_days, num_spat), np.nan)  # [days x spatial]

    best = {"score": -np.inf}
    summary_rows = []

    for s_idx, km in enumerate(SPATIAL_LIST_KM):
        fname = resolve_swot_file(mode, km)
        fpath = os.path.join(BASE_DIR, fname)
        if not os.path.exists(fpath):
            print(f"  [{km} km] NOT found: {fname}")
            continue
        print(f"  [{km} km] {fname}")

        df_swot = pd.read_csv(fpath)
        t_swot = pd.to_datetime(
            df_swot["time_datenum"] - 719529, unit="D", origin="unix"
        )
        t_swot_sec = (t_swot - base_time).dt.total_seconds().values

        # Common fine grid (6 h) for this spatial file
        dt_fine = 6 * 3600
        t_fine = np.arange(
            (t_swot.min() - base_time).total_seconds(),
            (t_swot.max() - base_time).total_seconds(),
            dt_fine,
        )

        # Raw SWOT interpolated onto the fine grid
        u_geo_fine = interp1d(
            t_swot_sec, df_swot["u_geo_re"], bounds_error=False, fill_value=np.nan
        )(t_fine)
        v_geo_fine = interp1d(
            t_swot_sec, df_swot["v_geo_re"], bounds_error=False, fill_value=np.nan
        )(t_fine)

        for d_idx, days in enumerate(TEMPORAL_RANGE):
            days = int(days)
            # Filtered SWOT (on the fine grid)
            u_swot_f = apply_gaussian_filter(u_geo_fine, days, dt_fine)
            v_swot_f = apply_gaussian_filter(v_geo_fine, days, dt_fine)

            best_pct_score = -np.inf
            for pct in PCT_RANGE:
                key_pct = round(float(pct), 4)
                u_a_nat, v_a_nat = adcp_filt_cache[(key_pct, days)]
                # Filtered ADCP -> fine grid (same as SWOT)
                u_a = interp1d(t_adcp_sec, u_a_nat, bounds_error=False)(t_fine)
                v_a = interp1d(t_adcp_sec, v_a_nat, bounds_error=False)(t_fine)

                mask = (
                    ~np.isnan(u_a)
                    & ~np.isnan(u_swot_f)
                    & ~np.isnan(v_a)
                    & ~np.isnan(v_swot_f)
                )
                if np.sum(mask) < MIN_POINTS:
                    continue

                # --- Anomalies: remove the COMMON-PERIOD mean ---
                ua = u_a[mask] - np.mean(u_a[mask])
                us = u_swot_f[mask] - np.mean(u_swot_f[mask])
                va = v_a[mask] - np.mean(v_a[mask])
                vs = v_swot_f[mask] - np.mean(v_swot_f[mask])

                if (
                    np.std(ua) == 0
                    or np.std(us) == 0
                    or np.std(va) == 0
                    or np.std(vs) == 0
                ):
                    continue

                r_u = np.corrcoef(ua, us)[0, 1]
                r_v = np.corrcoef(va, vs)[0, 1]
                score = (r_u + r_v) / 2.0

                if score > best_pct_score:
                    best_pct_score = score

                if score > best["score"]:
                    rmse_u = np.sqrt(np.mean((ua - us) ** 2)) * TO_CMS
                    rmse_v = np.sqrt(np.mean((va - vs) ** 2)) * TO_CMS

                    # Full series on the fine grid, NaN outside the common period
                    u_adcp_full = np.full_like(u_a, np.nan)
                    u_swot_full = np.full_like(u_a, np.nan)
                    v_adcp_full = np.full_like(u_a, np.nan)
                    v_swot_full = np.full_like(u_a, np.nan)
                    u_adcp_full[mask] = ua * TO_CMS
                    u_swot_full[mask] = us * TO_CMS
                    v_adcp_full[mask] = va * TO_CMS
                    v_swot_full[mask] = vs * TO_CMS

                    t_dn = to_matlab_datenum(t_fine, base_time)
                    valid_dn = t_dn[mask]

                    best = {
                        "score": score,
                        "t": t_dn.reshape(-1, 1),
                        "u_adcp": u_adcp_full.reshape(-1, 1),
                        "u_swot": u_swot_full.reshape(-1, 1),
                        "v_adcp": v_adcp_full.reshape(-1, 1),
                        "v_swot": v_swot_full.reshape(-1, 1),
                        "r_u": float(r_u),
                        "r_v": float(r_v),
                        "rmse_u": float(rmse_u),
                        "rmse_v": float(rmse_v),
                        "km": float(km),
                        "days": float(days),
                        "pct": float(pct),
                        "xlim_min": float(valid_dn.min()),
                        "xlim_max": float(valid_dn.max()),
                    }

            sens_matrix[d_idx, s_idx] = (
                best_pct_score if best_pct_score > -np.inf else np.nan
            )

        # best per file (for the summary CSV)
        col = sens_matrix[:, s_idx]
        if np.any(np.isfinite(col)):
            d_best = int(TEMPORAL_RANGE[np.nanargmax(col)])
            summary_rows.append(
                {
                    "Mode": mode,
                    "km": km,
                    "Best_day": d_best,
                    "Best_score": np.nanmax(col),
                }
            )

    if best["score"] == -np.inf:
        raise RuntimeError(
            f"No results generated for mode '{mode}'. "
            f"Check the paths in {BASE_DIR}."
        )

    print(
        f"  -> OPTIMUM {mode}: {best['km']:.0f} km | {best['days']:.0f} days | "
        f"{best['pct']*100:.0f}% Baro | score={best['score']:.3f} "
        f"(r_u={best['r_u']:.3f}, r_v={best['r_v']:.3f})"
    )

    best.pop("score", None)
    return {"sens_matrix": sens_matrix, "best": best}, summary_rows


results = {}
all_summary = []
for mode in ["near", "n9"]:
    results[mode], rows = process_mode(mode)
    all_summary.extend(rows)

# ==============================================================================
# 6. EXPORT -> .mat for MATLAB + summary CSV
# ==============================================================================
mat_dict = {
    "spatial_list": np.asarray(SPATIAL_LIST_KM, dtype=float).reshape(-1, 1),
    "temporal_range": np.asarray(TEMPORAL_RANGE, dtype=float).reshape(-1, 1),
    "near": results["near"],
    "n9": results["n9"],
}
savemat(OUTPUT_MAT, mat_dict, long_field_names=True)
print(f"\nSaved: {OUTPUT_MAT}")

pd.DataFrame(all_summary).to_csv(OUTPUT_CSV, index=False)
print(f"Saved: {OUTPUT_CSV}")
print("\nDone. Run the MATLAB plotting script with this .mat on the path.")
