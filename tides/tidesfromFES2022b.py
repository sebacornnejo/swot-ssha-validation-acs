#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Tidal Elevation Prediction using FES2022b (pyTMD)

Description:
    Generates high-resolution (1-minute) tidal elevation time series for 
    selected in-situ monitoring stations over the Argentine Continental Shelf.
    The temporal span covers the SWOT Cal/Val phase with a 7-day padding 
    at both boundaries. This padding ensures sufficient data coverage to 
    prevent boundary/edge effects during subsequent time-series filtering 
    and spectral analyses.
"""

import os
import sys
import lzma
import shutil
import numpy as np
import pandas as pd
import pyTMD.compute

# 1. Configuration and Directories
BASE_DIR  = '/Users/sebacornnejo/Google Drive/Mi unidad/DoctoradoDCAO/SWOT/tide_models'
OCEAN_DIR = os.path.join(BASE_DIR, 'fes2022b', 'ocean_tide')

stations = {
    'lander':         {'lon': -66.67000, 'lat': -52.14000},
    'puerto_deseado': {'lon': -65.81700, 'lat': -47.75000},
    'san_matias':     {'lon': -63.48030, 'lat': -41.18220},
    'hidra_norte':    {'lon': -68.22095, 'lat': -52.82061},
    'carina':         {'lon': -67.21960, 'lat': -52.75718},
    'vega_pleyade':   {'lon': -67.74600, 'lat': -53.29900},
}

# 2. Decompress FES2022b netCDF archives if required
if not os.path.isdir(OCEAN_DIR):
    sys.exit(f"ERROR: Model directory not found: {OCEAN_DIR}")
    
for filename in os.listdir(OCEAN_DIR):
    if filename.endswith('.nc.xz'):
        src_file = os.path.join(OCEAN_DIR, filename)
        dst_file = os.path.join(OCEAN_DIR, filename[:-3]) # Remove .xz extension
        if not os.path.exists(dst_file):
            print(f"Decompressing model file: {filename} ...")
            with lzma.open(src_file, 'rb') as f_in, open(dst_file, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

# 3. Generate continuous 1-minute time vector with +/- 1-week padding
# Original SWOT Cal/Val limits: 2023-03-29 03:00 UTC to 2023-07-09 11:00 UTC
START_DATE = pd.Timestamp('2023-03-22 03:00:00', tz='UTC')   # - 7 days padding
END_DATE   = pd.Timestamp('2023-07-16 11:00:00', tz='UTC')   # + 7 days padding

time_vector = pd.date_range(start=START_DATE, end=END_DATE, freq='1min')
print(f"Time vector generated: {len(time_vector):,} timestamps ({START_DATE} to {END_DATE})")

# 4. Compute elapsed time in seconds relative to the FES epoch (1992-01-01 UTC)
fes_epoch  = pd.Timestamp(1992, 1, 1, tz='UTC')
delta_time = (time_vector - fes_epoch).total_seconds().values

# 5. Loop over stations and compute tidal elevations
print("\n------------------------------------------------------------")
print("Initiating tidal predictions for in-situ stations")
print("------------------------------------------------------------\n")

for st_key, coords in stations.items():
    display_name = st_key.replace('_', ' ').title()
    print(f"Processing station: {display_name} ({coords['lat']:.5f} N, {coords['lon']:.5f} E)")
    
    try:
        # Compute tidal elevations using pyTMD
        z = pyTMD.compute.tide_elevations(
            coords['lon'], coords['lat'], delta_time,
            DIRECTORY   = BASE_DIR,
            MODEL       = 'FES2022',
            EPSG        = 4326,
            EPOCH       = (1992, 1, 1, 0, 0, 0),
            TYPE        = 'drift',
            TIME        = 'UTC',
            METHOD      = 'spline',
            CROP        = True,
            FILL_VALUE  = np.nan,
        )

        # Structure the output into a DataFrame
        df_out = pd.DataFrame({
            'Time_UTC':    time_vector.tz_convert('UTC').tz_localize(None),
            'Tide_Elev_m': z.flatten(),
            'Latitude':    coords['lat'],
            'Longitude':   coords['lon'],
        })

        # Export to Excel
        output_filename = f'Tide_Series_{st_key}_1min_FES2022b.xlsx'
        df_out.to_excel(output_filename, index=False)
        
        # Calculate dynamic range safely (ignoring NaNs)
        z_min = np.nanmin(z)
        z_max = np.nanmax(z)
        
        print(f"  -> Exported: {output_filename} | Tidal Range: [{z_min:.3f} m to {z_max:.3f} m]\n")

    except Exception as e:
        print(f"  -> ERROR computing tides for {display_name}: {e}\n")

print("================== PROCESSING COMPLETED ==================")