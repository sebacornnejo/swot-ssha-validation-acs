"""
Tidal Harmonic Constituent Extraction using GOT5.5

Description:
    Extracts tidal amplitude and phase for the primary harmonic constituents
    from the Global Ocean Tide (GOT5.5) model. This script automatically 
    processes multiple in-situ monitoring stations over the Argentine 
    Continental Shelf. 
    
    Since the total tide measured by bottom pressure sensors or tide gauges 
    includes both the ocean tide and the solid earth load tide, this routine 
    extracts both components, sums them in the complex plane, and outputs 
    the total geocentric tidal constituents (amplitude in meters, phase in degrees).
"""

import numpy as np
import pandas as pd
from pathlib import Path
import pyTMD.io as ptio

# 1. Configuration and Directories
BASE_DIR = Path("tide_models/GOT5.5")
OCEAN_DIR = BASE_DIR / "ocean_tides"
LOAD_DIR  = BASE_DIR / "load_tides"

# Define in-situ stations (Standard longitudes, converted internally)
stations = {
    'Carina':   {'lon': -67.21960, 'lat': -52.75718},
    'Lander':   {'lon': -66.67000, 'lat': -52.14000},
    'PDeseado': {'lon': -65.91420, 'lat': -47.75360},
    'Vega':     {'lon': -67.74600, 'lat': -53.29900},
    'sbe26SM':  {'lon': -63.77680, 'lat': -41.18190}
}

# 2. Validate model files
ocean_files = sorted(OCEAN_DIR.glob("*.nc"))
load_files  = sorted(LOAD_DIR.glob("*.nc"))

if not ocean_files or not load_files:
    raise FileNotFoundError(f"Model NetCDF files not found in {OCEAN_DIR} or {LOAD_DIR}")

print("------------------------------------------------------------")
print("Initiating GOT5.5 Harmonic Constituent Extraction")
print("------------------------------------------------------------\n")

# 3. Iterative Processing per Station
for st_name, coords in stations.items():
    print(f"Processing station: {st_name} ({coords['lat']:.5f} N, {coords['lon']:.5f} E)")
    
    # The GOT model utilizes longitude in the [0, 360) range
    lon_got = np.array([360.0 + coords['lon']] if coords['lon'] < 0 else [coords['lon']], dtype=float)
    lat_got = np.array([coords['lat']], dtype=float)
    
    try:
        # 3.1 Extract Ocean Tide constituents
        amp_ocean, ph_ocean, cons_ocean = ptio.GOT.extract_constants(
            lon_got, lat_got, ocean_files,
            grid="netcdf",
            method="spline",      # Smooth interpolation
            extrapolate=True,     # Nearest-neighbor for coastal boundary mismatch
            cutoff=np.inf,        # Unrestricted extrapolation distance
            scale=0.01,           # Conversion from cm to meters
            compressed=False
        )
        
        # 3.2 Extract Load Tide constituents
        amp_load, ph_load, cons_load = ptio.GOT.extract_constants(
            lon_got, lat_got, load_files,
            grid="netcdf",
            method="spline",
            extrapolate=True,
            cutoff=np.inf,
            scale=0.01,
            compressed=False
        )
        
        # Sanity check: Ensure constituent arrays align perfectly
        if cons_ocean != cons_load:
            raise ValueError(f"Constituent mismatch between Ocean and Load tides for {st_name}")
        
        # 3.3 Complex addition of Ocean and Load tides
        # Complex representation: Z = A * exp(-i * phase)
        c_ocean = amp_ocean[0] * np.exp(-1j * np.deg2rad(ph_ocean[0]))
        c_load  = amp_load[0]  * np.exp(-1j * np.deg2rad(ph_load[0]))
        c_total = c_ocean + c_load
        
        # Retrieve total amplitude (m) and phase (degrees)
        amp_total = np.abs(c_total)
        ph_total  = np.mod(-np.rad2deg(np.angle(c_total)), 360.0)
        
        # 4. Tabular Export
        df_out = pd.DataFrame({
            "constituent": cons_ocean,
            "amplitude":   amp_total,
            "phase_deg":   ph_total
        })
        
        output_filename = f"tides_constituents_{st_name}_GOT5_5.xlsx"
        df_out.to_excel(output_filename, index=False)
        
        print(f"  -> Successfully exported: {output_filename} ({len(df_out)} constituents)\n")
        
    except Exception as e:
        print(f"  -> ERROR processing {st_name}: {e}\n")

print("================== EXTRACTION COMPLETED ==================")