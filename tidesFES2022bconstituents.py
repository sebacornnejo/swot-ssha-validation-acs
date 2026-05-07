"""
FES2022b Tidal Constituents Extraction

Description:
    Extracts the amplitude and phase of tidal harmonic constituents from 
    the FES2022b global tide model for a defined set of in-situ monitoring 
    stations. The script utilizes bilinear interpolation to estimate the 
    values at the exact station coordinates and exports the results to 
    individual Excel files for subsequent validation and root-sum-square 
    (RSS) misfit analysis.
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import pyTMD.io as ptio

def main():
    # 1. Configuration and Model Directories
    base_dir = Path("tide_models/fes2022b")
    
    # Recursively search for all netCDF files
    model_files = sorted(base_dir.rglob("*.nc"))
    if not model_files:
        sys.exit(f"[ERROR] No .nc files found in directory: {base_dir}")
    print(f"Located {len(model_files)} netCDF files in {base_dir}")

    # 2. Station Definitions
    # Longitudes are natively mapped to [0, 360) by pyTMD, inputs can be [-180, 180]
    stations = {
        'Carina':   {'lon': -67.21960,  'lat': -52.75718},
        'Vega':     {'lon': -67.74600,  'lat': -53.29900},
        'Lander':   {'lon': -66.67000,  'lat': -52.14000},
        'PDeseado': {'lon': -65.899999, 'lat': -47.761325},
        'sbe26SM':  {'lon': -63.77680,  'lat': -41.18190}
    }

    # 3. Read Constituent List
    # Retrieves the list of available tidal constituents (e.g., 'M2', 'S2', 'O1')
    print("Reading tidal constituents list from model files...")
    constituents_list = ptio.FES.read_constants(
        model_files,
        type="z",
        version="FES2022",
        compressed=False
    ).fields
    
    print(f"Processing {len(constituents_list)} constituents: {constituents_list[:5]}...")

    # 4. Iterate over stations and extract harmonic data
    print("\n------------------------------------------------------------")
    print("Extracting Harmonic Constituents")
    print("------------------------------------------------------------\n")
    
    for st_name, coords in stations.items():
        lon_array = np.array([coords['lon']], dtype=float)
        lat_array = np.array([coords['lat']], dtype=float)
        
        print(f"Processing station: {st_name} ({coords['lat']:.5f} N, {coords['lon']:.5f} E)")
        
        # Extract amplitude and phase using bilinear spatial interpolation
        # Note: scale=0.01 converts native cm amplitudes to meters for elevation ('z')
        amp, ph = ptio.FES.extract_constants(
            lon_array, lat_array, model_files,
            type="z",
            version="FES2022",
            method="bilinear", 
            compressed=False,
            scale=0.01 
        )
        
        if amp.size == 0 or ph.size == 0:
            print(f"  -> [WARNING] Interpolation returned empty arrays for {st_name}. Check grid coverage.")
            continue

        # 5. Format and Export to Excel
        df_export = pd.DataFrame({
            "constituent": constituents_list,
            "amplitude":   amp[0],             # Amplitude in meters
            "phase_deg":   np.mod(ph[0], 360)  # Phase in degrees wrapped to [0, 360)
        })
        
        output_filename = f"tides_constituents_{st_name}_FES2022.xlsx"
        df_export.to_excel(output_filename, index=False)
        print(f"  -> Exported: {output_filename} ({len(df_export)} constituents)\n")

    print("================== PROCESSING COMPLETED ==================")

if __name__ == "__main__":
    main()