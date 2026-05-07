"""
Tidal Constituent Extraction using TPXO10-atlas-v2 (pyTMD)

Description:
    Extracts harmonic constituents (amplitude and phase) from the 
    TPXO10-atlas-v2 tidal model for selected in-situ monitoring stations 
    over the Argentine Continental Shelf. The extracted constituents are 
    exported to Excel files for subsequent comparative analysis (RSS misfit).
"""

import numpy as np
import pandas as pd
import pyTMD

# 1. Configuration and Directories
MODEL_DIR = "tide_models/"

# 2. Station Definitions
# Dictionary containing the target stations and their coordinates.
# Note: Keys are formatted to exactly match the expected filenames in MATLAB.
stations = {
    'Carina':   {'lon': -67.21960, 'lat': -52.75718},
    'Vega':     {'lon': -67.74600, 'lat': -53.29900}, # Vega-Pleyade
    'Lander':   {'lon': -66.67000, 'lat': -52.14000},
    'PDeseado': {'lon': -65.91420, 'lat': -47.75360}, # P. Deseado
    'sbe26SM':  {'lon': -63.77680, 'lat': -41.18190}  # SBE26 San Matias
}

# 3. Load Database and Initialize Model
print("Loading TPXO10-atlas-v2 tidal database...")
try:
    db = pyTMD.io.load_database()
    m  = pyTMD.io.model(directory=MODEL_DIR)
    
    # Select the elevation model for TPXO10-atlas-v2
    pytmd_model = m.elevation("TPXO10-atlas-v2")
except Exception as e:
    print(f"ERROR: Failed to load the tidal model. Check MODEL_DIR. Details: {e}")
    exit(1)

print("\n------------------------------------------------------------")
print("Extracting Harmonic Constituents")
print("------------------------------------------------------------\n")

# 4. Iterative Processing per Station
for st_key, coords in stations.items():
    print(f"Processing station: {st_key} ({coords['lat']:.5f} N, {coords['lon']:.5f} E)")
    
    # pyTMD requires coordinates as numpy arrays
    lon_arr = np.array([coords['lon']], dtype=float)
    lat_arr = np.array([coords['lat']], dtype=float)
    
    try:
        # 5. Extract amplitudes, phases, and constituent names
        # Uses linear interpolation and extrapolates to coastal nodes if necessary
        amp, ph, cons = pytmd_model.extract_constants(
            lon_arr,
            lat_arr,
            type=pytmd_model.type,
            crop=True,
            buffer=5,
            method="linear",
            extrapolate=True,
            cutoff=np.inf,
            append_node=False,
        )
        
        # 6. Structure output into a DataFrame
        df_out = pd.DataFrame({
            "constituent": cons,
            "amplitude":   amp[0],   # Extract scalar from the spatial array
            "phase_deg":   ph[0],    
        })
        
        # 7. Export to Excel
        output_filename = f"tides_constituents_{st_key}_TPXO10v2.xlsx"
        df_out.to_excel(output_filename, index=False)
        
        print(f"  -> Successfully exported: {output_filename}")
        
    except Exception as e:
        print(f"  -> ERROR extracting constituents for {st_key}: {e}")

print("\n================== EXTRACTION COMPLETED ==================\n")