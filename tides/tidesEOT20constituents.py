"""
EOT20 Tidal Constituent Extraction

Description:
    Extracts tidal harmonic constituents (amplitude and phase) from the 
    EOT20 global ocean tide model for a specific set of in-situ monitoring 
    stations over the Argentine Continental Shelf. The EOT20 model 
    automatically combines the ocean tide and load tide components.
"""

import numpy as np
import pandas as pd
import pyTMD

# 1. Load pyTMD database and EOT20 model
print("Loading EOT20 tidal model database...")
db = pyTMD.io.load_database()
m  = pyTMD.io.model(directory="tide_models/")
eot = m.elevation("EOT20")

# 2. Define target stations
# Note: The keys used here match the expected output filenames.
stations = {
    'Carina':   {'lon': -67.21960, 'lat': -52.75718},
    'Vega':     {'lon': -67.74600, 'lat': -53.29900},
    'Lander':   {'lon': -66.67000, 'lat': -52.14000},
    'PDeseado': {'lon': -65.91420, 'lat': -47.75360},
    'sbe26SM':  {'lon': -63.77680, 'lat': -41.18190}
}

print("\n------------------------------------------------------------")
print("Initiating harmonic constituent extraction")
print("------------------------------------------------------------\n")

# 3. Iterate over stations to extract and export constituents
for st_name, coords in stations.items():
    # pyTMD requires longitudes in the [0, 360) degree convention.
    # Convert from [-180, 180) to [0, 360)
    lon_360 = 360.0 + coords['lon'] if coords['lon'] < 0 else coords['lon']
    
    lon_arr = np.array([lon_360], dtype=float)
    lat_arr = np.array([coords['lat']], dtype=float)
    
    print(f"Processing station: {st_name} ({coords['lat']:.5f} N, {coords['lon']:.5f} E)")
    
    try:
        # Extract amplitude (meters), phase (degrees), and constituent names
        amp, ph, cons = eot.extract_constants(
            lon_arr, lat_arr,
            type=eot.type,
            crop=True,
            buffer=5,
            method="bilinear",
            extrapolate=True,
            cutoff=np.inf
        )
        
        # Construct DataFrame
        df_out = pd.DataFrame({
            "constituent": cons,
            "amplitude":   amp[0],   # Amplitude in meters
            "phase_deg":   ph[0],    # Phase in degrees
        })
        
        # Export to Excel
        output_filename = f"tides_constituents_{st_name}_EOT20.xlsx"
        df_out.to_excel(output_filename, index=False)
        print(f"  -> Exported: {output_filename}\n")
        
    except Exception as e:
        print(f"  -> ERROR extracting constituents for {st_name}: {e}\n")

print("================== EXTRACTION COMPLETED ==================")