#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
In-Situ Tidal Harmonic Analysis using pyFES

Description:
    Performs harmonic analysis on in-situ Sea Surface Height (SSH) time series 
    to extract tidal constituents. The script iterates over multiple stations 
    across the Argentine Continental Shelf, estimates the amplitude and phase 
    for a predefined set of constituents using pyFES, and exports the results 
    to Excel files for validation against global and regional tidal models.
"""

import os
import numpy as np
import pandas as pd
import pyfes

def extract_tidal_components(dates, ssh_array, constituents=None):
    """
    Performs harmonic analysis using pyFES and returns a DataFrame 
    with the estimated tidal constituents.
    
    Parameters:
    -----------
    dates : pandas.Series or numpy.ndarray
        Datetime array in 'datetime64[us]' format.
    ssh_array : numpy.ndarray
        Sea surface height array (or pressure).
    constituents : list of str, optional
        List of specific tidal constituents to extract.
        
    Returns:
    --------
    pandas.DataFrame
        Contains 'constituent' (name), 'amplitude' (same units as input), 
        and 'phase_deg' (phase in degrees referred to UTC).
    """
    # 1. Initialize WaveTable
    wt = pyfes.WaveTable(constituents)

    # 2. Compute leap seconds and nodal/astronomical modulations
    leap = pyfes.get_leap_seconds(dates.values)
    f, vu = wt.compute_nodal_modulations(dates.values, leap)

    # 3. Perform harmonic analysis to obtain the complex vector
    w = wt.harmonic_analysis(ssh_array, f, vu)

    # 4. Extract constituent names, amplitudes, and phases
    names = wt.keys() if hasattr(wt, 'keys') else list(wt)  
    amplitudes = np.abs(w)
    phases_deg = np.angle(w, deg=True)

    # 5. Build and return DataFrame
    df = pd.DataFrame({
        'constituent': names,
        'amplitude': amplitudes,
        'phase_deg': phases_deg
    })

    return df

if __name__ == "__main__":
    
    # Predefined list of tidal constituents to extract
    TARGET_CONSTITUENTS = [
        'Mm', 'Mf', 'Mtm', 'Msqm', '2Q1', 'Sigma1', 'Q1', 'Rho1', 'O1', 'MP1',
        'M11', 'M12', 'M13', 'Chi1', 'Pi1', 'P1', 'S1', 'K1', 'Psi1', 'Phi1',
        'Theta1', 'J1', 'OO1', 'MNS2', 'Eps2', '2N2', 'Mu2', '2MS2', 'N2',
        'Nu2', 'M2', 'MKS2', 'Lambda2', 'L2', '2MN2', 'T2', 'S2', 'R2', 'K2',
        'MSN2', 'Eta2', '2SM2', 'MO3', '2MK3', 'M3', 'MK3', 'N4', 'MN4', 'M4',
        'SN4', 'MS4', 'MK4', 'S4', 'SK4', 'R4', '2MN6', 'M6', 'MSN6', '2MS6',
        '2MK6', '2SM6', 'MSK6', 'S6', 'M8', 'MSf', 'Ssa', 'Sa'
    ]

    # Dictionary defining the input CSVs and target Excel outputs for each station.
    # Note: Ensure the 'input_file' paths match your local directory structure.
    stations = {
        'Carina': {
            'input_file': './serieCarina.csv',
            'output_file': 'tidal_pyFES_components_Carina.xlsx'
        },
        'Vega': {
            'input_file': './serieVega.csv',
            'output_file': 'tidal_pyFES_components_Vega.xlsx'
        },
        'Lander': {
            'input_file': './serieLander.csv',
            'output_file': 'tidal_pyFES_components_Lander.xlsx'
        },
        'PDeseado': {
            'input_file': './seriePuertoDeseado.csv',
            'output_file': 'tidal_pyFES_components_PDeseado.xlsx'
        },
        'sbe26SM': {
            'input_file': './serieSanMatias.csv',
            'output_file': 'tidal_pyFES_components_SanMatias.xlsx'
        }
    }

    print("------------------------------------------------------------")
    print("Initiating pyFES Harmonic Analysis for In-Situ Stations")
    print("------------------------------------------------------------\n")

    for st_key, st_info in stations.items():
        input_csv = st_info['input_file']
        output_xlsx = st_info['output_file']
        
        print(f"Processing station: {st_key}")
        
        if not os.path.exists(input_csv):
            print(f"  -> ERROR: Input file not found ({input_csv}). Skipping...\n")
            continue
            
        try:
            # Load dataset
            df_in = pd.read_csv(input_csv)
            
            # Extract SSH array and format datetime
            # Adjust column names ('SSH_T', 'Fecha') if they differ in other CSVs
            ssh_array = df_in['SSH_T'].values
            dates_datetime = pd.to_datetime(df_in['Fecha']).astype('datetime64[us]')
            
            # Extract components
            components_df = extract_tidal_components(
                dates_datetime,
                ssh_array,
                constituents=TARGET_CONSTITUENTS
            )
            
            # Export to Excel
            with pd.ExcelWriter(output_xlsx) as writer:
                components_df.to_excel(writer, sheet_name='Complete', index=False)
                
            print(f"  -> Successfully extracted {len(components_df)} constituents.")
            print(f"  -> Saved to: {output_xlsx}\n")
            
        except Exception as e:
            print(f"  -> ERROR processing {st_key}: {e}\n")

    print("================== ANALYSIS COMPLETED ==================")