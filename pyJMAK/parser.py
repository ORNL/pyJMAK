import re
from typing import Dict
import numpy as np
import os
"""
    Parse ABAQUS-formatted parameter/property tables into a nested dict. Unified for all alloy types
    Returns: tables[type][label] -> list of rows (each row is list)
"""

"""
def parse_abaqus_tables(inp_filename: str) -> Dict[str, dict]:
   
    tables = {}
    current_table = None
    current_label = None

    with open(inp_filename, 'r') as f:
        for raw in f:
            line = raw.strip()
            if line.startswith('*PARAMETER TABLE') or line.startswith('*PROPERTY TABLE'):
                type_match = re.search(r'TYPE\s*=\s*"?(.*?)"?(,|$)', line)
                label_match = re.search(r'LABEL\s*=\s*"?(.*?)"?(,|$)', line)
                current_table = type_match.group(1) if type_match else None
                current_label = label_match.group(1) if label_match else current_table
                tables.setdefault(current_table, {})
                tables[current_table][current_label] = []
            elif line.startswith('*') or line.startswith('**') or line == '':
                continue
            elif current_table:
                vals = [v.strip() for v in line.split(',')]
                if current_table == "ABQ_PHASE_TRANS_Transformations":
                    if all(v == '' for v in vals):
                        continue
                    # new row if first column is quoted name
                    if len(vals) > 0 and vals[0].startswith('"') and vals[0].endswith('"'):
                        parsed = []
                        for i, v in enumerate(vals):
                            v_clean = v.strip().strip('"')
                            if v_clean == '':
                                parsed.append(None)
                            elif i in [0, 7, 8, 9, 12, 13, 14, 15, 16, 17]:
                                parsed.append(v_clean)
                            else:
                                try:
                                    v_num = float(v_clean)
                                    if v_num == int(v_num):
                                        v_num = int(v_num)
                                    parsed.append(v_num)
                                except:
                                    parsed.append(v_clean if v_clean else None)
                        tables[current_table][current_label].append(parsed)
                    else:
                        # continuation line fill-in (if present)
                        if len(tables[current_table][current_label]) > 0:
                            last_row = tables[current_table][current_label][-1]
                            for i, v in enumerate(vals):
                                if i < len(last_row) and last_row[i] is None and v.strip():
                                    v_clean = v.strip().strip('"')
                                    if i in [0, 7, 8, 9, 12, 13, 14, 15, 16, 17]:
                                        last_row[i] = v_clean
                                    else:
                                        try:
                                            v_num = float(v_clean)
                                            if v_num == int(v_num): v_num = int(v_num)
                                            last_row[i] = v_num
                                        except:
                                            last_row[i] = v_clean
                else:
                    if len(vals) == 0 or all(v == '' for v in vals):
                        continue
                    parsed = []
                    for v in vals:
                        v_clean = v.strip().strip('"')
                        if v_clean == '':
                            continue
                        try:
                            v_num = float(v_clean)
                            if v_num == int(v_num):
                                v_num = int(v_num)
                            parsed.append(v_num)
                        except:
                            parsed.append(v_clean)
                    if parsed:
                        tables[current_table][current_label].append(parsed)
    return tables
"""



def parse_abaqus_tables(inp_filename):
    #Parse ABAQUS input file tables - unified for all alloy types
    tables = {}
    current_table = None
    current_label = None
    
    with open(inp_filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('*PARAMETER TABLE') or line.startswith('*PROPERTY TABLE'):
                type_match = re.search(r'TYPE\s*=\s*"?(.*?)"?(,|$)', line)
                label_match = re.search(r'LABEL\s*=\s*"?(.*?)"?(,|$)', line)
                current_table = type_match.group(1) if type_match else None
                current_label = label_match.group(1) if label_match else current_table
                tables.setdefault(current_table, {})
                tables[current_table][current_label] = []
            elif line.startswith('*') or line.startswith('**') or line == '':
                continue
            elif current_table:
                # Split by comma and keep empty fields
                vals = [v.strip() for v in line.split(',')]
                # For transformation, allow multi-line entries
                if current_table == "ABQ_PHASE_TRANS_Transformations":
                    if all(v == '' for v in vals):
                        continue
                    # New transformation row: must have a quoted name
                    if len(vals) > 0 and vals[0].startswith('"') and vals[0].endswith('"'):
                        parsed = []
                        for i, v in enumerate(vals):
                            v_clean = v.strip().strip('"')
                            if v_clean == '':
                                parsed.append(None)
                            elif i in [0, 7, 8, 9, 12, 13, 14, 15, 16, 17]:  # String columns
                                parsed.append(v_clean)
                            else:
                                try:
                                    v_num = float(v_clean)
                                    if v_num == int(v_num):
                                        v_num = int(v_num)
                                    parsed.append(v_num)
                                except:
                                    parsed.append(v_clean if v_clean else None)
                        tables[current_table][current_label].append(parsed)
                    else:
                        # Continuation line (fill missing fields)
                        if len(tables[current_table][current_label]) > 0:
                            last_row = tables[current_table][current_label][-1]
                            for i, v in enumerate(vals):
                                if i < len(last_row) and last_row[i] is None and v.strip():
                                    v_clean = v.strip().strip('"')
                                    if i in [0, 7, 8, 9, 12, 13, 14, 15, 16, 17]:
                                        last_row[i] = v_clean
                                    else:
                                        try:
                                            v_num = float(v_clean)
                                            if v_num == int(v_num):
                                                v_num = int(v_num)
                                            last_row[i] = v_num
                                        except:
                                            last_row[i] = v_clean
                else:
                    # Regular table parsing
                    if len(vals) == 0 or all(v == '' for v in vals):
                        continue
                    parsed = []
                    for v in vals:
                        v_clean = v.strip().strip('"')
                        if v_clean == '':
                            continue
                        try:
                            v_num = float(v_clean)
                            if v_num == int(v_num):
                                v_num = int(v_num)
                            parsed.append(v_num)
                        except:
                            parsed.append(v_clean)
                    if len(parsed) > 0:
                        tables[current_table][current_label].append(parsed)
    return tables 



def detect_alloy_type(tables):
    """
    Automatically detect alloy type based on phase names in input file.
    Returns: 'steel', 'ti64', 'generic', or specific alloy name
    """
    if "ABQ_PHASE_TRANS_SolidPhases" not in tables:
        return 'generic'
    
    solidphases = tables["ABQ_PHASE_TRANS_SolidPhases"]["ABQ_PHASE_TRANS_SolidPhases"]
    ph_names = [row[0].lower() for row in solidphases]
    
    # Steel detection
    steel_phases = ['austenite', 'ferrite', 'pearlite', 'bainite', 'martensite']
    if any(phase in ph_names for phase in steel_phases):
        return 'steel'
    
    # Ti-6Al-4V detection  
    ti64_phases = ['beta', 'alphaw', 'alphagb', 'alphaprime']
    if any(phase.replace('_', '').replace('-', '') in ph_names for phase in ti64_phases):
        return 'ti64'
    
    # Check for other common alloy indicators
    if 'gamma' in ph_names or 'delta' in ph_names:
        return 'nickel_superalloy'
        
    return 'generic'

def resample_temperature_data(time_orig, T_C_orig, dt_target):
    """Resample temperature data to uniform time step using interpolation."""
    # Create new uniform time array
    t_start = time_orig[0]
    t_end = time_orig[-1]
    time_new = np.arange(t_start, t_end + dt_target, dt_target)
    
    # Interpolate temperature to new time points
    T_C_new = np.interp(time_new, time_orig, T_C_orig)
    
    print(f"Resampled data: {len(time_orig)} -> {len(time_new)} points")
    print(f"Original time step range: {np.min(np.diff(time_orig)):.6f} - {np.max(np.diff(time_orig)):.6f} s")
    print(f"New uniform time step: {dt_target:.6f} s")
    
    return time_new, T_C_new

def loadTempFile(TEMP_FILE, alloy_type):
    Time_Scale = 5.0  # Time scale multiplier for auto-resampling
     # Check if temperature file is specified
    if TEMP_FILE is None:
        print("ERROR: No temperature file specified!")
        print("Please set the TEMP_FILE variable at the top of this script.")
        print("Example: TEMP_FILE = 'Temp-Steel.txt'")
        print("\nThe program will now exit.")
        exit(1)
    
    # Check if temperature file exists
    if not os.path.exists(TEMP_FILE):
        print(f"ERROR: Temperature file '{TEMP_FILE}' not found!")
        print("Please ensure the file exists in the current directory.")
        print("\nThe program will now exit.")
        exit(1)
    
    # Load temperature data
    try:
        try:
            #First try comma delimiter
            data = np.loadtxt(TEMP_FILE, delimiter=",")        
        except Exception as e:
            #If fails, try tab delimiter       
            data = np.loadtxt(TEMP_FILE, delimiter="\t")
        time_orig = data[:, 0]
        T_C_orig = data[:, 1]
    except Exception as e:
        print(f"ERROR: Failed to load temperature file '{TEMP_FILE}': {e}")
        print("Please ensure the file contains valid comma-separated data.")
        print("\nThe program will now exit.")
        exit(1)
    
    # Handle resampling (V7 style logic for Ti-6Al-4V)
    if alloy_type == 'ti64':
        # resampling logic
        if not np.allclose(np.diff(time_orig), np.mean(np.diff(time_orig)), rtol=0.2):
            mean_dt = np.mean(np.diff(time_orig))
            auto_dt = Time_Scale * mean_dt
            print(f"Time step is highly variable. Resampling to {Time_Scale}x mean time step: {auto_dt:.4f} seconds.")
            time, T_C = resample_temperature_data(time_orig, T_C_orig, auto_dt)
        else:
            time, T_C = time_orig, T_C_orig
    else:
        # Generic logic for other alloys
        time, T_C = time_orig, T_C_orig
    
    print(f"Loaded {len(time)} time points, T range: {T_C.min():.1f} to {T_C.max():.1f} °C")
    time = time.astype(np.float64, copy=False)
    T_C = T_C.astype(np.float64, copy=False)

    return time, T_C, time_orig, T_C_orig


