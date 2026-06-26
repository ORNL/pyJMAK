import matplotlib.pyplot as plt
import numpy as np




# ==================== TTT DIAGRAM PLOTTING FUNCTIONS ====================
def extract_ttt_data_from_tables(tables, alloy_type):
    """Extract TTT diagram data from parsed tables for plotting."""
    ttt_data = {}
    
    if alloy_type == 'steel':
        # Steel transformations
        transformations = {
            'AtoF': 'ABQ_PHASE_TRANS_Trans_AtoF_TTT_Diagram',
            'AtoP': 'ABQ_PHASE_TRANS_Trans_AtoP_TTT_Diagram', 
            'AtoB': 'ABQ_PHASE_TRANS_Trans_AtoB_TTT_Diagram',
            'AtoM': 'ABQ_PHASE_TRANS_Trans_AtoM_TTT_Diagram'
        }
        
        # Extract diffusleft transformations
        for trans_name, label in transformations.items():
            if trans_name == 'AtoM':
                # Martensitic transformation
                ttt_table = tables.get("ABQ_PHASE_TRANS_Martensitic_TTT_Diagram", {}).get(label, [])
                if len(ttt_table) > 0:
                    ttt_data[trans_name] = {
                        'start': np.array([ttt_table[0][0]]),  # Ms
                        'end': np.array([ttt_table[0][1]]),    # Mf
                        'temp': np.array([ttt_table[0][0]])    # Ms temperature
                    }
            else:
                # Diffusional transformations
                ttt_table = tables.get("ABQ_PHASE_TRANS_Diffusional_TTT_Diagram", {}).get(label, [])
                if len(ttt_table) > 0:
                    # Convert list of lists to numpy arrays
                    ttt_array = np.array(ttt_table)
                    ttt_data[trans_name] = {
                        'start': ttt_array[:, 0],  # Start times
                        'end': ttt_array[:, 1],    # End times  
                        'temp': ttt_array[:, 2]    # Temperatures
                    }
    
    elif alloy_type == 'ti64':
        # Ti-6Al-4V transformations
        transformations = {
            'BtoAw': 'ABQ_PHASE_TRANS_Trans_BtoAw_TTT_Diagram',
            'BtoAgb': 'ABQ_PHASE_TRANS_Trans_BtoAgb_TTT_Diagram',
            'BtoM': 'ABQ_PHASE_TRANS_Trans_BtoM_TTT_Diagram',
            'MtoBAw': 'ABQ_PHASE_TRANS_Trans_MtoBAw_TTT_Diagram'
        }
        
        for trans_name, label in transformations.items():
            if trans_name == 'BtoM':
                # Martensitic transformation
                ttt_table = tables.get("ABQ_PHASE_TRANS_Martensitic_TTT_Diagram", {}).get(label, [])
                if len(ttt_table) > 0:
                    ttt_data[trans_name] = {
                        'start': np.array([ttt_table[0][0]]),  # Ms
                        'end': np.array([ttt_table[0][1]]),    # Mf
                        'temp': np.array([ttt_table[0][0]])    # Ms temperature
                    }
            else:
                # Diffusional transformations
                ttt_table = tables.get("ABQ_PHASE_TRANS_Diffusional_TTT_Diagram", {}).get(label, [])
                if len(ttt_table) > 0:
                    # Convert list of lists to numpy arrays
                    ttt_array = np.array(ttt_table)
                    ttt_data[trans_name] = {
                        'start': ttt_array[:, 0],  # Start times
                        'end': ttt_array[:, 1],    # End times
                        'temp': ttt_array[:, 2]    # Temperatures
                    }
    
    return ttt_data



def create_standalone_ttt_diagram(tables, alloy_type):
    """Create a standalone TTT diagram plot."""
    ttt_data = extract_ttt_data_from_tables(tables, alloy_type)
    
    if not ttt_data:
        print("Warning: No TTT diagram data found - cannot create TTT plot")
        return None
    
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Set up the time range for the plot
    time_min = 0.1  # seconds
    time_max = 10000  # seconds
    
    # Define colors for different transformations
    if alloy_type == 'steel':
        colors = {
            'AtoF': 'blue',      # Ferrite
            'AtoP': 'green',     # Pearlite  
            'AtoB': 'red',       # Bainite
            'AtoM': 'black'      # Martensite
        }
        labels = {
            'AtoF': 'Ferrite',
            'AtoP': 'Pearlite',
            'AtoB': 'Bainite', 
            'AtoM': 'Martensite'
        }
    elif alloy_type == 'ti64':
        colors = {
            'BtoAw': 'orange',   # AlphaW
            'BtoAgb': 'green',   # AlphaGB
            'BtoM': 'red',       # Martensite
            'MtoBAw': 'purple'   # Reverse transformation
        }
        labels = {
            'BtoAw': 'AlphaW',
            'BtoAgb': 'AlphaGB',
            'BtoM': 'Martensite',
            'MtoBAw': 'M→Beta'
        }
    else:
        # Generic colors
        colors = {key: f'C{i}' for i, key in enumerate(ttt_data.keys())}
        labels = {key: key for key in ttt_data.keys()}
    
    # Plot each transformation
    for trans_name, data in ttt_data.items():
        if data is None or 'start' not in data or data['start'] is None:
            continue
            
        color = colors.get(trans_name, 'gray')
        label = labels.get(trans_name, trans_name)
        
        if trans_name.endswith('M'):  # Martensitic transformation
            # Plot horizontal lines for Ms and Mf
            if 'start' in data and len(data['start']) > 0:
                ms_temp = data['start'][0]
                ax.axhline(y=ms_temp, color=color, linestyle='-', alpha=0.8, 
                          label=f'{label} Start (Ms)', linewidth=2)
            
            if 'end' in data and len(data['end']) > 0:
                mf_temp = data['end'][0]
                ax.axhline(y=mf_temp, color=color, linestyle='--', alpha=0.8,
                          label=f'{label} End (Mf)', linewidth=2)
        else:
            # Diffusional transformations - plot start and end curves
            if 'start' in data and len(data['start']) > 0:
                ax.plot(data['start'], data['temp'], color=color, linestyle='-', 
                       alpha=0.8, label=f'{label} Start', linewidth=2)
            
            if 'end' in data and len(data['end']) > 0:
                ax.plot(data['end'], data['temp'], color=color, linestyle='--',
                       alpha=0.8, label=f'{label} End', linewidth=2)
    
    # Add labels and title
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Temperature (°C)')
    ax.set_title(f'Time-Temperature-Transformation (TTT) Diagram - {alloy_type.upper()}')
    ax.grid(True, alpha=0.3)
    
    # Set log scale for time axis
    ax.set_xscale('log')
    
    # Add legend
    ax.legend()
    
    # Set axis limits
    ax.set_xlim(time_min, time_max)
    
    # Set temperature limits based on data
    all_temps = []
    for data in ttt_data.values():
        if data and 'temp' in data and data['temp'] is not None:
            all_temps.extend(data['temp'])
    
    if all_temps:
        temp_min = min(all_temps) - 150
        temp_max = max(all_temps) + 50
        ax.set_ylim(temp_min, temp_max)
    
    plt.tight_layout()
    return fig


