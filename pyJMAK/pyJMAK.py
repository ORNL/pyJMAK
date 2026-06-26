import numpy as np
import os
import matplotlib.pyplot as plt
from . plotter import *
from . kinetics import compute_RLS_state, group_transformations, process_transformations
from . parser import *
from warnings import warn

class jmak:

    def __init__(self, VERBOSE=False, VERBOSE_DEBUG = False):
        self.fileNameTemp = None
        self.fileNamePhase = None
        self.alloyType = None
        self.tables = None
        self.ph_names = None
        self.F = None
        self.RLS_state = None
        self.time= None
        self.T_C = None
        self.time_orig= None
        self.T_orig = None
        self.compareWithExp = False
        self.VERBOSE = VERBOSE
        self.VERBOSE_DEBUG = VERBOSE_DEBUG
        self.gammaFactor = 1.0
        self.exp_data_files = None
        

    def loadInputFiles(self, filename_temp, filename_phase): 
        self.fileNameTemp = filename_temp
        self.fileNamePhase = filename_phase  
        self.tables = parse_abaqus_tables(self.fileNamePhase)
        self.alloy_type = detect_alloy_type(self.tables)
        self.time, self.T_C, self.time_orig, self.T_orig = loadTempFile(self.fileNameTemp, self.alloy_type)
        # if you want to reset time and temperature to original values, no interpolation
        #self.T_C = self.T_orig
        #self.time = self.time_orig

    def setVerbose(self, VERBOSE):
        self.VERBOSE = VERBOSE

    def setTimeTemp(self, time, temperature):
        self.time = time
        self.time_orig = time
        self.T_C = temperature
        self.T_orig = temperature

    def setExperimentalComparison(self, compareExp, exp_data_files=None):
        self.compareWithExp = compareExp
        if compareExp:
            self.exp_data_files = exp_data_files        
    
    def setAlloyType(self, alloy_type):
        self.alloy_type = alloy_type
    
    def setMaterialTables(self, tables):
        self.tables = tables

    def run_phase_trans_sim(self):
        
        time = self.time
        VERBOSE = self.VERBOSE
        T_C = self.T_C
        tables = self.tables        
        GAMMA_FACTOR = self.gammaFactor
        VERBOSE_DEBUG = self.VERBOSE_DEBUG  # Set to True for detailed debug output

        # Extract phases
        solidphases = tables["ABQ_PHASE_TRANS_SolidPhases"]["ABQ_PHASE_TRANS_SolidPhases"]
        ph_names = [row[0] for row in solidphases]
        ph_init = [row[1] for row in solidphases]
        nph = len(ph_names)
        
        
               

        #Assumption: Initially, raw material is present
        

        # Sanitize time array: remove duplicate or effectively-zero increments which
        # cause divide-by-zero inside `np.gradient` (produces RuntimeWarning).
        time = np.asarray(time, dtype=np.float64)
        T_C = np.asarray(T_C, dtype=np.float64)
        if len(time) > 1:
            dt = np.diff(time)
            # consider increments effectively zero if below this threshold
            tiny = 1e-15 * max(1.0, np.max(np.abs(time)))
            dup_mask = np.concatenate(([True], np.abs(dt) > tiny))
            if not np.all(dup_mask):
                #warn_count = len(time) - np.count_nonzero(dup_mask)
                #warn(f"Detected {warn_count} duplicate/zero-time steps; removing duplicates to avoid divide-by-zero in gradient.")
                time = time[dup_mask]                
                T_C = T_C[dup_mask]
        dTdt = np.gradient(T_C, time)
        self.T_C = T_C
        self.time = time
        F = np.zeros((len(time), nph), dtype=np.float64) #Intial phases when material is raw should be 0.
        F[0,:] = ph_init #Initial phase fractions from input file (upon solidification)

        # --- Get solidus and liquidus temperatures from input file ---
        melting_temp_table = tables.get("ABQ_PHASE_TRANS_MeltingTemperature", {}).get("ABQ_PHASE_TRANS_MeltingTemperature", [[980.0, 1000.0]])
        if len(melting_temp_table) > 0 and len(melting_temp_table[0]) >= 2:
            Tsol = melting_temp_table[0][0]
            Tliq = melting_temp_table[0][1]
        else:
            # Default values for Ti-6Al-4V
            Tsol, Tliq = 1608.0, 1640.0
        if VERBOSE:
            print(f"Using Solidus: {Tsol}°C, Liquidus: {Tliq}°C")

        # --- Initialize RLS state tracking ---    
        RLS_state = np.zeros(len(time))
        # Determine which phase solidifies first (if any)
        solidifying_idx = None
        for i, v in enumerate(ph_init):            
            if v == 1:
                solidifying_idx = i 
                break
        if solidifying_idx is not None: 
            has_melted = False #Assume material has not melted initially
        else: 
            has_melted = True #If no single solidifying phase, assume material has melted initially

        # Initialize RLS state
        RLS_state[0], has_melted = compute_RLS_state(T_C[0], Tsol, Tliq, has_melted)
        if VERBOSE: print(f"Initial RLS state: {RLS_state[0]:.2f}")
        
        #Initialize phase fractions based on initial RLS state
        if RLS_state[0] <=0.0: #Raw material OR LIQUID
            #If Raw material, set ALL phase fraction to zero
            F[0, :] = 0.0

        elif RLS_state[0] > 0.0 and has_melted:  # Liquid or solid after melting
            # Beta should be prominent if no significant alpha phases
            #Assumption: Other phases are alpha in Ti64
            
            #Assumption: One and only one phase solidifies initially and initial RSL_state = raw material
            if solidifying_idx: 
                total_other = sum(F[0, j] for j in range(nph) if j != solidifying_idx )
                if total_other < 0.1:
                    F[0, solidifying_idx] = max(F[0, solidifying_idx], 0.9)
            else:
                F[0, :] = ph_init  # Default to input fractions if no single solidifying phase found

        # --- Get all transformation rules ---
        ttype = "ABQ_PHASE_TRANS_Transformations"
        trans_rules = tables[ttype][ttype]
        #print(trans_rules)
        competing_transformations = group_transformations(trans_rules, ph_names)
        #print(f"Competing transformations grouped by parent phase: {competing_transformations}")

        # Initialize kinetic state variables
        km_active = False
        km_fp0 = None

        for i in range(1, len(time)):
            
            F[i,:] = F[i-1,:]
            T = T_C[i-1]
            dt = time[i] - time[i-1]
            rate = dTdt[i-1]

            # --- Update RLS state ---
            RLS_state[i], has_melted = compute_RLS_state(T, Tsol, Tliq, has_melted)
            is_liquid = -0.1 < RLS_state[i] <0.8
            was_liquid = -0.1 < RLS_state[i-1] <0.8 if i > 0 else False
            is_solidifying = (was_liquid and not is_liquid and has_melted)
            is_raw = RLS_state[i] < -0.1

            if VERBOSE_DEBUG:
                print(f"\nStep {i}: T={T:.2f}°C, dT/dt={rate:.2f}°C/s, time={time[i]:.4f}s, RLS={RLS_state[i]:.2f}, is_liquid={is_liquid}, is_solidifying={is_solidifying}, is_raw={is_raw}, has_melted={has_melted}")
                print(f"Initial Phase Fractions: " + ", ".join(f"{ph_names[j]}={F[i,j]:.4f}" for j in range(nph)))  

                    
            # LIQUID STATE LOGIC: All solid phases = 0 when above liquidus (material is liquid)
            if is_liquid or is_raw:
                # Set all solid phases to 0 (material is liquid)
                F[i,:] = 0.0
                if VERBOSE_DEBUG: print(f"Liquid or Raw: Phase fractions at the end: " + ", ".join(f"{ph_names[j]}={F[i,j]:.4f}" for j in range(nph)))
                continue
            
            # SOLIDIFICATION LOGIC: Material solidifies as phase 1 first, then transforms during cooling
            elif is_solidifying:                
                try:
                    F[i,:] = ph_init  
                    if VERBOSE_DEBUG:print("Solidifying: Phase fractions at the end: " + ", ".join(f"{ph_names[j]}={F[i,j]:.4f}" for j in range(nph)))
                    continue #Reset to initial phase fractions upon solidification          
                       
                except ValueError:
                    print("Warning:Please specify which phase will form upon solidification - cannot reset during solidification")
                    
            if VERBOSE_DEBUG:
                print(f"Phase fractions: " + ", ".join(f"{ph_names[j]}={F[i,j]:.4f}" for j in range(nph)))  

            km_active, km_fp0 = process_transformations(competing_transformations,  F, i, T, rate, dt, tables, km_active, km_fp0, GAMMA_FACTOR, tolerance=1e-15,VERBOSE = VERBOSE_DEBUG)

            """""
            # Start transformation processing for diffusional transformations
            process_diffusional_transformations(competing_transformations, F, i, T, rate, dt, tables, tolerance=1e-15, VERBOSE=VERBOSE_DEBUG)

            # Handle martensitic transformations separately
            for t_rule in trans_rules:
                if len(t_rule) < 12: continue
                trans_type = t_rule[7].strip().upper() if len(t_rule) > 7 and t_rule[7] else "DIFFUSIONAL"
                if trans_type == "MARTENSITIC":
                    # Process martensitic transformation using existing logic
                    km_active, km_fp0 = process_martensitic_transformation([t_rule], ph_names, F, i, T, rate, dt, tables, km_active, km_fp0, GAMMA_FACTOR, VERBOSE=VERBOSE_DEBUG)
            """""
            if VERBOSE_DEBUG:
                print(f"Phase fractions after JMAK/KM processing: " + ", ".join(f"{ph_names[j]}={F[i,j]:.4f}" for j in range(nph)))  

            F[i, :] = F[i, :].astype(np.float64, copy=False)
            # --- Phase conservation (all phases) ---
            total =  np.float64(np.sum(F[i,:]))
            if VERBOSE_DEBUG: print(f" Total phase sum before normalization: {total:.6f}")
            if abs(total - 1.0) > 1e-15:
                F[i,:] = F[i,:] / max(total, 1e-15)
                if VERBOSE_DEBUG:
                    print(f"Phase fractions after final normalization: " + ", ".join(f"{ph_names[j]}={F[i,j]:.4f}" for j in range(nph)))  
            if VERBOSE_DEBUG:print(f"End  Phase fractions: " + ", ".join(f"{ph_names[j]}={F[i,j]:.4f}" for j in range(nph)))  

        if VERBOSE:
            # Debug: Print transformation rules as parsed
            print("\n=== PARSED TRANSFORMATION RULES ===")
            trans_rules = tables["ABQ_PHASE_TRANS_Transformations"]["ABQ_PHASE_TRANS_Transformations"]
            for i, rule in enumerate(trans_rules):
                if len(rule) > 8:
                    name = rule[0] if rule[0] else "Unknown"
                    trate_cond = rule[1] if len(rule) > 1 else "N/A"
                    trate_lo = rule[2] if len(rule) > 2 else "N/A"
                    trate_hi = rule[3] if len(rule) > 3 else "N/A"
                    t_cond = rule[4] if len(rule) > 4 else "N/A"
                    t_lo = rule[5] if len(rule) > 5 else "N/A"
                    t_hi = rule[6] if len(rule) > 6 else "N/A"
                    trans_type = rule[7] if len(rule) > 7 else "N/A"
                    print(f"Rule {i} ({name}): TRate={trate_cond}[{trate_lo},{trate_hi}], Temp={t_cond}[{t_lo},{t_hi}], Type={trans_type}")
                else:
                    print(f"Rule {i}: {rule}")
            print("="*50)
        
        self.ph_names = ph_names
        self.F = F
        self.RLS_state = RLS_state

    def plot_results(self, filename=None, saveFig=False):        
        
        ph_names = self.ph_names
        F = self.F
        VERBOSE = self.VERBOSE
        if filename is None:
            filename = os.path.splitext(self.fileNameTemp)[0] + '_phaseFractions.png'

        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 9))

        # Temperature profile
        ax1.plot(self.time_orig, self.T_orig, 'r-', linewidth=2, label='Temperature')
        ax1.set_xlabel("Time (s)")
        ax1.set_ylabel("Temperature (°C)")
        ax1.set_title("Temperature Profile")
        ax1.grid(True, alpha=0.3)
        ax1.legend()

        # Phase evolution plot (model)
        model_colors = {}
        for i, ph in enumerate(ph_names):
            line, = ax2.plot(self.time, F[:, i], label=f"{ph} (Model)", linewidth=2)
            model_colors[ph] = line.get_color()  # store the color assigned by matplotlib

        ax2.set_xlabel("Time (s)")
        ax2.set_ylabel("Phase fraction")
        ax2.set_title("Phase Evolution")

        # RLS state plot
        ax3.plot(self.time, self.RLS_state, 'g-', linewidth=2, label='RLS State')
        ax3.axhline(y=-1, color='brown', linestyle='--', alpha=0.7, label='Raw (-1)')
        ax3.axhline(y=0, color='orange', linestyle='--', alpha=0.7, label='Liquid (0)')
        ax3.axhline(y=1, color='blue', linestyle='--', alpha=0.7, label='Solid (1)')
        ax3.set_xlabel("Time (s)")
        ax3.set_ylabel("RLS State")
        ax3.set_title("Raw-Liquid-Solid State")
        ax3.set_ylim(-1.1, 1.1)
        ax3.legend()
        ax3.grid(True, alpha=0.3)

        plt.tight_layout()

        # Experimental data overlay
        if self.compareWithExp:
            if not hasattr(self, "exp_data_files"):
                raise AttributeError("No experimental data files mapping (self.exp_data_files) found.")

            markers = ["o", "^", "s", "d", "x", "P", "v"]

            for idx, (phase_name, filepath) in enumerate(self.exp_data_files.items()):
                if phase_name not in self.ph_names:
                    print(f"Warning: Phase '{phase_name}' not found in model output; skipping.")
                    continue
                try:
                    exp_data = np.loadtxt(filepath)
                    exp_time, exp_frac = exp_data[:, 0], exp_data[:, 1]
                    ax2.plot(
                        exp_time, exp_frac,
                        marker=markers[idx % len(markers)],
                        markersize=0,
                        linestyle='--',
                        linewidth = 2.5,
                        label=f"{phase_name} (Exp)",
                        color=model_colors.get(phase_name, 'black')  # match model color, fallback to black
                    )

                    if self.VERBOSE:
                        print(f"Loaded experimental data for {phase_name} from '{filepath}'")

                        ph_idx = self.ph_names.index(phase_name) 
                        if ph_idx is not None:
                            final_model = F[-1, ph_idx]
                            final_exp = exp_frac[-1]
                            error_aw = abs(final_model - final_exp)
                            print(f"{phase_name} Final: Model={final_model:.3f}, Exp={final_exp:.3f}, Error={error_aw:.3f}")


                except FileNotFoundError:
                    print(f"Warning: File '{filepath}' not found for phase '{phase_name}'")
                except Exception as e:
                    print(f"Warning: Could not load data for '{phase_name}': {e}")
            ax2.set_title("Phase Transformation Evolution - Model vs Given Data")

        ax2.legend(loc='center left')
        ax2.grid(True, alpha=0.3)
        plt.show()
        if saveFig: 
            fig.savefig(filename, format='png', dpi=200)
            plt.close(fig)
            print(f"Plot saved as '{filename}'")
        if VERBOSE:  
            print("\n=== PHASE EVOLUTION SUMMARY ===")
            for i, ph in enumerate(ph_names):
                change = F[-1, i] - F[0, i]
                if abs(change) > 0.001:
                    print(f"{ph}: {F[0, i]:.3f} → {F[-1, i]:.3f} (Δ{change:+.3f})")

    def plot_TTT_diagram(self, saveFig=False):        
        
        ttt_fig = create_standalone_ttt_diagram(self.tables, self.alloy_type)
        if ttt_fig:
            if saveFig:
                ttt_filename = os.path.splitext(self.fileNameTemp)[0] + '_tttDiagram.png'
                ttt_fig.savefig(ttt_filename, dpi=300, bbox_inches='tight')
                print(f"TTT diagram saved as '{ttt_filename}'")

            plt.figure(ttt_fig.number)  # Make TTT figure active for display    
            plt.show()

