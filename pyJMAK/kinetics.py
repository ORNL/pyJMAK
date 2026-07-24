import numpy as np

def interp_table(T, arr, T_col=-1, y_col=0):
    """Interpolate values from a table."""
    arr = np.array(arr)
    if arr.shape[1] < abs(T_col): 
        return arr[0, y_col]
    Ts = arr[:, T_col]
    ys = arr[:, y_col]
    return np.interp(T, Ts, ys)

def group_transformations(trans_rules, ph_names):
    """Group transformations by parent phase for competing transformation logic."""
    competing_transformations = {}
    for t_rule in trans_rules:
        if len(t_rule) < 12: continue
        nparent = int(t_rule[10]) if len(t_rule) > 10 and t_rule[10] is not None else 1
        nchild = int(t_rule[11]) if len(t_rule) > 11 and t_rule[11] is not None else 1
        parent_names = []
        child_names = []
        for j in range(nparent):
            idx = 12 + j
            if len(t_rule) > idx and t_rule[idx] is not None and t_rule[idx].strip():
                parent_names.append(t_rule[idx])
        for j in range(nchild):
            idx = 15 + j
            if len(t_rule) > idx and t_rule[idx] is not None and t_rule[idx].strip():
                child_names.append(t_rule[idx])
        parent_idx = [ph_names.index(name) for name in parent_names if name in ph_names]
        child_idx = [ph_names.index(name) for name in child_names if name in ph_names]
        if len(parent_idx)==0 or len(child_idx)==0:
            continue

        trans_name = t_rule[0] if len(t_rule) > 0 and t_rule[0] else "Unknown"
        
        # Group competing transformations by parent phase
        if len(parent_names) > 0:
            parent_key = parent_names[0]  # Use first parent as key
            if parent_key not in competing_transformations:
                competing_transformations[parent_key] = []
            competing_transformations[parent_key].append({
                'trans_name': trans_name,
                'parent_idx': parent_idx,
                'child_idx': child_idx,
                'parent_names': parent_names,
                'child_names': child_names,
                'Trate_cond': t_rule[1] if len(t_rule) > 1 and t_rule[1] else "NONE",
                'Trate_lo': t_rule[2] if len(t_rule) > 2 and t_rule[2] is not None else -1e6,
                'Trate_hi': t_rule[3] if len(t_rule) > 3 and t_rule[3] is not None else 1e6,
                'T_cond': t_rule[4] if len(t_rule) > 4 and t_rule[4] else "NONE",
                'T_lo': t_rule[5] if len(t_rule) > 5 and t_rule[5] is not None else -1e6,
                'T_hi': t_rule[6] if len(t_rule) > 6 and t_rule[6] is not None else 1e6,
                'trans_type': t_rule[7].strip().upper() if len(t_rule) > 7 and t_rule[7] else "DIFFUSIONAL",
                'input_option': t_rule[8].strip().upper() if len(t_rule) > 8 and t_rule[8] else "JMA",
                'reversible': str(t_rule[9]).strip().upper() == "YES" if len(t_rule) > 9 and t_rule[9] else False
            })
    return competing_transformations
            
# ==================== KINETIC MODELS ====================        
def additivity_jmak_total(f_prev,f_total, k, n, dt, f_eq, forward=True, VERBOSE=False):
    """
    JMAK transformation operating directly on the child (transformed) phase fraction.
    Handles changing equilibrium fraction (f_eq) at each time step.
    Uses the additivity rule for non-isothermal transformations.
    """
    if f_eq <= 0.0:
        return 0.0
    
    if forward:
        if f_prev > (f_total*f_eq):
            if VERBOSE: print("additivity_jmak: f_prev > f_eq in additivity_jmak, returning f_eq")
            return f_eq
        try:
            if f_prev <= (f_total*f_eq):
                arg = 1.0 - (f_prev / (f_total*f_eq)  )              
                arg = min(max(arg, 1e-15), 1.0)
                xi = np.float64(((-1.0 / k) * np.log(arg)) ** (1.0 / n))
            else:
                xi = 0.0
            if VERBOSE: print(f"xi: {xi}, dt: {dt}, k: {k}, n: {n}, f_prev: {f_prev}, f_eq: {f_eq}, f_total: {f_total}")
            exp_arg = -k * ((xi + dt) ** n)
            exp_val = np.exp(exp_arg) #if exp_arg > -700 else 0.0
            f_new = f_total*f_eq * (1.0 - exp_val)
            if VERBOSE: print(f"f_new: {f_new}, f_prev: {f_prev}, f_eq: {f_eq}, f_total: {f_total}")
            if f_new <= f_prev:
                return f_prev
            return min(max(f_new, f_prev), f_eq)
        except:
            return f_prev
    

def km_fraction(T, Ms, gamma):
    """Koistinen-Marburger fraction calculation."""
    return 1.0 - np.exp(-gamma * (Ms - T))

def backward_euler_martensitic(f_prev, f_eq, dt, t_char):
    """
    Backward Euler method for martensitic transformations.
    Used for M-A transformations which are temperature-driven, not time-dependent.
    
    Solves: df/dt = -(f - f_eq) / t_char
    Solution: f_new = (f_prev + dt * f_eq / t_char) / (1 + dt / t_char)
    """
    if t_char <= 0 or dt <= 0:
        return f_prev
    rate_param = dt / t_char
    f_new = (f_prev + rate_param * f_eq) / (1.0 + rate_param)
    return max(min(f_new, 1.0), 0.0)

# ==================== RLS STATE FUNCTIONS ====================
def compute_RLS_state(T_C, Tsol, Tliq, has_melted):
    if T_C > Tliq:
        has_melted = True
        return 0.0, has_melted #Melted previously, liquid
    elif T_C < Tsol:
        if has_melted:
            return 1.0, has_melted #Melted previously, solid
        else:
            return -1.0, has_melted #Not melted previously, raw
    else:
        if has_melted:
            solid_fraction = (Tliq - T_C) / (Tliq - Tsol)
            return solid_fraction, has_melted
        else:
            liquid_fraction = (T_C - Tsol) / (Tliq - Tsol)
            return -1.0 + liquid_fraction, has_melted



def process_competing_diffusional_transformations(active_transformations, F, i, dt, tolerance, transformation_rates, VERBOSE):
    parent_phase_idx = active_transformations[0]['parent_idx'] 
    if VERBOSE: print(f"parent_phase_idx: {parent_phase_idx}, active transformations: {len(active_transformations)}")           
    f_p_prev = sum(F[i-1, p_idx] for p_idx in active_transformations[0]['parent_idx'])         
    F[i, parent_phase_idx] = np.float64(F[i-1, parent_phase_idx])
    f_p = np.float64(F[i, parent_phase_idx])
    if VERBOSE: print(f"Current parent phase fraction: {f_p}, previous: {f_p_prev}") 
    for trans in active_transformations:                
        if VERBOSE: print(f"Active transformations: {trans['trans_name']}")
        if VERBOSE: print(f"Child phase name : {trans['child_names']}")      

        f_eq_parent = trans['f_eq']
        f_eq_child = 1.0 - f_eq_parent

        f_c_prev = sum(F[i-1, c_idx] for c_idx in trans['child_idx'])
        f_total = f_p_prev + f_c_prev
                    
        # Calculate potential transformation rate
        if f_p_prev > ((f_total*f_eq_parent) + tolerance):                     
            if VERBOSE: print(f"Parent > f_eq: trans_name: {trans['trans_name']},f_total: {f_total}, f_eq_parent: {f_eq_parent}, f_tot*f_eq: {f_total*f_eq_parent}, f_prev_parent: {f_p_prev}")

            f_new = np.float64(additivity_jmak_total(f_c_prev, f_total, trans['k'], trans['n'], dt, f_eq_child, forward=True, VERBOSE=VERBOSE))
            delta_child = np.float64(f_new - f_c_prev)   
            delta_parent = delta_child        #matches with Abaqus                                         
            #delta_parent = np.float64(min(delta_child, f_p - f_eq_parent)) #does not match with Abaqus, but more logical
            transformation_rates.append(delta_parent)
            if VERBOSE: print(f"f_c_prev: {f_c_prev},  f_new: {f_new}, delta_child: {delta_child}, delta_parent: {delta_parent}, f_eq_parent: {f_eq_parent}, current_parent: {f_p}")
        
        elif f_p_prev < ((f_total*f_eq_parent) + tolerance) and trans['reversible']:
            if VERBOSE: print(f"Parent < f_eq: trans_name: {trans['trans_name']},f_eq_child: {f_eq_child}, f_c_prev: {f_c_prev}, f_p: {f_p}, f_prev_parent: {f_p_prev}")

            f_new_parent = additivity_jmak_total(f_p_prev, f_total, trans['k'], trans['n'], dt, f_eq_parent, forward=True, VERBOSE=VERBOSE)                        
            #delta_parent = np.float64(min(f_p_prev - f_new_parent, f_p_prev -f_eq_parent))
            delta_parent = min(f_p_prev-f_new_parent, 0.0)
            delta_child = delta_parent
            transformation_rates.append(delta_parent)

            #F[i,p] = min(f_new_parent, f_eq)
            #F[i,c] = max(f_c_prev + (f_p - F[i,p]), 0.0)
        else:
            transformation_rates.append(0.0)

    if VERBOSE: print(f"Transformation rates before normalization: {transformation_rates}, current_parent: {f_p}")
    # Normalize transformation rates to ensure mass conservation
    
    
    # Apply transformations proportionally
    for j, trans in enumerate(active_transformations):
        if abs(transformation_rates[j]) > 0.:                     
            
            f_c_sum = sum(F[i-1, c_idx] for c_idx in trans['child_idx'])
            f_p_sum = sum(F[i-1, p_idx] for p_idx in trans['parent_idx'])
            delta_parent_sum = transformation_rates[j] 
            if VERBOSE: print(f"Within trans:{trans['trans_name']}, delta_parent_sum: {delta_parent_sum}, f_c_sum: {f_c_sum}, f_p_sum: {f_p_sum}")
            for c_idx in trans['child_idx']:
                delta_child = (F[i-1, c_idx]/f_c_sum)*delta_parent_sum if (f_c_sum > 0 and len(trans['child_idx'])>1) else delta_parent_sum
                F[i, c_idx] = np.float64(F[i, c_idx] + delta_child)
                if VERBOSE: print(f"child phase {c_idx}: F[i, c_idx]: {F[i, c_idx]}, delta_child: {delta_child}")

            
            for p_idx in trans['parent_idx']:
                delta_parent = (F[i-1, p_idx]/f_p_sum)*delta_parent_sum if (f_p_sum > 0 and len(trans['child_idx'])>1) else delta_parent_sum
                F[i, p_idx] = np.float64(F[i, p_idx] - delta_parent)
                if VERBOSE: print(f"parent phase {p_idx}: F[i, p_idx]: {F[i, p_idx]}, delta_parent: {delta_parent}")
                                    
            if VERBOSE: print(f"After applying transformation: F[i, c_idx]: {F[i, c_idx]}, F[i, p_idx]: {F[i, p_idx]}")
    if VERBOSE: print(f"Sum of phases after applying transformations: {np.sum(F[i,:])}")



def process_transformations(competing_transformations,  F, i, T, rate, dt, tables, km_active, km_fp0, GAMMA_FACTOR, tolerance=1e-15,VERBOSE = False):
    """Process all transformations for a given time step."""

    for parent_key, trans_group in competing_transformations.items():
        if VERBOSE: print(f"Processing competing transformations for parent phase: {parent_key} ")
        # Evaluate each transformation in the group
        active_diff_transformations = []

        for trans_info in trans_group:
            trans_name = trans_info['trans_name']
            parent_idx = trans_info['parent_idx']
            child_idx = trans_info['child_idx']
            parent_names = trans_info['parent_names']
            child_names = trans_info['child_names']
            Trate_cond = trans_info['Trate_cond']
            Trate_lo = trans_info['Trate_lo']
            Trate_hi = trans_info['Trate_hi']
            T_cond = trans_info['T_cond']
            T_lo = trans_info['T_lo']
            T_hi = trans_info['T_hi']
            trans_type = trans_info['trans_type']
            input_option = trans_info['input_option']
            reversible = trans_info['reversible']

            # Check if this transformation is active at current T and rate
            triggered = True
            no_trigger_reasons = []
                
            
                # Temperature Rate Conditions (only apply if not NONE)
            if Trate_cond != "NONE":
                if Trate_cond == "INTERVAL" and not (Trate_lo < rate < Trate_hi): 
                    triggered = False
                    no_trigger_reasons.append(f"dT/dt={rate:.1f} not in [{Trate_lo},{Trate_hi}]")
                elif Trate_cond == "MIN":
                    # MIN condition: rate must be > threshold
                    # E.g., MIN[-410] means cooling rate > -410°C/s (slower cooling than threshold)
                    min_bound = Trate_lo if Trate_lo is not None and Trate_lo != -1e6 else Trate_hi
                    if not (rate > min_bound):
                        triggered = False
                        no_trigger_reasons.append(f"dT/dt={rate:.1f} < min threshold {min_bound} (too fast cooling)")
                elif Trate_cond == "MAX":
                    # MAX condition: rate must be < threshold  
                    # E.g., MAX[-410] means cooling rate < -410°C/s (faster cooling than threshold)
                    max_bound = Trate_hi if Trate_hi is not None and Trate_hi != 1e6 else Trate_lo
                    if not (rate < max_bound):
                        triggered = False
                        no_trigger_reasons.append(f"dT/dt={rate:.1f} > max threshold {max_bound} (too slow cooling)")
            
            # Temperature Conditions (only apply if not NONE)
            if T_cond != "NONE":
                if T_cond == "INTERVAL" and not (T_lo < T < T_hi): 
                    triggered = False
                    no_trigger_reasons.append(f"T={T:.1f} not in [{T_lo},{T_hi}]")
                elif T_cond == "MIN" and not (T >= T_lo): 
                    triggered = False
                    no_trigger_reasons.append(f"T={T:.1f} < min {T_lo}")
                elif T_cond == "MAX":
                    # For MAX condition, check against whichever bound is specified
                    max_bound = T_hi if T_hi is not None and T_hi != 1e6 else T_lo
                    if not (T <= max_bound):
                        triggered = False
                        no_trigger_reasons.append(f"T={T:.1f} > max {max_bound}")            

            if VERBOSE: print(f"Reason for not triggering {trans_name} transformation = {no_trigger_reasons}" )

            if T_cond == "NONE" and Trate_cond == "NONE":
                triggered = False  # No conditions specified, do not trigger
                ValueError(f"Error: Transformation rule {trans_name} must have at least one condition (T or dT/dt) specified.")

            if not triggered: 
                km_active = False
                continue
            
                # --- KINETICS: JMAK/TTT (diffusional) ---
            if trans_type == "MARTENSITIC" and input_option in ("KM", "TTT"):

                if input_option == "TTT":
                    # Get Ms, Mf from TTT diagram
                    ttt_diag_label = f"ABQ_PHASE_TRANS_Trans_{trans_name}_TTT_Diagram"
                    ttt_diag = tables.get("ABQ_PHASE_TRANS_Martensitic_TTT_Diagram", {}).get(ttt_diag_label, [])
                    if len(ttt_diag) > 0 and len(ttt_diag[0]) >= 2:
                        Ms, Mf = ttt_diag[0][0], ttt_diag[0][1]
                    else:
                        Ms, Mf = 400.0, 200.0  # Fallback values
                    # Get fMs, fMf from TTT constants
                    ttt_const_label = f"ABQ_PHASE_TRANS_Trans_{trans_name}_TTT_Constants"
                    ttt_const = tables.get("ABQ_PHASE_TRANS_Martensitic_TTT_Constants", {}).get(ttt_const_label, [])
                    if len(ttt_const) > 0 and len(ttt_const[0]) >= 2:
                        fMs, fMf = ttt_const[0][0], ttt_const[0][1]
                    else:
                        fMs, fMf = 0.0, 0.99  # Fallback values
                    # Calculate gamma
                    try:
                        gamma = -np.log((1-fMs)/(1-fMf)) / (Mf - Ms)
                        # Apply gamma scaling for TTT-based transformation
                        gamma = gamma * GAMMA_FACTOR
                    except Exception:
                        gamma = 0.005 * GAMMA_FACTOR  # Fallback with scaling
                    # Use Ms from TTT diagram
                else:  # KM option
                    km_label = f"ABQ_PHASE_TRANS_Trans_{trans_name}_KM_Coefficients"
                    gamma = tables["ABQ_PHASE_TRANS_Martensitic_KM_Coefficients"].get(km_label, [[0.005]])[0][0]
                    # Apply gamma scaling
                    gamma = gamma * GAMMA_FACTOR
                    # Ms from input file (T_hi)
                    Ms = T_hi
                
                # f_pr lookup from martensitic table if available, using km_fp0 as independent variable
                pr_label = f"ABQ_PHASE_TRANS_Trans_{trans_name}_ParentRetainedFrac"
                pr_table = tables.get("ABQ_PHASE_TRANS_Martensitic_ParentRetainedFrac", {}).get(pr_label, [])
                if not pr_table or len(pr_table) == 0:                
                    raise RuntimeError("ERROR: You must include a property table of type 'ABQ_PHASE_TRANS_Martensitic_ParentRetainedFrac' for martensitic transformations.")
                    
                    

                if not km_active:
                    km_fp0 = F[i-1,parent_idx[0]] + F[i-1,child_idx[0]]  # total parent before transformation
                    km_active = True
                
                # Interpolate f_pr using km_fp0 as independent variable (second column), retrieving retained proportion (first column)
                if len(pr_table) > 0:
                    f_pr = interp_table(km_fp0, pr_table, T_col=1, y_col=0)
                else:
                    f_pr = 0.0  # Default if no table
                
                if km_fp0 > f_pr + 1e-8:
                    # Direct KM fraction calculation
                    mart_frac = km_fraction(T, Ms, gamma)
                    
                    # Calculate maximum possible martensite fraction (limited by available parent)
                    max_mart_frac = km_fp0 - f_pr
                    
                    # Apply the transformation
                    fc_prev = F[i-1,child_idx[0]]
                    fc_new = min(mart_frac * max_mart_frac, max_mart_frac)  # Limit by available parent
                    fc_new = max(fc_new, fc_prev)  # Can only increase during cooling
                    
                    delta_fc = fc_new - fc_prev
                    actual_delta = min(delta_fc, F[i-1,parent_idx[0]])  # do not consume more parent than available
                    
                    F[i,child_idx[0]] = fc_prev + actual_delta
                    F[i,parent_idx[0]] = F[i-1,parent_idx[0]] - actual_delta
                else:
                    km_active = False

                
                # Reverse martensite on heating (if reversible)
                if reversible and (T > Ms) and (T < T_lo) and F[i,child_idx[0]] > 1e-8:
                    fc_prev = F[i-1,child_idx[0]]
                    parent_prev = F[i-1,parent_idx[0]]
                    dt = time[i] - time[i-1] if 'time' in locals() else 1.0 # type: ignore
                    t_char = 2.0  # You may want to set this to a material-specific value
                    fc_new = backward_euler_martensitic(fc_prev, 0.0, dt, t_char)
                    # Limit martensite fraction to available parent
                    fc_new = max(min(fc_new, km_fp0 - f_pr), 0.0)
                    delta_fc = fc_new - fc_prev
                    F[i,child_idx[0]] = fc_new
                    F[i,parent_idx[0]] = parent_prev - delta_fc
    
                
            elif trans_type == "DIFFUSIONAL" and input_option in ("JMA", "TTT"):
                # Get equilibrium fraction and kinetics
                eq_label = f"ABQ_PHASE_TRANS_Trans_{trans_name}_ParentEquiFrac"
                eq_table = tables["ABQ_PHASE_TRANS_Diffusional_ParentEquiFrac"].get(eq_label, [])
                f_eq = np.float64(interp_table(T, eq_table, T_col=1, y_col=0)) if len(eq_table) > 0 else 0.0
                
                # Get kinetic parameters - Steel uses TTT diagrams primarily
                if input_option == "JMA":
                    jma_label = f"ABQ_PHASE_TRANS_Trans_{trans_name}_JMA_Coefficients"
                    jma_table = tables["ABQ_PHASE_TRANS_Diffusional_JMA_Coefficients"].get(jma_label, [])
                    if len(jma_table) == 0:
                        continue
                    k, n_ = interp_table(T, jma_table, T_col=2, y_col=0), interp_table(T, jma_table, T_col=2, y_col=1)
                else:  # TTT option - preferred for steel
                    ttt_diag_label = f"ABQ_PHASE_TRANS_Trans_{trans_name}_TTT_Diagram"
                    ttt_diag = tables.get("ABQ_PHASE_TRANS_Diffusional_TTT_Diagram", {}).get(ttt_diag_label, [])
                    if len(ttt_diag) == 0:
                        continue
                    ttt_const_label = f"ABQ_PHASE_TRANS_Trans_{trans_name}_TTT_Constants"
                    ttt_const = tables.get("ABQ_PHASE_TRANS_Diffusional_TTT_Constants", {}).get(ttt_const_label, [[0.01, 0.5]])
                    f_start, f_end = ttt_const[0][0], ttt_const[0][1]
                    t_start = interp_table(T, ttt_diag, T_col=2, y_col=0)
                    t_end = interp_table(T, ttt_diag, T_col=2, y_col=1)
                    
                    try:
                        n_ = np.log(np.log(1-f_start)/np.log(1-f_end)) / np.log(t_start/t_end)
                        k = -np.log(1-f_end) / (t_end**n_)
                    except:
                        continue
                
                # Add to active transformations
                active_diff_transformations.append({
                    'trans_name': trans_name,
                    'parent_idx': parent_idx,
                    'child_idx': child_idx,
                    'parent_names': parent_names,
                    'child_names': child_names,
                    'f_eq': f_eq,
                    'k': k,
                    'n': n_,
                    'reversible': reversible,
                    'trans_type': trans_type                        
                })
    
    
        transformation_rates=[]
        # Now process active diffusional competing transformations
        if len(active_diff_transformations) > 0: 
            process_competing_diffusional_transformations(active_diff_transformations, F, i, dt, tolerance, transformation_rates, VERBOSE)
                              
    return km_active, km_fp0
