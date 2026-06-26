"""
Example file for Ti-Al-4V to run pyJMAK code for phase-transformation 
"""
import pyJMAK
import os

MATERIAL_INPUT_FILE = "abq_phase_trans_input-Ti64.inp"  # Set this to your abaqus-formatted material-specific input file: contains info about phases, transformation rules
TEMP_FILE = "temperature_center_node_2908_truncated.txt"   # Set this to your temperature file, contains 2 columns: time and temperature
COMPARE_WITH_EXPERIMENTAL = True
exp_data_files = {
        "Beta": "beta_AQ.txt", #Format: "phase_name": "file_name"
        "AlphaW": "alpha_w_AQ.txt",
        "Alphaprime": "alpha_prime_AQ.txt"
}

###############################################################################
#       change dir to current directory
###############################################################################
abspath = os.path.abspath(__file__)
path = os.path.dirname(abspath)
os.chdir(path)
###############################################################################

#Start
model = pyJMAK.jmak()

# Load the input file
model.loadInputFiles(TEMP_FILE, MATERIAL_INPUT_FILE)

# Main function computes the evolution of phase fractions over time
model.run_phase_trans_sim() 

# Indicate if the results are to be compared with the experimental/literature results
model.setExperimentalComparison(COMPARE_WITH_EXPERIMENTAL, exp_data_files) 

# Plot evolutions of temperature, phase fractions and RSL-state over time
model.plot_results(saveFig=True)
