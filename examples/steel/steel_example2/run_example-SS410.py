"""
Example file for SS410 to run pyJMAK code for phase-transformation 
"""
import pyJMAK
import os

MATERIAL_INPUT_FILE = "abq_phase_trans_input-SS410.inp"  # Set this to your abaqus-formatted material-specific input file: contains info about phases, transformation rules
TEMP_FILE = "time_temp.csv"   # Set this to your temperature file, contains 2 columns: time and temperature


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

# Plot evolutions of temperature, phase fractions and RSL-state over time
model.plot_results(saveFig=True)