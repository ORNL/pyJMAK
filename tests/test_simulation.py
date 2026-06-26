import numpy as np
import pyJMAK

def test_run_minimal():
    # Minimal "tables" with one phase to ensure simulation returns expected shape
    tables = {
        "ABQ_PHASE_TRANS_SolidPhases": {
            "ABQ_PHASE_TRANS_SolidPhases": [["Beta", 1.0]]
        },
        "ABQ_PHASE_TRANS_MeltingTemperature": {
            "ABQ_PHASE_TRANS_MeltingTemperature": [[980.0, 1000.0]]
        },
        "ABQ_PHASE_TRANS_Transformations": {
            "ABQ_PHASE_TRANS_Transformations": []
        }
    }
    time = np.linspace(0, 10, 11)
    T_C = np.linspace(20, 20, 11)
    model = pyJMAK.jmak(verbose=False)
    model.setAlloyType('ti64')
    model.setTimeTemp(time, T_C)
    model.setMaterialTables(tables)
    model.run_phase_trans_sim()
    assert len(model.ph_names) == 1
    assert model.F.shape == (len(time), 1)
    assert model.RLS_state.shape[0] == len(time)
