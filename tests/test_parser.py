import tempfile
from pyJMAK.parser import parse_abaqus_tables

def test_parse_basic_table():
    sample = '''*PARAMETER TABLE, TYPE="ABQ_PHASE_TRANS_SolidPhases", LABEL="ABQ_PHASE_TRANS_SolidPhases"
    "Beta", 1.0
    "AlphaW", 0.0
    '''
    tf = tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.inp')
    tf.write(sample)
    tf.flush()
    tf.close()
    tables = parse_abaqus_tables(tf.name)
    assert "ABQ_PHASE_TRANS_SolidPhases" in tables
    assert "ABQ_PHASE_TRANS_SolidPhases" in tables["ABQ_PHASE_TRANS_SolidPhases"]
    rows = tables["ABQ_PHASE_TRANS_SolidPhases"]["ABQ_PHASE_TRANS_SolidPhases"]
    assert rows[0][0] == "Beta"
