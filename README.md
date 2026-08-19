# PyJMAK — An Open-Source Python Toolkit for Modeling Solid-State Metallurgical Phase Transformations 

A framework validated for Ti-4Al-4V and general steel. Applicable to heat treatment, with initial phase fractions provided

A lightweight Python toolkit for simulating solid-state phase transformations using common kinetic models:

- **JMAK** (Johnson–Mehl–Avrami–Kolmogorov) additivity
- **Koistinen–Marburger (KM)** model for martensitic transformation
- Simple handling of multiple martensitic phases

It supports parsing Abaqus-like material definition tables and producing time-resolved phase fraction predictions.  
Includes both a **Python API** and a **command-line interface (CLI)**.

---

## ✨ Features
- Parse Abaqus-formatted phase transformation tables (`parse_abaqus_tables`)
- Run phase transformation simulations (`run_phase_trans_sim`)
- JMAK and KM kinetic models
- Compare simulated results with experimental data
- Command-line tool: `jmak-run`

---

## 📦 Installation

### MacOS / Linux
```bash
# Create and activate a virtual environment
python3 -m venv myenv
source myenv/bin/activate

# Upgrade packaging tools
pip install --upgrade pip setuptools wheel

# Install jmak in editable mode
pip install -e .

# Optionally, install test dependencies
pip install -e '.[test]' pytest
```

### Windows (PowerShell)
```powershell
# Create and activate a virtual environment
python -m venv myenv
myenv\Scripts\Activate.ps1

# Upgrade packaging tools
pip install --upgrade pip setuptools wheel

# Install jmak in editable mode
pip install -e .

# Optionally, install test dependencies
pip install -e ".[test]" pytest
```

---

## ▶️ Usage

### Python API
```python
python3 examples/Ti6Al4V/Ti64_example2/run_example-Ti64.py 
--------------------
import pyJMAK
model = pyJMAK.jmak()
# Load the input file
model.loadInputFiles(TEMP_FILE, MATERIAL_INPUT_FILE)
# Main function computes the evolution of phase fractions over time
model.run_phase_trans_sim() 

```

### CLI
```bash
# Basic run
jmak-run -i ../examples/steel/steel_example1/abq_phase_trans_input-Steel.inp -t ../examples/steel/steel_example1/Temp-1.txt 
```

---

## 🧪 Running Tests
```bash
pytest -q
```

---

## 📂 Input File Guidelines
- The **material input file** follows Abaqus material definition syntax for phase transformation data.
- Refer to the **Abaqus documentation** for the correct formatting of the phase-specific tables.
- When comparing with experimental results:
  - **Phase names must match exactly** (case-sensitive) between the experimental data and the material input file.

---

## 🛠 Development 
Clone this repository and install in development mode:
```bash
git clone https://github.com/yourusername/pyjmak.git (update)
cd pyjmak
pip install -e '.[test]'
```


