# DRIDmetric - Example

This README walks you through computing the DRID metric for an Aβ42 trajectory using the Python class provided in this project. The example uses the two files in `data/`:

* `Abeta42.xtc` — GROMACS trajectory
* `Abeta42.tpr` — GROMACS run input/topology

The DRID calculation is performed by the `DRID` class defined in `calculate_DRID_metric.py`. A helper script `calc_DRID.py` shows how to use it.

> Directory layout:
>
> ```text
> example/
> ├─ calc_DRID.py
> ├─ data/
> │   ├─ Abeta42.xtc
> │   └─ Abeta42.tpr
> ```

---

## 1) Install

### Python

Use Python ≥ 3.9 (tested with 3.12).

### Create a virtual environment

```bash
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate
```

### Install packckage + dependencies

```bash
pip install -e git+http://github.com/MoSchaeffler/DRIDmetric.git#egg=DRIDmetric
```

Dependencies:

* `numpy`
* `tqdm`
* `MDAnalysis`

The dependencies should be intalled automatically alongside the package

---

## 2) Running the example

The main way to use this package is via the provided `calc_DRID.py` script. 
Edit: 
- the file to point to your trajectory and topology 
- the centroid selection
- atom selection for reference atoms (usually "protein" i.e. atoms not included in centroids) 

Selections follow [MDAnalysis selection syntax](https://userguide.mdanalysis.org/stable/selections.html).

Example (for the Aβ42 test data in `data/`):

```python
from DRIDmetric import DRID

# Input files
top = "data/Abeta42.tpr"
traj = "data/Abeta42.xtc"

# Atom/centroid selections (MDAnalysis syntax)
sel_cent = "(name CA and resid 28) or (name CA and resid 23) or (name CA and resid 1) or (name CA and resid 42) or (name CA and resid 19) or (name CA and resid 34)"
sel_atoms = "protein"

# Output name
outname = "results/DRID_abeta42"

# Run DRID
calc = DRID(top, traj, sel_atoms, sel_cent)
calc.run(outname)
```

This will produce a NumPy `.npy` file (`results/DRID_abeta42.npy`) containing the framewise DRID metric.

**Note:**
The centroid selction is specific to each molecule. Either make an educated choice e.g. choose CA of residues that are crucial to the process under study, or choose CA of equidistantly spaced residues e.g.  "(name CA and resid 1) or (name CA and resid 6) or (name CA and resid 11) ...". While a more educated selection can improve results, consequent analysis of the dynamics of the system with for example the freeEnergyCalculation package has been shown to be relatively robust against the selection of centroids (https://doi.org/10.1039/D4CC02856B).

---

## 3) Interpreting the output

* Output is a NumPy array of shape `(N_frames, N_centroids, 3)`.

  * Axis 0: trajectory frames
  * Axis 1: selected centroids
  * Axis 2: first three DRID moments `(μ, ν, ξ)`

Load it back in Python with:

```python
import numpy as np
DRID_data = np.load("results/DRID_abeta42.npy")
print(DRID_data.shape)
```

---

## 4) Citation

If you use this code/data, please cite:

* The DRIDmetric repository
* The Aβ42 paper: [https://doi.org/10.1039/D4CC02856B](https://doi.org/10.1039/D4CC02856B)

---

*Happy analyzing!*

