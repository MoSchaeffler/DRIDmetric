# DRIDmetric

**DRIDmetric** is a Python package to calculate the first three moments of the **Distribution of Reciprocal Interatomic Distances (DRID)** for molecular dynamics (MD) trajectories.  

The DRID representation is a **structure-preserving dimensionality reduction method** that captures essential kinetic and structural information, making it highly suitable for clustering, free energy landscape analysis, and kinetic modeling.

---

## Features

- Compute **frame-wise DRID vectors** from MD trajectories  
- Based on **reciprocal interatomic distances** to emphasize short-range contacts  
- Extract **first three moments** of distance distributions for selected centroids:
  - Mean (μ)  
  - Variance (ν)  
  - Skewness (ξ)  
- Simple integration with [**MDAnalysis**](https://www.mdanalysis.org/)  
- Outputs results as NumPy arrays for direct use in clustering or ML pipelines  

---

## Mathematical Definition

Given:
- a set of **m centroids** \( C = \{c_i\} \),  
- a set of **N reference atoms** \( A = \{a_j\} \), excluding covalently bound neighbors,  

the distribution of **reciprocal distances** is defined by:  

\[
d_{ij} = \| c_i - a_j \|, \quad \text{with reciprocal } \frac{1}{d_{ij}}
\]

The **first three moments** for centroid \(i\) are:

\[
\mu_i = \frac{1}{N-1-n_b^i} \sum_j \frac{1}{d_{ij}}
\]

\[
\nu_i = \left[ \frac{1}{N-1-n_b^i} \sum_j \frac{1}{(d_{ij}-\mu_i)^2} \right]^{1/2}
\]

\[
\xi_i = \left[ \frac{1}{N-1-n_b^i} \sum_j \frac{1}{(d_{ij}-\mu_i)^3} \right]^{1/3}
\]

where \( n_b^i \) is the number of covalent bonds for centroid \(c_i\).  

The DRID vector of a frame is then a **3m-dimensional vector**:  

\[
(\mu_1,\nu_1,\xi_1,\;\ldots,\;\mu_m,\nu_m,\xi_m)
\]

The **distance metric** between two conformations \(j\) and \(k\) in DRID space is:  

\[
s_{jk} = \frac{1}{3m} \sum_{i=1}^m \left[ (\mu_i^j-\mu_i^k)^2 + (\nu_i^j-\nu_i^k)^2 + (\xi_i^j-\xi_i^k)^2 \right]^{1/2}
\]

---

## Installation

Install directly within your environment

```bash
pip install git+http://github.com/MoSchaeffler/DRIDmetric.git#egg=DRIDmetric
```

or clone and install in editable mode

```bash
git clone https://github.com/MoSchaeffler/DRIDmetric.git
cd DRIDmetric
pip install -e .
```

This installs:
- `numpy`
- `tqdm`
- `MDAnalysis`

## `DRID` class

```python
from DRIDmetric import DRID

DRID(
    top: str,
    traj: str,
    atom_selection: str,
    centroid_selection: str,
)
```

**Parameters**

- `top` *(str)*  
  Topology file (e.g., `.gro`, `.tpr`, `.pdb`, …).

- `traj` *(str)*  
  Trajectory file (e.g., `.xtc`, `.trr`, `.dcd`, …).

- `centroid_selection` *(str)*  
  MDAnalysis **selection string** defining centroid atoms, e.g.:  
  `"(name CA and resid 28) or (name CA and resid 23) or (name CA and resid 1) ..."`

- `atom_selection` *(str)*  
  MDAnalysis **selection string** for the reference atom group (atoms compared **to** the centroids), e.g.:  
  `"protein"`

> Notes  
> • Covalently bonded neighbors of each centroid are automatically excluded from its reference set.  
> • The class iterates frames of `traj`, computing \((\mu,\nu,\xi)\) per centroid and saving an array of shape `(n_frames, n_centroids, 3)`.

#### Methods

- `run(outname: str = "DRID") -> np.array`  
  Computes DRID for all frames and saves a NumPy file `outname.npy`.


## Usage

**Minimal Example**
```python
from DRIDmetric import DRID

# Input files
top = "/path/to/your/topology.tpr"
traj = "/path/to/your/trajectory.xtc"

# MDAnalysis selections
sel_atoms = "protein"
sel_cent  = "(name CA and resid 1) or (name CA and resid 11) or (name CA and resid 21)"

drid = DRID(top, traj, sel_atoms, sel_cent)
drid.run("DRID_output")
```

Please refer to the **example** for detailed usage instructions

## Reference

If you use this package, please cite:

- Schäffler, Wales, & Strodel (2024) *Chem. Commun.*  
  “The energy landscape of Aβ42: a funnel to disorder for the monomer becomes a folding funnel for self-assembly.”

Background on DRID:

- Zhou & Caflisch (2012) *J. Chem. Theory Comput.*  
- Chakraborty, Straub & Thirumalai (2023) *Sci. Adv.*