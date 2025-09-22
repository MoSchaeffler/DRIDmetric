#!/usr/bin/python3

from DRIDmetric import DRID

###############################
############ Input ############
###############################

"""
define topology and trajectory path
"""

top  = "./data/Abeta42.tpr"
traj = "./data/Abeta42.xtc"

"""
define selection for centroids and reference atoms in the MDAnalysis selection style
    - https://userguide.mdanalysis.org/stable/selections.html

The selection below was used for Abeta42 in https://doi.org/10.1039/D4CC02856B
"""

sel_cent = "(name CA and resid 28) or (name CA and resid 23) or (name CA and resid 1) or (name CA and resid 42) or (name CA and resid 19) or (name CA and resid 34)"
sel_atoms = "protein"

"""
Define name of the output file containing the framewise DRID metric of the system
saved as .npy array
"""

outname = "DRID_Abeta42"

###############################
############  Run  ############
###############################

df = DRID(top, traj, sel_atoms, sel_cent)
df.run(outname)
