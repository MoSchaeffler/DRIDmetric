#!/usr/bin/python3


import numpy as np
from tqdm import tqdm
import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array as darray


class DRID(object):
    """Calculate the firs three moments of the 
    Distribution of Reciprocal Interatomic Distances (DRID) 
    for a given Molecular Dynamics trajectory

    Parameters
    ----------
    top:
        Topology file (GRO/TPR/PDB/etc.)
    
    traj:
        Trajectory file (XTC/TRR/DATA/etc.)

    centroid_selection:
        string with MDAnalysis selection syntax for the centroids
        e.g. "(name CA and resid 28) or (name CA and resid 23) or (name CA and resid 1) ..."
        
    atom_selection:
        string with MDAnalysis selection syntax for the reference group
        e.g. "protein"
    """

    def __init__(self,
                 top:str,
                 traj:str,
                 atom_selection:str,
                 centroid_selection:str
                 ) -> None:

        # create universe
        self.u = mda.Universe(top, traj)
        self.Nf = len(self.u.trajectory)

        # get disjoined atom/centroid groups
        self.centroids = self.u.select_atoms(centroid_selection)
        self.Nc = len(self.centroids)

        print("Number of centroids selected: {}".format(self.Nc))

        self.atom_selection = atom_selection

        group = self.u.select_atoms(atom_selection)
        self.atoms = group - self.centroids
        self.Na = len(self.atoms)

        print("Number of atoms selected: {}".format(self.Na))

    def run(self,
            outname:str="DRID"
    ) -> None:
        """Method for calculating the DRID metric for each frame of the input trajecory
        saves the output as numpy array

        Parameters
        ----------
        outname:
            name of the output file
        """

        print("Start Calculating DRID.")

        drid = np.zeros((self.Nf, self.Nc, 3))

        for f in tqdm(range(self.Nf), desc="Processing frames"):

            self.u.trajectory[f]
            # box = self.u.dimensions

            mu, nu, xi = np.zeros(self.Nc), np.zeros(self.Nc), np.zeros(self.Nc)

            for i in range(self.Nc):

                # centroid
                cent = self.centroids[i]

                # get reference group for centroid excluding bonded neigbours
                group = self.u.select_atoms(self.atom_selection)
                atoms = group - self.centroids

                bonds = cent.get_connections("bonds")
                for b in bonds:
                    atoms = atoms - b.atoms

                # calculate moments

                mu[i], nu[i], xi[i] = self.moments(cent.position, atoms.positions)

            drid[f, :, 0] = mu
            drid[f, :, 1] = nu
            drid[f, :, 2] = xi

        np.save(outname, drid)
        print("DRID calculation complete and saved.")


    def moments(self, 
                c_pos:np.ndarray,
                atom_pos:np.ndarray
    ) -> tuple[float, float, float]:
        """Method for calculating the first three moments of the DRID

        Parameters
        ----------
        c_pos:
            Position vector of the centroid
        atom_pos:
            (n,3) dimensional array of refernce atom position vectors

        Returns
        -------
        Tuple:
            Float values of first three moments
        """

        dij = darray(c_pos, atom_pos)
        s = np.sum(1 / dij)

        mu = 1 / (len(atom_pos) - 1) * s

        s1 = np.sum(1 / (dij - mu) ** 2)
        s2 = np.sum(1 / (dij - mu) ** 3)

        nu = (1 / (len(atom_pos) - 1) * s1) ** (1 / 2)
        xi = (1 / (len(atom_pos) - 1) * s2) ** (1 / 3)

        return mu, nu, xi
