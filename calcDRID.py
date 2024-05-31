#!/usr/bin/python3

import numpy as np

# MD Analysis for loading trajectories
import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array as darray


class DRID(object):

    def __init__(self, top, traj, atom_selection, centroid_selection):

        # create universe
        self.u = mda.Universe(top, traj)

        self.Nf = len(self.u.trajectory)

        print("Trajectory length {}".format(self.Nf))

        # get disjoined atom/centroid groups
        self.centroids = self.u.select_atoms(centroid_selection)
        self.Nc = len(self.centroids)

        print("Number of centroids selected: {}".format(self.Nc))

        self.atom_selection = atom_selection

        group = self.u.select_atoms(atom_selection)
        self.atoms = group - self.centroids
        self.Na = len(self.atoms)

        print("Number of atoms selected: {}".format(self.Na))

    def run(self, outname="DRID"):

        print("Start Calculating DRID")

        drid = np.zeros((self.Nf, self.Nc, 3))

        for f in range(self.Nf):

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

            # write to coord file
            drid[f, :, 0] = mu
            drid[f, :, 1] = nu
            drid[f, :, 2] = xi

            try:
                if f % int(self.Nf / 10) == 1:
                    print(int(f / self.Nf * 100), '%')
            except:
                continue

        np.save(outname + ".npy", drid)

    # first moment of DRID
    def moments(self, c_pos, atom_pos):

        dij = darray(c_pos, atom_pos)
        s = np.sum(1 / dij)

        mu = 1 / (len(atom_pos) - 1) * s

        s1 = np.sum(1 / (dij - mu) ** 2)
        s2 = np.sum(1 / (dij - mu) ** 3)

        nu = (1 / (len(atom_pos) - 1) * s1) ** (1 / 2)
        xi = (1 / (len(atom_pos) - 1) * s2) ** (1 / 3)

        return mu, nu, xi
