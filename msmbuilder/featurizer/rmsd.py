# Authors: Matthew Harrigan <matthew.harrigan@outlook.com>
#          Robert McGibbon <rmcgibbo@gmail.com>
#          Brooke Husic <brookehusic@gmail.com>
# Copyright (c) 2016, Stanford University and the Authors
# All rights reserved.

import warnings

import mdtraj as md

from msmbuilder import libdistance
from .base import Featurizer


class RMSDFeaturizer(Featurizer):
    """Featurizer based on RMSD to one or more reference structures.

    This featurizer inputs a trajectory to be analyzed ('traj') and a
    reference trajectory ('ref') and outputs the RMSD of each frame of
    traj with respect to each frame in ref. The output is a numpy array
    with n_rows = traj.n_frames and n_columns = ref.n_frames.

    Parameters
    ----------
    reference_traj : md.Trajectory
        The reference conformations to superpose each frame with respect to
    atom_indices : np.ndarray, shape=(n_atoms,), dtype=int
        The indices of the atoms to superpose and compute the distances with.
        If not specified, all atoms are used.
    trj0
        Deprecated. Please use reference_traj.
    """

    def __init__(self, reference_traj=None, atom_indices=None, trj0=None):
        if trj0 is not None:
            warnings.warn("trj0 is deprecated. Please use reference_traj",
                          DeprecationWarning)
            reference_traj = trj0
        else:
            if reference_traj is None:
                raise ValueError("Please specify a reference trajectory")

        self.atom_indices = atom_indices
        if self.atom_indices is not None:
            self.sliced_reference_traj = reference_traj.atom_slice(
                self.atom_indices)
        else:
            self.sliced_reference_traj = reference_traj

    def partial_transform(self, traj):
        """Featurize an MD trajectory into a vector space via distance
        after superposition

        Parameters
        ----------
        traj : mdtraj.Trajectory
            A molecular dynamics trajectory to featurize.

        Returns
        -------
        features : np.ndarray, shape=(n_frames, n_ref_frames)
            The RMSD value of each frame of the input trajectory to be
            featurized versus each frame in the reference trajectory. The
            number of features is the number of reference frames.

        See Also
        --------
        transform : simultaneously featurize a collection of MD trajectories
        """
        if self.atom_indices is not None:
            sliced_traj = traj.atom_slice(self.atom_indices)
        else:
            sliced_traj = traj
        result = libdistance.cdist(
            sliced_traj, self.sliced_reference_traj, 'rmsd'
        )
        return result
