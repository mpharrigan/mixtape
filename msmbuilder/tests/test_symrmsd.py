from msmbuilder.libdistance import cdist_symrmsd, dist_symrmsd, dist, cdist
import mdtraj as md
import numpy as np
from random import shuffle
import itertools


def get_shuffle_traj():
    """Make a trajectory with two particles at each point x=10, 20, 30, 40
    with some random noise. The second one has the atom groups permuted.
    """
    top = md.Topology()
    c = top.add_chain()
    r = top.add_residue('HET', c)
    for _ in range(8):
        top.add_atom('C', md.element.carbon, r)

    x1 = []
    x2 = []
    perminds = list(range(4))
    shuffle(perminds)
    for i, j in zip(range(4), perminds):
        x1 += [[10 * (i + 1), 0, 0]]
        x1 += [[10 * (i + 1), 0, 0]]

        x2 += [[10 * (j + 1), 0, 0]]
        x2 += [[10 * (j + 1), 0, 0]]

    xyz1 = np.asarray([x1] * 7, dtype=np.float32)
    xyz2 = np.asarray([x2] * 7, dtype=np.float32)
    xyz1 += np.random.rand(*xyz1.shape)
    xyz2 += np.random.rand(*xyz2.shape)

    traj1 = md.Trajectory(xyz1, top)
    traj2 = md.Trajectory(xyz2, top)
    return traj1, traj2


def test_normal():
    # Make sure we *need* permutation by making sure rmsd is big without
    for attempts in range(10):
        traj1, traj2 = get_shuffle_traj()
        traj1.center_coordinates()
        traj2.center_coordinates()
        dists = dist(traj1, traj2[0], 'rmsd')
        assert dists.ndim == 1
        if np.any(dists > 5):
            return True
        print("attempt", attempts)
    assert False


def test_dist():
    traj1, traj2 = get_shuffle_traj()
    groups = list(zip(range(0, 7, 2), range(1, 8, 2)))
    print("Groups", groups)
    assert len(groups) == 4

    perms = list(itertools.permutations(range(4)))
    assert len(perms) == 4 * 3 * 2

    dists = dist_symrmsd(traj1, traj2[0], np.asarray(groups), np.asarray(perms))
    assert dists.ndim == 1
    assert np.all(dists < 5)


def test_cdist():
    traj1, traj2 = get_shuffle_traj()
    groups = list(zip(range(0, 7, 2), range(1, 8, 2)))
    print("Groups", groups)
    assert len(groups) == 4

    perms = list(itertools.permutations(range(4)))
    assert len(perms) == 4 * 3 * 2

    dists = cdist_symrmsd(traj1, traj2, np.asarray(groups), np.asarray(perms))
    assert dists.ndim == 2
    assert dists.shape == (len(traj1), len(traj2))
    assert np.all(dists < 5)
