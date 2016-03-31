# Author: Robert McGibbon <rmcgibbo@gmail.com>
# Contributors: Matthew Harrigan <matthew.harrigan@outlook.com>
# Copyright (c) 2015, Stanford University
# All rights reserved.

from __future__ import absolute_import, print_function, division

import errno
import getpass
import glob
import os
import re
import socket
import sys
import warnings
from collections import Sequence
from datetime import datetime
from os.path import join, exists, expanduser

import mdtraj as md
import numpy as np
import tables
from mdtraj.core.trajectory import _parse_topology

from . import version

_PYTABLES_DISABLE_COMPRESSION = tables.Filters(complevel=0)

__all__ = ['dataset']


def dataset(path, mode='r', fmt=None, verbose=False, **kwargs):
    """Open a dataset object

    MSMBuilder supports several dataset 'formats' for storing
    lists of sequences on disk.

    This function can also be used as a context manager.

    Parameters
    ----------
    path : str
        The path to the dataset on the filesystem
    mode : {'r', 'w', 'a'}
        Open a dataset for reading, writing, or appending. Note that
        some formats only support a subset of these modes.
    fmt : {'dir-npy', 'hdf5', 'mdtraj'}
        The format of the data on disk

        ``dir-npy``
            A directory of binary numpy files, one file per sequence

        ``hdf5``
            A single hdf5 file with each sequence as an array node

        ``mdtraj``
            A read-only set of trajectory files that can be loaded
            with mdtraj

    verbose : bool
        Whether to print information about the dataset

    """

    if mode == 'r' and fmt is None:
        fmt = _guess_format(path)
    elif mode in 'wa' and fmt is None:
        raise ValueError('mode="%s", but no fmt. fmt=%s' % (mode, fmt))

    if fmt == 'dir-npy':
        return NumpyDirDataset(path, mode=mode, verbose=verbose, **kwargs)
    elif fmt == 'mdtraj':
        return MDTrajDataset(path, mode=mode, verbose=verbose, **kwargs)
    elif fmt == 'hdf5':
        return HDF5Dataset(path, mode=mode, verbose=verbose)
    elif fmt.endswith("-union"):
        raise ValueError("union datasets have been removed. "
                         "Please use msmbuilder.featurizer.FeatureUnion")
    else:
        raise NotImplementedError("Unknown format fmt='%s'" % fmt)


def _guess_format(path):
    """Guess the format of a dataset based on its filename / filenames.
    """
    if os.path.isdir(path):
        return 'dir-npy'

    if path.endswith('.h5') or path.endswith('.hdf5'):
        # TODO: Check for mdtraj .h5 file
        return 'hdf5'

    # TODO: What about a list of trajectories, e.g. from command line nargs='+'
    return 'mdtraj'


class _BaseDataset(Sequence):
    _PROVENANCE_TEMPLATE = '''MSMBuilder Dataset:
  MSMBuilder:\t{version}
  Command:\t{cmdline}
  Path:\t\t{path}
  Username:\t{user}
  Hostname:\t{hostname}
  Date:\t\t{date}
  Comments:\t\t{comments}
'''
    _PREV_TEMPLATE = '''
== Derived from ==
{previous}
'''

    def __init__(self, path, mode='r', verbose=False):
        self.path = path
        self.mode = mode
        self.verbose = verbose

        if mode not in ('r', 'w', 'a'):
            raise ValueError('mode must be one of "r", "w", "a"')
        if mode in 'wa':
            if mode == 'w' and exists(path):
                raise ValueError('File exists: %s' % path)
            try:
                os.makedirs(path)
            except OSError:
                pass
            self._write_provenance()

    def create_derived(self, out_path, comments='', fmt=None):
        if fmt is None:
            out_dataset = self.__class__(out_path, mode='w',
                                         verbose=self.verbose)
        else:
            out_dataset = dataset(out_path, mode='w', verbose=self.verbose,
                                  fmt=fmt)
        out_dataset._write_provenance(previous=self.provenance,
                                      comments=comments)
        return out_dataset

    def apply(self, fn):
        for key in self.keys():
            yield fn(self.get(key))

    def _build_provenance(self, previous=None, comments=''):
        val = self._PROVENANCE_TEMPLATE.format(
            version=version.full_version,
            cmdline=' '.join(sys.argv),
            user=getpass.getuser(),
            hostname=socket.gethostname(),
            path=self.path,
            comments=comments,
            date=datetime.now().strftime("%B %d, %Y %I:%M %p"))
        if previous:
            val += self._PREV_TEMPLATE.format(previous=previous)
        return val

    def fit_with(self, estimator):
        """Call the fit method of the estimator on this dataset

        Parameters
        ----------
        estimator : BaseEstimator
            estimator.fit will be called on this dataset.

        Returns
        -------
        estimator
            The fit estimator.
        """
        estimator.fit(self)
        return estimator

    def transform_with(self, estimator, out_ds, fmt=None):
        """Call the partial_transform method of the estimator on this dataset

        Parameters
        ----------
        estimator : object with ``partial_fit`` method
            This object will be used to transform this dataset into a new
            dataset. The estimator should be fitted prior to calling
            this method.
        out_ds : str or Dataset
            This dataset will be transformed and saved into out_ds. If
            out_ds is a path, a new dataset will be created at that path.
        fmt : str
            The type of dataset to create if out_ds is a string.

        Returns
        -------
        out_ds : Dataset
            The tranformed dataset.
        """
        if isinstance(out_ds, str):
            out_ds = self.create_derived(out_ds, fmt=fmt)
        elif isinstance(out_ds, _BaseDataset):
            err = "Dataset must be opened in write mode."
            assert out_ds.mode in ('w', 'a'), err
        else:
            err = "Please specify a dataset path or an existing dataset."
            raise ValueError(err)
        for key in self.keys():
            out_ds[key] = estimator.partial_transform(self[key])
        return out_ds

    def fit_transform_with(self, estimator, out_ds, fmt=None):
        """Create a new dataset with the given estimator.

        The estimator will be fit by this dataset, and then each trajectory
        will be transformed by the estimator.

        Parameters
        ----------
        estimator : BaseEstimator
            This object will be fit and used to transform this dataset
            into a new dataset.
        out_ds : str or Dataset
            This dataset will be transformed and saved into out_ds. If
            out_ds is a path, a new dataset will be created at that path.
        fmt : str
            The type of dataset to create if out_ds is a string.

        Returns
        -------
        out_ds : Dataset
            The transformed dataset.

        Examples
        --------
        diheds = dataset("diheds")
        tica = diheds.fit_transform_with(tICA(), 'tica')
        kmeans = tica.fit_transform_with(KMeans(), 'kmeans')
        msm = kmeans.fit_with(MarkovStateModel())
        """
        self.fit_with(estimator)
        return self.transform_with(estimator, out_ds, fmt=fmt)

    @property
    def provenance(self):
        raise NotImplementedError('implemented in subclass')

    def _write_provenance(self, previous=None, comments=''):
        raise NotImplementedError('implemented in subclass')

    def __len__(self):
        return sum(1 for _ in self.keys())

    def __getitem__(self, i):
        if isinstance(i, slice):
            err = "Please index datasets with explicit indices or ds[:]."
            assert i.start is None and i.step is None and i.stop is None, err
            return [self.get(i) for i in self.keys()]
        return self.get(i)

    def __setitem__(self, i, x):
        return self.set(i, x)

    def __iter__(self):
        for key in self.keys():
            yield self.get(key)

    def keys(self):
        # keys()[i], get(i) and set(i, x) should all follow
        # the same ordering convention for the indices / items.
        raise NotImplementedError('implemeneted in subclass')

    def items(self):
        for key in self.keys():
            yield (key, self.get(key))

    def get(self, i):
        raise NotImplementedError('implemeneted in subclass')

    def set(self, i, x):
        raise NotImplementedError('implemeneted in subclass')

    def close(self):
        pass

    def flush(self):
        pass

    def __enter__(self):
        return self

    def __exit__(self, *exc_info):
        self.close()


class NumpyDirDataset(_BaseDataset):
    """A dataset consisting of a directory full of binary numpy files

    Parameters
    ----------
    path : str
    mode : {'r', 'w', 'a'}
        Read, write, or append. If mode is set to 'a' or 'w',
        duplicate keys will be overwritten.
    verbose : bool, default=False
        Whether to print saving and loading status to stdout
    zero_pad : bool, default=True
        Whether to pad filenames with zeros so lexicographical sort works.
        If set to false, you must do a "version sort" or "natural sort" when
        interacting with the files from outside a dataset object. For example,
        use `ls -v` in bash. Note that the dataset object will always sort
        files and directories correctly by casting their filename to an
        integer.
    strict_filenames : bool, default=False
        If this and zero_pad is true, only read from files that have been
        zero-padded. Exclude files and directories with insufficient zero
        padding. If zero_pad is false, this must be false. If this is false,
        cast filenames to integers as best you can.

    Examples
    --------
    with NumpyDirDataset("path", 'w') as ds:
        ds[4] = np.arange(30)
    """

    provenance_fn = 'PROVENANCE.txt'

    def __init__(self, path, mode='r', verbose=False,
                 zero_pad=True, strict_filenames=None):
        super(NumpyDirDataset, self).__init__(path, mode=mode, verbose=verbose)

        if (not zero_pad) and strict_filenames:
            raise ValueError("If you're not writing zero-padded filenames, "
                             "it's probably a bad idea to set it so you won't "
                             "be able to read them! strict_filenames cannot "
                             "be true while zero_pad is false")

        # Starting in v3.4, we allow an arbitrary depth of nested integer
        # indices (for example project/run/clone in foldingathome)
        # Use self.inter_* variables to configure this

        if zero_pad:
            # We write out padded to 8 places. If the number does not fit into
            # this width, it will spill over and not truncate (good)
            # We formulate the regex to accept 8 or more digits
            self.item_fmt = "{index:08d}.npy"
            self.inter_fmt = "{index:05d}"
        else:
            self.item_fmt = "{index:d}.npy"
            self.inter_fmt = "{index:d}"

        if strict_filenames:
            self.item_re = re.compile(r'(\d{8}\d*).npy')
            self.inter_re = re.compile(r'(\d{5}\d*)')
        else:
            self.item_re = re.compile(r'(\d+).npy')
            self.inter_re = re.compile(r'(\d+)')

    def _get_filename(self, index, make_dirs=False):
        # Handle tuples
        try:
            filename = join(self.path, *[self.inter_fmt.format(index=i)
                                         for i in index[:-1]])
            i = index[-1]

            if make_dirs:
                try:
                    os.makedirs(filename)
                except OSError as e:
                    if e.errno != errno.EEXIST:
                        raise
        except TypeError:
            filename = self.path
            i = index

        filename = join(filename, self.item_fmt.format(index=i))
        return filename

    def get(self, i, mmap=False):
        mmap_mode = 'r' if mmap else None

        filename = self._get_filename(i)
        if self.verbose:
            print('[NumpydirDataset] loading %s' % filename)
        try:
            return np.load(filename, mmap_mode)
        except IOError as e:
            raise IndexError(e)

    def set(self, i, x):
        if self.mode not in 'wa':
            raise IOError('Dataset not opened for writing')

        filename = self._get_filename(i, make_dirs=True)
        if self.verbose:
            print('[NumpydirDataset] saving %s' % filename)
        return np.save(filename, x)

    def _keys_recurse(self, curdir, inds):
        for fn in os.listdir(curdir):
            full_fn = join(curdir, fn)
            if fn == self.provenance_fn:
                # Prevent warnings for provenance file
                continue

            if os.path.isdir(full_fn):
                match = self.inter_re.match(fn)
                if match is not None:
                    newinds = inds + (int(match.group(1)),)
                    newdir = full_fn
                    # wish we had yield from...
                    for ret in self._keys_recurse(newdir, newinds):
                        yield ret
                else:
                    warnings.warn("Unknown directory found in '{}': '{}'"
                                  .format(curdir, fn))
                continue

            # What if there is a file with the same index as a directory?
            # I think this is ok
            match = self.item_re.match(fn)
            if match is not None:
                key = int(match.group(1))
                if len(inds) > 0:
                    yield inds + (key,)
                else:
                    yield key
            else:
                warnings.warn("Unknown file found in '{}': '{}'"
                              .format(curdir, fn))

    def keys(self):
        root = os.path.expanduser(self.path)
        # wish we had yield from...
        for key in sorted(self._keys_recurse(root, tuple()),
                          key=_tuples_and_ints):
            yield key

    @property
    def provenance(self):
        try:
            with open(join(self.path, self.provenance_fn), 'r') as f:
                return f.read()
        except IOError:
            return 'No available provenance'

    def _write_provenance(self, previous=None, comments=''):
        with open(join(self.path, self.provenance_fn), 'w') as f:
            p = self._build_provenance(previous=previous, comments=comments)
            f.write(p)


class HDF5Dataset(_BaseDataset):
    _ITEM_FORMAT = 'arr_%d'
    _ITEM_RE = re.compile('arr_(\d+)')

    def __init__(self, path, mode='r', verbose=False):
        if mode not in ('r', 'w'):
            raise ValueError('mode must be one of "r", "w"')
        if mode == 'w':
            if exists(path):
                raise ValueError('File exists: %s' % path)

        self._handle = tables.open_file(path, mode=mode,
                                        filters=_PYTABLES_DISABLE_COMPRESSION)
        self.path = path
        self.mode = mode
        self.verbose = verbose

        if mode == 'w':
            self._write_provenance()

    def __getstate__(self):
        # pickle does not like to pickle the pytables handle, so...
        # self.flush()
        return {'path': self.path, 'mode': self.mode, 'verbose': self.verbose}

    def __setstate__(self, state):
        self.path = state['path']
        self.mode = state['mode']
        self.verbose = state['verbose']
        self._handle = tables.open_file(self.path, mode=self.mode,
                                        filters=_PYTABLES_DISABLE_COMPRESSION)

    def get(self, i, mmap=False):
        return self._handle.get_node('/', self._ITEM_FORMAT % i)[:]

    def keys(self):
        nodes = self._handle.list_nodes('/')
        for node in sorted(nodes, key=lambda x: _keynat(x.name)):
            match = self._ITEM_RE.match(node.name)
            if match:
                yield int(match.group(1))

    def set(self, i, x):
        if 'w' not in self.mode:
            raise IOError('Dataset not opened for writing')

        try:
            self._handle.create_carray('/', self._ITEM_FORMAT % i, obj=x)
        except tables.exceptions.NodeError:
            self._handle.remove_node('/', self._ITEM_FORMAT % i)
            self.set(i, x)

    @property
    def provenance(self):
        try:
            return self._handle.root._v_attrs['provenance']
        except KeyError:
            return 'No available provenance'

    def _write_provenance(self, previous=None, comments=''):
        p = self._build_provenance(previous=previous, comments=comments)
        self._handle.root._v_attrs['provenance'] = p

    def close(self):
        if hasattr(self, '_handle'):
            self._handle.close()

    def flush(self):
        self._handle.flush()

    def __del__(self):
        self.close()


class MDTrajDataset(_BaseDataset):
    _PROVENANCE_TEMPLATE = '''MDTraj dataset:
  path:\t\t{path}
  topology:\t{topology}
  stride:\t{stride}
  atom_indices\t{atom_indices}
'''

    def __init__(self, path, mode='r', topology=None, stride=1,
                 atom_indices=None, verbose=False):
        if mode != 'r':
            raise ValueError('mode must be "r"')
        self.path = path
        self.topology = topology
        self.stride = stride
        self.atom_indices = atom_indices
        self.verbose = verbose

        if isinstance(path, list):
            self.glob_matches = [expanduser(fn) for fn in path]
        else:
            self.glob_matches = sorted(glob.glob(expanduser(path)),
                                       key=_keynat)

        if topology is None:
            self._topology = None
        else:
            self._topology = _parse_topology(os.path.expanduser(topology))

    def get(self, i):
        if self.verbose:
            print('[MDTraj dataset] loading %s' % self.filename(i))

        if self._topology is None:
            t = md.load(self.filename(i), stride=self.stride,
                        atom_indices=self.atom_indices)
        else:
            t = md.load(self.filename(i), stride=self.stride,
                        atom_indices=self.atom_indices, top=self._topology)
        return t

    def filename(self, i):
        return self.glob_matches[i]

    def iterload(self, i, chunk):
        if self.verbose:
            print('[MDTraj dataset] iterloading %s' % self.filename(i))

        if self._topology is None:
            return md.iterload(
                self.filename(i), chunk=chunk, stride=self.stride,
                atom_indices=self.atom_indices)
        else:
            return md.iterload(
                self.filename(i), chunk=chunk, stride=self.stride,
                atom_indices=self.atom_indices, top=self._topology)

    def keys(self):
        return iter(range(len(self.glob_matches)))

    @property
    def provenance(self):
        return self._PROVENANCE_TEMPLATE.format(
            path=self.path, topology=self.topology,
            atom_indices=self.atom_indices, stride=self.stride)


def _tuples_and_ints(x):
    try:
        return tuple(x)
    except TypeError:
        return (x,)


def _dim_match(arr):
    if arr.ndim == 1:
        return arr[:, np.newaxis]
    return arr


def _keynat(string):
    """A natural sort helper function for sort() and sorted()
    without using regular expression.
    """
    r = []
    for c in string:
        if c.isdigit():
            if r and isinstance(r[-1], int):
                r[-1] = r[-1] * 10 + int(c)
            else:
                r.append(int(c))
        else:
            r.append(9 + ord(c))
    return r
