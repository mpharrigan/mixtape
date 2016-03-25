from __future__ import print_function, absolute_import

import sys
from os.path import join as pjoin

import mdtraj
import numpy as np
from setuptools import setup, Extension, find_packages

try:
    sys.dont_write_bytecode = True
    sys.path.insert(0, '.')
    from basesetup import (write_version_py,
                           CompilerDetection,
                           check_dependencies)
finally:
    sys.dont_write_bytecode = False

if '--debug' in sys.argv:
    sys.argv.remove('--debug')
    DEBUG = True
else:
    DEBUG = False
if '--disable-openmp' in sys.argv:
    sys.argv.remove('--disable-openmp')
    DISABLE_OPENMP = True
else:
    DISABLE_OPENMP = False

try:
    import Cython
    from Cython.Distutils import build_ext

    if Cython.__version__ < '0.18':
        raise ImportError()
except ImportError:
    print('Cython version 0.18 or later is required. '
          'Try "conda install cython"')
    sys.exit(1)

mdtraj_capi = mdtraj.capi()

# #########################
VERSION = '3.4.0.dev0'
ISRELEASED = False
__version__ = VERSION
# #########################

CLASSIFIERS = """\
Intended Audience :: Science/Research
Intended Audience :: Developers
License :: OSI Approved :: GNU Lesser General Public License v2 or later (LGPLv2+)
Programming Language :: C++
Programming Language :: Python
Development Status :: 5 - Production/Stable
Topic :: Software Development
Topic :: Scientific/Engineering
Operating System :: POSIX
Operating System :: Unix
Operating System :: MacOS
Programming Language :: Python :: 2
Programming Language :: Python :: 2.7
Programming Language :: Python :: 3
Programming Language :: Python :: 3.4
Programming Language :: Python :: 3.5
"""

if any(cmd in sys.argv for cmd in ('install', 'build', 'develop')):
    check_dependencies((
        ('numpy',),
        ('scipy',),
        ('pandas',),
        ('six',),
        ('mdtraj',),
        ('sklearn', 'scikit-learn'),
        ('numpydoc',),
        ('tables', 'pytables'),
    ))

# Where to find extensions
MSMDIR = 'msmbuilder/msm/'
HMMDIR = 'msmbuilder/hmm/'
CLUSTERDIR = 'msmbuilder/cluster/'

# Configure compiler with debug and openmp settings
compiler = CompilerDetection(DISABLE_OPENMP)
with open('msmbuilder/src/config.pxi', 'w') as f:
    f.write('\n'.join([
        'DEF DEBUG = {debug}'.format(debug=DEBUG),
        'DEF OPENMP = {openmp}'.format(openmp=compiler.openmp_enabled),
    ]))

extensions = []
extensions.append(
    Extension('msmbuilder.example_datasets._muller',
              sources=['msmbuilder/example_datasets/_muller.pyx'],
              include_dirs=[np.get_include()]))

extensions.append(
    Extension('msmbuilder.msm._markovstatemodel',
              sources=[pjoin(MSMDIR, '_markovstatemodel.pyx'),
                       pjoin(MSMDIR, 'src/transmat_mle_prinz.c')],
              include_dirs=[pjoin(MSMDIR, 'src'), np.get_include()]))

extensions.append(
    Extension('msmbuilder.tests.test_cyblas',
              sources=['msmbuilder/tests/test_cyblas.pyx'],
              include_dirs=['msmbuilder/src', np.get_include()]))

extensions.append(
    Extension('msmbuilder.msm._ratematrix',
              sources=[pjoin(MSMDIR, '_ratematrix.pyx')],
              language='c++',
              extra_compile_args=compiler.compiler_args_openmp,
              libraries=compiler.compiler_libraries_openmp,
              include_dirs=['msmbuilder/src', np.get_include()]))

extensions.append(
    Extension('msmbuilder.decomposition._speigh',
              sources=['msmbuilder/decomposition/_speigh.pyx'],
              language='c++',
              extra_compile_args=compiler.compiler_args_openmp,
              libraries=compiler.compiler_libraries_openmp,
              include_dirs=['msmbuilder/src', np.get_include()]))

extensions.append(
    Extension('msmbuilder.msm._metzner_mcmc_fast',
              sources=[pjoin(MSMDIR, '_metzner_mcmc_fast.pyx'),
                       pjoin(MSMDIR, 'src/metzner_mcmc.c')],
              libraries=compiler.compiler_libraries_openmp,
              extra_compile_args=compiler.compiler_args_openmp,
              include_dirs=[pjoin(MSMDIR, 'src'), np.get_include()]))

extensions.append(
    Extension('msmbuilder.libdistance',
              language='c++',
              sources=['msmbuilder/libdistance/libdistance.pyx'],
              # msvc needs to be told "libtheobald", gcc wants just "theobald"
              libraries=['{}theobald'.format('lib' if compiler.msvc else '')],
              include_dirs=["msmbuilder/libdistance/src",
                            mdtraj_capi['include_dir'],
                            np.get_include()],
              library_dirs=[mdtraj_capi['lib_dir']],
              ))

extensions.append(
    Extension('msmbuilder.cluster._kmedoids',
              language='c++',
              sources=[pjoin(CLUSTERDIR, '_kmedoids.pyx'),
                       pjoin(CLUSTERDIR, 'src', 'kmedoids.cc')],
              include_dirs=[np.get_include()]))

extensions.append(
    Extension('msmbuilder.hmm.gaussian',
              language='c++',
              sources=[pjoin(HMMDIR, 'gaussian.pyx'),
                       pjoin(HMMDIR, 'src/GaussianHMMFitter.cpp')],
              libraries=compiler.compiler_libraries_openmp,
              extra_compile_args=(compiler.compiler_args_sse3
                                  + compiler.compiler_args_openmp),
              include_dirs=[np.get_include(),
                            HMMDIR,
                            pjoin(HMMDIR, 'src/include/'),
                            pjoin(HMMDIR, 'src/')]))

extensions.append(
    Extension('msmbuilder.hmm.vonmises',
              language='c++',
              sources=[pjoin(HMMDIR, 'vonmises.pyx'),
                       pjoin(HMMDIR, 'src/VonMisesHMMFitter.cpp'),
                       pjoin(HMMDIR, 'cephes/i0.c'),
                       pjoin(HMMDIR, 'cephes/chbevl.c')],
              libraries=compiler.compiler_libraries_openmp,
              extra_compile_args=(compiler.compiler_args_sse3
                                  + compiler.compiler_args_openmp),
              include_dirs=[np.get_include(),
                            HMMDIR,
                            pjoin(HMMDIR, 'src/include/'),
                            pjoin(HMMDIR, 'src/'),
                            pjoin(HMMDIR, 'cephes/')]))

write_version_py(VERSION, ISRELEASED, filename='msmbuilder/version.py')

setup(
    name='msmbuilder',
    author='Robert McGibbon',
    author_email='rmcgibbo@gmail.com',
    description="MSMBuilder: Statistical models for biomolecular dynamics",
    version=VERSION,
    url='https://github.com/msmbuilder/msmbuilder',
    platforms=['Linux', 'Mac OS-X', 'Unix'],
    classifiers=CLASSIFIERS.splitlines(),
    packages=find_packages('msmbuilder'),
    package_data={'msmbuilder.tests': ['workflows/*']},
    entry_points={
        'console_scripts':
            ['msmb = msmbuilder.scripts.msmb:main']
    },
    zip_safe=False,
    ext_modules=extensions,
    cmdclass={'build_ext': build_ext}
)
