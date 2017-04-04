.. highlight:: bash
.. _tutorial:

Tutorial
========

#.
   .. figure:: _static/fspeptide.png
       :align: right

       Fs Peptide

   The command-line tool is ``msmb``. Check that it's installed correctly and
   read the help options::

    msmb -h

   This tutorial assumes you are using Python 3 on Linux. MSMBuilder works
   on many more platforms than this (python 2, mac, windows, ...) but you
   won't be able to follow the tutorial verbatim.

#. In your own research, you probably have a large molecular dynamics
   dataset that you wish to analyze. For this tutorial, we will perform a
   quick analysis on a simple system: the Fs peptide. Download the sample
   trajectories from `figshare <https://figshare.com/articles/Fs_MD_Trajectories/1030363>`_::

    mkdir fs-trajectories
    cd fs-trajectories/
    wget https://ndownloader.figshare.com/articles/1030363/versions/1
    # Figshare doesn't give you a good filename, but rest assured, this is a zip
    unzip 1
    rm 1

   Now you have some trajectory files to work with as well as the OpenMM
   scripts that were used to generate them::

    $ ls -C
    100-fs-peptide-400K.pdb  trajectory-5.xtc   trajectory-13.xtc  trajectory-21.xtc
    fs-peptide.pdb           trajectory-6.xtc   trajectory-14.xtc  trajectory-22.xtc
    post-process.py          trajectory-7.xtc   trajectory-15.xtc  trajectory-23.xtc
    simulate_fs.py           trajectory-8.xtc   trajectory-16.xtc  trajectory-24.xtc
    trajectory-1.xtc         trajectory-9.xtc   trajectory-17.xtc  trajectory-25.xtc
    trajectory-2.xtc         trajectory-10.xtc  trajectory-18.xtc  trajectory-26.xtc
    trajectory-3.xtc         trajectory-11.xtc  trajectory-19.xtc  trajectory-27.xtc
    trajectory-4.xtc         trajectory-12.xtc  trajectory-20.xtc  trajectory-28.xtc

   At this point, you should have a directory named ``fs-trajectories`` with 28 ``xtc``
   trajectories in it and a ``pdb`` topology file named ``fs-peptide.pdb``.
   Feel free to load one of these trajectories in VMD to get a sense of
   what they look like.


#. We'll use MSMBuilder to create a set of sample python scripts to sketch
   out our project::

    msmb TemplateProject -h # read the possible options
    msmb TemplateProject

   This will generate a hierarchy of small scripts that will handle the
   boilerplate and organization inherent in MSM analysis.::

    .
    ├── analysis
    │   ├── dihedrals
    │   │   ├── featurize-plot.py
    │   │   ├── featurize.py
    │   │   └── tica
    │   │       ├── cluster
    │   │       │   ├── cluster-plot.py
    │   │       │   ├── cluster.py
    │   │       │   ├── msm
    │   │       │   │   ├── microstate-plot.py
    │   │       │   │   ├── microstate.py
    │   │       │   │   ├── microstate-traj.py
    │   │       │   │   ├── timescales-plot.py
    │   │       │   │   └── timescales.py
    │   │       │   ├── sample-clusters-plot.py
    │   │       │   └── sample-clusters.py
    │   │       ├── tica-plot.py
    │   │       ├── tica.py
    │   │       ├── tica-sample-coordinate-plot.py
    │   │       └── tica-sample-coordinate.py
    │   ├── gather-metadata-plot.py
    │   ├── gather-metadata.py
    │   ├── top.pdb -> ../top.pdb
    │   └── trajs -> ../trajs
    ├── msmb-test-install.py
    └── README.md


   Each subsequent step in the MSM construction pipeline is a subdirectory.

   If you're comfortable with version control, the ``TemplateProject`` command
   also writes sample ``.gitignore`` files so you can put your scripts under
   version control (Optional!)::

    git init
    git add .
    git commit -m "Initial commit"

#. When using the ``TemplateProject`` command and/or when following this
   tutorial, we're a little picky about filenames (if you're a power user
   or are using msmbuilder as a library, you can use whatever filenames you
   want). In particular, we require a topology named ``top.pdb`` and trajectories
   in a folder named ``trajs``. Right now, your topology is named ``fs-trajectories/fs-peptide.pdb``
   and your trajectories are in a folder named ``fs-trajectories``. Let's use symlinks to
   help us out::

    ln -s fs-trajectories/fs-peptide.pdb top.pdb
    ln -s fs-trajectories/ trajs

   If your topology is a different format (e.g. Amber ``prmtop``), name it
   ``top.prmtop`` and use ``msmb TemplateProject --topology_fext=prmtop``
   on the command line.::

    $ ls
    total 16K
    drwxr-xr-x 5   4.0K Apr  4 15:48 analysis
    drwxr-xr-x 2   4.0K Apr  4 15:46 fs-trajectories
    -rw-r--r-- 1   1.1K Apr  4 15:48 README.md
    -rw-r--r-- 1   2.6K Apr  4 15:48 msmb-test-install.py
    lrwxrwxrwx 1     30 Apr  4 16:00 top.pdb -> fs-trajectories/fs-peptide.pdb
    lrwxrwxrwx 1     16 Apr  4 16:00 trajs -> fs-trajectories/



#. Begin our analysis::

    cd analysis/

#. We first generate a list of our trajectories and associated metadata.
   Examine the ``gather-metadata.py`` script. Note that it uses the ``xtc``
   filenames to generate an integer key for each trajectory. The script
   also extracts the length of each trajectory and stores the xtc filename::

    python gather-metadata.py
    python gather-metadata-plot.py


   .. figure:: _static/lengths-hist.png
       :align: right

       Sometimes you'll have many different length-ed trajectories and
       this histogram will be interesting. All of our trajectories are 500 ns
       though.

   **Important:** For your own projects, you *will* have to modify this short
   script (``gather-metadata.py``). Make sure you set ``step_ps`` to the length
   of each frame (in picoseconds). From the fs-peptide description on figshare,
   we can read that frames were saved every 50ps. Note that the example
   script uses ``NumberedRunsParser`` because our trajectories are all numbered
   with integers. You can use the ``msmbuilder.io.GenericParser`` for more
   control (it uses regular expressions).

   The plot script contains several example functions of computing statistics
   on your dataset including aggregate length. It will also generate an ``html``
   rendering of the table of metadata. **Exercise:** modify the
   script to genereate ``png`` images intead of ``pdf`` vector graphics for
   the plots. Editor's note: use ``pdf`` for preperation of manuscripts
   because you can
   `infinitely resize <https://en.wikipedia.org/wiki/Vector_graphics>`_ your plots.

#. We'll start reducing the dimensionality of our dataset by transforming
   the raw Cartesian coordinates into biophysical "features". We'll use
   dihedral angles. The templated project also includes subfolders ``landmarks``
   and ``rmsd`` for alternative approaches, but we'll ignore those for now::

    cd dihedrals/

#. Examine the ``featurize.py`` script. Note that it loops through our trajectories
   using the convenience function ``itertrajs`` (which only ever holds one
   trajectory in RAM) and calls ``DihedralFeaturizer.partial_transform()``
   on each. Read more about :ref:`featurizers<featurization>` and MSMBuilder
   :ref:`API patterns<apipatterns>`. Run the scripts::

    python featurize.py
    python featurize-plot.py

   The plots will show you a box and whisker plot of each feature value. This
   is not very useful, but we wanted to make sure you can plot something
   for each step.
   **Exercise:** include `chi1 and chi2 <https://swissmodel.expasy.org/course/text/chapter3.htm>`_
   angles in addition to the default `phi and psi <https://en.wikipedia.org/wiki/Ramachandran_plot>`_
   angles.

#. Dihedrals are too numerous to be interpretable. We can use :ref:`tica<decomposition>`
   to learn a small number of "kinetic coordinates" from our data::

    cd tica/

#. Examine ``tica.py``. Note that it loads the feature trajectories, learns
   a model from them by calling ``fit()`` and then transforms the feature trajectories
   into "tica trajectories" by calling ``partial_transform()``
   on each (see :ref:`api patterns<apipatterns>`). The MSMBuilder API *does not*
   keep track of units. Our data was saved every 50 ps (Editor's note: this is
   way too frequent for a "real" simulation). The template script for learning
   our tica model sets the ``lag_time`` parameter to ``10``. This means 10 steps
   in our data. This translates to 500 ps here. Let's use something a little
   longer like 5 ns (= 100 steps). Edit the ``lag_time`` parameter to 100 and
   learn the model::

    vim tica.py # edit lag_time -> 100
    python tica.py
    python tica-plot.py

#.
   .. figure:: _static/tica-heatmap.png
       :align: right

       tICA heatmaps provide a convenient 2d projection of your data
       onto which you can overlay more interesting info.

   The tICA plotting script makes a 2d histogram of our data. Note the apparent
   free energy well on the left of the figure. We might suspect that this is
   the folded state and the x-axis is an unfolding coordinate. We'll use
   this tica heatmap as a background for our further plots. tICA is extremely
   useful at taking hundreds of dihedral angles (for example) and distilling it
   into a handful of coordinates that we can plot.

#. We can sample configurations along a tIC to inspect what that tIC "means".
   Another common strategy for interpreting tICs is to inspect prominent
   (most non-zero) coefficients corresponding to particular features (dihedrals).
   A common tactic is to color residues based on their tIC loading. Example
   scripts to set up VMD for this will be included in a later release. Here,
   we simply draw configurations along a tIC direction::

    python tica-sample-coordinate.py
    python tica-sample-coordinate-plot.py

#.
   .. figure:: _static/tica-movie.gif
        :align: right

        The first tIC is roughly a folding coordinate.

   This produces a trajectory of conformations, saved as ``tica-dimension-0.xtc``.
   **Exercise:** Save the conformations as a ``dcd`` trajectory instead.
   You can load this trajectory in VMD and inspect the particular tIC::

    vmd top.pdb tica-dimension-0.xtc

   Align the structures and apply some "smoothing".
   **Exercise:** Sample the second tIC. Note that it probably isn't an
   interesting coordinate in this case.

#. We can group conformations that interconvert rapidly by using off-the-shelf
   clustering algorithms on our kinetic coordinates (tICs)::

    cd cluster/

   By default, we generate 500 clusters using a form of KMeans. Read more
   about :ref:`clustering<cluster>`.
   **Exercise:** try a different number of clusters or a different
   clustering algorithm. Run the clustering scripts::

    python cluster.py
    python cluster-plot.py

   Note that the tIC heatmap provides a convenient space onto which we project
   our cluster centers.

#. With our states defined, we count the transitions between them. An :ref:`MSM<msm>`
   is simply states and rates. First we make a "microstate" MSM consisting
   of many, small states::

    cd msm/

#.

   .. figure:: _static/msm-microstates.png
        :align: right

        The microstate centers are shown as circles on the tIC heatmap.
        They are sized according to state population. They are colored
        according to the first dynamical eigenvector. The slowest
        processes is a transition from red states to blue.

   The MSM lag-time is a parameter that cannot be optimized using :ref:`gmrq<gmrq>`.
   You can use the ``timescales.py`` script to check how the model timescales
   would react to changing the lag-time. We'll just use a lag-time of 5 ns.
   Remember from above that we have to keep track of units. 5 ns is 100 steps.
   Edit ``microstate.py`` and set ``lag_time = 100``::

    vim microstate.py # edit lag_time
    python microstate.py
    python microstate-plot.py

#. Generate a sample trajectory from the MSM::

    python microstate-traj.py

   By default, each frame will be 1 lag-time unit. Here, that is 5 ns.
   **Exercise:** Use the ``n_steps`` and ``stride`` parameter to sample
   a 200 frame movie with 50 ns steps.
   You can load the trajectory in VMD and watch the Fs-peptide stochastically
   fold and unfold::

    vmd top.pdb msm-traj.xtc
