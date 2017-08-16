Installation instructions for PySheaf using a conda environment
===============================================================

PySheaf is written for python 2.7. To install and test PySheaf we recommend using a virtual environment. The instructions below will walk you through the following steps:

1. Download and install miniconda (or anaconda if you prefer)
2. Create a :py:mod:`conda` environment
3. Download and install the PySheaf repository
4. Run unit tests

Download and install miniconda:
-------------------------------

*If you already have an anaconda or miniconda distribution you may skip this part.*

Anaconda and Miniconda are source distributions of Python used by the mathematics and scientific communities that simplify package management and version control.  Both use the :py:mod:`conda` package management system and both provide environments to isolate python packages and their dependencies. Anaconda includes :py:mod:`conda`, :py:mod:`conda-build`, :py:mod:`python` and over 150 python packages. Miniconda is a smaller version that only includes :py:mod:`conda`, :py:mod:`python`, and the packages they depend on. For our purposes we will install miniconda to create a minimal environment adequate for testing and using pysheaf.
 
To download miniconda go to `<http://www.continuum.io/downloads>`_ and chose your operating system. Both Anaconda and Miniconda install multiple Python versions. The python version you choose will be used by default but in a virtual environment we will specify the version. Installation instructions are described on the site.

You can test that your installation was successful by entering::

  conda list

on the command line. To view the python versions available enter::

  conda search python

Create a :py:mod:`conda` environment
------------------------------------

The :py:mod:`conda` environment will expose only the version of python and the python packages that you choose. On the command line enter::

  conda create --name pysheafenv python=2.7

Conda environments are stored within your miniconda (anaconda) directory in the envs folder. To see a list of your environments enter::

  conda info -e

To remove the pysheafenv environment enter::

  conda remove --name pysheafenv --all

To use the environment you must activate it. On Windows::

  activate pysheafenv

On OS X or Linux::

  source activate pysheafenv

This will add (pysheafenv) to your command prompt.

To deactivate the environment on Windows::

  deactivate

On OS X or Linux::

  source deactivate

Download and install the PySheaf repository
-------------------------------------------

Activate pysheafenv as described above. To see what packages were installed enter::

  conda list

Change to the directory where you wish PySheaf to be installed.

The :py:mod:`pip` package was installed automatically with python 2.7. We will need the :py:mod:`pip` package manager to install PySheaf directly from github.

Because we wish to test the installation using the :py:mod:`unittest` module we will use the -e for editable option. This will also permit us to review the code from the current directory::

  pip install -e git+https://github.com/kb1dds/pysheaf.git#egg=pysheaf

Depending on your system (Windows in particular) this may or may not install all of the dependencies. If you get an error message that something is missing you can install the missing dependencies using :py:mod:`conda`.  PySheaf requires :py:mod:`numpy`, :py:mod:`networkx`, :py:mod:`scipy`, and :py:mod:`deap`.
To install these using :py:mod:`conda`::

  conda install numpy networkx scipy
  conda install -c conda-forge deap

**Note:** If you only wish to install the PySheaf package for purposes of imports and do not wish to run the tests or review the source code then use the following for pip install::

  pip install git+https://github.com/kb1dds/pysheaf.git

If you look at the installed packages for this environment::

  conda list

you will see PySheaf has been added to the packages. (If your default python installation is python 2.7 and you do not wish to use a virtual environment you may try this pip command directly without the use of :py:mod:`conda`. Success will depend on your system's installed packages and configuration.)

Run unit tests
--------------

Assuming you have installed the package into PySheaf you will now change to that directory and run the tests::

  cd src/pysheaf
  python -m unittest discover

If everything is working you will get a bunch of warnings and a line that looks something like::

  Ran 57 tests in 0.070s  
  OK
