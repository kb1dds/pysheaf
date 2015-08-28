# pysheaf
Python Cellular Sheaf Library

This is part of the [AU Sheaves Simplex Activity](http://www.american.edu/cas/darpasheaves/).

The pip install follows [this example](https://pip.pypa.io/en/latest/reference/pip_install.html#git).

    pip install -e git+https://github.com/phreed/pysheaf.git@master#egg=pysheaf

This works well with a prior installation of Miniconda.
Then an appropriate environment is constructed into which pysheaf is added as a package.

    conda create -n sheaves jupyter numpy scipy matplotlib networkx
    activate sheaves
    pip install -e git+https://github.com/phreed/pysheaf.git@master#egg=pysheaf

