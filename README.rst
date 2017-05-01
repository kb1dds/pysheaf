PySheaf: Sheaf-theoretic toolbox
================================

This repository consists of Python 2.7 libraries for manipulating cell complexes and data structures on them:

1. Local and relative homology of abstract simplicial complexes [0]_

2. Sheaves of sets or vector spaces on cell complexes [1]_, [2]_, [3]_, [4]_, [5]_, [6]_

Documentation:
--------------

Full documentation is at `<http://kb1dds.github.io/pysheaf/>`_

Download:
---------

You can install by cloning this repo (there is no PyPI package).  For Linux, you can do::

  $ git clone https://github.com/kb1dds/pysheaf.git
  $ pip install pysheaf

See the `documentation <http://kb1dds.github.io/pysheaf/install.html>`_ for full details!

Usage:
------

The general plan of usage is

1. First (usually on paper!) lay out the cell complex that will serve as the base for your sheaf.  *Label each cell with a unique index, starting from zero.*  

2. Determine all of the stalks over each cell, and the restriction maps from lower dimension to higher.  Restriction maps can be a mixture of `numpy` matrices or instances of `LinearMorphism`, `SetMorphism`, or `SheafMorphism`.
   
3. Construct a `Sheaf` instance from a list of `SheafCell` instances, all at once using the information from the previous two steps.
   
4. Analyze the resulting sheaf:
   
   a. If you have data, you can build `Section` instances whose list of `SectionCell` instances refer to your sheaf.
      
   b. You can compute cohomology of the `Sheaf` instance as well.

Have a look at the `pysheaf/tests` folder for some examples!  

Notes:
------
The pysheaf.py module should be able to duplicate the functionality of simplicialHomology.py, albeit somewhat more slowly because the base space constructions are more general and less efficient.

Finally, this code is under active development, so not everything works as it should.  Stay away from `Sheaf.localSectional()` and `Cell.localPairComplex()` as these don't currently work correctly!  If you find anything that you can correct, feel free to send me suggestions!

| Thanks!
| Michael Robinson
| American University
| kb1dds@gmail.com
| michaelr@american.edu

.. [0] Alan Hatcher, "Algebraic Topology", Cambridge, 2002, https://www.math.cornell.edu/~hatcher/AT/ATpage.html

.. [1] http://www.drmichaelrobinson.net/sheaftutorial/index.html

.. [2] https://www.youtube.com/user/drmichaelrobinson

.. [3] Cliff Joslyn, Emilie Hogan, Michael Robinson, "Towards a topological framework for integrating semantic information sources," Semantic Technologies for Intelligence, Defense, and Security (STIDS), 2014. http://ceur-ws.org/Vol-1304/STIDS2014_P2_JoslynEtAl.pdf

.. [4] Michael Robinson, "Sheaves are the canonical datastructure for information integration," Information Fusion, 36 (2017), pp. 208-224. (preprint version is http://arxiv.org/abs/1603.01446)

.. [5] Michael Robinson, "Sheaf and cosheaf methods for analyzing multi-model systems," http://arXiv.org/abs/1604.04647

.. [6] Michael Robinson, *Topological Signal Processing*, Springer, 2014.
