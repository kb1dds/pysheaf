PySheaf: Sheaf-theoretic toolbox
================================

This repository consists of Python 3.6 libraries for manipulating cell complexes and sheaves of sets or vector spaces on cell complexes [1]_, [2]_, [3]_, [4]_, [5]_, [6]_

Documentation:
--------------

Full (very out-of-date) documentation for PySheaf verson v0.xx is at `<http://kb1dds.github.io/pysheaf/>`_

Right now, the best strategy is to look at the example code!

Download:
---------

You can install by cloning this repo (there is no PyPI package).  For Linux, you can do::

  $ git clone https://github.com/kb1dds/pysheaf.git
  $ pip install pysheaf

See the `documentation <http://kb1dds.github.io/pysheaf/install.html>`_ for full details!  (Though it's a bit out of date...)

Usage:
------

The general plan of usage is

1. First (usually on paper!) lay out the cell complex that will serve as the base for your sheaf.  *Give each cell a unique label.*  

2. Determine all of the stalks over each cell, and the restriction maps.  Restriction maps can be a mixture of `numpy` matrices or arbitrary single-input Python function objects.
   
3. Construct a `Sheaf` instance and add each of your cells as `Cell` instances with the `Sheaf.AddCell` method.  Make sure to use your unique label for each `Cell`, because that is how PySheaf identifies them! Once you've done that, create each restriction as a `Coface` instance and add it to the sheaf using the `Sheaf.AddCoface` method.  The `Sheaf.AddCoface` method will connect the listed `Cell`s based on their labels.  `Cell`s and `Coface`s can be added later if you want, and they can be added in any order provided any `Coface` refers to `Cell`s that already exist.

4. Install some data into the sheaf by way of an `Assignment` to some of the `Cell`s.  

5. Analyze the sheaf and its data:
  a. You can compute consistency radius with `Sheaf.ComputeConsistencyRadius()`
  b. You can improve the consistency radius by extending or altering the values of the assignment with `Sheaf.FuseAssignment()`.  This will only alter Cells whose `Cell.mOptimizationCell` attribute is `True`.  You can also change the optimization algorithm if you want.
  c. You can find all star open sets whose local consistency is less than a desired bound using `Sheaf.CellIndexesLessThanConsistencyThreshold()`.

Have a look at the example code for some ideas!  

This code is under active development, so not everything works as it should.  If you find anything that you can correct, feel free to send me suggestions!

| Thanks!
| Michael Robinson
| American University
| kb1dds@gmail.com
| michaelr@american.edu

.. [1] http://www.drmichaelrobinson.net/sheaftutorial/index.html

.. [2] https://www.youtube.com/user/drmichaelrobinson

.. [3] Cliff Joslyn, Emilie Hogan, Michael Robinson, "Towards a topological framework for integrating semantic information sources," Semantic Technologies for Intelligence, Defense, and Security (STIDS), 2014. http://ceur-ws.org/Vol-1304/STIDS2014_P2_JoslynEtAl.pdf

.. [4] Michael Robinson, "Sheaves are the canonical datastructure for information integration," Information Fusion, 36 (2017), pp. 208-224. (preprint version is http://arxiv.org/abs/1603.01446)

.. [5] Michael Robinson, "Sheaf and cosheaf methods for analyzing multi-model systems," http://arXiv.org/abs/1604.04647

.. [6] Michael Robinson, *Topological Signal Processing*, Springer, 2014.
