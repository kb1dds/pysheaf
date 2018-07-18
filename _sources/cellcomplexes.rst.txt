The CellComplex type
====================

The :py:class:`CellComplex` class consists of a list of :py:class:`Cell` instances and methods that manipulate the complex as a whole.  It is also the base class for the :py:class:`Sheaf` class.  The indices into the list of :py:class:`Cell` instances are used throughout PySheaf, and are the usual way to refer to individual :py:class:`Cell` instances when they are in context in a :py:class:`CellComplex`.  Because the indices are necessary to construct the :py:class:`Cofaces` as well, it is usually necessary to determine the necessary cells ahead of time, and then build the :py:class:`CellComplex` instance all at once. 

.. py:class:: CellComplex

   .. py:attribute:: cells

      List of :py:class:`Cell` instances that form this cell.  Indices into this list are used throughout PySheaf, and they are generally expected to be static once created.

   .. py:method:: homology(k,subcomplex=None,compactSupport=False,tol=1e-5)

      Compute the degree `k` homology of the :py:class:`CellComplex`.  If you want relative homology, the `subcomplex` field specifies a list of indices into :py:class:`CellComplex.cells` for the relative subcomplex.  If you want compactly supported homology (if you don't know what that means, you don't) then set `compactSupport=True`.  The `tol` argument sets the tolerance below which a singular value is said to be zero, and thus is to be considered part of the kernel. This returns a :py:class:`numpy.ndarray` whose columns are the generators for homology.

   .. py:method:: boundary(k,subcomplex=None,compactSupport=False)

      Compute the degree `k` boundary map of the :py:class:`CellComplex`, returning it as a :py:class:`numpy.ndarray`.   If you want relative homology, the `subcomplex` field specifies a list of indices into :py:class:`CellComplex.cells` for the relative subcomplex.  If you want compactly supported homology (if you don't know what that means, you don't) then set `compactSupport=True`.  

   .. py:method:: components(cells=[])
		  
      Compute connected components of the :py:class:`CellComplex`.   The optional argument `cells` specifies list of permissible indices into :py:attr:`CellComplex.cells`.  Returns a list of lists of indices into :py:attr:`CellComplex.cells`.
   
   .. py:method:: star(cells)

      Compute the star of a list of `cells` (specified as a list of indices into :py:attr:`CellComplex.cells`) in the topology of the :py:class:`CellComplex`.  Returns a list of indices into :py:attr:`CellComplex.cells`.

   .. py:method:: closure(cells)

      Compute the closure of a list of `cells` (specified as a list of indices into :py:attr:`CellComplex.cells`) in the topology of the :py:class:`CellComplex`.  Returns a list of indices into :py:attr:`CellComplex.cells`.

   .. py:method:: skeleton(k)

      Compute the dimension `k` skeleton of the :py:class:`CellComplex`.  Returns a list of indices into :py:attr:`CellComplex.cells`.

   .. py:method::  cofaces(c,cells=[])

      Return a generator that iterates over over :py:class:`Coface` instances (of *all* dimensions) for a cell whose index in :py:attr:`CellComplex.cells` is `c`.  The optional argument specifies a list of indices into :py:attr:`CellComplex.cells` that are permitted to be traversed.

      .. warning :: duplicate :py:class:`Coface` instances are possible!

Constructing :py:class:`CellComplex` instances
----------------------------------------------

A :py:class:`Cell` represents a single topological disk that is present in a given :py:class:`CellComplex`.  Mostly, it contains references to other cells by their respective indices into :py:attr:`CellComplex.cells`.

.. py:class:: Cell

   Base class representing a topological disk of a definite dimension.

   .. py:attribute:: dimension

      The dimension of the disk that this :py:class:`Cell` represents.  The actual points of the disk are *not* represented, merely its dimension.  (Note: this is *not* the dimension of the stalk over the cell in a :py:class:`SheafCell`)

   .. py:attribute:: compactClosure

      Flag that specifies if the topological closure of the :py:class:`Cell` in the :py:class:`CellComplex` is compact.  Usually this should be `True`, as only those cells with compact closure are included in a homology calculation.  Roughly speaking, those cells that have "missing" boundaries do not have compact closure.

   .. py:attribute:: name

      An optional name for the :py:class:`Cell`, which is generally not used by PySheaf.

   .. py:attribute:: cofaces

      A list of :py:class:`Coface` instances, specifying each coface of this cell.  It is assumed that this coface points to a strictly higher-dimensional cell, and you will encounter endless loops if this assumption is violated.  It is *not* assumed that the cofaces are all *one* dimension higher, though.  It is not necessary to specify a transitive closure -- all cofaces -- as this can be determined by the containing :py:class:`CellComplex` as needed using :py:meth:`CellComplex.cofaces()`.

The :py:class:`Coface` class specifies a single coface relation, in the context of a :py:class:`CellComplex`.

.. py:class:: Coface
   
   Class representing a coface relation between two :py:class:`Cell` instances.  The lower-dimension cell is implied to be the one holding this instance as its :py:attr:`Cell.cofaces` attribute, so this class *only* refers to the higher-dimension cell.

   .. py:attribute:: index

      The index of the higher-dimension cell in the containing :py:class:`CellComplex`.  

   .. py:attribute:: orientation

      The orientation of this coface relation, usually either +1 or -1.

:py:class:`CellComplex` instances are best built all at once.  So for instance, a cell complex consisting of four vertices, named `A`, `B`, `C`, `D`, five edges `AB`, `AC`, `BC`, `BD`, `CD`, and one triangle `ABC` is constructed thusly::

      pysheaf.CellComplex([pysheaf.Cell(dimension=0,
                                        compactClosure=True,
					name='A',
					cofaces=[pysheaf.Coface(index=4,orientation=1),   # Index 4 = 'AB'
                                                 pysheaf.Coface(index=5,orientation=1)]), # Index 5 = 'AC'
                           pysheaf.Cell(dimension=0,
                                        compactClosure=True,
					name='B',
					cofaces=[pysheaf.Coface(index=4,orientation=-1),  # Index 4 = 'AB'
				                 pysheaf.Coface(index=6,orientation=1),   # Index 6 = 'BC'
						 pysheaf.Coface(index=7,orientation=1)]), # Index 7 = 'BD'
 			   pysheaf.Cell(dimension=0,
			                compactClosure=True,
                                        name='C',
   					cofaces=[pysheaf.Coface(index=5,orientation=-1),  # Index 5 = 'AC'
			                         pysheaf.Coface(index=6,orientation=-1),  # Index 6 = 'BC'
						 pysheaf.Coface(index=8,orientation=1)]), # Index 8 = 'CD'
			   pysheaf.Cell(dimension=0,
                                        compactClosure=True,
					name='D',
					cofaces=[pysheaf.Coface(index=7,orientation=-1),  # Index 7 = 'BD'
				                 pysheaf.Coface(index=8,orientation=-1)]),# Index 4 = 'CD'
			   pysheaf.Cell(dimension=1,
			                compactClosure=True,
					name='AB',
					cofaces=[pysheaf.Coface(index=9,orientation=1)]), # Index 9 = 'ABC'
			   pysheaf.Cell(dimension=1,
                                        compactClosure=True,
					name='AC',
					cofaces=[pysheaf.Coface(index=9,orientation=-1)]),# Index 9 = 'ABC' 
			   pysheaf.Cell(dimension=1,
                                        compactClosure=True,
					name='BC',
					cofaces=[pysheaf.Coface(index=9,orientation=1)]), # Index 9 = 'ABC'
			   pysheaf.Cell(dimension=1,
                                        compactClosure=True,
					name='BD',
					cofaces=[]),
			   pysheaf.Cell(dimension=1,
                                        compactClosure=True,
					name='CD',
					cofaces=[]),
			   pysheaf.Cell(dimension=2,
                                        compactClosure=True,
                                        name='ABC',
 					cofaces=[])])

Subclasses of :py:class:`CellComplex`
-------------------------------------

Since certain kinds of cell complex are specified combinatorially, PySheaf provides subclasses of :py:class:`CellComplex` whose constructors are a bit less verbose.

.. py:class:: AbstractSimplicialComplex(CellComplex)

   An abstract simplicial complex is defined as a list of lists.  Each list specifies a top dimensional simplex (a *toplex*).

   .. py:method:: __init__(toplexes,maxdim=None)

      It is only necessary to pass the constructor a generating set of `toplexes`, as all other simplices will be generated as :py:class:`Cell` instances as appropriate.  Because of this, an :py:class:`AbstractSimplicialComplex` should not be constructed with high-dimensional simplices!  To avoid problems, you may need to set the `maxdim` to be the largest dimension you want to have constructed.  Simplices are sorted from greatest dimension to least in the resulting :py:attr:`AbstractSimplicialComplex.cells` list.

.. py:class:: DirectedGraph(CellComplex)

   A directed graph consists of a list of ordered pairs of vertices.  Strictly speaking, this is a directed, weighted *multi* graph, since duplicate edges are allowed.
   
   .. py:method:: __init__(graph,vertex_capacity=-1)
		  
        Create a cell complex from a directed graph, where `graph` is a list of pairs (src,dest) or triples (src,dest,capacity) of numbers representing vertices.  The optional `vertex_capacity` sets all the weights to the same value.  Orientations of the resulting :py:class:`Coface` instances are taken from the ordering of the vertices in the tuples in `graph`.
	
        The vertex labeled `None` in any of these tuples is an external connection.  In the resulting :py:class:`DirectedGraph.cells`, the cells are indexed as follows:
        1. First all of the edges (in the order given),
        2. then all vertices (in the order they are given; not by their numerical values)
