The :py:class:`Sheaf` type
==========================

The :py:class:`Sheaf` class derives from :py:class:`CellComplex` to describe its base space.  The additional attributes and methods of :py:class:`Sheaf` describe the "vertical fiber" structure of the sheaf.

.. py:class:: Sheaf(CellComplex)

   The base sheaf class in PySheaf consists of a list :py:attr:`cell` of :py:class:`SheafCell` instances, describing the cells and the stalks over them.  Restriction maps are built into :py:class:`SheafCoface` instances that are stored with each :py:class:`SheafCell`.

   .. py:method:: coboundary(k,compactSupport=False)

      Compute the degree `k` coboundary map of the :py:class:`Sheaf`, returning it as a :py:class:`numpy.ndarray`.   If you want compactly supported cohomology (if you don't know what that means, you don't) then set `compactSupport=True`.

   .. py:method:: cohomology(k,compactSupport=False)

      Compute the degree `k` sheaf cohomology for the :py:class:`Sheaf`.   If you want compactly supported cohomology (if you don't know what that means, you don't) then set `compactSupport=True`.  This returns a :py:class:`numpy.ndarray` whose columns are the generators for cohomology.
		  
   .. py:method:: cobetti(k)

      The dimension of the degree `k` sheaf cohomology space.

      .. warning :: This is not the dimension of the :py:meth:`CellComplex.homology()`, which is inherited into :py:class:`Sheaf`!

.. py:class:: SheafCell(Cell)

   .. py:attribute:: stalkDim

      The dimension of the stalk over this cell.

      .. warning :: This has nothing to do with :py:attr:`SheafCell.dimension` inherited from :py:attr:`Cell`.

   .. py:attribute:: metric

      The metric on the stalk, used to measure distances between values in that stalk.  This must by a function object that takes two arguments, each of which is a :py:class:`numpy.ndarray` and produces a numerical value.  By default, it is the Euclidean distance given by :py:func:`numpy.linalg.norm()`.
		     
   .. py:attribute:: cofaces

      Like :py:class:`Cell`, this should be a list of coface relations, but unlike :py:class:`Cell`, they must be :py:class:`SheafCoface` instances!

   .. py:attribute:: bounds
		     
      The portion of the stalk in which legal values live.  This is used by the optimizer in :py:meth:`Sheaf.fuseAssignment` to set bounds on which values are used.  This should either be `None` (in which the stalk represents the entire vector space) or a list of `(min,max)` pairs of length :py:attr:`SheafCell.stalkDim`.  Use `None` for any bound you wish to ignore. 

.. py:class:: SheafCoface(Coface)

   A coface relation in a sheaf on a cell complex is built just like a :py:class:`Coface` in a :py:class:`CellComplex`, but with the addition of a :py:attr:`restriction`.
	      
   .. py:attribute:: restriction

      A restriction in a sheaf generally needs to be a morphism in some category.   It must be an object that supports *composition*, namely a class that implements a multiplication operator.  In most examples, this is one of :py:class:`SetMorphism` (for functions between sets), :py:class:`LinearMorphism` (for linear maps), or a :py:class:`SheafMorphism` (for a :py:class:`Sheaf` of sheaves).  Note that if you construct a :py:class:`SheafCoface` by passing a :py:class:`numpy.ndarray`, PySheaf will construct a :py:class:`LinearMorphism` restriction automatically.

Morphisms: :py:class:`SetMorphism`, :py:class:`LinearMorphism`, and :py:class:`SheafMorphism`
---------------------------------------------------------------------------------------------

Since cellular sheaves are special functors from the face category of a cell complex to some other *data category*, PySheaf supports three kinds of *data categories*:

1. Sets,
2. Finite dimensional vector spaces (via :py:mod:`numpy`), and
3. Sheaves (of some of other type).

The restrictions in a given :py:class:`SheafCoface` instance are therefore of a corresponding morphism class:

1. :py:class:`SetMorphism`,
2. :py:class:`LinearMorphism`, and
3. :py:class:`SheafMorphism`.

The simplest of these is :py:class:`SetMorphism`.

.. py:class:: SetMorphism
	      
   This represents a *set* morphism, otherwise known as a *function* between sets.  This is implemented by a single attribute

   .. py:attribute:: fcn

   which is a function object taking one argument.

   :py:class:`SetMorphism` objects support a multiplication operator, which composes their respective :py:attr:`fcn` attributes to form a new function object.  They also support call semantics, so you can simply call a :py:class:`SetMorphism` object as a function to access its :py:attr:`fcn` attribute.

Namely, if you say::

  foo = pysheaf.SetMorphism( lambda x : x**2 )
  bar = pysheaf.SetMorphism( lambda y : 3*y )

then::

  foo(3)

returns 9 and::

  baaz = foo * bar
  baaz(1)

is also 9.

A :py:class:`Sheaf` with only :py:class:`SetMorphism` restrictions does not allow you to compute :py:meth:`Sheaf.cohomology()`.  For that, you need linearity, which is implemented by the following subclass of :py:class:`SetMorphism`.

.. py:class:: LinearMorphism(SetMorphism)

   This implements a linear map, encoded as a :py:class:`numpy.ndarray`.  Since it subclasses :py:class:`SetMorphism`, it inherits composition (as multiplication, which is of course *also* matrix multiplication) and call semantics.  It also stores the matrix as a new attribute

   .. py:attribute:: matrix

   as you might expect.

When constructing a :py:class:`SheafCoface`, if you pass an :py:class:`numpy.ndarray` as the `restriction` argument, PySheaf will automatically create :py:class:`LinearMorphism` objects as the restriction.

The final kind of morphism that is supported is :py:class:`SheafMorphism`, which supports composition by implementing a multiplication operator, but *not* call semantics since sheaf morphisms are *not* functions.  (It is true that they induce functions of various kinds, but PySheaf refrains from demanding that any particular kind of induced maps be computed by default.)

.. py:class:: SheafMorphism

   A sheaf morphism from one sheaf to another consists of lists :py:attr:`destinations` and :py:attr:`maps` which correspond to the :py:attr:`Sheaf.cells` list in the domain sheaf.

   .. py:attribute:: destinations

      List of indices into the codomain sheaf's :py:attr:`Sheaf.cells` list for each component map of the sheaf morphism.  Entries correspond to the :py:attr:`Sheaf.cells` list in the domain sheaf.

   .. py:attribute:: maps

      List of component maps (each of which may be any morphism class, like :py:class:`SetMorphism`, :py:class:`LinearMorphism`, or even :py:class:`SetMorphism`) corresponding to the :py:attr:`Sheaf.cells` list in the domain sheaf.
	 
Constructing :py:class:`Sheaf` instances
----------------------------------------

:py:class:`Sheaf` objects are constructed in essentially the same way as :py:class:`CellComplex` objects.  Determining the indices for the :py:attr:`Sheaf.cells` list is crucial, as each :py:attr:`SheafCoface.index` will refer into it.  Changes to the base space -- the inherited structure from :py:class:`CellComplex` -- are not easy and will generally involve many updates.  Additionally, each :py:attr:`SheafCoface.restriction` ought to be known beforehand, though these can be changed at run time if needed.
