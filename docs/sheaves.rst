The :py:class:`Sheaf` class
===========================

The :py:class:`Sheaf` class derives from :py:class:`CellComplex` to describe its base space.  The additional attributes and methods of :py:class:`Sheaf` describe the "vertical fiber" structure of the sheaf.

.. py:class:: Sheaf(CellComplex)

   The base sheaf class in PySheaf consists of a list :py:attr:`cell` of :py:class:`SheafCell` instances, describing the cells and the stalks over them.  Restriction maps are built into :py:class:`SheafCoface` instances that are stored with each :py:class:`SheafCell`.

   .. py:method:: coboundary(k,compactSupport=False)

      Compute the degree `k` coboundary map of the :py:class:`Sheaf`, returning it as a :py:class:`numpy.ndarray`.   If you want compactly supported cohomology (if you don't know what that means, you don't) then set `compactSupport=True`.

   .. py:method:: cohomology(k,compactSupport=False)

      Compute the degree `k` sheaf cohomology for the :py:class:`Sheaf`.   If you want compactly supported cohomology (if you don't know what that means, you don't) then set `compactSupport=True`.  This returns a :py:class:`numpy.ndarray` whose columns are the generators for cohomology.
		  
   .. py:method:: betti(k)

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

.. py:class:: SheafCoface(Coface)

   A coface relation in a sheaf on a cell complex is built just like a :py:class:`Coface` in a :py:class:`CellComplex`, but with the addition of a :py:attr:`restriction`.
	      
   .. py:attribute:: restriction

      A restriction in a sheaf generally needs to be a morphism in some category.   It must be an object that supports *composition*, namely a class that implements a multiplication operator.  In most examples, this is a :py:class:`numpy.ndarray`.

Constructing :py:class:`Sheaf` instances
----------------------------------------

:py:class:`Sheaf` objects are constructed in essentially the same way as :py:class:`CellComplex` objects.  Determining the indices for the :py:attr:`Sheaf.cells` list is crucial, as each :py:attr:`SheafCoface.index` will refer into it.  Changes to the base space -- the inherited structure from :py:class:`CellComplex` -- are not easy and will generally involve many updates.  Additionally, each :py:attr:`SheafCoface.restriction` ought to be known beforehand, though these can be changed at run time if needed.
