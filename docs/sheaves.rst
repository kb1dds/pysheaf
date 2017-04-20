The :py:class:`Sheaf` type
====================

The :py:class:`Sheaf` class derives from :py:class:`CellComplex` to describe its base space, the additional attributes and methods describe the "vertical fiber" structure of the sheaf.

.. py:class:: Sheaf(CellComplex)

   The base sheaf class in PySheaf consists of a list :py:attr:`cell` of :py:class:`SheafCell` instances, describing the cells and the stalks over them.  Restriction maps are built into :py:class:`SheafCoface` instances that are stored with each :py:class:`SheafCell`.


   .. py:method:: coboundary(k)

   .. py:method:: cohomology(k)
		  
   .. py:method:: betti(k)

   The dimension of the :math:`k` sheaf cohomology space.

   .. warning :: This is not the dimension of the :py:meth:`CellComplex.homology()`!
		  
.. py:class:: SheafCell(Cell)

.. py:class:: SheafCoface(Coface)
	      
Constructing :py:class:`Sheaf` instances
----------------------------------------


Morphisms: :py:class:`SetMorphism` and :py:class:`LinearMorphism`
-----------------------------------------------------------------

