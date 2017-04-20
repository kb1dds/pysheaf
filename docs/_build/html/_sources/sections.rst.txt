The :py:class:`Section` class
=============================

A :py:class:`Section` represents a local section, which are instances of data stored in sheaves.  They are constructed with a :py:class:`Sheaf` instance in mind, and therefore contain indices that refer into the :py:attr:`Sheaf.cells` list.  A given :py:class:`Sheaf` might have several :py:class:`Section` instances.


.. py:class:: Section

   .. py:attribute:: sectionCells

   The list of :py:class:`SectionCell` instances corresponding to this local section.  Although mathematically local sections are not multi-valued, it is possible that duplicates are present.

   .. py:method:: extend()

.. py:class:: SectionCell


Data fusion with :py:class:`Section` and :py:class:`Sheaf` instances
--------------------------------------------------------------------

.. py:method:: Sheaf.consistencyRadius()


.. py:method:: Sheaf.fuseAssignment()
