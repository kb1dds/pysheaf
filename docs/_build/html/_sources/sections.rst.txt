The :py:class:`Section` class
=============================

A :py:class:`Section` represents a local section, which is an instance of data stored in a sheaf.  A :py:class:`Section` is constructed with a :py:class:`Sheaf` instance in mind, and therefore contains indices that refer into the :py:attr:`Sheaf.cells` list.  A given :py:class:`Sheaf` might have several :py:class:`Section` instances.

.. py:class:: Section

   .. py:attribute:: sectionCells

      The list of :py:class:`SectionCell` instances corresponding to this local section.  Although mathematically local sections are not multi-valued, it is possible that duplicates are present as there are no checks for this.

   .. py:method:: support()

      List the :py:attr:`Sheaf.cells` indices in the support of this :py:class:`Section`
   
   .. py:method:: extend(sheaf,cell,value=None,tol=1e-5)
      
      Extend this :py:class:`Section` to another cell whose index is `cell` in the :py:attr:`Sheaf.cells` list of the `sheaf` and returns `True` if this can be done consistently according to the tolerance `tol`.  You can optionally specify a `value` from the stalk over that cell; in this case the method can be used to test if this is a consistent choice or not.

.. py:class:: SectionCell

   A single value from the stalk over a given cell, consisting of the

   .. py:attribute:: value

   itself, which ought to be of the appropriate type (can be passed to the functions :py:attr:`SheafCell.metric` and/or :py:attr:`SheafCoface.restriction`).  One also needs to specify

   .. py:attribute:: support

   which records the :py:attr:`Sheaf.cells` index of the cell whose stalk from which this value was taken. 


Data fusion with :py:class:`Section` and :py:class:`Sheaf` instances
--------------------------------------------------------------------

Given some data in a :py:class:`Section` for a :py:class:`Sheaf` you can measure the overall consistency radius with

.. py:method:: Sheaf.consistencyRadius(assignment, tol=1e-5)

   where `assignment` is the :py:class:`Section` to be tested.  The optional `tol` specifies the tolerance for consistency, to be used in conjunction with each :py:attr:`SheafCell.metric` in the :py:class:`Sheaf`.

Similarly, if you want the nearest global section to your data, you can call
	       
.. py:method:: Sheaf.fuseAssignment(assignment, tol=1e-5)

   which returns a new :py:class:`Section` instance that is the global section nearest to your `assignment`.  In this case, the tolerance `tol` is passed to :py:func:`scipy.optimize.minimize`.
