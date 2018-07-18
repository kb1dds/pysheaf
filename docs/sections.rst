The :py:class:`Assignment` class
=============================

A :py:class:`Assignment` represents a local section, which is an instance of data stored in a sheaf.  A :py:class:`Assignment` is constructed with a :py:class:`Sheaf` instance in mind, and therefore contains indices that refer into the :py:attr:`Sheaf.cells` list.  A given :py:class:`Sheaf` might have several :py:class:`Assignment` instances.

.. py:class:: Assignment

   .. py:attribute:: assignmentCells

      The list of :py:class:`AssignmentCell` instances corresponding to this assignment.  Although mathematically assignments and/or sections are not multi-valued, it is possible that duplicates are present as there are no checks for this.

   .. py:method:: support()

      List the :py:attr:`Sheaf.cells` indices in the support of this :py:class:`Assignment`
   
   .. py:method:: extend(sheaf,cell,value=None,tol=1e-5)
      
      Extend this :py:class:`Assignment` to another cell whose index is `cell` in the :py:attr:`Sheaf.cells` list of the `sheaf` and returns `True` if this can be done consistently according to the tolerance `tol`.  You can optionally specify a `value` from the stalk over that cell; in this case the method can be used to test if this is a consistent choice or not.

   .. py:method:: maximalExtend(sheaf,multiassign=False,tol=1e-5)

      Extend this :py:class:`Assignment` to a maximal assignment that's non-conflicting (if `multiassign=False`) or one in which multiple values can be given to a given cell (if `multiassign=True`).

.. py:class:: AssignmentCell

   A single value from the stalk over a given cell, consisting of the

   .. py:attribute:: value

   itself, which ought to be of the appropriate type (can be passed to the functions :py:attr:`SheafCell.metric` and/or :py:attr:`SheafCoface.restriction`).  One also needs to specify

   .. py:attribute:: support

   which records the :py:attr:`Sheaf.cells` index of the cell whose stalk from which this value was taken. 


Data fusion with :py:class:`Assignment` and :py:class:`Sheaf` instances
--------------------------------------------------------------------

Given some data in a :py:class:`Assignment` for a :py:class:`Sheaf` you can measure the overall consistency radius with

.. py:method:: Sheaf.consistencyRadius(assignment, testSupport=None, tol=1e-5)

   where `assignment` is the :py:class:`Assignment` to be tested.  The optional `tol` specifies the tolerance for consistency, to be used in conjunction with each :py:attr:`SheafCell.metric` in the :py:class:`Sheaf`.

   The optional `testSupport` parameter is a list of cell indices on which consistency is to be assessed.  If it is listed as `None`, then the entire base space is to be tested.

   .. warning :: It is not assumed that `testSupport` is an open set in the topology of the underlying base space.  :py:meth:`Sheaf.consistencyRadius` automatically extends the `assignment` to be supported on the star over `testSupport`. 

   .. warning :: Consistency radius is measured using the cells specified in the :py:class:`Assignment` and all cells that are specified as :py:attr:`Sheaf.Cofaces`.  Preimages through :py:attr:`Sheaf.restriction` maps are not computed, so values on faces are not tested.

If you want to extend your :py:class:`Assignment` to be supported on all cells of its :py:class:`Sheaf`, leaving the existing :py:class:`Assignment` unchanged, while extending it as much as possible.  This is done either via :py:meth:`Assignment.maximalExtend` or the following:

.. py:method:: Sheaf.minimizeConsistencyRadius(assignment, activeCells=None, testSupport=None, method='nelder-mead', ord = np.inf, options={}, tol=1e-5)

   This constructs a new :py:class:`Assignment` instance based on an existing `assignment`.  The `activeCells` is the set of cells whose values are to be changed.  If `activeCells` is `None`, all cells outside the support of the assignment will be changed, but nothing in the support of the assignment will be changed.

   As in other methods, `testSupport` is the set of cells over which consistency radius is measured.

   Currently, any optimization supported by `scipy.optimize.minimize` is supported as a `method` oprtion, and `tol` is the passed to :py:func:`scipy.optimize.minimize`.

On the other hand, if you want the nearest global section to your data, you can call
	       
.. py:method:: Sheaf.fuseAssignment(assignment, activeCells=None, testSupport=None, method='SLSQP', options={}, tol=1e-5)

   which returns a new :py:class:`Assignment` instance that is the global section nearest to your `assignment`.  In this case, the tolerance `tol` is passed to :py:func:`scipy.optimize.minimize`.

   As in :py:meth:`Sheaf.consistencyRadius`, the `testSupport` specifies a list of cells under which consistency is measured.

   The `activeCells` argument is a list of cells whose stalks are allowed to be changed by the fusion process.  If passed as `None`, all values in the stalks over all cells may be changed.

   The method is a string, specifying the optimizer method to be used.  There are currently three optimizers implemented, 'KernelProj', 'SLSQP', and 'GA':

   'KernelProj': Uses kernel projection for sheaves of vector spaces.  In this case, every restriction map must be given as a :py:class:`LinearMorphism`, and :py:meth:`Sheaf.isLinear` must return `True`.  Kernel projection is usually the fastest and most accurate method if it is available.
   'SLSQP': This algorithm is :py:meth:`scipy.optimize.minimize` default for bounded optimization
   'GA': This genetic algorithm was implemented using DEAP for optimizations over nondifferentiable functions.  For this algorithm, it takes `options`: a dictionary to store changes to parameters, the keys must be identical to the current parameters.
   
    1. `initial_pop_size` - the number of individuals in the starting population
    2. `mutation_rate` - the proportion of the offspring (newly created individuals each round) that are from mutations rather than mating
    3. `num_generations` - the number of iterations that the genetic algorithm runs
    4. `num_ele_Hallfame` - the number of top individuals that should be reported in the hall of fame (hof)
 
Covers of :py:class:`CellComplex` instances based on a :py:class:`Assignment` of a :py:class:`Sheaf`
-------------------------------------------------------------------------------------------------

One can restrict attention to portions of a :py:class:`Assignment` instance.  This allows its consistency with its :py:class:`Sheaf` to be assessed locally.  By restricting the consistency testing to a list `testSupport` of cells in :py:meth:`Sheaf.consistencyRadius`, callers can examine consistency on a part of the base space of the :py:class:`Sheaf`.

Given a `threshold` for consistency, one can compute cover in which each list consists of cells that are consistent to that `threshold`.  This is computed by

.. py:method:: Sheaf.consistentCollection(self,assignment,threshold,testSupport=None,tol=1e-5)

   which relies on the aforementioned `threshold` and a :py:class:`Assignment` instance `assignment`.

   As in :py:meth:`Sheaf.consistencyRadius`, `testSupport` specifies what portion of the base :py:class:`CellComplex` is being analyzed.
