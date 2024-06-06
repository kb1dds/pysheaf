"""
# Python 3.6 Sheaf theoretic toolbox

This library provides tools for manipulating sheaves of pseudometric spaces on partially ordered sets.  The primarly class is a `Sheaf`, which subclasses `NetworkX.DiGraph`; the idea is that the `DiGraph` specifies the Hasse diagram for the partially ordered set. 

Most of the details you need to get started are in the `pysheaf` module.  Start there first.  For specific sheaves that havew linear maps as restrictions, you might consider looking at `dataTools`.  Finally, there are several detailed examples to explore.

## General usage instructions

1. First (usually on paper!) lay out the cell complex that will serve as the base for your sheaf.  *Give each cell a unique label.*  

2. Determine all of the stalks over each cell, and the restriction maps.  Restriction maps can be a mixture of `numpy` matrices or arbitrary single-input Python function objects.
   
3. Construct a `Sheaf` instance and add each of your cells as `Cell` instances with the `Sheaf.AddCell` method.  Make sure to use your unique label for each `Cell`, because that is how PySheaf identifies them! Once you've done that, create each restriction as a `Coface` instance and add it to the sheaf using the `Sheaf.AddCoface` method.  The `Sheaf.AddCoface` method will connect the listed `Cell`s based on their labels.  `Cell`s and `Coface`s can be added later if you want, and they can be added in any order provided any `Coface` refers to `Cell`s that already exist.

4. Install some data into the sheaf by way of an `Assignment` to some of the `Cell`s.  

5. Analyze the sheaf and its data:
  a. You can compute consistency radius with `Sheaf.ComputeConsistencyRadius()`
  b. You can improve the consistency radius by extending or altering the values of the assignment with `Sheaf.FuseAssignment()`.  This will only alter Cells whose `Cell.mOptimizationCell` attribute is `True`.  You can also change the optimization algorithm if you want.
  c. You can find all star open sets whose local consistency is less than a desired bound using `Sheaf.CellIndexesLessThanConsistencyThreshold()`.

"""

from .pysheaf import *
