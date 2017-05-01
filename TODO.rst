PySheaf master tasks list:
==========================

1. Implement maximal local sections algorithm as a method in :py:class:`Sheaf`.  See `<https://arxiv.org/abs/1612.00397>`_

2. Finish implementing and testing local homology for cell complexes:
   
 a. :py:meth:`CellComplex.inducedMapLocalHomology()`
    
 b. Clean up :py:meth:`CellComplex.localPairComplex()`

3. Finish implementing and testing :py:class:`AmbiguitySheaf`

5. Replace list semantics throughout with dictionaries.  Especially where random accesses are important, for instance in
   
 a. :py:class:`CellComplex`
    
 b. :py:class:`Sheaf`

6. Generalize :py:func:`inducedMap` to handle categories other than finite dimensional vector spaces
