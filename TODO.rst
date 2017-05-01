PySheaf master tasks list:
==========================

1. Implement maximal local sections algorithm as a method in `Sheaf`.  See `<https://arxiv.org/abs/1612.00397>`_

2. Finish implementing and testing local homology for cell complexes:
   
   a. `CellComplex.inducedMapLocalHomology()`
    
   b. Clean up `CellComplex.localPairComplex()`

3. Finish implementing and testing `AmbiguitySheaf`

4. Replace list semantics throughout with dictionaries.  Especially where random accesses are important, for instance in
   
   a. `CellComplex`
      
   b. `Sheaf`

5. Generalize `inducedMap` to handle categories other than finite dimensional vector spaces
