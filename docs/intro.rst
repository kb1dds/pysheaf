Introduction
============

PySheaf is a Python module that implements cellular sheaves with a view towards computing useful invariants such as cohomology, consistency radius, and induced maps.  It is based on the idea that cell complexes are topological spaces, and that sheaves on cell complexes have at least that structure with some additional information as well.  The design follows the principles set out in the following book, which might be helpful to consult:

Michael Robinson, *Topological Signal Processing*, Springer, 2014.

Although PySheaf can compute (co)homology of sheaves and cell complexes, the primary workflow for sheaves is

1. Construct the :py:class:`Sheaf` instance, which involves defining lists of :py:class:`SheafCell` and :py:class:`SheafCoface` instances.  Presently, the :py:class:`Sheaf` instance will remain fixed once constructed.  Therefore, make sure to have all stalks and restrictions defined at this point!
2. If you want to compute cohomology of the :py:class:`Sheaf`, you can do so using the :py:meth:`Sheaf.cohomology()` or :py:meth:`Sheaf.cobetti()`
3. Construct various :py:class:`Section` instances for the data you have.  Key point: a :py:class:`Section` refers to a *local* section.  Its :py:class:`SectionCell` members refer to the :py:class:`Sheaf` you've just finished building, so you must built the :py:class:`Sheaf` first!
4. Process the data using the :py:class:`Sheaf` and :py:class:`Section` instances:
   
   a. You can measure the consistency using :py:meth:`Sheaf.consistencyRadius()`
   b. You can fit a nearest global section using :py:meth:`Sheaf.fuseAssignment()`
   c. You may also use :py:meth:`Section.extend()`

5. If you want to relate your :py:class:`Sheaf` to others, you may construct a :py:class:`SheafMorphism` instance, which incidentally may *also* be used as a restriction morphism if you want to build a :py:class:`Sheaf` *of* sheaves!

For cell complexes, you can do steps (1) and (2).  Instead of a :py:class:`Sheaf`, you define a :py:class:`CellComplex` built from lists of :py:class:`Cell` and :py:class:`Coface` instances.  
   
PySheaf is very much under active development and exploration, since there aren't well-established ways of computing using sheaves to process data.  Expect some rough edges, but feel free to suggest improvements!
