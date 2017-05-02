Introduction
============

PySheaf is a Python module that implements cellular sheaves with a view towards computing useful invariants such as cohomology, consistency radius, and induced maps.  It is based on the idea that cell complexes are topological spaces, and that sheaves on cell complexes have at least that structure with some additional information as well.  The design follows the principles set out in the following book, which might be helpful to consult:

Michael Robinson, *Topological Signal Processing*, Springer, 2014.

The primary workflow for sheaves is

1. Construct the :py:class:`Sheaf` instance, which involves defining lists of :py:class:`SheafCell` and :py:class:`SheafCoface` instances.  Presently, the :py:class:`Sheaf` instance will remain fixed once constructed.  Therefore, make sure to have all stalks and restrictions defined at this point!
2. If you want to compute cohomology of the :py:class:`Sheaf`, you can do so using the :py:meth:`Sheaf.cohomology()` or :py:meth:`Sheaf.betti()`

For cell complexes, you can do both steps as well.  Instead of a :py:class:`Sheaf`, you define a :py:class:`CellComplex` built from lists of :py:class:`Cell` and :py:class:`Coface` instances.  
   
PySheaf is very much under active development and exploration, since there aren't well-established ways of computing using sheaves to process data.  Expect some rough edges, but feel free to suggest improvements!  Have a look at the latest version at the PySheaf GitHub repository: `<https://github.com/kb1dds/pysheaf>`_
