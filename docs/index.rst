.. PySheaf documentation master file, created by
   sphinx-quickstart on Tue Apr 11 15:52:28 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to PySheaf's documentation!
===================================

The PySheaf library is an experimental platform for computational and applied sheaf theory.   Have a look at the latest version at the PySheaf GitHub repository: `<https://github.com/kb1dds/pysheaf>`_

Sheaves are managed by subclassing from :py:class:`CellComplex`, which is a little awkward as :py:class:`SheafCell` is subclassed from :py:class:`Cell`, etc.  However, the idea is that everything that can be done to the base space of a sheaf is inherited from its underlying :py:class:`CellComplex` structure.

Sections of sheaves have their own :py:class:`Section` class.  Truly, the data structure does not self-check for consistency, so really they are more of an "assignment" or "serration" than a section.  However, they do know how to extend themselves.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   intro
   install
   cellcomplexes
   sheaves
   sections


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
