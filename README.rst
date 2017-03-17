Distribution-based OTU calling
==============================

*dbotu3* is a new implementation of Sarah Preheim's dbOTU_ algorithm.  The
scope is narrower, the numerical comparisons are faster, and the interface is
more user-friendly.

.. _dbOTU: http://aem.asm.org/content/79/21/6593.long

Read the documentation_ for:

- a guide to getting started,
- an explanation of the algorithm, and
- the API reference.

.. _documentation: http://dbotu3.readthedocs.io/en/latest/

You can also read our new manuscript_ for more technical details about the
algorithm.  The Alm Lab website_ also has a short page with information.

.. _manuscript: http://dx.doi.org/10.1101/076927
.. _website: http://almlab.mit.edu/dbotu3.html

Installation
------------

dbotu3 is on PyPi_ and can be installed with pip.

.. _PyPi: https://pypi.python.org/pypi/dbotu

Requirements
------------

- Numpy, SciPy, BioPython_, Pandas_
- Levenshtein_

.. _BioPython: http://biopython.org
.. _Pandas: http://pandas.pydata.org
.. _Levenshtein: https://pypi.python.org/pypi/python-Levenshtein

Version history
---------------

- 1.1: Corrected error where sequence IDs that could be read as integers would not be found in the table
- 1.2: Python 2 compatibility, tox test framework, warnings for improperly-formatted sequence count tables
- 1.2.1: Added setup requirements

To-do
-----

- Better coverage for unit tests

Author
------

Scott Olesen / *swo at alum.mit.edu*
