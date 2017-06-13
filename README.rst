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

You can also read our new paper_ for more technical details about the
algorithm.  The Alm Lab website_ also has a short page with information.

.. _paper: https://doi.org/10.1371/journal.pone.0176335
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
- 1.3.0: Improved OTU file header. Split the log file into a debug and progress log.
- 1.4.0: Made an improvement to the Levenshtein-based genetic dissimilarity metric.
- 1.4.1: Account for pandas API change to ``MultiIndex``
- 1.5.0: Added the restart and rep seq scripts

To-do
-----

- Testing for the restart scripts
- Better coverage for unit tests

Citation
--------

If you use dbOTU3 in a scientific paper, we ask that you cite the
original dbOTU publication (Preheim *et al*.) or the dbOTU3 publication:

Preheim *et al*. Distribution-Based Clustering: Using Ecology To Refine the
Operational Taxonomic Unit. *Appl Environ Microbiol* (2013) doi:10.1128/AEM.00342-13.

Olesen SW, Duvallet C, and Alm EJ. dbOTU3: A new implementation of
distribution-based OTU calling. *PLoS ONE* (2017) doi:10.1371/journal.pone.0176335.

Author
------

If you find a bug or have a request for a new feature, open an issue_.

.. _issue: https://github.com/swo/dbotu3/issues

Scott Olesen / *swo at alum.mit.edu*
