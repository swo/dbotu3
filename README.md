# Distribution-based OTU calling

*dbotu3* is a new implementation of Sarah Preheim's
[dbOTU](http://aem.asm.org/content/79/21/6593.long) algorithm.
The scope is narrower, the numerical comparisons are faster, and the interface
is more user-friendly.

Read the [documentation](http://dbotu3.readthedocs.io/en/latest/) for:

- a guide to getting started,
- an explanation of the algorithm, and
- the API reference.

You can also read our new [manuscript](http://dx.doi.org/10.1101/076927) for
more technical details about the algorithm.
The Alm Lab [website](http://almlab.mit.edu/dbotu3.html) also has a short page
with information.

## Installation

dbotu3 is on [PyPi](https://pypi.python.org/pypi/dbotu) and can be installed with pip.

## Requirements

- Numpy, SciPy, [BioPython](http://biopython.org), [Pandas](http://pandas.pydata.org)
- [Levenshtein](https://pypi.python.org/pypi/python-Levenshtein)

## Version history

- 1.1: Corrected error where sequence IDs that could be read as integers would not be found in the table

## To-do

- Better coverage for unit tests

## Author

Scott Olesen / *swo at alum.mit.edu*
