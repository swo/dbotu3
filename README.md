# Distribution-based OTU calling

A new implementation of Sarah Preheim's
[dbOTU](http://aem.asm.org/content/79/21/6593.long) algorithm.
The scope is narrower, the numerical comparisons are faster, and the interface
may be more user-friendly.

Read the [documentation](http://dbotu3.readthedocs.io/en/latest/) for:

- a guide to getting started,
- an explanation of the algorithm, and
- the API reference.

## Requirements
- Python 3
- Numpy, SciPy, BioPython, Pandas
- Levenshtein

## To-do

- Benchmark the quality of the output against gold standards and previous algorithm.
- Benchmark the speed against the previous algorithm.
- Figure out how to better integrate this into existing pipelines. Maybe there should be different kinds of output, like a table showing which sequence got assigned to which OTU.
- Better coverage for unit tests

## Authors

* **Scott Olesen** - *swo at mit.edu*
