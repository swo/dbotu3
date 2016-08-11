# Distribution-based OTU calling

A modification of Sarah Preheim's [dbOTU](http://aem.asm.org/content/79/21/6593.long) caller.
The scope is narrower, the numerical comparisons are faster, and the interface may be more
user-friendly.

Read the [documentation](http://dbotu2.readthedocs.io/en/latest/).

## Getting started

### Prerequisites

- A table of sequence counts. The first column is sequence IDs; the rest of the column headers are sample names. Each cell is the number of times that sequence appears in that sample.
- A fasta file containing the sequences to be processed into OTUs. The sequences should *not* be aligned.
- Python 3
- Numpy, SciPy, BioPython

### Installing

I'm not sure the package is fully functional. You can set it up in "development" mode with:

    python3 setup.py develop

After that, you should be able to `import dbotu` and get everything. If you want to "unlink"
this development version, you can

    python3 setup.py develop --uninstall

### Example workflow

First, gather your data (your cleaned sequences in a fasta and your table of sequence counts as described above).

Second, decide on a maximum genetic distance. The documentation explains this. As a rule of thumb,  (maximum acceptable number of errors) x (length of kmer). The default k is 8.

Then, feed the sequences and the table of sequence counts by samples into this program:

    dbotu.py my-sequence-table.txt my-fasta.fasta X -o my-otu-table.txt

Finally, chimera-check the OTUs, possibly with the script in `tools/`.

## Running the tests

The testing framework is [py.test](http://docs.pytest.org/en/latest/). The tests (in `test/`) can be run with `make test` in the top directory.

## Contributing

(placeholder)

## To-do

- Benchmark the quality of the output against gold standards and previous algorithm.
- Benchmark the speed against the previous algorithm.
- Figure out how to better integrate this into existing pipelines. Maybe there should be different kinds of output, like a table showing which sequence got assigned to which OTU.
- Improve speed of k-mer computations

## Authors

* **Scott Olesen** - *swo at mit.edu*
