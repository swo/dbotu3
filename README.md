# Distribution-based OTU calling

A modification of Sarah Preheim's [dbOTU](http://aem.asm.org/content/79/21/6593.long) caller.
The scope is narrower, the numerical comparisons are faster, and the interface may be more
user-friendly.

## Scope of the project

The original pipeline was:

1. Process 16S data up to dereplicated, provenanced sequences.
2. Align those reads. Using the alignment, make a phylogenetic tree and a "distance matrix" showing the genetic distance between sequences.
3. Feed the distance matrix and the table of sequence counts into the algorithm proper, which groups the sequences into OTUs.

This project is aiming to replace step 3. In outline, step 3 meant:

1. Make the most abundant sequence an OTU.
2. For each sequence (in order of decreasing abundance), find the set of OTUs that meet "abundance" and "genetic" cutoffs. The abundance cutoff requires that the candidate sequence be some fold smaller than the OTU (e.g., so that it can be considered sequencing error). The genetic cutoff requires that the candidate sequence be sufficiently similar to the OTU.  
3. If no OTUs meet these two criteria, make this sequence into an OTU.
4. If OTUs do meet these criteria, then, starting with the most closely-genetically-related OTU, check if this sequence is distributed differently among the samples than that OTU. If the distributions are sufficiently similar, merge this sequence into that OTU and go on to the next sequence.
5. If this candidate sequence does not have a distribution across sample sufficiently similar to an existing OTU, then make this sequence a new OTU.
6. Move on to the next candidate sequence.

The original implementation took a genetic distance matrix (a Jukes-Cantor distance
computed using FastTree) as input. When the number of sequences to be clustered
became too large, the entire matrix would not fit into memory. In this implementation,
the genetic dissimilarities between sequences and OTUs is computed using k-mers.

In the original algorithm, two sequences were considered differently
distributed based on chi-square test. Because many of the comparisons involved
sequences with small numbers of counts, the p-value from the chi-square test
had to be computed empirically, which was expensive. In this software, I use
a likelihood ratio test, which I think performs very well.

More details about the math of DBC are in the `doc` folder, which is currently a markdown
file intended to be turned into a pdf with

    pandoc --to latex --output math.pdf math.md

## Getting started

### Prerequisites

- A table of sequence counts. The first column is sequence IDs; the rest of the column headers are sample names. Each cell is the number of times that sequence appears in that sample.
- A fasta file containing the sequences to be processed into OTUs. The sequences should *not* be aligned.
- Python 3
- Numpy, SciPy

### Installing

(placeholder)

### Example workflow

Decide on the maximum genetic distance you're willing to accept. This is encoded in a squared
Euclidean distance between k-mer profiles (i.e., Equation 1 in [here](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2674673/)). As as a rule-of-thumb, you can say (maximum acceptable number of errors) x (length of kmer).
The default k is 8.

Feed the sequences and the table of sequence counts by samples into this program:

    dbotu.py my-sequence-table.txt my-fasta.fasta X -o my-otu-table.txt

Chimera-check the OTUs, possibly with the script in `tools/`.

## Running the tests

The testing framework is [py.test](http://docs.pytest.org/en/latest/). The tests (in `test/`) can be run with `make test` in the top directory.

## Contributing

(placeholder)

## To-do

- Benchmark the quality of the output against gold standards and previous algorithm.
- Benchmark the speed against the previous algorithm.
- Figure out how to better integrate this into existing pipelines. Maybe there should be different kinds of output, like a table showing which sequence got assigned to which OTU.

## Authors

* **Scott Olesen** - *swo at mit.edu*
