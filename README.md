# Distribution-based OTU calling

A rework of Sarah Preheim's [dbOTU](http://aem.asm.org/content/79/21/6593.long) caller.

## Scope of the project

The original pipeline was:

1. Process 16S data up to dereplicated, provenanced sequences.
2. Align those reads. Using the alignment, make a phylogenetic tree and a "distance matrix" showing the genetic distance between sequences.
3. Feed the distance matrix and the table of sequence counts into the algorithm proper, which groups the sequences into OTUs.

This project is aiming to replace step 3. In outline, step 3 meant:

1. Make the most abundant sequence an OTU.
2. For each sequence (in order of decreasing abundance), find the set of OTUs that is above some abundance cutoff with respect to this sequence. If no OTUs exceed that cutoff, make this sequence a new OTU and go to the next sequence.
3. Otherwise, find the set of OTUs among that those satisfied the abundance criterion that also satisfy a genetic distance cutoff. If no OTUs are within that cutoff, make this sequence a new OTU and go to the next sequence.
4. Otherwise, starting with the most closely-genetically-related OTU, check if this sequence is distributed differently among the samples than that OTU. If not, merge this sequence into that OTU and go on to the next sequence.
5. Otherwise, make this sequence a new OTU and move on.

In the original algorithm, two sequences were considered differently
distributed based on $\chi^2$ test. Because many of the comparisons involved
sequences with small numbers of counts, the $p$-value from the $\chi^2$ test
had to be computed empirically, which was expensive. In this software, I use
a likelihood ratio test, which I think performs very well.

More details about the math of DBC are in the `doc` folder, which is currently a markdown
file intended to be turned into a pdf with

    pandoc --to latex --output math.pdf math.md

## Getting started

### Prerequisites

- A matrix describing genetic distances between your sequences. This software is currently configured to accept the input that comes from [FastTree](http://www.microbesonline.org/fasttree/)'s `-makematrix` option.
- Python 3
- Numpy, SciPy

### Installing

(placeholder)

## Running the tests

The testing framework is [py.test](http://docs.pytest.org/en/latest/). The tests (in `test/`) can be run with `make test` in the top directory.

## Contributing

## Authors

* **Scott Olesen** - *swo at mit.edu*

## License

(placeholder)
