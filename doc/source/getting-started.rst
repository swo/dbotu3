===============
Getting started
===============

Installation
============

Sorry, we're still in development, so there isn't a flashy way to install yet.
I recommend you clone the github repo::

    cd directory_for_code
    git clone https://github.com/swo/dbotu2.git

Then set it up in "development" mode with::

    python3 setup.py develop

After that, you should be able to ``import dbotu`` and get everything. If you want to "unlink"
this development version, you can::

    python3 setup.py develop --uninstall

Getting your data in shape
==========================

To run this software, you will need:

- A table of sequence counts. The first column is sequence IDs; the rest of the
  column headers are sample names. Each cell is the number of times that
  sequence appears in that sample.
- A fasta file containing the sequences to be processed into OTUs. The
  sequences should *not* be aligned.

The table of sequence counts will be read into memory, but the fasta file
will be indexed. As per the algorithm, the sequences will be processed in
order of decreasing abundance, but neither the table nor the fasta file need
to be in any particular order. The software will throw an error if there are
sequence IDs in the table that are not in the fasta.

Deciding on parameters
======================

You'll need to pick values for a few parameters:

- The *abundance cutoff*. The original paper suggests using 10 to create OTUs
  that account just for sequencing error and 0 to create OTUs that merge
  ecological populations. The default is 10.
- The *genetic cutoff*. As described in :ref:`genetic_section`, the genetic
  cutoff is a :math:`k`-mer distance. As a rule of thumb, you can use the product
  of the maximum number of acceptable differences between sequences times the
  word size. For example, if you have 300 nucleotide amplicons, you want your
  OTUs to be no larger than 90% (i.e., sequences no more than 10% different from
  the OTU's representative sequence), and you are using the default word size of 8,
  then you would want a genetic cutoff :math:`300 \times 0.1 \times 8 = 240`.
  There is no default value.
- The *distribution cutoff*. This is the :math:`p`-value from the statistical
  test of distribution described in :ref:`distribution_section`. The default
  value is :math:`0.0005` (as suggested in the original publication).
- The *word size*. The default value is 8. I would change this value only if
  you validate a different size as described in :ref:`evaluating-genetic-section`. 

Engage!
=======

Push the button::

    dbotu.py my-sequence-table.txt my-fasta.fasta D -o my-otu-table.txt

If you wanted to change some of the default parameters, you can find the
relevant options using::

    dbotu.py --help

Evaluate
========

If you can validate your choices for parameters, do so. You might also want
to chimera-check the OTUs, possibly with the script in ``tools/``.
