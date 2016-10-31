===============
Getting started
===============

Installation
============

You will probably want to use `dbotu.py` as a command-line script, so I
recommend you clone the github repo::

    cd directory_for_code
    git clone https://github.com/swo/dbotu3.git

If you want to be able to import the code, you can set it up in "development"
mode with::

    python3 setup.py develop

After that, you should be able to ``import dbotu`` and get everything. If you
want to "unlink" this development version, you can::

    python3 setup.py develop --uninstall

Getting your data in shape
==========================

To run this software, you will need:

- A table of sequence counts. The first column is sequence IDs; the rest of the
  column headers are sample names. Each cell is the number of times that
  sequence appears in that sample.
- A fasta file containing the sequences to be processed into OTUs. The
  sequences should *not* be aligned, but they *should* be trimmed to the
  same length if they are single-end reads. (Paired-end reads should have been
  merged.)

The table of sequence counts will be read into memory, but the fasta file
will be indexed. As per the algorithm, the sequences will be processed in
order of decreasing abundance, but neither the table nor the fasta file need
to be in any particular order. The software will throw an error if there are
sequence IDs in the table that are not in the fasta.

Deciding on parameters
======================

You'll need to pick values for a few parameters:

- The *abundance cutoff*. The original paper suggests using :math:`10.0` to create OTUs
  that account just for sequencing error and :math:`0.0` to create OTUs that merge
  ecological populations. The default is :math:`10.0`.
- The *genetic cutoff*. The original paper suggests using :math:`0.10` as a cutoff
  of genetic dissimilarity. This means that sequences that are :math:`10\%` different
  will definitely not be put into the same OTU. The default is :math:`0.10`.
- The *distribution cutoff*. This is the :math:`p`-value from the statistical
  test of distribution described in :ref:`distribution_section`. The default
  value is :math:`0.0005` (as suggested in the original publication), although some
  testing has suggested that smaller :math:`p`-values might be more sensible.

Engage!
=======

Push the button::

    dbotu.py my-sequence-table.txt my-fasta.fasta -o my-otu-table.txt

If you wanted to change some of the default parameters, you can find the
relevant options using::

    dbotu.py --help

The `--output` option specifies where the resulting OTU table should go. The
`--membership` option specifies that a QIIME-style membership file should be
created (one line for each OTU; the representative sequence ID is the first
field, all member sequence IDs are tab-separated after that). The `--log`
option give some verbose information about exactly what tests were run for
which sequences.

Evaluate
========

If you can validate your choices for parameters, do so. You might also want
to chimera-check the OTUs, possibly with the script in ``tools/``.
