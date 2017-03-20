===============
Getting started
===============

Installation
============

You can download dbOTU3 using ``pip install dbotu``. dbOTU3 should be compatible
with Python 2 and 3. Installing with ``pip`` will add ``dbotu.py`` to your path.

Getting your data in shape
==========================

To run this software, you will need:

- A table of sequence counts. **This table has a very specific format**: the
  file is tab-separated; the first column is sequence IDs; the rest of the
  column headers are sample names. Each cell is the number of times that
  sequence appears in that sample. If you use convert a BIOM_ file to a TSV
  using ``biom convert --to-tsv``, you will need to remove the first line
  (i.e., ``# Constructed from biom file``).
- A fasta file containing the sequences to be processed into OTUs. The
  sequences should *not* be aligned. Single-end reads should be trimmed to
  the same length. Paired-end reads should be merged but not trimmed.

.. _BIOM: http://biom-format.org/

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

The ``--output`` option specifies where the resulting OTU table should go. The
``--membership`` option specifies that a QIIME-style membership file should be
created (one line for each OTU; the representative sequence ID is the first
field, all member sequence IDs are tab-separated after that). The ``--log``
option give some verbose information about exactly what tests were run for
which sequences.

Evaluate
========

If you can validate your choices for parameters, do so. You might also want to
chimera-check the OTUs, possibly with a script like my `UCHIME chimera checker
<https://github.com/swo/uchime-chimera-check>`_.

If you want to get into the specifics of what the algorithm did, you can read
the log file (produced by using the ``--log`` option). The log file has 5 types
of lines, all of which are tab-separated:

- Lines like ``A abundance_check B C`` show that the abundance criterion was
  applied to sequence ``A`` and all existing OTUs, of which OTUs ``B`` and
  ``C`` passed. Thus, the genetic similarity of ``A`` will be tested against
  ``B`` and ``C``. If no fields follow ``abundance_check``, no OTUs passed the
  abundance criterion.
- Lines with ``genetic_check`` are like ``abundance_check``: the sequence in
  the first field was sufficiently genetically similar to the OTUs after the
  ``genetic_check`` field to qualify for a distrubtion test.
- Lines like ``A distribution_check B 0.001`` mean that the distribution of
  sequence ``A`` was compared with that of OTU ``B`` and the distribution
  criterion returned a :math:`p`-value of 0.001.
- Lines like ``A new_otu`` show that sequence ``A`` was made into a new OTU.
- Lines like ``A new_otu B`` show that sequence ``A`` was merged into OTU ``B``.
