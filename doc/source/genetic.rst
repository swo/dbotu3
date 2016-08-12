.. _genetic_section:

===========================
The genetic criterion
===========================

A critical piece of the dbOTU algorithm is determining which sequences
are sufficiently genetically dissimilar that they belong in different
OTUs regardless of their distribution across samples.

Review of previous approaches
=============================

Jukes-Cantor distance
---------------------

In dbOTU1, this genetic criterion was articulated
in terms of the ``makematrix`` output from
FastTree_.
First, the sequences were aligned to a common
16S alignment using PyNAST_. FastTree was then used to compute the
distances. As far as I can tell, FastTree simply computes the
Jukes-Cantor_ distance :math:`-\frac{3}{4} \log (1 - \frac{4d}{3})`,
where :math:`d` is the proportion of positions that differ (i.e., the
number of mismatches in the aligned sequences divided by the length of
the aligned sequences, ignoring gaps).

.. _Jukes-Cantor: https://en.wikipedia.org/wiki/Models_of_DNA_evolution#JC69_model_.28Jukes_and_Cantor.2C_1969.29.5B1.5D
.. _FastTree: http://www.microbesonline.org/fasttree/
.. _PyNAST: http://biocore.github.io/pynast/

In fact, the genetic criterion was evaluated using the
minimum of the aligned and unaligned Jukes-Cantor distances. This was a
weird hack: sometimes the alignment led to a greater
distance between two sequences than would be computed by just comparing
the unaligned sequences.

Hamming distance
----------------

In dbOTU2, the sequences were aligned with
PyNAST_, but the distances between the sequences were computed one at a
time to ease memory needs. Sequences were read one at a time from the
input fasta, and OTUs were stored in memory so that a sequence was only
compared to the OTUs to which it might potentially belong. In this
implementation, the dissimilarity between OTUs was just the proportion :math:`d`
of positions that differ (ignoring positions in which either or both
sequences had a gap).

Like in dbOTU1, the actual distance criterion was checked against the
minimum of the aligned and unaligned distances.

Evaluation of previous approaches
=================================

Aside from leading to occasionally confusing results, the alignment step
is moderately slow. Taking the minimum of two distances is an unpleasant
hack.

However, in general, the dissimilarities between sequences measured in different
ways are well-correlated. For example, it was not critical to use the
Jukes-Cantor distance rather than the simple Hamming distance.

This implementation
===================

New approach used in this implementation
----------------------------------------

This implementation uses a :math:`k`-mer-based distance to evaluate the
genetic criterion. A :math:`k`-mer-based method does not require
aligning the sequences, which eliminates the alignment step and the
confusing results requiring the minimum of the two distances. If
:math:`k`-mer compositions are stored in a sparse way, the memory
requirements for the sequence-based and :math:`k`-mer-based approaches
are similar.

The disadvantage to using a :math:`k`-mer-based is that the distance
metrics are not as familiar: it is easy to articulate that two sequences
should be at most 90% similar to be considered potentially part of the
same OTU.

I found that, for small genetic distances, :math:`k`-mer distances
correlate very well with the Jukes-Cantor distance and Hamming distance.

Selecting a distance criterion
------------------------------

In this implementation, the distance between two sequences is defined as

.. math::

   d_\text{kmer} \equiv \sum_i \left( c_{1i} - c_{2i} \right) ^2,

where :math:`c_{1i}` is the number of times the :math:`i`-th
:math:`k`-mer appears in sequence 1. In general, I found that the
:math:`k`-mer distance relates to the proportion :math:`d` of mismatched 
positions by

.. math::

   d_\text{kmer} \approx k L d,

where :math:`k` is the length of the :math:`k`-mer and :math:`L` is the
length of the sequence.

This approximation is easy to derive if you assume that the sequence is
long compared to the :math:`k`-mer, there are few differences between
the two sequences, the differences are separated from one another by a
number of positions greater than :math:`k`, and that each :math:`k`-mer
either does not appear in a sequence or appears only once. In this
regime, each different position leads to :math:`k` of the :math:`k`-mers being
different, changing from zero to one count or vice versa.

.. _evaluating-genetic-section:

Evaluating a distance criterion
-------------------------------

It may be worth testing this relationship for your own data. I verified
this rule of thumb by selecting a set of sequences from a dataset,
computing their dissimilarities using a variety of methods, and showing
that, for small genetic dissimilarities, all the results correlate well.
I regressed :math:`d_\text{kmer}` against :math:`d` to confirm that, in
the appropriate regime, the slope was approximately :math:`kL`.
