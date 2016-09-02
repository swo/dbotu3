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

This implementation uses the Levenshtein edit distance to measure the
dissimilarity between two sequences. A sequence is disqualified from
being merged into an OTU if

.. math::

   \frac{E}{\tfrac{1}{2}(\ell_\text{seq} + \ell_\text{OTU})}

is greater than some threshold, where :math:`E` is the Levenshtein
edit distance between the sequence and the OTU, :math:`\ell_\text{seq}`
is the length of the candidate sequence, and :math:`\ell_\text{OTU}` is
the length of the OTU.

.. _evaluating-genetic-section:

Evaluating a distance criterion
-------------------------------

I found that this dissmilarity correlates with the dissimilarity computed
by a pairwise alignment just as well as the original implementation's
dissimilarity metric (i.e., the minimum of the aligned and unaligned
sequence dissimilarity).
It may be worth testing this relationship for your own data.
