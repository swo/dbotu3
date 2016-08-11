==================
Overview of dbOTUs
==================

dbOTU-calling algorithm
=======================

The original pipeline, described in Preheim *et al.* [#preheim]_ was:

1. Process 16S data up to dereplicated, provenanced sequences.
2. Align those reads. Using the alignment, make a phylogenetic tree and a "distance matrix" showing the genetic distance between sequences.
3. Feed the distance matrix and the table of sequence counts into the algorithm proper, which groups the sequences into OTUs.

In outline, step 3 meant:

1. Make the most abundant sequence an OTU.
2. For each sequence (in order of decreasing abundance), find the set of OTUs that meet "abundance" and "genetic" cutoffs. The abundance cutoff requires that the candidate sequence be some fold smaller than the OTU (e.g., so that it can be considered sequencing error). The genetic cutoff requires that the candidate sequence be sufficiently similar to the OTU.  
3. If no OTUs meet these two criteria, make this sequence into an OTU.
4. If OTUs do meet these criteria, then, starting with the most closely-genetically-related OTU, check if this sequence is distributed differently among the samples than that OTU. If the distributions are sufficiently similar, merge this sequence into that OTU and go on to the next sequence.
5. If this candidate sequence does not have a distribution across sample sufficiently similar to an existing OTU, then make this sequence a new OTU.
6. Move on to the next candidate sequence.

Previous implementations
========================

The original implementation (`github <https://github.com/spacocha/Distribution-based-clustering>`_),
in Perl and shell scripts,
took a genetic distance matrix (a Jukes-Cantor distance
computed using FastTree) as input. When the number of sequences to be clustered
became too large, the entire matrix would not fit into memory.

The second implementation (`github <https://github.com/spacocha/dbOTUcaller>`_),
in Python interfaced with R using ``r2py``, took a set of aligned sequences as
input and computed the Hamming distance between these sequences as necessary.
This reduced the memory required (since it was no longer an entire matrix of all
pairwise distances).

In these implementations, two sequences were considered differently
distributed based on a :math:`\chi^2` test. Because many of the comparisons involved
sequences with small numbers of counts, the :math:`p`-value from the :math:`\chi^2` test
had to be computed empirically, which was expensive.

This implementation
===================

This implementation aims to improve speed and ease of use. Rather than aligning
sequences, this implementation uses a :math:`k`-mer-based distance to evaluate
the genetic dissimilarity criterion. The merit of this choice is discussed in
:ref:`distribution_section`.

Rather than using an empirical :math:`\chi^2` test, this implementation uses a
likelihood ratio test. The merit of this choice is discussed elsewhere in
:ref:`genetic_section`.


.. [#preheim] Preheim *et al.*, Distribution-Based Clustering: Using Ecology To
   Refine the Operational Taxonomic Unit. *Appl Environ Microbiol* (2013)
   doi:`10.1128/AEM.00342-13 <http://dx.doi.org/10.1128/AEM.00342-13>`_.
