=====================================
The algorithm and its implementations
=====================================

To help keep things straight, I will distinguish between the dbOTU algorithm
and its three implementations dbOTU1, dbOTU2, and dbOTU3.

I will also define the word *provenience* (from the `archaeological term <https://en.wikipedia.org/wiki/Provenance>`_
referring to the physical location where an artifact was found) to mean the
information about how many times each unique DNA sequence appeared in each
sample in the sequencing library. (QIIME calls this the "`OTU mapping <http://qiime.org/scripts/merge_otu_maps.html>`_",
which I find confusing because "mapping" refers to many things.)

dbOTU algorithm
===============

Motivation
----------

The algorithm aims to separate genetically-similar sequences that appear to be
ecologically distinct (or, conversely, to join less-genetically-similar
sequences that appear to be ecologically identical). For example, if two sequences
differ by only one nucleotide and you had no provenience data, you would
probably put those sequences into the same OTU. However, if the two sequences
never appeared together in the same sample, you would probably conclude that
that one nucleotide difference corresponds to two distinct groups of organisms,
one which lives in one group of samples, the other living in the other.

Conversely, if two sequences had a few nucleotides different, you might, without
provenience data, place them into different OTUs. However, if the two sequences
appeared in the same ratio in all samples (e.g., sequence 2 was always almost
exactly ten times less abundant than sequence 1), you would probably conclude
that the second sequence was either sequencing error or a member of the same
ecological population as the first sequence.

As a cheesy example, consider my last name, Olesen. You might think this is just
an Ellis Island mispelling of the more common name Olsen. However, the name
Olsen is present among the abundant Norwegians and the rarer Danes, while
Olesen is only present among the Danes. This provenience data would lead you
to (correctly) conclude that Olesen is a distinctly Danish name that is the
the result of differences between the Norwegian and Danish languages and
orthography. (The Swedish equivalent, Olsson, is so "genetically" different that
you probably do not need provenience data to know that it has a different
"ecology".)

Mechanics
---------

The original pipeline, described in Preheim *et al.* [#preheim]_ was:

1. Process 16S data up to dereplicated, provenienced sequences.
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

The implementations vary in terms of:

* The exact input files they required
* How they evaluated the genetic (i.e., sequence similarity) criterion
* How they evaluated the distribution (i.e., ecological similarity) criterion
* The details of the software

dbOTU1
------

The original implementation (`dbOTU1 <https://github.com/spacocha/Distribution-based-clustering>`_),
coded in Perl and shell scripts,
took a genetic distance matrix (a Jukes-Cantor distance
computed using FastTree_) as input.
using that distance matrix.

In this implementation, the genetic criterion was evaluated using the
minimum of the aligned and unaligned Jukes-Cantor distances. This was a
weird hack: sometimes the alignment, made using NAST [#nast]_ (actually the
PyNAST_ implementation), led to a greater
distance between two sequences than would be computed by just comparing
the unaligned sequences.

.. _FastTree: http://www.microbesonline.org/fasttree/
.. _PyNAST: http://biocore.github.io/pynast/

In this implementation the distribution criterion was evaluated using
the |chisq-test|_ function in R_,
called in a separate process from a Perl script.
Many of the comparisons involved
sequences with small numbers of counts, for which the asymptotic (i.e., commonly-used)
calculation of the :math:`p`-value of a :math:`\chi^2` test is not accurate. This implementation
therefore used a simulated :math:`p`-value, available through the R
commands ``simulate.p.value`` option. This empirical calculation
required many simulated contingency tables, which was expensive.

.. _R: https://www.r-project.org/about.html
.. |chisq-test| replace:: ``chisq.test``
.. _chisq-test: https://stat.ethz.ch/R-manual/R-devel/library/stats/html/chisq.test.html

dbOTU2
------

The second implementation (`dbOTU2 <https://github.com/spacocha/dbOTUcaller>`_),
coded in Python 2 and interfaced with R using `r2py <http://rpy2.bitbucket.org/>`_,
took a set of aligned sequences as
input and computed the Hamming distance between these sequences as necessary.
This reduced the memory required (since it was no longer an entire matrix of all
pairwise distances).

Like the first implementation, this one used the minimum of the aligned and
unaligned sequences.

Like the first implementation, this one used R's ``chisq.test``, but this time
called via ``r2py`` from the Python script. This removed the need for hacky
temporary files, but it was still slow and required R and Python to talk nicely
to one another.

This implementation also allowed for the distribution criterion to be articulated
in terms of the `Jensen-Shannon divergence <https://en.wikipedia.org/wiki/Jensen%E2%80%93Shannon_divergence>`_
(JSD). The JSD had some advantages over the :math:`\chi^2` test but suffered some
of the same weaknesses, as will be reviewed in :ref:`distribution_section`.


This implementation
===================

This implementation, dbOTU3, aims to improve speed and ease of use. It is written
in pure Python 3 and aims to `do one thing <https://en.wikipedia.org/wiki/Unix_philosophy#Do_One_Thing_and_Do_It_Well>`_,
namely to turn sequence and provenience data into OTUs.

Rather than requiring aligned sequences, this implementation uses a Levenshtein
edit distance as an approximation for the aligned sequence dissimilarity.
The merit of this choice is discussed in :ref:`genetic_section`.

Rather than using an empirical :math:`\chi^2` test, this implementation uses a
likelihood ratio test. The merit of this choice is discussed in
:ref:`distribution_section`.

A more thorough comparison of the implementations and an evaluation of the
accuracy and speed of this new implementation is in a separate technical
manuscript [#dbotu3]_ (although note the :ref:`genetic_caveat`).


.. [#preheim] Preheim *et al.* Distribution-Based Clustering: Using Ecology To
   Refine the Operational Taxonomic Unit. *Appl Environ Microbiol* (2013)
   doi:`10.1128/AEM.00342-13 <http://dx.doi.org/10.1128/AEM.00342-13>`_.

.. [#nast] DeSantis *et al.* NAST: a multiple sequence alignment server for
   comparative analysis of 16S rRNA genes. *Nucleic Acids Res* (2006)
   doi:`10.1093/nar/gkl244 <https://dx.doi.org/10.1093/nar/gkl244>`_.

.. [#dbotu3] SW Olesen, C Duvallet, EJ Alm. dbOTU3: A new implementation of
   distribution-based OTU calling. *bioRxiv* (2016) doi:`10.1101/076927 <http://dx.doi.org/10.1101/076927>`_.
