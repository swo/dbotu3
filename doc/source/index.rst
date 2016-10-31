.. dbotu documentation master file, created by
   sphinx-quickstart on Thu Aug 11 11:26:23 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Distribution-based OTU calling
==============================

This software, dbOTU3, is the third major implementation of the
distribution-based OTU-calling algorithm formulated by Preheim *et al.*
[#preheim]_, an extremely accurate algorithm for grouping DNA sequences from
microbial communities into OTUs for ecological or biomedical research.

Unlike most OTU-calling approaches, which group sequences based only on the
similarities of the sequences themselves, this algorithm also uses information
about the distribution of sequences across samples. This allows dbOTU to
distinguish ecologically-distinct but sequence-similar organisms or
populations.

This documentation includes a guide to getting started, description of the
algorithm, a comparison of the methodologies of the different implementations,
and an API reference.

.. toctree::
   :maxdepth: 2

   getting-started
   dbotu
   genetic
   distribution
   development
   api


Indices and tables
==================

* :ref:`genindex`
* :ref:`search`

.. [#preheim] Preheim *et al.* Distribution-Based Clustering: Using Ecology To
   Refine the Operational Taxonomic Unit. *Appl Environ Microbiol* (2013)
   doi:`10.1128/AEM.00342-13 <http://dx.doi.org/10.1128/AEM.00342-13>`_.
