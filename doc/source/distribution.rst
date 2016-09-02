.. _distribution_section:

==========================
The distribution criterion
==========================

The mathematically critical part of the dbOTU algorithm is deciding whether or
not two OTUs are distributed identically.

Review of previous approaches
=============================

:math:`\chi^2` test
-------------------

The original paper and implementations dbOTU1 and dbOTU2 both articulate
the distribution criterion in terms of the :math:`\chi^2` test. The idea
is that Pearson's :math:`\chi^2` `test of independence <https://en.wikipedia.org/wiki/Pearson%27s_chi-squared_test#Test_of_independence>`_
evaluates whether the values of two outcomes (i.e., the distributions of
the two OTUs) are statistically independent.

The :math:`\chi^2` test runs into two problems. First, when counts are small, the
asymptotic (and easy-to-compute) :math:`\chi^2` statistic is not
accurate, so the :math:`p`-value needs to be simulated empirically,
which is computationally expensive. Second, when counts are large,
the :math:`\chi^2` test tends to be oversensitive: the sorts of
variations that we don't find unusual in microbiome data are
construed as true differences by the statistic.

Jensen-Shannon divergence
-------------------------

The original publication suggests that the distribution criterion can
be articulated as a cutoff on the JSD between two the distribution of
two OTUs (rather than on the :math:`p`-value from the :math:`\chi^2` test).

A JSD cutoff works well for large counts: it accords better with the sorts of
differences we would consider meaningful for microbiome data. This avoids
the oversensitivity problem that the :math:`\chi^2` test exhibits for
large numbers of counts.

However, because the JSD works with proportions and not counts, it
tends to be overly sensitive when the counts for one OTU are small.
For example, if one OTU has only one count in one sample, the JSD
treats that the same as if that OTU has a million counts in only that
sample.

Other approaches
----------------

Other distribution criteria have been informally evaluated. For example,
correlation coefficients (e.g., Kendall :math:`\tau`) have similar
strengths and weaknesses as the JSD: they perform well when the
number of counts is small but tend to be overly sensitive when the
counts are small.

This implementation
===================

I found that, for a few example cases, the asymptotic
`likelihood ratio test <https://en.wikipedia.org/wiki/Likelihood-ratio_test>`_,
which is quick to compute, gives results similar to the simulated
:math:`\chi^2` test. In the few examples I tested, the :math:`p`-value
for the likelihood ratio test were around two- to five-fold greater than the
:math:`p`-value from the simulated :math:`\chi^2` test. I consider this
a push in the right direction, since the simulated :math:`\chi^2` test
tended to be too sensitive.

Theoretical motivation
----------------------

The :math:`\chi^2` test measures whether two outcomes are independent.
In other words, it asks, given that the sequence has some total number
of counts, what is the probability that those counts would have been
distributed across samples in the same way that the potential parent OTU's
counts were distributed?

This formulation of the question is theoretically unpleasant because
this is not how sequence counts are distributed. It is not that a sequence
gets some total number of counts and those counts get distributed across
samples. Instead, each sample gets a total number of counts (based on
its fractional representation in the library and the depth of sequencing),
and that sample's counts get distributed among the sequences in that sample.

The likelihood ratio test evaluates a more natural description of
the question of ecologically identical distributions: allowing for some
Poisson error, is it more likely that the relative abundance of the
candidate sequence in each sample is proportional to the relative abundance
of the OTU in each sample (where the constant of proportionality is the
same for all samples)?

For one special case, this likelihood ratio test also gives a more sensible
answer that the :math:`\chi^2` test. If a sequence and an OTU are present
in only one sample, a distribution test should not disqualify them from
being merged: the data are perfectly concordant with the model that they
are distributed "the same". The likelihood ratio test gives :math:`p = 1`
in this case, but the :math:`\chi^2` test can give really low :math:`p`-values,
often below the default threshold!

Formulation of the null and alternative hypotheses
--------------------------------------------------

First, some definitions. There are :math:`N` samples. In sample :math:`i`, the number of reads
assigned to the more abundant OTU is :math:`a_i`; the candidate sequence
has :math:`b_i`. Define also :math:`A = \sum_{i=1}^N a_i` and similarly
:math:`B`.

The alternative hypothesis is that the OTU and sequence are distributed
differently, that is, that each of the :math:`a_i` and :math:`b_i` are
all drawn from different random variables. Technical replicates from
sequence data seem to be well-modeled by Poisson random variables [1]_,
so I formulate this hypothesis as

.. math::

   H_1: a_i \sim \mathrm{Poisson}(\lambda_{\mathrm{a}i}) \text{ and } b_i \sim \mathrm{Poisson}(\lambda_{\mathrm{b}i}),

where there are no constraints on the relationships between the Poisson
parameters.

The null model asserts that there is some relationship between the sequence
and the OTU, specifically that candidate sequence is distributed "the same as" the
OTU, just rescaled to some lower abundance. I articulate this as

.. math::

   H_0: a_i \sim \mathrm{Poisson}(\lambda_i) \text{ and } b_i \sim \mathrm{Poisson}(\sigma \lambda_i),

that is, that the :math:`a_i` are all free to be distributed
differently from one another, but each :math:`b_i` is constrained to be
distributed according a Poisson random variable that has the same mean
as the random variable corresponding to :math:`a_i`, just rescaled by some common scaling
factor :math:`\sigma`. (Because the sequence is less abundant than the OTU, we
expect that :math:`0 < \sigma < 1`.)

Maximum likelihood of models
----------------------------

The likelihood ratio test will require the maximum likelihood with
respect to those parameters. Not surprisingly, this requirement implies
that :math:`\lambda_{\mathrm{a}i} = a_i` and
:math:`\lambda_{\mathrm{b}i} = b_i`, i.e., that the best estimates for
the Poisson parameters are just the single number that distribution
produces. For the null hypothesis, we find that
:math:`\sigma = \frac{B}{A}`, i.e., the scaling factor is just the ratio
of the total counts for the sequence and OTU, and

.. math::

   \lambda_i = \frac{A}{A + B}(a_i + b_i).

Inserting these variables shows that the log likelihoods are:

.. math::

   \begin{aligned}
   \mathcal{L}_0 &= -(A + B) + \sum_i \left( a_i \ln a_i + b_i \ln b_i - \ln a_i! - \ln b_i! \right) \\
   \mathcal{L}_1 &= -(A + B) + A \ln A + B \ln B - (A + B) \ln (A + B) \\
     &\quad + \sum_i \left[ (a_i + b_i) \ln (a_i + b_i) - \ln a_i! - \ln b_i! \right]
   \end{aligned}

The difference between the logs can be conveniently written in terms of
a helper function:

.. math::

   f(\boldsymbol{x}) \equiv \sum_i x_i \ln x_i - \left( \sum_i x_i \right) \ln \left( \sum_i x_i \right)

so that

.. math::

   \mathcal{L}_1 - \mathcal{L}_0 = f(\boldsymbol{a} + \boldsymbol{b}) - f(\boldsymbol{a}) - f(\boldsymbol{b}).

The statistic
-------------

The likelihood ratio test uses the statistic
:math:`\Lambda = -2 \left( \mathcal{L}_1 - \mathcal{L}_0 \right)`, which
is distributed according to a :math:`\chi^2` distribution with
:math:`(N - 1)` degrees of freedom. (This is the difference in the
number of parameters in the two models: the alternative has :math:`2N`,
i.e., one for the OTU and the sequence in each sample, and the null has :math:`N + 1`, one
for each sample and the scaling factor :math:`\sigma`.) The cumulative
distribution function of :math:`\chi^2` at :math:`\Lambda` is easy to
compute.
