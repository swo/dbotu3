from __future__ import division
import pytest
from dbotu import *

import numpy as np, pandas as pd
import scipy.stats, scipy.optimize


def test_init():
    otu = OTU('id', 'ACGTACGT', [1, 2, 3, 4])
    assert otu.abundance == 1 + 2 + 3 + 4

def test_equal():
    otu1 = OTU('id', 'ACGT', [1, 2, 3])
    otu2 = OTU('id', 'ACGT', [1, 2, 3])
    assert otu1 == otu2


class TestAbsorb:
    def test1(self):
        otu1 = OTU('1', 'ACGT', [1, 2, 3, 4])
        otu2 = OTU('2', 'ACGT', [1, 2, 3, 4])
        otu1.absorb(otu2)
        assert all(otu1.counts == np.array([2, 4, 6, 8]))
        assert all(otu2.counts == np.array([1, 2, 3, 4]))


class TestDistanceTo:
    def test1(self):
        otu1 = OTU('1', 'AAGT', [])
        otu2 = OTU('2', 'ACGT', [])
        assert otu1.distance_to(otu2) == 1.0 / 4


class TestD:
    def test1(self):
        '''check that the difference of log likelihoods is what I computed by hand'''
        # pick some example count data
        # "a" is the count data for one OTU; "b" is the other
        a = np.array([1, 23, 2, 3])
        A = np.sum(a)
        b = np.array([10, 1797, 20, 5])
        B = np.sum(b)
        n = len(a)

        # (negative of) the log likelihood for the two models
        # 1: there are separate Poisson means for each OTU and each sample (there are 2n parameters in par)
        # 2: the means for a and b are related by a scaling factor (par[0]) (there are n + 1 parameters)
        neg_ll1 = lambda par: -sum(scipy.stats.poisson.logpmf(a, par[0: n] ** 2)) - sum(scipy.stats.poisson.logpmf(b, par[n:] ** 2))
        neg_ll2 = lambda par: -sum(scipy.stats.poisson.logpmf(a, par[1:] ** 2)) - sum(scipy.stats.poisson.logpmf(b, (par[0] * par[1:]) ** 2))

        # guess that the parameters for model 1 is just the counts for each
        par1 = np.sqrt(np.concatenate((a, b)))
        res1 = scipy.optimize.minimize(neg_ll1, x0=par1, method='Nelder-Mead')

        # guess that the parameters for model 2 are: scaled B / A and then the sum of counts for each sample scaled by A/A+B
        par2 = np.sqrt(np.insert(A / (A + B) * (a + b), 0, B / A))
        res2 = scipy.optimize.minimize(neg_ll2, x0=par2, method='Nelder-Mead')

        # using the results of those two maximizations, compute the D value numerically
        diff = -(res2.fun - res1.fun)
        numerical_D = -2.0 * diff

        # compare that numerical value against the analytical expression in DBCaller
        precision = 6
        assert round(numerical_D, precision) == round(OTU._D(a, b), precision)


class TestDistributionTestPval:
    def test1(self):
        a = np.array([1, 23, 2, 3])
        b = np.array([10, 1797, 20, 5])
        pval = OTU._distribution_test_pval(a, b)
        assert round(pval, 7) == 8.5e-5
