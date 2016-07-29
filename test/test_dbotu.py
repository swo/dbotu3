import pytest
from dbotu import *

import numpy as np, pandas as pd
import scipy.stats, scipy.optimize

@pytest.fixture
def caller():
    names = ['seq1', 'seq2', 'seq3']
    matrix = pd.DataFrame(np.array([[0.0, 0.1, 0.2], [0.1, 0.0, 0.3], [0.2, 0.3, 0.0]]), index=names, columns=names)
    print(matrix)
    table = pd.DataFrame(np.array([[0, 10, 20], [0, 1, 2], [10, 0, 0]]), index=['seq1', 'seq2', 'seq3'], columns=['sample1', 'sample2', 'sample3'])
    return DBCaller(table, matrix, 0.1, 0.0, 0.001)


class TestInit:
    def test_seq_list(self, caller):
        '''they should be in abundance order'''
        assert (caller.seqs == ['seq1', 'seq3', 'seq2']).all()

    def test_seq_abunds(self, caller):
        assert (caller.seq_abunds == [30, 10, 3]).all()

    def test_fail_on_asymmetric(self):
        names = ['seq1', 'seq2']
        matrix = pd.DataFrame(np.array([[1, 2], [3, 4]]), index=names, columns=names)
        table = pd.DataFrame(np.array([[0, 0], [0, 0]]), index=names, columns=['sample1', 'sample2'])
        with pytest.raises(RuntimeError):
            DBCaller(table, matrix, 0.1, 0.0, 0.001)

    def test_fail_on_nonsquare(self):
        matrix = pd.DataFrame(np.array([[0, 0, 0], [0, 0, 0]]), index=['a', 'b'], columns=['a', 'b', 'c'])
        table = pd.DataFrame(np.array([[0, 0], [0, 0]]), index=['a', 'b'], columns=['sample1', 'sample2'])
        with pytest.raises(ValueError):
            DBCaller(table, matrix, 0.1, 0.0, 0.001)


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
        assert round(numerical_D, precision) == round(DBCaller._D(a, b), precision)
