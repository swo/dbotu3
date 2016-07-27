#!/usr/bin/env python3

'''
confirm that the differences in log likelihood is what I computed by hand
'''

import numpy as np
import scipy.optimize, scipy.stats

x1s = [1, 23, 2, 3]
x2s = [10, 1797, 20, 5]
xs = x1s + x2s
A = sum(x1s)
B = sum(x2s)
n = len(x1s)

def ll1(par):
    print("ll1:", par)
    print("ll1:", x1s, x2s)
    l1s = par[0: n]
    l2s = par[n:]
    assert len(l1s) == n
    assert len(l2s) == n

    ll = sum(scipy.stats.poisson.logpmf(x1s, l1s ** 2)) + sum(scipy.stats.poisson.logpmf(x2s, l2s ** 2))

    return ll

def ll2(par):
    print("ll2:", par)
    sigma = par[0]
    ls = par[1:]
    assert len(ls) == n

    ll = sum(scipy.stats.poisson.logpmf(x1s, ls ** 2)) + sum(scipy.stats.poisson.logpmf(x2s, (sigma * ls) ** 2))

    return ll

def max_ll(ll_fun, par0):
    min_f = lambda x: -ll_fun(x)
    res = scipy.optimize.minimize(min_f, x0=par0, method='Nelder-Mead', options={'maxiter': 1e4})

    if not res.success:
        raise RuntimeError("optimization failed: {}".format(res.message))

    return -res.fun, res.x

def expected_diff_ll():
    return A * np.log(A) + B * np.log(B) - (A + B) * np.log(A + B) + sum([(a + b) * np.log(a + b) - a * np.log(a) - b * np.log(b) for a, b in zip(x1s, x2s)])

max_ll1, fopt1 = max_ll(ll1, np.sqrt(xs))
max_ll2, fopt2 = max_ll(ll2, np.sqrt([B / A] + [(A / (A + B)) * (a + b) for a, b in zip(x1s, x2s)]))

print("fopt1:", np.sqrt(fopt1))
print("fopt2:", np.sqrt(fopt2))

diff = max_ll2 - max_ll1
print("computed diff:", diff)
print("expected diff:", expected_diff_ll())

D = -2.0 * diff
p = scipy.stats.chi2.sf(D, df=(n - 1))
q = scipy.stats.chi2.cdf(D, df=(n - 1))
print("D:", D)
print("p:", p)
print("q:", q)
