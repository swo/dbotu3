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
