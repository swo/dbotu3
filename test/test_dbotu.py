import pytest
from dbotu import *

import numpy as np, pandas as pd, io
import scipy.stats, scipy.optimize
from Bio import SeqIO

@pytest.fixture
def caller():
    fasta_content = io.StringIO("\n".join(['>seq1', 'AAAAAAA', '>seq2', 'AAATAAA', '>seq3', 'AATTTAAA']))
    fasta_content.seek(0)
    log_fh = io.StringIO()
    records = SeqIO.parse(fasta_content, 'fasta')
    table = pd.DataFrame(np.array([[0, 10, 20], [0, 1, 2], [10, 0, 0]]), index=['seq1', 'seq2', 'seq3'], columns=['sample1', 'sample2', 'sample3'])
    return DBCaller(table, records, word_size=3, max_dist=0.1, min_fold=0.0, threshold_pval=0.001, log=log_fh)

def test_seq_abunds(caller):
    assert all(caller.seq_abunds == pd.Series([30, 10, 3], index=['seq1', 'seq3', 'seq2']))

def test_empty_otus(caller):
    assert caller.otus == []

def test_process_one(caller):
    caller._process_record(next(caller.records))
    assert caller.otus[0] == OTU('seq1', 'AAAAAAA', [0, 10, 20], 3)

def test_process_two(caller):
    for i in range(2):
        caller._process_record(next(caller.records))

    assert caller.otus[0] == OTU('seq1', 'AAAAAAA', [0, 10, 20], 3)
    assert caller.otus[1] == OTU('seq3', 'AATTTAA', [10, 0, 0], 3)
