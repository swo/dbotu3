import pytest
from dbotu import *

import numpy as np, pandas as pd
import scipy.stats, scipy.optimize
from Bio import SeqIO
import os.path, io

input_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', 'data', 'input')
fasta_fn = os.path.join(input_dir, 'seq.fa')
table_fn = os.path.join(input_dir, 'counts.txt')

@pytest.fixture
def caller():
    records = SeqIO.index(fasta_fn, 'fasta')
    table = read_sequence_table(table_fn)
    return DBCaller(table, records, max_dist=0.10, min_fold=10.0, threshold_pval=0.0005)

def test_init(caller):
    assert all(caller.seq_abunds[0:3] == pd.Series([2523, 1831, 1793], index=['seq9', 'seq0', 'seq2']))

def test_process_one(caller):
    caller._process_record(caller.seq_abunds.index[0])
    assert caller.otus[0] == OTU('seq9', 'GGAATATTGGTCAATGGAGGCAACTCTGAACCAGCCATGCCGCGTGCAGGATGACGGCCCTATGGGTTGTAAACTGCTTTTGTACCAGAGAAAACCCGAGTACGTGTACTCGGTTGATAGTATGGTAAGAATAAGCATCGGCTAACTTCGTGCCAGCAGCCGCGGTAAGACGAAGGATGCAAGCGTTATCCGGATTCATTGGGTTTAAAGGGTGCGTAGGCGGACCTGTAAGTCAGTGGTGAAATCTCTTTGCTTAACAAAGAAACTGCCATTGATACTGCAAGTCTAGAGTATAGATGACGTTGGCGGAATATGACATGTAGTGGTGAAATACTTAGATATGTCATAGAACACCGATTGCGAAGGCAGC', [947, 749, 274, 162, 391])

def test_process_until_first_merge(caller):
    for i in range(17):
        caller._process_record(caller.seq_abunds.index[i])

    seq9_counts = np.array([947, 749, 274, 162, 391])
    seq1537_counts = np.array([38, 25, 8, 5, 10])
    new_counts = seq9_counts + seq1537_counts

    assert caller.membership['seq9'] == ['seq9', 'seq1537']
    assert caller.otus[0] == OTU('seq9', 'GGAATATTGGTCAATGGAGGCAACTCTGAACCAGCCATGCCGCGTGCAGGATGACGGCCCTATGGGTTGTAAACTGCTTTTGTACCAGAGAAAACCCGAGTACGTGTACTCGGTTGATAGTATGGTAAGAATAAGCATCGGCTAACTTCGTGCCAGCAGCCGCGGTAAGACGAAGGATGCAAGCGTTATCCGGATTCATTGGGTTTAAAGGGTGCGTAGGCGGACCTGTAAGTCAGTGGTGAAATCTCTTTGCTTAACAAAGAAACTGCCATTGATACTGCAAGTCTAGAGTATAGATGACGTTGGCGGAATATGACATGTAGTGGTGAAATACTTAGATATGTCATAGAACACCGATTGCGAAGGCAGC', new_counts)

def test_full_process():
    records = SeqIO.index(fasta_fn, 'fasta')
    table = read_sequence_table(table_fn)
    test_log_fh = io.StringIO()
    caller = DBCaller(table, records, max_dist=0.10, min_fold=10.0, threshold_pval=0.0005, log=test_log_fh)
    caller.generate_otu_table()

    test_log_content = test_log_fh.getvalue()

    test_otu_fh = io.StringIO()
    caller.write_otu_table(test_otu_fh)
    test_otu_content = test_otu_fh.getvalue()

    test_membership_fh = io.StringIO()
    caller.write_membership(test_membership_fh)
    test_membership_content = test_membership_fh.getvalue()

    output_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', 'data', 'output')
    otu_fn, log_fn, membership_fn = [os.path.join(output_dir, fn) for fn in ['otu.txt', 'log.txt', 'membership.txt']]

    with open(otu_fn) as f:
        otu_content = f.read()

    with open(log_fn) as f:
        log_content = f.read()

    with open(membership_fn) as f:
        membership_content = f.read()

    assert otu_content == test_otu_content
    assert log_content == test_log_content
    assert membership_content == test_membership_content
