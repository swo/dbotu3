import pytest
from dbotu import *

import numpy as np, pandas as pd
import scipy.stats, scipy.optimize
from Bio import SeqIO
import os.path, itertools
# For Python 2 compatibility
try:
    from StringIO import StringIO
except:
    from io import StringIO

def named_stringio(contents, name='stringio'):
    x = StringIO(contents)
    x.name = name
    return x


input_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', 'data', 'input')
fasta_fn = os.path.join(input_dir, 'seq.fa')
table_fn = os.path.join(input_dir, 'counts.txt')

@pytest.fixture
def caller():
    records = SeqIO.index(fasta_fn, 'fasta')
    table = read_sequence_table(table_fn)
    return DBCaller(table, records, max_dist=0.10, min_fold=10.0, threshold_pval=0.0005)

def test_init(caller):
    assert all(caller.seq_abunds[0:3] == pd.Series([11184, 8896, 8723], index=['seq106', 'seq53', 'seq86']))

def test_process_one(caller):
    caller._process_record(caller.seq_abunds.index[0])
    assert caller.otus[0] == OTU('seq106', 'TGCTGCCTCCCGTAGGAGTCTGGGCCGTGTCTCAGTCCCAGTGTGGCCGGTCACCCTCTCAGGTCGGCTACGCATCGTCGCCTTGGTGGGCCGTTACCCCGCCAACCAGCTAATGCGCCATAAGTCCATCCTCTACCAGTGCTTTGAGGCACTTTTAATACGGTCACCATGCAGTGTCCGTACCTAT', [2484, 2151, 2341, 3474, 724, 10])

def test_process_until_first_merge(caller):
    member_id = 'seq32'
    otu_id = 'seq84'
    otu_seq = 'TGCTGCCTCCCGTAGGAGTTTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTAGGCCGTTACCCCACCAACTAGCTAATCAGACGCGGGTCCATCTCATACCACCGGAGTTTTTCACACTGTGCCATGCAGCACTGTGCGCTTATGCGGTAT'
    member_counts = np.array([26, 28, 40, 5, 0, 0])
    otu_counts = np.array([462, 398, 622, 17, 96, 1])
    merged_counts = member_counts + otu_counts

    for i in itertools.count():
        caller._process_record(caller.seq_abunds.index[i])

        if caller.seq_abunds.index[i] == member_id:
            break

    assert caller.membership[otu_id] == [otu_id, member_id]
    print([o.name for o in caller.otus])
    assert [o for o in caller.otus if o.name == otu_id][0] == OTU(otu_id, otu_seq, merged_counts)

def test_full_process():
    records = SeqIO.index(fasta_fn, 'fasta')
    table = read_sequence_table(table_fn)

    table_fh = open(table_fn)

    test_otu_fh = named_stringio('', 'otu_filename')
    test_log_fh = named_stringio('', 'log_filename')
    test_membership_fh = named_stringio('', 'membership_filename')

    call_otus(table_fh, fasta_fn, test_otu_fh, gen_crit=0.10, abund_crit=10.0, pval_crit=0.0005, log=test_log_fh, membership=test_membership_fh)

    test_log_content = test_log_fh.getvalue()
    test_otu_content = test_otu_fh.getvalue()
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
    assert membership_content == test_membership_content

    # I don't compare log_content because sometimes the last decimal place in
    # the numerical values is different. I could write a function that compares
    # the values in the log.
    #
    # assert log_content == test_log_content

improper_input_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'improper_data')
csv_table_fn = os.path.join(improper_input_dir, 'counts.csv')
biom_table_fn = os.path.join(improper_input_dir, 'counts_biom.txt')

def test_warn_csv_format():
    with pytest.warns(RuntimeWarning):
        read_sequence_table(csv_table_fn)

def test_warn_biom_format():
    with pytest.warns(RuntimeWarning):
        read_sequence_table(biom_table_fn)
