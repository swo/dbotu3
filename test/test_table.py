import pytest
from dbotu import read_sequence_table
from Bio import SeqIO
import os
try:
    from StringIO import StringIO
except:
    from io import StringIO



@pytest.fixture
def read_table():
    table_fh = StringIO('\n'.join(['\tS1\tS2\tS3', '1\t0\t1\t5', '2\t10\t25\t5', '3\t4\t0\t0', '4\t0\t0\t1']))
    return read_sequence_table(table_fh)

def test_shape(read_table):
    assert read_table.shape == (4,3)

def test_index(read_table):
    assert read_table.index[0] == str(1)

def test_seqIDs(read_table):
    table = read_table
    fasta_fh = StringIO('\n'.join(['>1', 'ATCG', '>2', 'TCGA', '>3', 'AGCT', '>4', 'TTCC']))
    fasta_fh.seek(0)
    records = SeqIO.to_dict(SeqIO.parse(fasta_fh, 'fasta'))
    missing_ids = [seq_id for seq_id in table.index if seq_id not in records]
    assert len(missing_ids) == 0

