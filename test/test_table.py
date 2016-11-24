import pytest
from dbotu import read_sequence_table
from Bio import SeqIO
import os

tablefn = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'table_test.counts')
fastafn = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'table_test.fasta')

@pytest.fixture
def read_table():
    return read_sequence_table(tablefn)

def test_shape(read_table):
    assert read_table.shape == (4,3)

def test_index(read_table):
    assert read_table.index[0] == str(1)

def test_seqIDs(read_table):
    table = read_table
    records = SeqIO.index(fastafn, 'fasta')
    missing_ids = [seq_id for seq_id in table.index if seq_id not in records]
    assert len(missing_ids) == 0


