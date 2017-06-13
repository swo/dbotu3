#!/usr/bin/env python3
#
# author: scott olesen <swo@alum.mit.edu>

from __future__ import print_function
import argparse, sys
import pandas as pd
from Bio import SeqIO
from dbotu import read_sequence_table

def extract_rep_seqs(otu_table, fasta, output):
    table = read_sequence_table(otu_table)
    seq_counts = table.sum(axis=1).sort_values(ascending=False)
    records = SeqIO.index(fasta, 'fasta')

    def output_records():
        for seq, counts in seq_counts.iteritems():
            record = records[seq]
            record.id = '{};size={}'.format(seq, counts)
            record.description = ''
            yield record

    SeqIO.write(output_records(), output, 'fasta')

if __name__ == '__main__':
    p = argparse.ArgumentParser(description='Extract dbOTU3 representative sequences')
    p.add_argument('otu_table', type=argparse.FileType('r'), help='OTU table output by dbOTU3')
    p.add_argument('fasta', help='sequence file (used as dbOTU3 input)')
    p.add_argument('--output', '-o', default=sys.stdout, type=argparse.FileType('w'), metavar='FILE', help='representative sequences fasta')
    args = p.parse_args()

    extract_rep_seqs(args.otu_table, args.fasta, args.output)
