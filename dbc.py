#!/usr/bin/env python3
#
# author: scott olesen <swo@mit.edu>

import argparse, sys
import pandas as pd

def read_otu_table(fn):
    return pd.read_table(fn, index_col=0, header=0)

def read_fasttree_matrix(fn):
    # FastTree put it in a funny format: the number of rows in its own line at the top,
    # then each row gets a label
    matrix = pd.read_table(fn, sep='\s+', skiprows=1, header=None, index_col=0)
    matrix.columns = matrix.index
    return matrix

def merge_seq(table, abundances, otu, merged_seq):
    table.loc[otu] += table.loc[merged_seq]
    table = table.drop(merged_seq)

    abundances[otu] += abundances[merged_seq]
    abundances = abundances.drop(merged_seq)

    return table, abundances

def abundance_test(x, y):
    pass


if __name__ == '__main__':
    p = argparse.ArgumentParser(description='')
    p.add_argument('matrix', type=argparse.FileType('r'), help='distance matrix from FastTree')
    p.add_argument('table', type=argparse.FileType('r'), help='sequence count table')
    p.add_argument('--abund', type=float, default=0.0, help='minimum fold difference for comparing two OTUs')
    p.add_argument('--dist', type=float, default=0.1, help='maximum genetic difference for comparing two OTUs')
    p.add_argument('--pval', type=float, default=0.001, help='minimum p-value for merging OTUs')
    p.add_argument('--output', '-o', default=sys.stdout, type=argparse.FileType('w'))
    p.add_argument('--save', default=None, help='save OTU table?')
    args = p.parse_args()

    assert args.abund >= 0.0

    # read in the sequences table
    seq_table = read_otu_table(args.table)

    # read in the distance matrix
    matrix = read_fasttree_matrix(args.matrix)

    # get a list of the names of the sequences in order of their (decreasing) abundance 
    seq_abunds = seq_table.sum(axis=1).sort_values(ascending=False)

    # initialize an OTU table
    otu_table = pd.DataFrame(columns=seq_table.columns)
    otu_abunds = pd.Series(index=otu_table.columns)

    for seq in seq_abunds.index:
        # which OTUs pass the genetic cutoff criteria?
        dists_to_otus = matrix[seq][otu_table.index]
        candidate_otus = dists_to_otus[dists_to_otus < args.dist].index

        # which of those OTUs also pass the abundance criterion?
        candidate_abunds = otu_abunds[candidate_otus]
        cutoff = seq_abunds[seq] * args.abund
        candidate_otus = candidate_abunds[candidate_abunds > cutoff].index

        # do any remaining candidates pass the distribution test?
        merged = False
        for otu in candidate_otus:
            abund_pval = abundance_test(otu_table.loc[otu], seq_table.loc[seq])
            if abund_pval > args.pval:
                # merge this sequence into that otu
                otu_table.loc[otu] += seq_table.loc[seq]
                otu_abunds[otu] += seq_abunds[seq]
                args.output.write('\t'.join([seq, 'match', otu]) + '\n')
                break

        if not merged:
            # this is its own OTU
            otu_table = otu_table.append(seq_table.loc[seq])
            otu_abunds = otu_abunds.append(pd.Series(seq_abunds[seq], index=[seq]))
            args.output.write('\t'.join([seq, 'otu', '']) + '\n')

        # make sure the otu abundances have remained sorted
        otu_abunds.sort_values(ascending=False, inplace=True)

    if args.save is not None:
        otu_table.to_csv(args.save, sep='\t')
