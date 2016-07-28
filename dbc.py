#!/usr/bin/env python3
#
# author: scott olesen <swo@mit.edu>

import argparse, sys
import pandas as pd, numpy as np
import scipy.stats

def read_otu_table(fn):
    return pd.read_table(fn, index_col=0, header=0).astype(int)

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

def D_helper(x):
    x = np.array(x)
    x = x[x > 0]
    return np.sum(x * np.log(x)) - (np.sum(x) * np.log(np.sum(x)))

def D(x, y):
    x = np.array(x)
    y = np.array(y)
    return -2.0 * (D_helper(x + y) - D_helper(x) - D_helper(y))

def abundance_test_pval(x, y):
    assert len(x) == len(y)
    df = len(x) - 1
    return scipy.stats.chi2.sf(D(x, y), df=df)


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
        # get a list of OTUs in order of increasing distance
        dists_to_otus = matrix[seq][otu_table.index]
        g_candidates = dists_to_otus[dists_to_otus < args.dist].sort_values(ascending=True).index
        print('genetic_check', seq, *g_candidates, sep='\t')

        # which of those OTUs also pass the abundance criterion?
        candidate_abunds = otu_abunds[g_candidates]
        cutoff = seq_abunds[seq] * args.abund
        ga_candidates = candidate_abunds[candidate_abunds > cutoff].index
        print('abundance_check', seq, *ga_candidates, sep='\t')

        # do any remaining candidates pass the distribution test?
        merged = False
        for otu in ga_candidates:
            test_pval = abundance_test_pval(otu_table.loc[otu], seq_table.loc[seq])
            print('distribution_check', seq, otu, test_pval, sep='\t')
            if test_pval > args.pval:
                # merge this sequence into that otu
                #args.output.write('> {} had dist: {}\n'.format(otu, otu_table.loc[otu]))
                otu_table.loc[otu] += seq_table.loc[seq]
                otu_abunds[otu] += seq_abunds[seq]
                print('match', seq, otu, sep='\t')
                args.output.write('\t'.join([seq, 'match', otu]) + '\n')
                #args.output.write('> {} has dist: {}\n'.format(seq, seq_table.loc[seq]))
                #args.output.write('> {} now has dist: {}\n'.format(otu, otu_table.loc[otu]))
                #args.output.write('> pval: {}\n'.format(test_pval))
                merged = True
                break

        if not merged:
            # this is its own OTU
            otu_table = otu_table.append(seq_table.loc[seq])
            otu_abunds = otu_abunds.append(pd.Series(seq_abunds[seq], index=[seq]))
            print('otu', seq, sep='\t')
            args.output.write('\t'.join([seq, 'otu', '']) + '\n')
            #args.output.write('> {} has dist: {}\n'.format(seq, seq_table.loc[seq]))

        # make sure the otu abundances have remained sorted
        otu_abunds.sort_values(ascending=False, inplace=True)

    if args.save is not None:
        otu_table.to_csv(args.save, sep='\t')
