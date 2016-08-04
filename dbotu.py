#!/usr/bin/env python3
#
# author: scott olesen <swo@mit.edu>

import argparse, sys
import pandas as pd, numpy as np
from Bio import SeqIO
import scipy.stats

import kmer

class DBCaller:
    '''
    Object for processing the sequence table and distance matrix into an OTU table.
    '''
    def __init__(self, seq_table, records, k, max_dist, min_fold, threshold_pval, log=None, verbose=False):
        '''
        seq_table: pandas.DataFrame
          Samples on the columns; sequences on the rows
        records: sequence of Bio.Seq
          unaligned input sequences
        k: int
          k for the k-mers
        max_dist: float
          kmer distance cutoff above which a sequence will not be merged into an OTU
        min_fold: float
          Multiply the sequence's abundance by this fold to get the minimum abundance
          of an OTU for merging
        threshold_pval: float
          P-value below which a sequence will not be merged into an OTU
        log: filehandle
          Log file reporting how sequences are merged or placed into their own OTUs.
          Tab-delimited, one line per sequence. If the sequence goes into its own
          OTU, the line is sequence-"otu". If the sequence goes into an existing OTU,
          the line is sequence-"match"-otu.
        verbose: bool
          If true, put extra lines into the log file. For every genetic threshold
          check, a line "genetic_check"-sequence-[OTUs that meet the genetic threshold].
          For every abundance check, a line "abundance_check"-sequence-[OTUs that meet].
          For every distribution test, a line "distribution_check"-sequence-OTU-pvalue.
        '''
        self.seq_table = seq_table
        self.records = records
        self.k = k
        self.max_dist = max_dist
        self.min_fold = min_fold
        self.threshold_pval = threshold_pval
        self.log = log
        self.verbose = verbose

        if self.verbose and self.log is None:
            raise RuntimeError("verbose option requires a log")

        # get a list of the names of the sequences in order of their (decreasing) abundance
        self.seq_abunds = self.seq_table.sum(axis=1).sort_values(ascending=False)

        # initialize the OTU information
        self.otu_table = pd.DataFrame(columns=self.seq_table.columns)
        self.otu_abunds = pd.Series(index=self.otu_table.columns)
        self.otu_kmers = {}

    def ga_candidates(self, record):
        '''OTUs that meet the genetic and abundance criteria'''
        #dists_to_otus = self.matrix[seq][self.otu_table.index]
        seq_kmers = kmer.kmer_dict(str(record.seq), self.k)
        dists_to_otus = pd.Series({otu: kmer.distance(seq_kmers, self.otu_kmers[otu]) for otu in self.otu_table.index})
        candidates1 = dists_to_otus[dists_to_otus < self.max_dist].sort_values(ascending=True).index

        for c in candidates1:
            print('genetic dist from', c, dists_to_otus[c], file=self.log)

        if self.verbose:
            print('genetic_check', *candidates1, sep='\t', file=self.log)

        candidate_abunds = self.otu_abunds[candidates1]
        cutoff = self.seq_abunds[record.id] * self.min_fold
        candidates2 = candidate_abunds[candidate_abunds > cutoff].index

        if self.verbose and len(candidates1) > 0:
            print('abundance_check', *candidates2, sep='\t', file=self.log)

        return candidates2

    def merge_seq_into_otu(self, seq_id, otu):
        '''Merge sequence into OTU, adjusting the OTU abundance'''
        self.otu_table.loc[otu] += self.seq_table.loc[seq_id]
        self.otu_abunds[otu] += self.seq_abunds[seq_id]

    def create_otu(self, record):
        '''
        Create a new OTU by starting a new line in the OTU table,
        creating an OTU abundance entry, and ensuring the the OTU
        abundances remain sorted.
        '''
        self.otu_table = self.otu_table.append(self.seq_table.loc[record.id])
        self.otu_abunds = self.otu_abunds.append(pd.Series(self.seq_abunds[record.id], index=[record.id]))
        self.otu_abunds.sort_values(ascending=False, inplace=True)
        self.otu_kmers[record.id] = kmer.kmer_dict(str(record.seq), self.k)

    @staticmethod
    def _D_helper(x):
        '''A helper function for the _D method'''
        x = np.array(x)
        x = x[x > 0]
        return np.sum(x * np.log(x)) - (np.sum(x) * np.log(np.sum(x)))

    @classmethod
    def _D(cls, x, y):
        '''
        Statistic for the likelihood ratio test. See docs for mathematical derivation.
        '''
        x = np.array(x)
        y = np.array(y)
        return -2.0 * (cls._D_helper(x + y) - cls._D_helper(x) - cls._D_helper(y))

    @classmethod
    def abundance_test_pval(cls, x, y):
        '''
        P-value from the likelihood ratio test comparing the distribution of the abundances
        of two taxa (x and y). See docs for explanation of the test.
        '''
        assert len(x) == len(y)
        df = len(x) - 1
        return scipy.stats.chi2.sf(cls._D(x, y), df=df)

    def _process_record(self, record):
        '''
        Process the next sequence: run the genetic, abundance, and distribution checks, either
        merging the sequence into an existing OTU or creating a new OTU.
        '''
        assert record.id in self.seq_table.index

        if self.verbose:
            print('seq', record.id, sep='\t', file=self.log)

        for otu in self.ga_candidates(record):
            test_pval = self.abundance_test_pval(self.otu_table.loc[otu], self.seq_table.loc[record.id])

            if self.verbose:
                print('distribution_check', otu, test_pval, sep='\t', file=args.log)

            if test_pval > self.threshold_pval:
                if self.log is not None:
                    print('match', record.id, otu, sep='\t', file=args.log)

                self.merge_seq_into_otu(record.id, otu)
                return otu

        # form own otu
        if self.log is not None:
            print('otu', record.id, sep='\t', file=args.log)

        self.create_otu(record)
        return record

    def generate_otu_table(self):
        '''
        Process all the input sequences to make an OTU table.

        returns: pandas.DataFrame
          OTU table (which can also be found at instance.otu_table)
        '''
        for record in self.records:
            self._process_record(record)

        return self.otu_table


def read_sequence_table(fn):
    '''
    Read in a table of sequences. Expect a header and the sequence IDs in the
    first column. Samples are on the columns.

    fn: filename (or handle)

    returns: pandas.DataFrame
    '''
    return pd.read_table(fn, index_col=0, header=0).astype(int)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description='', )
    p.add_argument('table', type=argparse.FileType('r'), help='sequence count table')
    p.add_argument('fasta', type=argparse.FileType('r'), help='sequences (unaligned)')
    p.add_argument('dist', type=float, help='maximum kmer difference for comparing two OTUs (recommended: kmer size * maximum number of acceptable mismatches)')
    p.add_argument('--word_length', '-k', type=int, default=8, help='kmer size (default: 8)')
    p.add_argument('--abund', '-a', type=float, default=10.0, help='minimum fold difference for comparing two OTUs (0=no abundance criterion; default 10.0)')
    p.add_argument('--pval', '-p', type=float, default=0.0005, help='minimum p-value for merging OTUs (default: 0.0005)')
    p.add_argument('--log', '-l', default=None, type=argparse.FileType('w'), help='log output')
    p.add_argument('--verbose', '-v', action='store_true', help='record checks in log?')
    p.add_argument('--output', '-o', default=sys.stdout, help='OTU table output (default: stdout)')
    args = p.parse_args()

    assert args.dist >= 0
    assert args.abund >= 0.0
    assert args.dist >= 0.0
    assert args.pval >= 0.0 and args.pval <= 1.0

    # read in the sequences table
    seq_table = read_sequence_table(args.table)

    # set up the input fasta records
    records = SeqIO.parse(args.fasta, 'fasta')

    # generate the caller object
    caller = DBCaller(seq_table, records, args.word_length, args.dist, args.abund, args.pval, args.log, args.verbose)
    caller.generate_otu_table()

    caller.otu_table.to_csv(args.output, sep='\t')
