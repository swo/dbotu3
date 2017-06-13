#!/usr/bin/env python3
#
# author: scott olesen <swo@alum.mit.edu>

from __future__ import print_function, division
import argparse, sys, warnings, os, datetime
import pandas as pd, numpy as np
import Levenshtein
from Bio import SeqIO
import scipy.stats

class OTU:
    '''
    Object for keeping track of an OTU's distribution and computing genetic distances
    '''
    def __init__(self, name, sequence, counts):
        '''
        name: str
          OTU ID
        sequence: str
          OTU's nucleotide sequence
        counts: numpy.Array
          length of sequence should be the number of samples
        '''
        # make this assertion so that lists of counts don't get concatenated
        self.name = name
        self.sequence = sequence
        self.counts = np.array(counts)

        self.abundance = sum(self.counts)

    def __eq__(self, other):
        return self.name == other.name and self.sequence == other.sequence and all(self.counts == other.counts)

    def __repr__(self):
        return "OTU(name={}, sequence={}, counts={})".format(repr(self.name), repr(self.sequence), repr(self.counts))

    def absorb(self, other):
        '''
        Add another OTU's counts to this one

        other: OTU

        returns: nothing
        '''
        self.counts += other.counts
        self.abundance += other.abundance

    def distance_to(self, other):
        '''
        Length-adjusted Levenshtein "distance" to other OTU

        other: OTU
          distance to this OTU

        returns: float
        '''
        ops = Levenshtein.editops(self.sequence, other.sequence)
        return len(ops) / (len(self.sequence) + len([o for o in ops if o[0] == 'delete']))

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
    def _distribution_test_pval(cls, x, y):
        '''
        P-value from the likelihood ratio test comparing the distribution of the abundances
        of two taxa (x and y). See docs for explanation of the test.
        '''
        assert len(x) == len(y)
        df = len(x) - 1
        return scipy.stats.chi2.sf(cls._D(x, y), df=df)

    def distribution_pval(self, other):
        '''
        P-value from the likelihood ratio test comparing the distribution of the abundances
        of two OTU objects. See docs for explanation of the test.

        other: OTU

        returns: float
        '''
        return self._distribution_test_pval(self.counts, other.counts)


class DBCaller:
    '''
    Object for processing the sequence table and distance matrix into an OTU table.
    '''
    def __init__(self, seq_table, records, max_dist, min_fold, threshold_pval, log=None, debug=None):
        '''
        seq_table: pandas.DataFrame
          Samples on the columns; sequences on the rows
        records: index of Bio.Seq
          Indexed, unaligned input sequences. This could come from BioPython's
          SeqIO.to_dict or SeqIO.index.
        max_dist: float
          genetic distance cutoff above which a sequence will not be merged into an OTU
        min_fold: float
          Multiply the sequence's abundance by this fold to get the minimum abundance
          of an OTU for merging
        threshold_pval: float
          P-value below which a sequence will not be merged into an OTU
        log: filehandle
          Log file reporting which sequences became OTUs and which were merged
        debug: filehandle
          Log file reporting the abundance, genetic, and distribution checks
        '''
        self.seq_table = seq_table
        self.records = records
        self.max_dist = max_dist
        self.min_fold = min_fold
        self.threshold_pval = threshold_pval
        self.progress_log = log
        self.debug_log = debug

        # get a list of the names of the sequences in order of their (decreasing) abundance
        self.seq_abunds = self.seq_table.sum(axis=1).sort_values(ascending=False)

        # check that all sequence IDs in the table are in the fasta
        missing_ids = [seq_id for seq_id in self.seq_abunds.index if seq_id not in self.records]
        if len(missing_ids) > 0:
            raise RuntimeError("{} sequence IDs found in the sequence table but not in the fasta: {}".format(len(missing_ids), missing_ids))

        # initialize OTU information
        self.membership = {}
        self.otus = []

    def _print_debug_log(self, *fields):
        '''
        Write fields to the debug log file (if present).

        returns: nothing
        '''
        if self.debug_log is not None:
            print(*fields, sep='\t', file=self.debug_log)

    def _print_progress_log(self, *fields):
        '''
        Write fields to progress log file (if present).

        returns: nothing
        '''
        if self.progress_log is not None:
            print(*fields, sep='\t', file=self.progress_log)

    def ga_matches(self, candidate):
        '''
        OTUs that meet the genetic and abundance criteria

        candidate: OTU
          sequence to evaluate

        returns: nothing
        '''

        # find abundance matches
        min_abundance = self.min_fold * candidate.abundance
        abundance_matches = [otu for otu in self.otus if otu.abundance > min_abundance]

        self._print_debug_log(candidate.name, 'abundance_check', *[otu.name for otu in abundance_matches])

        if len(abundance_matches) == 0:
            return []
        else:
            # find genetic matches (in order of decreasing genetic distance)
            matches_distances = [(otu.distance_to(candidate), otu) for otu in abundance_matches]
            matches_distances.sort(key=lambda x: (x[0], -x[1].abundance, x[1].name))
            matches = [otu for dist, otu in matches_distances if dist < self.max_dist]

            self._print_debug_log(candidate.name, 'genetic_check', *[otu.name for otu in matches])

            return matches

    def _process_record(self, record_id):
        '''
        Process the next sequence: run the genetic, abundance, and distribution checks, either
        merging the sequence into an existing OTU or creating a new OTU.

        record_id: str
          ID of the sequence to be processed. Should match fasta and seq table.

        returns: nothing
        '''
        assert record_id in self.seq_table.index
        record = self.records[record_id]

        candidate = OTU(record.id, str(record.seq), self.seq_table.loc[record.id])

        merged = False
        for otu in self.ga_matches(candidate):
            test_pval = candidate.distribution_pval(otu)

            self._print_debug_log(candidate.name, 'distribution_check', otu.name, test_pval)

            if test_pval > self.threshold_pval:
                self._merge_sequence(candidate, otu)
                merged = True
                break

        if not merged:
            self._make_otu(candidate)

    def _merge_sequence(self, member, otu):
        '''
        Merge member into OTU, making log entries if applicable.

        member, otu: OTU
          objects to be merged

        returns: nothing
        '''
        otu.absorb(member)
        self.membership[otu.name].append(member.name)

        self._print_progress_log(member.name, otu.name)
        self._print_debug_log(member.name, 'merged_into', otu.name)

    def _make_otu(self, otu):
        '''
        Make the sequence into its own OTU, recording log entries if applicable.

        otu: OTU
          sequence object to be designated as an OTU

        returns: nothing
        '''
        self.otus.append(otu)
        self.membership[otu.name] = [otu.name]

        self._print_progress_log(otu.name)
        self._print_debug_log(otu.name, 'new_otu')

    def run(self):
        '''
        Process all the input sequences in order of their abundance.

        returns: nothing
        '''

        for record_id in self.seq_abunds.index:
            self._process_record(record_id)

    def otu_table(self):
        '''
        Generate OTU table.

        returns: pandas.DataFrame
        '''
        sorted_otus = sorted(self.otus, key=lambda otu: otu.abundance, reverse=True)

        otu_table = pd.DataFrame([otu.counts for otu in sorted_otus], index=[otu.name for otu in sorted_otus])
        otu_table.columns = self.seq_table.columns

        return otu_table

    def write_otu_table(self, output):
        '''
        Write the QIIME-style OTU table to a file.

        output: filehandle

        returns: nothing
        '''

        self.otu_table().to_csv(output, sep='\t', index_label='OTU_ID')

    def write_membership(self, output):
        '''
        Write the QIIME-style OTU mapping information to a file.

        output: filehandle

        returns: nothing
        '''

        sorted_otus = sorted(self.otus, key=lambda otu: otu.abundance, reverse=True)
        for otu in sorted_otus:
            print(otu.name, *self.membership[otu.name], sep='\t', file=output)


def read_sequence_table(fn):
    '''
    Read in a table of sequences. The table must be tab-separated with exactly
    one header line of a field naming the sequences (e.g., "OTU", "OTU_ID",
    "seq", etc.) followed by tab-separated sample names. Sequence names are the
    first field of the following rows. The cells in the table are the counts of
    that sequence in that sample.

    fn: filename (or handle)

    returns: pandas.DataFrame
    '''

    # to ensure that the first column is read as a string,
    # the whole table must initially read as strings
    df = pd.read_table(fn, dtype={0: str}, header=0) # read in all columns as strings

    # BIOM format things will complain here
    if type(df.index) is pd.MultiIndex:
        warnings.warn('Table was parsed with unusual indexes. Does this table comply with the tab-separated format?', RuntimeWarning)

    # TSV formats will complain here
    if df.shape[1] == 1:
        warnings.warn('Only one tab-separated column detected. Does this table comply with the tab-separated format?', RuntimeWarning)

    df.index = df.iloc[:,0] # specify the index as the first column
    df = df.iloc[:,1:].astype(int) # cast all data columns as integers
    return df

def call_otus(seq_table_fh, fasta_fn, output_fh, gen_crit, abund_crit, pval_crit, log=None, membership=None, debug=None):
    '''
    Read in input files, call OTUs, and return output.

    seq_table_fh: filehandle
      sequence count table, tab-separated
    fasta_fn: str
      sequences fasta filename
    output_fh: filehandle
      place to write main output OTU table
    gen_crit, abund_crit, pval_crit: float
      threshold values for genetic criterion, abundance criterion, and distribution criterion (pvalue)
    log, membership, debug: filehandles
      places to write supplementary output
    '''

    # ensure valid argument values
    assert gen_crit >= 0
    assert abund_crit >= 0.0
    assert pval_crit >= 0.0 and pval_crit <= 1.0

    # read in the sequences table
    seq_table = read_sequence_table(seq_table_fh)

    # set up the input fasta records
    records = SeqIO.index(fasta_fn, 'fasta')

    # generate the caller object
    caller = DBCaller(seq_table, records, gen_crit, abund_crit, pval_crit, log, debug)

    # write the setup values to the log file (if present)
    caller._print_progress_log('---')
    caller._print_progress_log('time_started', datetime.datetime.now())
    caller._print_progress_log('genetic_criterion_threshold', gen_crit)
    caller._print_progress_log('abundance_criterion_threshold', abund_crit)
    caller._print_progress_log('distribution_criterion_threshold', pval_crit)
    caller._print_progress_log('sequence_table_filename', os.path.realpath(seq_table_fh.name))
    caller._print_progress_log('fasta_filename', os.path.realpath(fasta_fn))
    caller._print_progress_log('otu_table_output_filename', os.path.realpath(output_fh.name))
    caller._print_progress_log('progress_log_output_filename', os.path.realpath(log.name))

    if membership is not None:
        caller._print_progress_log('membership_output_filename', os.path.realpath(membership.name))

    if debug is not None:
        caller._print_progress_log('debug_log_output_filename', os.path.realpath(debug.name))

    caller._print_progress_log('---')

    # run it!
    caller.run()
    caller.write_otu_table(output_fh)

    if membership is not None:
        caller.write_membership(membership)

if __name__ == '__main__':
    p = argparse.ArgumentParser(description='dbOTU3: Distribution-based OTU calling')
    p.add_argument('table', type=argparse.FileType('r'), help='sequence count table')
    p.add_argument('fasta', help='sequences (trimmed if unpaired, merged if paired-end; always unaligned)')

    g = p.add_argument_group(title='criteria')
    g.add_argument('--dist', '-d', type=float, default=0.1, metavar='D', help='maximum genetic dissimilarity between sequences; more dissimilar sequence pairs do not pass the genetic criterion (default: 0.1)')
    g.add_argument('--abund', '-a', type=float, default=10.0, metavar='A', help='minimum fold difference for comparing two OTUs (0=no abundance criterion; default 10.0)')
    g.add_argument('--pval', '-p', type=float, default=0.0005, metavar='P', help='minimum p-value for merging OTUs (default: 0.0005)')

    g = p.add_argument_group(title='output options')
    g.add_argument('--output', '-o', default=sys.stdout, type=argparse.FileType('w'), metavar='FILE', help='OTU table output (default: stdout)')
    g.add_argument('--membership', '-m', default=None, type=argparse.FileType('w'), metavar='FILE', help='QIIME-style OTU mapping file output')
    g.add_argument('--log', '-l', default=None, type=argparse.FileType('w'), metavar='FILE', help='progress log output')
    g.add_argument('--debug', default=None, type=argparse.FileType('w'), metavar='FILE', help='debug log output')
    args = p.parse_args()

    call_otus(args.table, args.fasta, args.output, args.dist, args.abund, args.pval, log=args.log, membership=args.membership, debug=args.debug)
