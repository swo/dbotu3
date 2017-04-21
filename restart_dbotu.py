#!/usr/bin/env python3
#
# author: scott olesen <swo@alum.mit.edu>

from __future__ import print_function

import argparse, yaml, os.path
from Bio import SeqIO
import dbotu

def restart_dbotu_run(log_fn):
    log_fn = os.path.realpath(log_fn)

    with open(log_fn) as f:
        header, progress = yaml.load_all(f)

    # assert logs are the same
    if not header['progress_log_output_filename'] == log_fn:
        raise RuntimeError('log filenames do not match')

    # extract header data
    gen_crit = header['genetic_criterion_threshold']
    abund_crit = header['abundance_criterion_threshold']
    pval_crit = header['distribution_criterion_threshold']
    seq_table_fh = open(header['sequence_table_filename'])
    fasta_fn = header['fasta_filename']
    output_fh = open(header['otu_table_output_filename'], 'w')
    log = open(log_fn, 'a')  # append to existing log

    if 'membership_output_filename' in header:
        membership_fh = open(header['membership_output_filename'], 'w')
    else:
        membership_fh = None

    if 'debug_log_output_filename' in header:
        debug = open(header['debug_log_output_filename'], 'a')  # append to existing debug log
    else:
        debug = None

    seq_table = dbotu.read_sequence_table(seq_table_fh)
    records = SeqIO.index(fasta_fn, 'fasta')
    caller = dbotu.DBCaller(seq_table, records, gen_crit, abund_crit, pval_crit, log, debug)

    # assert that the already-processed sequences are in a block at the top of the caller's to-do list
    # swo> sort the progress keys by their abundance to maintain order
    # swo> or, don't use a dictionary for the log, since the ordering does matter
    already_processed_seq_idx = sorted([caller.seq_abunds.index.get_loc(seq) for seq in progress.keys()])
    n_already_processed = len(already_processed_seq_idx)
    if not already_processed_seq_idx == list(range(n_already_processed)):
        raise RuntimeError('{} sequences were already processed, but they are not 0 through {}, instead {}'.format(n_already_processed, n_already_processed-1, already_processed_seq_idx))

    # make all the otus
    for seq_id, otu_id in progress.items():
        if seq_id == otu_id:
            assert otu_id in caller.seq_table.index
            record = caller.records[otu_id]
            otu = dbotu.OTU(record.id, str(record.seq), caller.seq_table.loc[record.id])
            caller._make_otu(otu)

    # merge all the sequences
    for seq_id, otu_id in progress.items():
        if seq_id != otu_id:
            # find the otu with that name
            otu = [o for o in caller.otus if o.name == otu_id][0]

            # make the member sequence
            assert seq_id in caller.seq_table.index
            record = caller.records[seq_id]
            member = dbotu.OTU(record.id, str(record.seq), caller.seq_table.loc[record.id])

            caller._merge_sequence(member, otu)

    # process the rest of the samples
    for record_id in caller.seq_abunds.index[n_already_processed:]:
        caller._process_record(record_id)

    # write the output
    caller.write_otu_table(output_fh)

    if membership_fh is not None:
        caller.write_membership(membership_fh)

if __name__ == '__main__':
    p = argparse.ArgumentParser(description='Restart a dbOTU3 run')
    p.add_argument('log', help='progress log from stopped run')
    args = p.parse_args()

    restart_dbotu_run(args.log)
