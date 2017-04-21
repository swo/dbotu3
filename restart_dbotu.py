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

    # create the caller. give it log file as None so that it doesn't write anything
    # to the log file yet
    caller = dbotu.DBCaller(seq_table, records, gen_crit, abund_crit, pval_crit, log=None, debug=debug)

    for elt_i, elt in enumerate(progress):
        if isinstance(elt, str):
            # this is a new otu
            seq_id = elt
            assert seq_id in caller.seq_table.index
            record = caller.records[seq_id]
            otu = dbotu.OTU(record.id, str(record.seq), caller.seq_table.loc[record.id])
            caller._make_otu(otu)
            print(elt_i, 'new otu', seq_id)
        elif isinstance(elt, list) and len(elt) == 2:
            # this is a merging
            seq_id, otu_id = elt

            # find the otu with that name
            otu = [o for o in caller.otus if o.name == otu_id][0]

            # make the member sequence
            assert seq_id in caller.seq_table.index
            record = caller.records[seq_id]
            member = dbotu.OTU(record.id, str(record.seq), caller.seq_table.loc[record.id])

            caller._merge_sequence(member, otu)
            print(elt_i, 'merge', seq_id, otu_id)
        else:
            raise RuntimeError('progress log item "{}" is not a string or 2-item list'.format(elt))

        # make sure that sequence we processed was the right one
        expected_id = caller.seq_abunds.index[elt_i]
        if not seq_id == expected_id:
            raise RuntimeError('the {}-th sequence already processed was "{}" but the data makes it look like it should be "{}"'.format(elt_i + 1, member_id, expected_id))

    # now that we are processing real data, the logging should start again
    # swo> or, read in the original log, then create a new one with a new header, updated with the start time
    caller.progress_log = log

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
