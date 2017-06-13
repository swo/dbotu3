#!/usr/bin/env python3
#
# author: scott olesen <swo@alum.mit.edu>

from __future__ import print_function

import argparse, os.path, datetime, shutil, textwrap
from Bio import SeqIO
import dbotu

REQUIRED_HEADER_KEYS = ['genetic_criterion_threshold', 'abundance_criterion_threshold',
        'distribution_criterion_threshold', 'sequence_table_filename',
        'fasta_filename', 'otu_table_output_filename', 'progress_log_output_filename']

OPTIONAL_HEADER_KEYS = ['membership_output_filename', 'debug_log_output_filename',
        'restarted_run', 'time_started', 'time_restarted']

HEADER_KEYS = REQUIRED_HEADER_KEYS + OPTIONAL_HEADER_KEYS

def backup_log(log_fn, force=False):
    target_fn = log_fn + '.bak'

    if os.path.exists(target_fn) and not force:
        raise RuntimeError('backup file destination {} already exists. delete it or re-run with --force'.format(target_fn))

    shutil.copy2(log_fn, target_fn)

def parse_log(log_fh):
    '''
    Parse a log file. Ensure that all the fields are as expected before
    proceeding.
    '''

    header = {}
    progress = []
    state = 'start'
    for line_i, line in enumerate(log_fh):
        line = line.rstrip()
        if state == 'start':
            # the first line we should see is the three dashes
            if line == '---':
                state = 'header'
            else:
                raise RuntimeError('malformed log file: expect "---" on first line')
        elif state == 'header':
            if line == '---':
                # leaving the header
                # check that header has the required keys
                missing_keys = [k for k in REQUIRED_HEADER_KEYS if k not in header]
                if len(missing_keys) > 0:
                    missing_keys_str = ', '.join('"' + k + '"' for k in missing_keys)
                    raise RuntimeError('malformed log file: header does not have required keys: ' + missing_keys_str)

                state = 'progress'
            else:
                fields = line.split('\t')
                if len(fields) == 2:
                    key, value = fields
                    # check that the key is one of the known header keys
                    if key in HEADER_KEYS:
                        header[key] = value
                    else:
                        raise RuntimeError('malformed log file: line {} (in header) has unknown key "{}"'.format(line_i + 1, key))
                else:
                    raise RuntimeError('malformed log file: line {} (in header) has {} (not 2) tab-separated fields'.format(line_i + 1, len(fields)))
        elif state == 'progress':
            fields = line.split('\t')
            if len(fields) in [1, 2]:
                progress.append((line_i + 1, fields))
            else:
                raise RuntimeError('malformed log file: line {} (in progress) has {} (not 1 or 2) tab-separated fields'.format(line_i + 1, len(fields)))

    return header, progress

def restart_dbotu_run(log_fn, drop_last=False):
    log_fn = os.path.realpath(log_fn)

    with open(log_fn) as f:
        header, progress = parse_log(f)

    # drop the last entry in the progress section
    if drop_last:
        progress.pop()

    # assert logs are the same
    if not header['progress_log_output_filename'] == log_fn:
        raise RuntimeError('log filenames do not match')

    # extract header data
    gen_crit = float(header['genetic_criterion_threshold'])
    abund_crit = float(header['abundance_criterion_threshold'])
    pval_crit = float(header['distribution_criterion_threshold'])
    seq_table_fh = open(header['sequence_table_filename'])
    fasta_fn = header['fasta_filename']
    output_fh = open(header['otu_table_output_filename'], 'w')

    if 'membership_output_filename' in header:
        membership_fh = open(header['membership_output_filename'], 'w')
    else:
        membership_fh = None

    if 'debug_log_output_filename' in header:
        debug = open(header['debug_log_output_filename'], 'w')
        print('restarted_run', file=debug)
    else:
        debug = None

    seq_table = dbotu.read_sequence_table(seq_table_fh)
    records = SeqIO.index(fasta_fn, 'fasta')

    # check that the progress items in the log are all known sequences
    for line_no, elt in progress:
        seqs_missing_from_table = [x for x in elt if x not in seq_table.index]
        seqs_missing_from_fasta = [x for x in elt if x not in records]
        if len(seqs_missing_from_table) > 0:
            raise RuntimeError('bad log file: seq(s) {} on line {} not in seq table'.format(seqs_missing_from_table, line_no))
        if len(seqs_missing_from_fasta) > 0:
            raise RuntimeError('bad log file: seq(s) {} on line {} not in fasta'.format(seqs_missing_from_fasta, line_no))

    # open the new log file
    log = open(log_fn, 'w')  # append to existing log

    # update the header and write it to the new log file
    header.update({'restarted_run': True, 'time_restarted': datetime.datetime.now()})
    print('---', file=log)
    for k, v in header.items():
        print(k, v, sep='\t', file=log)
    print('---', file=log)

    # create the caller. give it log file as None so that it doesn't write anything
    # to the log file yet
    caller = dbotu.DBCaller(seq_table, records, gen_crit, abund_crit, pval_crit, log, debug)

    for progress_i, progress_elt in enumerate(progress):
        line_no, elt = progress_elt
        if len(elt) == 1:
            # this is a new otu
            seq_id = elt[0]
            assert seq_id in caller.seq_table.index
            record = caller.records[seq_id]
            otu = dbotu.OTU(record.id, str(record.seq), caller.seq_table.loc[record.id])
            caller._make_otu(otu)
        elif len(elt) == 2:
            # this is a merging
            seq_id, otu_id = elt

            # find the otu with that name
            otu = [o for o in caller.otus if o.name == otu_id][0]

            # make the member sequence
            assert seq_id in caller.seq_table.index
            record = caller.records[seq_id]
            member = dbotu.OTU(record.id, str(record.seq), caller.seq_table.loc[record.id])

            caller._merge_sequence(member, otu)
        else:
            raise RuntimeError('progress log item on line {} is not a string or 2-item list'.format(line_no))

        # make sure that sequence we processed was the right one
        expected_id = caller.seq_abunds.index[progress_i]
        if not seq_id == expected_id:
            raise RuntimeError('the {}-th sequence already processed was "{}" but the data makes it look like it should be "{}"'.format(progress_i + 1, member_id, expected_id))

    # process the rest of the samples
    for record_id in caller.seq_abunds.index[progress_i + 1:]:
        caller._process_record(record_id)

    # write the output
    caller.write_otu_table(output_fh)

    if membership_fh is not None:
        caller.write_membership(membership_fh)

if __name__ == '__main__':
    epilog = textwrap.dedent('''This script reads the log file and writes a new
    file at that location, overwriting the old file. The default behavior is to
    also back up the original log file by copying it to a new file with the
    additional extension .bak. If a file with that .bak name already exists,
    this script will quit. If --force is on, then the .bak will be overwritten.
    If --no_backup is on, then the original log file will not be backed up. The
    --drop_last option aims to sanitize log files produced by dbOTU runs that
    were interrupted (e.g., jobs on a compute cluster that got killed) by
    dropping the last entry in the log file (i.e., the last new OTU creation or
    sequence-into-OTU merger). If the log file was truncated in the middle of
    a line after the header, then this should be equivalent to dropping the
    last line of the log file.''')

    p = argparse.ArgumentParser(description='Restart a dbOTU3 run', epilog=epilog)
    p.add_argument('log', help='progress log from stopped run')
    p.add_argument('--force', '-f', action='store_true', help='force overwrite of backup log file?')
    p.add_argument('--no_backup', '-b', action='store_false', dest='backup', help='do not create backup?')
    p.add_argument('--drop_last', '-d', action='store_true', help='drop last log progress item (new OTU or seq assignment)?')
    args = p.parse_args()

    if args.force and not args.backup:
        raise RuntimeError("incompatible command-line arguments: can't force overwrite if not writing a backup")

    if args.backup:
        backup_log(args.log, args.force)

    restart_dbotu_run(args.log, drop_last=args.drop_last)
