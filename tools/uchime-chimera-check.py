#!/usr/bin/env python3

'''
Progressively check for chimeric sequences in the style of UPARSE-OTU.

UPARSE-OTU reads in sequences from a database and classifies them as
(among other things) chimeric or not. The database is initialized with
the most abundant sequence. If the next-most-abundant sequence is not
classified as a match with an existing OTU or a chimera formed from
existing OTUs, it is added to the OTU database.

This script mimics this behavior by calling UPARSE-REF once per
sequence in the input. If the sequence is classified as "other",
then it is added to the OTU database used for classifying later
sequences. The results of each iteration of UPARSE-REF are output.

This script calls usearch, which must be installed and on the PATH.

author: Scott Olesen <swo@mit.edu>
'''

from Bio import SeqIO
import argparse, subprocess, tempfile, os.path, sys

def results(otus):
    database_otus = [next(otus)]

    for query in otus:
        tmpdir = tempfile.TemporaryDirectory()
        db_fn = os.path.join(tmpdir.name, 'db.fa')
        query_fn = os.path.join(tmpdir.name, 'query.fa')
        up_fn = os.path.join(tmpdir.name, 'out.up')

        SeqIO.write(database_otus, db_fn, 'fasta')
        SeqIO.write([query], query_fn, 'fasta')
        result = subprocess.run(['usearch', '-uparse_ref', query_fn, '-db', \
            db_fn, '-strand', 'plus', '-uparseout', up_fn], \
            check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

        with open(up_fn) as f:
            contents = f.read().rstrip()

        tmpdir.cleanup()

        query_label, classification, hit_pid, mp_pid, otu_label = contents.split("\t")
        assert classification in ['perfect', 'good', 'noisy', 'chimera', 'other']

        if classification == 'other':
            database_otus.append(query)

        yield contents


if __name__ == '__main__':
    p = argparse.ArgumentParser('progressive chimera checking a la UPARSE-OTU')
    p.add_argument('otus', type=argparse.FileType('r'), help='fasta of otus in abundance order')
    p.add_argument('--output', '-o', type=argparse.FileType('w'), default=sys.stdout, help='uparse output')
    args = p.parse_args()

    otus = SeqIO.parse(args.otus, 'fasta')
    for line in results(otus):
        print(line, file=args.output)
