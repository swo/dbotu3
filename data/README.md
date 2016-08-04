# Reproducing the output
From the sequence count data (`input/counts.txt`), input sequences
(`input/seq.fa`), and k-mer distance 320 (which I computed as average sequence
length 400 * acceptable error rate 0.1 * k-mer size 8), you should be able
to reproduce the contents of `output/` using:

    dbotu.py data/input/counts.txt data/input/seq.fa 320 --verbose --log data/output/log --output data/output/otu.txt

# The "gold standard"
The contents of `gold-standard/` are the results of the original algorithm.
