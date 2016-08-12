# Reproducing the output
From the sequence count data (`input/counts.txt`), input sequences
(`input/seq.fa`), and k-mer distance 320 (which I computed as average sequence
length 400 * acceptable error rate 0.1 * k-mer size 8), you should be able
to reproduce the contents of `output/` using:

    dbotu.py input/counts.txt input/seq.fa 320 --log output/log.txt --output output/otu.txt --membership output/membership.txt

# The "gold standard"
The contents of `gold-standard/` are the results of the original algorithm.
