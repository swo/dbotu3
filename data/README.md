# Reproducing the output
From the sequence count data (`input/counts.txt`) and input sequences
(`input/seq.fa`), you should be able
to reproduce the contents of `output/` using:

    dbotu.py input/counts.txt input/seq.fa --output output/otu.txt --membership output/membership.txt --log output/log.txt 

# The "gold standard"
The contents of `gold-standard/` are the results of the original algorithm.
