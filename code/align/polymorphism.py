# calculate polymorphism rate between reference genomes in 100-bp
# windows across each chromosome

import sys
from misc import read_fasta

headers, seqs = read_fasta.read_fasta(sys.argv[1])
a = dict(zip(headers, seqs))
nseqs = len(seqs)
nsites = len(seqs[0])
