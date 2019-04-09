import gzip
import numpy as np
from typing import Tuple, List


def read_fasta(fn: str, gz: bool = False) -> Tuple[
        List[str], np.array]:
    '''
    Read the provided fasta file, returning the
    headers (lines startin with >) and sequences
    '''

    headers = []
    seqs = []

    f = None
    if gz:
        f = gzip.open(fn, 'rb')
    else:
        f = open(fn, 'r')
    line = f.readline()

    while line[0] != '>':
        line = f.readline()

    while True:
        h = line[:-1]
        s = []
        line = f.readline()
        while line != '' and line[0] != '>':
            s += list(line[:-1])
            line = f.readline()
        headers.append(h)
        seqs.append(s)
        if line == '':
            break
    f.close()

    return headers, np.asarray(seqs)
