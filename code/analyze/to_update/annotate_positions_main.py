# columns:
# position in strain
# position in cer ref
# gene
# in ORF?

import sys
import os
import gzip
from annotate_positions import (get_genes, get_orfs, write_annotated_file)
import global_params as gp
from align import align_helpers

# ======
# get strains
# ======

i = int(sys.argv[1])
s = align_helpers.get_strains(align_helpers.flatten(gp.non_ref_dirs.values()))
strain, d = s[i]

# ======
# get genes on each chromosome
# ======

genes_by_chrm = {}
for chrm in gp.chrms:
    fn = gp.analysis_out_dir_absolute + gp.master_ref + '_chr' + chrm + \
         '_genes.txt'
    genes_by_chrm[chrm] = get_genes(fn)

# ======
# loop through all strains and chromosomes, generating annotated
# position file for each
# ======

coord_dir = gp.analysis_out_dir_absolute + 'coordinates/'
if not os.path.exists(coord_dir + 'annotated'):
    os.makedirs(coord_dir + 'annotated')

for chrm in gp.chrms:

    print(strain, chrm)

    fn = strain + '_to_' + gp.master_ref + '_chr' + chrm + '.txt.gz'

    fn_orfs = d + 'orfs/' + strain + '_chr' + chrm + \
        '_orfs' + gp.fasta_suffix
    orfs = get_orfs(fn_orfs)

    fn_out = coord_dir + 'annotated/' + fn
    coords = [float(line)
              for line in gzip.open(coord_dir + fn, 'rb').readlines()]
    write_annotated_file(coords, genes_by_chrm[chrm], orfs, fn_out)
