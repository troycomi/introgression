import os

# d = '/tigress/AKEY/akey_vol2/aclark4/nobackup/100_genomes/genomes_gb/orfs/'
d = '../../data/CBS432/orfs/'
fns = os.listdir(d)
for fn in fns:
    new_fn = fn[len('S288c_CBS432_'):]
    os.system('mv ' + d + fn + ' ' + d + new_fn)
