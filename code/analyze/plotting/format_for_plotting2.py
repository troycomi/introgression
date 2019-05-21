import sys
from analyze.to_update import gene_predictions
import global_params as gp


tag = sys.argv[1]
suffix = ''
if len(sys.argv == 3):
    suffix = sys.argv[2]

fn = gp.analysis_out_dir_absolute + tag + '/' + \
     'introgressed_blocks_par' + suffix + '_' + tag + '_summary_plus.txt'
region_summary = gene_predictions.read_region_summary(fn)

sep = '\t'

# ======
# for plot: lengths of all introgressed regions
# ======

# one table for each tag
# strain chrm region_length

lengths_all = []
for region in region_summary:
    length = int(region_summary[region]['end']) - \
             int(region_summary[region]['start']) + 1
