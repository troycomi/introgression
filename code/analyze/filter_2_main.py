# two levels of filtering:
# 1. remove regions that don't look confidently introgressed at all,
#    based on fraction gaps/masked, number of matches to S288c and not S288c
#    --> _filtered1
# 2. remove regions that we can't confidently pin on a specific reference,
#    based on whether it matches similarly to other reference(s)
#    --> _filtered2

# do second level of filtering here, based on previously selected
# thresholds

import sys
from analyze import predict
from analyze.filter_helpers import (write_filtered_line,
                                    passes_filters2)
import global_params as gp
from misc import read_table
from misc.region_reader import Region_Reader


def main() -> None:
    '''
    Perform second stage of filtering
    Input files:
    -blocks_{species}_filtered1.txt
    regions/{species}.fa.gz
    regions/{species}.pkl

    Output files:
    -blocks_{species}_filtered2.txt
    -blocks_{species}_filtered2intermediate.txt
    '''
    args = predict.process_predict_args(sys.argv[2:])
    threshold = float(sys.argv[1])
    out_dir = gp.analysis_out_dir_absolute + args['tag']

    for species_from in args['known_states'][1:]:

        print(species_from)

        region_summary, fields = read_table.read_table_rows(
            f'{out_dir}/blocks_{species_from}_{args["tag"]}_filtered1.txt',
            '\t')

        fields2i = fields + ['alternative_states', 'alternative_ids',
                             'alternative_P_counts']
        fields2 = fields

        with open(f'{out_dir}/blocks_{species_from}_{args["tag"]}'
                  '_filtered2intermediate.txt', 'w') as f_out2i, \
                open(f'{out_dir}/blocks_{species_from}_{args["tag"]}'
                     '_filtered2.txt', 'w') as f_out2, \
                Region_Reader(f'{out_dir}/regions/{species_from}.fa.gz',
                              as_fa=True) as region_reader:

            f_out2i.write('\t'.join(fields2i) + '\n')
            f_out2.write('\t'.join(fields2) + '\n')

            for region_id, header, seqs in \
                    region_reader.yield_fa(region_summary.keys()):
                region = region_summary[region_id]

                seqs = seqs[:-1]

                # filtering stage 2: things that we're confident in calling
                # introgressed from one species specifically
                p, alt_states, alt_ids, alt_P_counts = passes_filters2(
                    region, seqs, threshold)
                region['alternative_states'] = ','.join(alt_states)
                region['alternative_ids'] = ','.join([str(x) for x in alt_ids])
                region['alternative_P_counts'] = ','.join(
                    [str(x) for x in alt_P_counts])
                write_filtered_line(f_out2i, region_id, region, fields2i)

                if p:
                    write_filtered_line(f_out2, region_id, region, fields2)


if __name__ == '__main__':
    main()
