# two levels of filtering:
# 1. remove regions that don't look confidently introgressed at all,
#    based on fraction gaps/masked, number of matches to S288c and not S288c
#    --> _filtered1
# 2. remove regions that we can't confidently pin on a specific reference,
#    based on whether it matches similarly to other reference(s)
#    --> _filtered2

# just do the first level here, then run filter_2_thresholds_main.py
# to choose filtering thresholds for next level


import sys
from analyze import predict
from analyze.filter_helpers import passes_filters1, write_filtered_line
import global_params as gp
from misc import read_table
from misc.region_reader import Region_Reader


def main():
    args = predict.process_predict_args(sys.argv[1:])
    out_dir = gp.analysis_out_dir_absolute + args['tag']

    for species_from in args['known_states'][1:]:

        print(species_from)

        region_summary, fields = read_table.read_table_rows(
            f'{out_dir}/blocks_{species_from}_{args["tag"]}_quality.txt',
            '\t')

        fields1i = fields + ['reason']
        fields1 = fields

        with open(f'{out_dir}/blocks_{species_from}_{args["tag"]}'
                  '_filtered1intermediate.txt', 'w') as f_out1i, \
                open(f'{out_dir}/blocks_{species_from}_{args["tag"]}'
                     '_filtered1.txt', 'w') as f_out1, \
                Region_Reader(f'{out_dir}/regions/{species_from}.fa.gz',
                              as_fa=True) as region_reader:

            f_out1i.write('\t'.join(fields1i) + '\n')
            f_out1.write('\t'.join(fields1) + '\n')

            for region_id, header, seqs in region_reader.yield_fa():
                region = region_summary[region_id]
                info_string = seqs[-1]
                seqs = seqs[:-1]

                # filtering stage 1: things that we're confident in calling not
                # S288c
                p, reason = passes_filters1(region, info_string)
                region['reason'] = reason
                write_filtered_line(f_out1i, region_id, region, fields1i)

                if p:
                    write_filtered_line(f_out1, region_id, region, fields1)


if __name__ == "__main__":
    main()
