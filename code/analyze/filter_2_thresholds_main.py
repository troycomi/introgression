# explore different thresholds for calling introgressions for specific
# strains

# specifically, try a range of thresholds, and for each one, calculate
# fraction of introgressions we've classified as 1 strain or every
# possible combination of strains

# then we'll make some plots in R to see if there's a sort of obvious
# place to draw the line

import sys
from analyze import predict
from analyze.filter_helpers import passes_filters2
import global_params as gp
from misc import read_table
from misc.region_reader import Region_Reader


thresholds = [.999, .995, .985, .975, .965, .955, .945,
              .935, .925, .915, .905, .89, .87, .86]
# thresholds = [.99, .98, .97, .96, .95, .94, .93, .92,
#               .91, .9, .88, .85, .82, .8, .75, .7, .6, .5]
# thresholds = [1]


def main():
    args = predict.process_predict_args(sys.argv[1:])
    out_dir = gp.analysis_out_dir_absolute + args['tag']

    open_mode = 'w'
    with open(f'{out_dir}/filter_2_thresholds_{args["tag"]}.txt', open_mode)\
            as writer:
        if open_mode == 'w':
            writer.write(
                'threshold\tpredicted_state\talternative_states\tcount\n')

        data_table = {}
        for species_from in args['known_states'][1:]:
            print(f'* {species_from}')

            region_summary, fields = read_table.read_table_rows(
                f'{out_dir}/blocks_{species_from}'
                f'_{args["tag"]}_filtered1.txt',
                '\t')

            with Region_Reader(f'{out_dir}/regions/{species_from}.fa.gz',
                               as_fa=True) as region_reader:
                for region_id, header, seqs in \
                        region_reader.yield_fa(region_summary.keys()):

                    region = region_summary[region_id]
                    seqs = seqs[:-1]

                    for threshold in thresholds:
                        _, alt_states, _, _ = \
                            passes_filters2(region, seqs, threshold)

                        record_data_hit(data_table,
                                        threshold,
                                        species_from,
                                        ','.join(sorted(alt_states)))

        for threshold in thresholds:
            for species in args['known_states'][1:]:
                d = data_table[threshold][species]
                for key in d.keys():
                    writer.write(f'{threshold}\t{species}\t{key}\t{d[key]}\n')


def record_data_hit(data_dict, threshold, species, key):
    '''
    adds an entry to the data table or increments if exists
    '''
    if threshold not in data_dict:
        data_dict[threshold] = {}

    if species not in data_dict[threshold]:
        data_dict[threshold][species] = {}

    if key not in data_dict[threshold][species]:
        data_dict[threshold][species][key] = 0

    data_dict[threshold][species][key] += 1


if __name__ == "__main__":
    main()
