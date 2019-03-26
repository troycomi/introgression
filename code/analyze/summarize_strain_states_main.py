import sys
import itertools
from analyze import predict
from collections import defaultdict
import global_params as gp
from misc import read_table


def main():
    args = predict.process_predict_args(sys.argv[1:])

    d = defaultdict(lambda: defaultdict(int))
    outdir = gp.analysis_out_dir_absolute + args['tag']
    states = args['known_states'][1:]
    for species_from in states:

        print(species_from)

        regions1, _ = read_table.read_table_rows(
            f'{outdir}/blocks_{species_from}_'
            f'{args["tag"]}_filtered1intermediate.txt', '\t')
        regions2, _ = read_table.read_table_rows(
            f'{outdir}/blocks_{species_from}_'
            f'{args["tag"]}_filtered2intermediate.txt', '\t')

        for region_id, region1 in regions1.items():

            strain = region1['strain']
            length = int(region1['end']) - int(region1['start']) + 1

            d[strain][f'num_regions_{species_from}'] += 1
            d[strain]['num_regions_total'] += 1
            d[strain][f'num_bases_{species_from}'] += length
            d[strain]['num_bases_total'] += length

            if regions1[region_id]['reason'] != '':
                continue

            d[strain][f'num_regions_{species_from}_filtered1'] += 1
            d[strain]['num_regions_total_filtered1'] += 1
            d[strain][f'num_bases_{species_from}_filtered1'] += length
            d[strain]['num_bases_total_filtered1'] += length

            alt_states = regions2[region_id]['alternative_states'].split(',')
            for species_from_alt in alt_states:
                d[strain][f'num_regions_{species_from_alt}'
                          '_filtered2_inclusive'] += 1
                d[strain][f'num_bases_{species_from_alt}'
                          '_filtered2_inclusive'] += length
                if species_from_alt == species_from:
                    d[strain]['num_regions_total_filtered2_inclusive'] += 1
                    d[strain]['num_bases_total_filtered2_inclusive'] += length

            if len(alt_states) == 1:
                d[strain][f'num_regions_{species_from}'
                          '_filtered2'] += 1
                d[strain]['num_regions_total_filtered2'] += 1
                d[strain][f'num_bases_{species_from}'
                          '_filtered2'] += length
                d[strain]['num_bases_total_filtered2'] += length

            else:
                d[strain]['num_bases_' +
                          '_or_'.join(sorted(alt_states)) +
                          '_filtered2i'] += length

            d[strain][f'num_bases_{len(alt_states)}_filtered2i'] += length

    with open(
            '/home/tcomi/projects/aclark4_introgression/100_genomes_info.txt',
            'r') as reader:
        strain_info = [line[:-1].split('\t') for line in reader]
    strain_info = {x[0].lower(): (x[5], x[3], x[4]) for x in strain_info}

    for strain in d.keys():
        d[strain]['population'] = strain_info[strain][0]
        d[strain]['geographic_origin'] = strain_info[strain][1]
        d[strain]['environmental_origin'] = strain_info[strain][2]

    fields = ['population', 'geographic_origin', 'environmental_origin'] +\
        [f'num_regions_{x}' for x in states] +\
        ['num_regions_total'] +\
        [f'num_regions_{x}_filtered1' for x in states] +\
        ['num_regions_total_filtered1'] +\
        [f'num_regions_{x}_filtered2' for x in states] +\
        ['num_regions_total_filtered2'] +\
        [f'num_regions_{x}_filtered2_inclusive' for x in states] +\
        ['num_regions_total_filtered2_inclusive'] +\
        [f'num_bases_{x}' for x in states] +\
        ['num_bases_total'] +\
        [f'num_bases_{x}_filtered1' for x in states] +\
        ['num_bases_total_filtered1'] +\
        [f'num_bases_{x}_filtered2' for x in states] +\
        ['num_bases_total_filtered2'] +\
        [f'num_bases_{x}_filtered2_inclusive' for x in states] +\
        ['num_bases_total_filtered2_inclusive']

    r = sorted(gp.alignment_ref_order[1:])
    for n in range(2, len(r)+1):
        for combo in itertools.combinations(r, n):
            fields += ['num_bases_' + '_or_'.join(combo) + '_filtered2i']
        fields += ['num_bases_' + str(n) + '_filtered2i']

    with open(f'{outdir}/state_counts_by_strain.txt', 'w') as writer:
        writer.write('strain\t' + '\t'.join(fields) + '\n')

        for strain in sorted(d.keys()):
            writer.write(f'{strain}\t' +
                         '\t'.join([str(d[strain][x]) for x in fields]) +
                         '\n')


if __name__ == '__main__':
    main()
