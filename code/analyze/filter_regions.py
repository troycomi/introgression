import global_params as gp
from misc import seq_functions
import numpy as np
from typing import List, Dict, TextIO, Tuple
import sys
from contextlib import ExitStack
from analyze import predict
from misc import read_table
from misc.region_reader import Region_Reader


def write_filtered_line(writer: TextIO,
                        region_id: str,
                        region: Dict,
                        fields: List) -> None:
    '''
    Write the region id and values in "region" dict to open file writer
    '''
    writer.write(f'{region_id}\t'
                 + '\t'.join([str(region[field])
                              for field in fields[1:]])
                 + '\n')


def filter_introgressed(region: Dict,
                        info: str,
                        reference_species: str) -> Tuple[bool, str]:
    '''
    filtering out things that we can't call introgressed in general
    with confidence (i.e. doesn't seem like a strong case against
    being S288c)
    Return true if the region passes the filter, or false with a string
    specifying which filter failed
    Tests:
    -fraction of gaps masked in reference > 0.5
    -fraction of gaps masked in predicted species > 0.5
    -number of matches to predicted > 7
    -number of matches to predicted > number matches to reference
    -divergence with predicted species
    '''

    predicted_species = region['predicted_species']

    aligned_length = (int(region['end']) - int(region['start']) + 1)

    # FILTER: fraction gaps + masked
    fraction_gaps_masked_threshold = .5
    # num_sites_nonmask_x is number of sites at which neither
    # reference x nor the test sequence is masked or has a gap or
    # unsequenced character
    fraction_gaps_masked_r = \
        1 - region['num_sites_nonmask_' + reference_species] / aligned_length
    fraction_gaps_masked_s = \
        1 - region['num_sites_nonmask_' + predicted_species] / aligned_length

    if fraction_gaps_masked_r > fraction_gaps_masked_threshold:
        return False, f'fraction gaps/masked in master = '\
            f'{fraction_gaps_masked_r}'
    if fraction_gaps_masked_s > fraction_gaps_masked_threshold:
        return False, f'fraction gaps/masked in predicted = '\
            f'{fraction_gaps_masked_s}'

    # FILTER: number sites analyzed by HMM that match predicted (P)
    # reference (C)
    count_P = info.count('P')
    count_C = info.count('C')
    number_match_only_threshold = 7
    if count_P < number_match_only_threshold:
        return False, f'count_P = {count_P}'
    if count_P <= count_C:
        return False, f'count_P = {count_P} and count_C = {count_C}'

    # FILTER: divergence with predicted reference and master reference
    # (S288c)
    id_predicted = float(region['match_nongap_' + predicted_species]) / \
        float(region['num_sites_nongap_' + predicted_species])
    id_master = float(region['match_nongap_' + reference_species]) / \
        float(region['num_sites_nongap_' + reference_species])

    if id_master >= id_predicted:
        return False, f'id with master = {id_master} '\
            f'and id with predicted = {id_predicted}'
    if id_master < .7:
        return False, f'id with master = {id_master}'

    return True, ''


def filter_ambiguous(region: Dict,
                     seqs: np.array,
                     threshold: float,
                     refs: List[str]) -> Tuple[bool,
                                               List[str],
                                               List[float],
                                               List[int]]:
    '''
    filter out things we can't assign to one species specifically;
    return the other reasonable alternatives if we're filtering
    it out
    Returns a tuple of:
    True if the region passes the filter
    A list of likely species for the region
    A list of fraction of matching sequence for each species
    A list of total matching sites
    Fails the filter if number of matches and fraction matching are >= more
    than one state for the region
    '''

    s = region['predicted_species']

    ids = {}
    P_counts = {}

    seqs = np.asarray(seqs)
    # skip any gap or unsequenced in ref or test
    # also skip if ref and test equal (later test ri == test but not ref)
    skip = np.any(
        (seqs[0] == gp.gap_symbol,
         seqs[0] == gp.unsequenced_symbol,
         seqs[-1] == gp.gap_symbol,
         seqs[-1] == gp.unsequenced_symbol,
         seqs[0] == seqs[-1]),
        axis=0)

    for ri, ref in enumerate(refs):
        if ri == 0:
            continue
        r_match, r_total = seq_functions.seq_id(seqs[-1], seqs[ri])
        if r_total != 0:
            ids[ref] = r_match / r_total
            P_counts[ref] = np.sum(
                np.logical_and(
                    np.logical_not(skip),
                    seqs[ri] == seqs[-1]))

    alts = {}
    for r in ids.keys():
        # TODO should threshold be the same for both?
        if ids[r] >= threshold * ids[s] and \
           P_counts[r] >= threshold * P_counts[s]:
            alts[r] = (ids[r], P_counts[r])

    alt_states = sorted(alts.keys(), key=lambda x: alts[x][0], reverse=True)
    alt_ids = [alts[state][0] for state in alt_states]
    alt_P_counts = [alts[state][1] for state in alt_states]

    if len(alts) > 1:
        return False, alt_states, alt_ids, alt_P_counts

    return True, alt_states, alt_ids, alt_P_counts


def main(thresholds=[]):
    '''
    Perform first step of filtering
    Input files:
    -blocks_{species}_quality.txt

    Output files:
    -blocks_{species}_filtered1intermediate.txt
    -blocks_{species}_filtered1.txt
    -regions/{species}.fa.gz
    -regions/{species}.pkl
    '''
    # thresholds = [.999, .995, .985, .975, .965, .955, .945,
    #               .935, .925, .915, .905, .89, .87, .86]
    args = predict.process_predict_args(sys.argv[2:])
    out_dir = gp.analysis_out_dir_absolute + args['tag']
    threshold = float(sys.argv[1])

    with ExitStack() as stack:
        if thresholds != []:
            threshold_writer = stack.enter_context(
                open(f'{out_dir}/filter_2_thresholds_{args["tag"]}.txt', 'w'))
            threshold_writer.write(
                'threshold\tpredicted_state\talternative_states\tcount\n')

        data_table = {}

        for species_from in args['known_states'][1:]:

            print(species_from)

            region_summary, fields = read_table.read_table_rows(
                f'{out_dir}/blocks_{species_from}_{args["tag"]}_quality.txt',
                '\t')

            fields1i = fields + ['reason']
            fields1 = fields
            fields2i = fields + ['alternative_states', 'alternative_ids',
                                 'alternative_P_counts']
            fields2 = fields

            with open(f'{out_dir}/blocks_{species_from}_{args["tag"]}'
                      '_filtered1intermediate.txt', 'w') as f_out1i, \
                    open(f'{out_dir}/blocks_{species_from}_{args["tag"]}'
                         '_filtered1.txt', 'w') as f_out1, \
                    open(f'{out_dir}/blocks_{species_from}_{args["tag"]}'
                         '_filtered2intermediate.txt', 'w') as f_out2i, \
                    open(f'{out_dir}/blocks_{species_from}_{args["tag"]}'
                         '_filtered2.txt', 'w') as f_out2, \
                    Region_Reader(f'{out_dir}/regions/{species_from}.fa.gz',
                                  as_fa=True) as region_reader:

                f_out1i.write('\t'.join(fields1i) + '\n')
                f_out1.write('\t'.join(fields1) + '\n')
                f_out2i.write('\t'.join(fields2i) + '\n')
                f_out2.write('\t'.join(fields2) + '\n')

                for region_id, header, seqs in region_reader.yield_fa():
                    region = region_summary[region_id]
                    info_string = seqs[-1]
                    seqs = seqs[:-1]

                    # filtering stage 1: things that we're confident in
                    # calling not S288c
                    p, reason = filter_introgressed(region,
                                                    info_string,
                                                    args['known_states'][0])
                    region['reason'] = reason
                    write_filtered_line(f_out1i, region_id, region, fields1i)

                    if p:
                        write_filtered_line(f_out1, region_id, region, fields1)

                        for thresh in thresholds:
                            _, alt_states, _, _ = \
                                filter_ambiguous(region, seqs, thresh,
                                                 args['known_states'])

                            record_data_hit(data_table,
                                            thresh,
                                            species_from,
                                            ','.join(sorted(alt_states)))

                        (p, alt_states,
                         alt_ids, alt_P_counts) = filter_ambiguous(
                            region, seqs, threshold, args['known_states'])
                        region['alternative_states'] = ','.join(alt_states)
                        region['alternative_ids'] = ','.join(
                            [str(x) for x in alt_ids])
                        region['alternative_P_counts'] = ','.join(
                            [str(x) for x in alt_P_counts])
                        write_filtered_line(f_out2i, region_id,
                                            region, fields2i)

                        if p:
                            write_filtered_line(f_out2, region_id,
                                                region, fields2)

        for thresh in thresholds:
            for species in args['known_states'][1:]:
                d = data_table[thresh][species]
                for key in d.keys():
                    threshold_writer.write(
                        f'{thresh}\t{species}\t{key}\t{d[key]}\n')


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
