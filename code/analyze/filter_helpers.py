import global_params as gp
from misc import seq_functions
import numpy as np
from typing import List, Dict, TextIO, Tuple


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


def passes_filters(region: Dict) -> bool:
    '''
    test if the supplied region satisfies:
    -Fraction of gaps and masked < 0.5
    -Number of matching > 7
    -Divergence < 0.7
    '''
    # fraction gaps + masked filter
    fraction_gaps_masked_threshold = .5
    fraction_gaps_masked = \
        (float(region['number_gaps']) +
         float(region['number_masked_non_gap'])) / \
        (int(region['end']) - int(region['start']) + 1)
    if fraction_gaps_masked > fraction_gaps_masked_threshold:
        return False

    # number sites match only par filter
    number_match_only_threshold = 7
    number_match_only = int(region['number_match_ref2_not_ref1'])
    if number_match_only < number_match_only_threshold:
        return False

    # divergence from cer filter (idea is that poor alignments will
    # result in much larger divergence than we'd expect)
    id_ref1_threshold = .7
    id_ref1 = float(region['number_match_ref1']) / \
        (float(region['aligned_length']) - float(region['number_gaps']))
    if id_ref1 < id_ref1_threshold:
        return False

    return True


def passes_filters1(region: Dict,
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


def passes_filters2(region: Dict,
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
