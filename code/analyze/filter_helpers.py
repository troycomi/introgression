import global_params as gp
from misc import seq_functions
import numpy as np


def write_filtered_line(f, region_id, region, fields):
    f.write(f'{region_id}\t'
            + '\t'.join([str(region[field])
                         for field in fields[1:]])
            + '\n')


def passes_filters(region):

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


def passes_filters1(region, info_string):
    # filtering out things that we can't call introgressed in general
    # with confidence (i.e. doesn't seem like a strong case against
    # being S288c)

    r = gp.alignment_ref_order[0]
    s = region['predicted_species']

    aligned_length = (int(region['end']) - int(region['start']) + 1)

    # FILTER: fraction gaps + masked
    fraction_gaps_masked_threshold = .5
    # num_sites_nonmask_x is number of sites at which neither
    # reference x nor the test sequence is masked or has a gap or
    # unsequenced character
    fraction_gaps_masked_r = \
        1 - float(region['num_sites_nonmask_' + r]) / aligned_length
    fraction_gaps_masked_s = \
        1 - float(region['num_sites_nonmask_' + s]) / aligned_length

    if fraction_gaps_masked_r > fraction_gaps_masked_threshold:
        return False, f'fraction gaps/masked in master = '\
            f'{fraction_gaps_masked_r}'
    if fraction_gaps_masked_s > fraction_gaps_masked_threshold:
        return False, f'fraction gaps/masked in predicted = '\
            f'{fraction_gaps_masked_s}'

    # FILTER: number sites analyzed by HMM that match predicted
    # reference
    count_P = info_string.count('P')
    count_C = info_string.count('C')
    number_match_only_threshold = 7
    if count_P < number_match_only_threshold:
        return False, f'count_P = {count_P}'
    if count_P <= count_C:
        return False, f'count_P = {count_P} and count_C = {count_C}'

    # FILTER: divergence with predicted reference and master reference
    # (S288c)
    id_predicted = float(region['match_nongap_' + s]) / \
        float(region['num_sites_nongap_' + s])
    id_master = float(region['match_nongap_' + r]) / \
        float(region['num_sites_nongap_' + r])

    if id_master >= id_predicted:
        return False, f'id with master = {id_master} '\
            f'and id with predicted = {id_predicted}'
    if id_master < .7:
        return False, f'id with master = {id_master}'

    return True, ''


def passes_filters2(region, seqs, threshold):
    # filter out things we can't assign to one species specifically;
    # also return the other reasonable alternatives if we're filtering
    # it out

    refs = gp.alignment_ref_order
    s = region['predicted_species']

    ids = {}
    totals = {}
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
            totals[ref] = r_total
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
