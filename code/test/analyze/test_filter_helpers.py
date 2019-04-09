from analyze import filter_helpers
from io import StringIO
import numpy as np
from misc import read_fasta
import os
import warnings
from pytest import approx


def test_write_filtered_line():
    # single value, first field is ignored
    output = StringIO()
    filter_helpers.write_filtered_line(output, 'r1', {'chr': 'I'}, ['', 'chr'])

    assert output.getvalue() == 'r1\tI\n'

    # no value
    output = StringIO()
    filter_helpers.write_filtered_line(output, 'r1', {}, [])

    assert output.getvalue() == 'r1\t\n'

    # two values
    output = StringIO()
    filter_helpers.write_filtered_line(output, 'r1',
                                       {'a': 'b', 'c': 'd'},
                                       ['', 'c', 'a'])

    assert output.getvalue() == 'r1\td\tb\n'


def test_passes_filters():
    # check gaps + number masked / end-start+1 > 0.5
    region = {'number_gaps': 1,
              'number_masked_non_gap': 0,
              'start': 0,
              'end': 1,
              'number_match_ref2_not_ref1': 0,
              'number_match_ref1': 0,
              'aligned_length': 0,
              }
    assert filter_helpers.passes_filters(region) is False
    region = {'number_gaps': 1,
              'number_masked_non_gap': 1,
              'start': 0,
              'end': 1,
              'number_match_ref2_not_ref1': 0,
              'number_match_ref1': 0,
              'aligned_length': 0,
              }
    assert filter_helpers.passes_filters(region) is False

    # check match only > 7
    region = {'number_gaps': 0,
              'number_masked_non_gap': 0,
              'start': 0,
              'end': 1,
              'number_match_ref2_not_ref1': 6,
              'number_match_ref1': 0,
              'aligned_length': 0,
              }
    assert filter_helpers.passes_filters(region) is False

    # check divergences (match_ref1 / aligned - gapped) < 0.7
    region = {'number_masked_non_gap': 0,
              'start': 0,
              'end': 1,
              'number_match_ref2_not_ref1': 7,
              'number_match_ref1': 6,
              'aligned_length': 11,
              'number_gaps': 1}
    assert filter_helpers.passes_filters(region) is False

    # passes
    region = {'number_gaps': 0,
              'number_masked_non_gap': 0,
              'start': 0,
              'end': 1,  # fraction gaps > 0.5
              'number_match_ref2_not_ref1': 7,  # >= 7
              'number_match_ref1': 7,  # div >= 0.7
              'aligned_length': 10,
              }
    assert filter_helpers.passes_filters(region) is True


def test_passes_filters1(mocker):
    # fail fraction gapped on reference
    region = {'predicted_species': 'pred',
              'start': 0,
              'end': 9,
              'num_sites_nonmask_ref': 4,
              'num_sites_nonmask_pred': 0,
              'match_nongap_pred': 0,
              'num_sites_nongap_pred': 0,
              'match_nongap_ref': 0,
              'num_sites_nongap_ref': 0,
              }

    assert filter_helpers.passes_filters1(region, '', 'ref') == \
        (False, 'fraction gaps/masked in master = 0.6')

    # fail fraction gapped on predicted
    region = {'predicted_species': 'pred',
              'start': 0,
              'end': 9,
              'num_sites_nonmask_ref': 5,
              'num_sites_nonmask_pred': 3,
              'match_nongap_pred': 0,
              'num_sites_nongap_pred': 0,
              'match_nongap_ref': 0,
              'num_sites_nongap_ref': 0,
              }

    assert filter_helpers.passes_filters1(region, '', 'ref') == \
        (False, 'fraction gaps/masked in predicted = 0.7')

    # fail match counts
    region = {'predicted_species': 'pred',
              'start': 0,
              'end': 9,
              'num_sites_nonmask_ref': 5,
              'num_sites_nonmask_pred': 5,
              'match_nongap_pred': 0,
              'num_sites_nongap_pred': 0,
              'match_nongap_ref': 0,
              'num_sites_nongap_ref': 0,
              }

    assert filter_helpers.passes_filters1(region, 'CP', 'ref') == \
        (False, 'count_P = 1')
    assert filter_helpers.passes_filters1(region,
                                          'CCCCCCCCPPPPPPP', 'ref') == \
        (False, 'count_P = 7 and count_C = 8')

    # fail divergence, master >= pred
    region = {'predicted_species': 'pred',
              'start': 0,
              'end': 9,
              'num_sites_nonmask_ref': 5,
              'num_sites_nonmask_pred': 5,
              'match_nongap_pred': 5,
              'num_sites_nongap_pred': 10,
              'match_nongap_ref': 6,
              'num_sites_nongap_ref': 10,
              }

    assert filter_helpers.passes_filters1(region, 'CPPPPPPP', 'ref') == \
        (False, 'id with master = 0.6 and id with predicted = 0.5')

    # fail divergence, master >= 0.7
    region = {'predicted_species': 'pred',
              'start': 0,
              'end': 9,
              'num_sites_nonmask_ref': 5,
              'num_sites_nonmask_pred': 5,
              'match_nongap_pred': 8,
              'num_sites_nongap_pred': 10,
              'match_nongap_ref': 6,
              'num_sites_nongap_ref': 10,
              }

    assert filter_helpers.passes_filters1(region, 'CPPPPPPP', 'ref') == \
        (False, 'id with master = 0.6')

    # passes
    region = {'predicted_species': 'pred',
              'start': 0,
              'end': 9,
              'num_sites_nonmask_ref': 5,
              'num_sites_nonmask_pred': 5,
              'match_nongap_pred': 8,
              'num_sites_nongap_pred': 10,
              'match_nongap_ref': 7,
              'num_sites_nongap_ref': 10,
              }

    assert filter_helpers.passes_filters1(region, 'CPPPPPPP', 'ref') == \
        (True, '')


def test_passes_filters2(mocker):
    mocker.patch('analyze.filter_helpers.gp.gap_symbol', '-')
    mocker.patch('analyze.filter_helpers.gp.unsequenced_symbol', 'n')

    region = {'predicted_species': '1',
              }
    seqs = [list('attatt'),  # reference
            list('aggcat'),  # 4 / 5, p = 2
            list('a--tta'),  # 2 / 4, p = 1
            list('nng---'),  # no matches, '3' not in outputs
            list('attatt'),  # 2 / 5, p = 0
            list('ag-tat')]  # test sequence

    threshold = 0
    filt, states, ids, p_count = filter_helpers.passes_filters2(
        region, seqs, threshold, ['ref', '1', '2', '3', '4'])
    assert filt is False
    assert states == ['1', '2', '4']
    assert ids == [0.8, 0.5, 0.4]
    assert p_count == [2, 1, 0]

    threshold = 0.1
    filt, states, ids, p_count = filter_helpers.passes_filters2(
        region, seqs, threshold, ['ref', '1', '2', '3', '4'])
    assert filt is False
    assert states == ['1', '2']
    assert ids == [0.8, 0.5]
    assert p_count == [2, 1]

    threshold = 0.9
    filt, states, ids, p_count = filter_helpers.passes_filters2(
        region, seqs, threshold, ['ref', '1', '2', '3', '4'])
    assert filt is True
    assert states == ['1']
    assert ids == [0.8]
    assert p_count == [2]


def test_passes_filters2_on_region(mocker):
    mocker.patch('analyze.filter_helpers.gp.gap_symbol', '-')
    mocker.patch('analyze.filter_helpers.gp.unsequenced_symbol', 'n')

    fa = os.path.join(os.path.split(__file__)[0], 'r10805.fa')

    if os.path.exists(fa):
        headers, seqs = read_fasta.read_fasta(fa, gz=False)
        seqs = seqs[:-1]
        p, alt_states, alt_ids, alt_P_counts = filter_helpers.passes_filters2(
            {'predicted_species': 'N_45'}, seqs, 0.1,
            ['S288c', 'CBS432', 'N_45', 'DBVPG6304', 'UWOPS91_917_1'])
        assert p is False
        assert alt_states == ['CBS432', 'N_45', 'UWOPS91_917_1', 'DBVPG6304']
        assert alt_ids == approx([0.9983805668016195, 0.994331983805668,
                                  0.9642857142857143, 0.9618506493506493])
        assert alt_P_counts == [145, 143, 128, 129]

        p, alt_states, alt_ids, alt_P_counts = filter_helpers.passes_filters2(
            {'predicted_species': 'N_45'}, seqs, 0.98,
            ['S288c', 'CBS432', 'N_45', 'DBVPG6304', 'UWOPS91_917_1'])
        assert p is False
        assert alt_states == ['CBS432', 'N_45']
        assert alt_ids == approx([0.9983805668016195, 0.994331983805668])
        assert alt_P_counts == [145, 143]

    else:
        warnings.warn('Unable to test with datafile r10805.fa')
