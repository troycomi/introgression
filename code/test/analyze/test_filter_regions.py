from analyze import filter_regions
from io import StringIO
from misc import read_fasta
import os
import warnings
from pytest import approx


def test_main_no_thresh(mocker, capsys):
    mocker.patch('sys.argv', ['', '0.1'])
    mocker.patch('analyze.filter_regions.predict.process_predict_args',
                 return_value={
                     'known_states': ['state1', 'state2'],
                     'tag': 'tag'
                 })
    mocker.patch('analyze.filter_regions.gp.analysis_out_dir_absolute',
                 '/dir')
    mocker.patch('analyze.filter_regions.read_table.read_table_rows',
                 return_value=({'r1': {}, 'r2': {'a': 1}}, ['regions']))
    mocked_file = mocker.patch('analyze.filter_regions.open')

    mock_read = mocker.patch('analyze.filter_regions.Region_Reader')
    mock_read().__enter__().yield_fa.return_value = iter([
        ('r1', ['> seq', '> info'], ['atcg', 'x..']),
        ('r2', ['> seq', '> info'], ['atcg', 'x..'])])

    mock_filter1 = mocker.patch('analyze.filter_regions.filter_introgressed',
                                side_effect=[(False, 'test'),  # r1
                                             (True, '')])  # r2
    mock_filter2 = mocker.patch(
        'analyze.filter_regions.filter_ambiguous',
        side_effect=[
            (True, ['1'], [0.8], [2]),
            (False, ['1', '2'], [0.8, 0.5], [2, 1, 0])
        ])
    mock_write = mocker.patch('analyze.filter_regions.write_filtered_line')

    filter_regions.main()

    captured = capsys.readouterr().out
    assert captured == 'state2\n'

    assert mock_read.call_count == 2  # called once during setup
    mock_read.assert_called_with('/dirtag/regions/state2.fa.gz', as_fa=True)

    assert mocked_file.call_args_list == [
        mocker.call('/dirtag/blocks_state2_tag_filtered1intermediate.txt',
                    'w'),
        mocker.call('/dirtag/blocks_state2_tag_filtered1.txt', 'w'),
        mocker.call('/dirtag/blocks_state2_tag_filtered2intermediate.txt',
                    'w'),
        mocker.call('/dirtag/blocks_state2_tag_filtered2.txt', 'w'),
    ]

    # just headers, capture others
    assert mocked_file().__enter__().write.call_args_list == [
        mocker.call('regions\treason\n'),
        mocker.call('regions\n'),
        mocker.call('regions\talternative_states\t'
                    'alternative_ids\talternative_P_counts\n'),
        mocker.call('regions\n'),
    ]

    assert mock_filter1.call_count == 2
    # seems like this references the object, which changes after call
    assert mock_filter1.call_args_list == [
        mocker.call({'reason': 'test'}, 'x..', 'state1'),
        mocker.call({'a': 1, 'reason': '',
                     'alternative_states': '1',
                     'alternative_ids': '0.8',
                     'alternative_P_counts': '2'}, 'x..', 'state1')
    ]

    assert mock_filter2.call_args_list == [
        mocker.call({'a': 1, 'reason': '',
                     'alternative_states': '1',
                     'alternative_ids': '0.8',
                     'alternative_P_counts': '2'},
                    ['atcg'], 0.1, ['state1', 'state2']),
    ]
    assert mock_write.call_args_list == [
        mocker.call(mocker.ANY, 'r1', {'reason': 'test'},
                    ['regions', 'reason']),
        mocker.call(mocker.ANY, 'r2', {'a': 1, 'reason': '',
                                       'alternative_states': '1',
                                       'alternative_ids': '0.8',
                                       'alternative_P_counts': '2'
                                       },
                    ['regions', 'reason']),
        mocker.call(mocker.ANY, 'r2', {'a': 1, 'reason': '',
                                       'alternative_states': '1',
                                       'alternative_ids': '0.8',
                                       'alternative_P_counts': '2'
                                       },
                    ['regions']),
        mocker.call(mocker.ANY, 'r2', {'a': 1, 'reason': '',
                                       'alternative_states': '1',
                                       'alternative_ids': '0.8',
                                       'alternative_P_counts': '2'
                                       },
                    ['regions', 'alternative_states', 'alternative_ids',
                     'alternative_P_counts']),
        mocker.call(mocker.ANY, 'r2', {'a': 1, 'reason': '',
                                       'alternative_states': '1',
                                       'alternative_ids': '0.8',
                                       'alternative_P_counts': '2'
                                       },
                    ['regions']),
    ]


def test_main(mocker, capsys):
    mocker.patch('sys.argv', ['', '0.1'])
    mocker.patch('analyze.filter_regions.predict.process_predict_args',
                 return_value={
                     'known_states': ['state1', 'state2'],
                     'tag': 'tag'
                 })
    mocker.patch('analyze.filter_regions.gp.analysis_out_dir_absolute',
                 '/dir')
    mocker.patch('analyze.filter_regions.read_table.read_table_rows',
                 return_value=({'r1': {}, 'r2': {'a': 1}}, ['regions']))
    mocked_file = mocker.patch('analyze.filter_regions.open')

    mock_read = mocker.patch('analyze.filter_regions.Region_Reader')
    mock_read().__enter__().yield_fa.return_value = iter([
        ('r1', ['> seq', '> info'], ['atcg', 'x..']),
        ('r2', ['> seq', '> info'], ['atcg', 'x..'])])

    mock_filter1 = mocker.patch('analyze.filter_regions.filter_introgressed',
                                side_effect=[(False, 'test'),  # r1
                                             (True, '')])  # r2
    mock_filter2 = mocker.patch(
        'analyze.filter_regions.filter_ambiguous',
        side_effect=[
            (False, ['1', '2'], [0.8, 0.5], [2, 1, 0]),
            (True, ['1'], [0.8], [2]),
            (True, ['1'], [0.8], [2]),
            (False, ['1', '2'], [0.8, 0.5], [2, 1, 0])
        ])
    mock_write = mocker.patch('analyze.filter_regions.write_filtered_line')

    filter_regions.main([0.99])

    captured = capsys.readouterr().out
    assert captured == 'state2\n'

    assert mock_read.call_count == 2  # called once during setup
    mock_read.assert_called_with('/dirtag/regions/state2.fa.gz', as_fa=True)

    assert mocked_file.call_count == 5
    assert mocked_file.call_args_list == [
        mocker.call('/dirtag/filter_2_thresholds_tag.txt', 'w'),
        mocker.call('/dirtag/blocks_state2_tag_filtered1intermediate.txt',
                    'w'),
        mocker.call('/dirtag/blocks_state2_tag_filtered1.txt', 'w'),
        mocker.call('/dirtag/blocks_state2_tag_filtered2intermediate.txt',
                    'w'),
        mocker.call('/dirtag/blocks_state2_tag_filtered2.txt', 'w'),
    ]

    # just headers, capture others
    assert mocked_file().__enter__().write.call_args_list == [
        mocker.call('threshold\tpredicted_state\talternative_states\tcount\n'),
        mocker.call('regions\treason\n'),
        mocker.call('regions\n'),
        mocker.call('regions\talternative_states\t'
                    'alternative_ids\talternative_P_counts\n'),
        mocker.call('regions\n'),
        mocker.call('0.99\tstate2\t1,2\t1\n')
    ]

    assert mock_filter1.call_count == 2
    # seems like this references the object, which changes after call
    assert mock_filter1.call_args_list == [
        mocker.call({'reason': 'test'}, 'x..', 'state1'),
        mocker.call({'a': 1, 'reason': '',
                     'alternative_states': '1',
                     'alternative_ids': '0.8',
                     'alternative_P_counts': '2'}, 'x..', 'state1')
    ]

    assert mock_filter2.call_args_list == [
        mocker.call({'a': 1, 'reason': '',
                     'alternative_states': '1',
                     'alternative_ids': '0.8',
                     'alternative_P_counts': '2'},
                    ['atcg'], 0.99, ['state1', 'state2']),
        mocker.call({'a': 1, 'reason': '',
                     'alternative_states': '1',
                     'alternative_ids': '0.8',
                     'alternative_P_counts': '2'},
                    ['atcg'], 0.1, ['state1', 'state2']),
    ]
    assert mock_write.call_args_list == [
        mocker.call(mocker.ANY, 'r1', {'reason': 'test'},
                    ['regions', 'reason']),
        mocker.call(mocker.ANY, 'r2', {'a': 1, 'reason': '',
                                       'alternative_states': '1',
                                       'alternative_ids': '0.8',
                                       'alternative_P_counts': '2'
                                       },
                    ['regions', 'reason']),
        mocker.call(mocker.ANY, 'r2', {'a': 1, 'reason': '',
                                       'alternative_states': '1',
                                       'alternative_ids': '0.8',
                                       'alternative_P_counts': '2'
                                       },
                    ['regions']),
        mocker.call(mocker.ANY, 'r2', {'a': 1, 'reason': '',
                                       'alternative_states': '1',
                                       'alternative_ids': '0.8',
                                       'alternative_P_counts': '2'
                                       },
                    ['regions', 'alternative_states', 'alternative_ids',
                     'alternative_P_counts']),
        mocker.call(mocker.ANY, 'r2', {'a': 1, 'reason': '',
                                       'alternative_states': '1',
                                       'alternative_ids': '0.8',
                                       'alternative_P_counts': '2'
                                       },
                    ['regions']),
    ]


def test_record_data_hit():
    dt = {}
    filter_regions.record_data_hit(dt, 0.9, 's1', 'k1')
    assert dt == {0.9: {'s1': {'k1': 1}}}
    filter_regions.record_data_hit(dt, 0.9, 's1', 'k1')
    filter_regions.record_data_hit(dt, 0.9, 's1', 'k1')
    assert dt == {0.9: {'s1': {'k1': 3}}}
    filter_regions.record_data_hit(dt, 0.9, 's1', 'k2')
    assert dt == {
        0.9: {
            's1': {'k1': 3, 'k2': 1}
        }
    }
    filter_regions.record_data_hit(dt, 0.9, 's2', 'k2')
    assert dt == {
        0.9: {
            's1': {'k1': 3, 'k2': 1},
            's2': {'k2': 1}
        }
    }
    filter_regions.record_data_hit(dt, 0.8, 's2', 'k2')
    assert dt == {
        0.9: {
            's1': {'k1': 3, 'k2': 1},
            's2': {'k2': 1}
        },
        0.8: {
            's2': {'k2': 1}
        }
    }
    filter_regions.record_data_hit(dt, 0.9, 's2', 'k2')
    assert dt == {
        0.9: {
            's1': {'k1': 3, 'k2': 1},
            's2': {'k2': 2}
        },
        0.8: {
            's2': {'k2': 1}
        }
    }


def test_write_filtered_line():
    # single value, first field is ignored
    output = StringIO()
    filter_regions.write_filtered_line(output, 'r1', {'chr': 'I'}, ['', 'chr'])

    assert output.getvalue() == 'r1\tI\n'

    # no value
    output = StringIO()
    filter_regions.write_filtered_line(output, 'r1', {}, [])

    assert output.getvalue() == 'r1\t\n'

    # two values
    output = StringIO()
    filter_regions.write_filtered_line(output, 'r1',
                                       {'a': 'b', 'c': 'd'},
                                       ['', 'c', 'a'])

    assert output.getvalue() == 'r1\td\tb\n'


def test_filter_introgressed(mocker):
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

    assert filter_regions.filter_introgressed(region, '', 'ref') == \
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

    assert filter_regions.filter_introgressed(region, '', 'ref') == \
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

    assert filter_regions.filter_introgressed(region, 'CP', 'ref') == \
        (False, 'count_P = 1')
    assert filter_regions.filter_introgressed(region,
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

    assert filter_regions.filter_introgressed(region, 'CPPPPPPP', 'ref') == \
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

    assert filter_regions.filter_introgressed(region, 'CPPPPPPP', 'ref') == \
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

    assert filter_regions.filter_introgressed(region, 'CPPPPPPP', 'ref') == \
        (True, '')


def test_filter_ambiguous(mocker):
    mocker.patch('analyze.filter_regions.gp.gap_symbol', '-')
    mocker.patch('analyze.filter_regions.gp.unsequenced_symbol', 'n')

    region = {'predicted_species': '1',
              }
    seqs = [list('attatt'),  # reference
            list('aggcat'),  # 4 / 5, p = 2
            list('a--tta'),  # 2 / 4, p = 1
            list('nng---'),  # no matches, '3' not in outputs
            list('attatt'),  # 2 / 5, p = 0
            list('ag-tat')]  # test sequence

    threshold = 0
    filt, states, ids, p_count = filter_regions.filter_ambiguous(
        region, seqs, threshold, ['ref', '1', '2', '3', '4'])
    assert filt is False
    assert states == ['1', '2', '4']
    assert ids == [0.8, 0.5, 0.4]
    assert p_count == [2, 1, 0]

    threshold = 0.1
    filt, states, ids, p_count = filter_regions.filter_ambiguous(
        region, seqs, threshold, ['ref', '1', '2', '3', '4'])
    assert filt is False
    assert states == ['1', '2']
    assert ids == [0.8, 0.5]
    assert p_count == [2, 1]

    threshold = 0.9
    filt, states, ids, p_count = filter_regions.filter_ambiguous(
        region, seqs, threshold, ['ref', '1', '2', '3', '4'])
    assert filt is True
    assert states == ['1']
    assert ids == [0.8]
    assert p_count == [2]


def test_filter_ambiguous_on_region(mocker):
    mocker.patch('analyze.filter_regions.gp.gap_symbol', '-')
    mocker.patch('analyze.filter_regions.gp.unsequenced_symbol', 'n')

    fa = os.path.join(os.path.split(__file__)[0], 'r10805.fa')

    if os.path.exists(fa):
        headers, seqs = read_fasta.read_fasta(fa, gz=False)
        seqs = seqs[:-1]
        p, alt_states, alt_ids, alt_P_counts = filter_regions.filter_ambiguous(
            {'predicted_species': 'N_45'}, seqs, 0.1,
            ['S288c', 'CBS432', 'N_45', 'DBVPG6304', 'UWOPS91_917_1'])
        assert p is False
        assert alt_states == ['CBS432', 'N_45', 'UWOPS91_917_1', 'DBVPG6304']
        assert alt_ids == approx([0.9983805668016195, 0.994331983805668,
                                  0.9642857142857143, 0.9618506493506493])
        assert alt_P_counts == [145, 143, 128, 129]

        p, alt_states, alt_ids, alt_P_counts = filter_regions.filter_ambiguous(
            {'predicted_species': 'N_45'}, seqs, 0.98,
            ['S288c', 'CBS432', 'N_45', 'DBVPG6304', 'UWOPS91_917_1'])
        assert p is False
        assert alt_states == ['CBS432', 'N_45']
        assert alt_ids == approx([0.9983805668016195, 0.994331983805668])
        assert alt_P_counts == [145, 143]

    else:
        warnings.warn('Unable to test with datafile r10805.fa')
