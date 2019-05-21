from analyze import filter_regions
import pytest
from io import StringIO
from misc import read_fasta
import os
import numpy as np
import warnings
from analyze.introgression_configuration import Configuration


@pytest.fixture
def filterer():
    config = Configuration()
    config.set('symbols',
               introgressed='int_{state}.txt',
               introgressed_intermediate='int_int_{state}.txt',
               ambiguous='amb_{state}.txt',
               ambiguous_intermediate='amb_int_{state}.txt',
               filter_sweep='sweep.txt',
               filter_threshold=0.1,
               regions='region_{state}.fa.gz',
               region_index='region_{state}.pkl',
               quality_blocks='block_{state}_quality.txt')
    config.add_config({
            'analysis_params':
            {'reference': {'name': 'ref'},
             'known_states': [
                 {'name': 'pred'},
                 {'name': 'pred2'},
             ],
             }
        })
    config.set('states')
    return filter_regions.Filterer(config)


class NoCloseStringIO(StringIO):
    def close(self):
        pass

    def super_close(self):
        super(StringIO).close(self)


def test_run_no_thresh_file(filterer, mocker):
    mocker.patch('analyze.filter_regions.read_table.read_table_rows',
                 return_value=({'r1': {
                     'predicted_species': 'pred',
                     'start': 0,
                     'end': 9,
                     'num_sites_nonmask_ref': 4,
                     'num_sites_nonmask_pred': 0,
                     'match_nongap_pred': 0,
                     'num_sites_nongap_pred': 0,
                     'match_nongap_ref': 0,
                     'num_sites_nongap_ref': 0,
                 }, 'r2': {
                     'predicted_species': 'pred',
                     'start': 0,
                     'end': 9,
                     'num_sites_nonmask_ref': 5,
                     'num_sites_nonmask_pred': 5,
                     'match_nongap_pred': 8,
                     'num_sites_nongap_pred': 10,
                     'match_nongap_ref': 7,
                     'num_sites_nongap_ref': 10,
                 }, 'r3': {
                     'predicted_species': 'pred',
                     'start': 0,
                     'end': 9,
                     'num_sites_nonmask_ref': 5,
                     'num_sites_nonmask_pred': 5,
                     'match_nongap_pred': 8,
                     'num_sites_nongap_pred': 10,
                     'match_nongap_ref': 7,
                     'num_sites_nongap_ref': 10,
                 }}, ['regions']))

    files = [NoCloseStringIO() for i in range(8)]
    mocked_file = mocker.patch('analyze.filter_regions.open',
                               side_effect=files)

    mock_read = mocker.patch('analyze.filter_regions.Region_Reader')
    mock_read().__enter__().yield_fa.return_value = iter([
        ('r1', ['> seq', '> info'], ['atcg', 'x..']),
        ('r2', ['> seq', '> info'], ['attatt', 'aggcat', 'attatt',
                                     'ag-tat', np.array(list('CPPPPPPP'))]),
        ('r3', ['> seq', '> info'], ['actata', 'attatt', 'nng---',
                                     'ag-tat', np.array(list('CPPPPPPP'))])])

    filterer.config.filter_sweep = None
    filterer.run([.9])

    assert mock_read.call_count == 3  # called once during setup
    mock_read.assert_called_with('region_pred2.fa.gz', as_fa=True)

    assert mocked_file.call_args_list == [
        mocker.call('int_pred.txt', 'w'),
        mocker.call('int_int_pred.txt', 'w'),
        mocker.call('amb_pred.txt', 'w'),
        mocker.call('amb_int_pred.txt', 'w'),
        mocker.call('int_pred2.txt', 'w'),
        mocker.call('int_int_pred2.txt', 'w'),
        mocker.call('amb_pred2.txt', 'w'),
        mocker.call('amb_int_pred2.txt', 'w'),
    ]

    assert files[0].getvalue() == 'regions\nr2\t\nr3\t\n'  # pass filter 1
    assert files[1].getvalue() == (
        'regions\treason\n'
        'r1\tfraction gaps/masked in master = 0.6\n'
        'r2\t\n'
        'r3\t\n'
    )
    assert files[2].getvalue() == 'regions\nr3\t\n'  # pass filter 2
    assert files[3].getvalue() == (
        'regions\talternative_states\talternative_ids\talternative_P_counts\n'
        'r2\tpred,pred2\t1.0,1.0\t0,0\n'
        'r3\tpred\t1.0\t0\n'
    )
    # files 4:8 are just headers


def test_run_no_thresh(filterer, mocker):
    mocker.patch('analyze.filter_regions.read_table.read_table_rows',
                 return_value=({'r1': {
                     'predicted_species': 'pred',
                     'start': 0,
                     'end': 9,
                     'num_sites_nonmask_ref': 4,
                     'num_sites_nonmask_pred': 0,
                     'match_nongap_pred': 0,
                     'num_sites_nongap_pred': 0,
                     'match_nongap_ref': 0,
                     'num_sites_nongap_ref': 0,
                 }, 'r2': {
                     'predicted_species': 'pred',
                     'start': 0,
                     'end': 9,
                     'num_sites_nonmask_ref': 5,
                     'num_sites_nonmask_pred': 5,
                     'match_nongap_pred': 8,
                     'num_sites_nongap_pred': 10,
                     'match_nongap_ref': 7,
                     'num_sites_nongap_ref': 10,
                 }, 'r3': {
                     'predicted_species': 'pred',
                     'start': 0,
                     'end': 9,
                     'num_sites_nonmask_ref': 5,
                     'num_sites_nonmask_pred': 5,
                     'match_nongap_pred': 8,
                     'num_sites_nongap_pred': 10,
                     'match_nongap_ref': 7,
                     'num_sites_nongap_ref': 10,
                 }}, ['regions']))

    files = [NoCloseStringIO() for i in range(8)]
    mocked_file = mocker.patch('analyze.filter_regions.open',
                               side_effect=files)

    mock_read = mocker.patch('analyze.filter_regions.Region_Reader')
    mock_read().__enter__().yield_fa.return_value = iter([
        ('r1', ['> seq', '> info'], ['atcg', np.array(list('x..'))]),
        ('r2', ['> seq', '> info'], ['attatt', 'aggcat', 'attatt',
                                     'ag-tat', np.array(list('CPPPPPPP'))]),
        ('r3', ['> seq', '> info'], ['actata', 'attatt', 'nng---',
                                     'ag-tat', np.array(list('CPPPPPPP'))])])

    mock_log = mocker.patch('analyze.filter_regions.log')

    filterer.run()

    assert mock_log.info.call_args_list == [
        mocker.call('pred'),
        mocker.call('pred2'),
    ]

    assert mock_read.call_count == 3  # called once during setup
    mock_read.assert_called_with('region_pred2.fa.gz', as_fa=True)

    assert mocked_file.call_args_list == [
        mocker.call('int_pred.txt', 'w'),
        mocker.call('int_int_pred.txt', 'w'),
        mocker.call('amb_pred.txt', 'w'),
        mocker.call('amb_int_pred.txt', 'w'),
        mocker.call('int_pred2.txt', 'w'),
        mocker.call('int_int_pred2.txt', 'w'),
        mocker.call('amb_pred2.txt', 'w'),
        mocker.call('amb_int_pred2.txt', 'w'),
    ]

    assert files[0].getvalue() == 'regions\nr2\t\nr3\t\n'  # pass filter 1
    assert files[1].getvalue() == (
        'regions\treason\n'
        'r1\tfraction gaps/masked in master = 0.6\n'
        'r2\t\n'
        'r3\t\n'
    )
    assert files[2].getvalue() == 'regions\nr3\t\n'  # pass filter 2
    assert files[3].getvalue() == (
        'regions\talternative_states\talternative_ids\talternative_P_counts\n'
        'r2\tpred,pred2\t1.0,1.0\t0,0\n'
        'r3\tpred\t1.0\t0\n'
    )


def test_run(filterer, mocker):
    mocker.patch('analyze.filter_regions.read_table.read_table_rows',
                 return_value=({'r1': {
                     'predicted_species': 'pred',
                     'start': 0,
                     'end': 9,
                     'num_sites_nonmask_ref': 4,
                     'num_sites_nonmask_pred': 0,
                     'match_nongap_pred': 0,
                     'num_sites_nongap_pred': 0,
                     'match_nongap_ref': 0,
                     'num_sites_nongap_ref': 0,
                 }, 'r2': {
                     'predicted_species': 'pred',
                     'start': 0,
                     'end': 9,
                     'num_sites_nonmask_ref': 5,
                     'num_sites_nonmask_pred': 5,
                     'match_nongap_pred': 8,
                     'num_sites_nongap_pred': 10,
                     'match_nongap_ref': 7,
                     'num_sites_nongap_ref': 10,
                 }, 'r3': {
                     'predicted_species': 'pred',
                     'start': 0,
                     'end': 9,
                     'num_sites_nonmask_ref': 5,
                     'num_sites_nonmask_pred': 5,
                     'match_nongap_pred': 8,
                     'num_sites_nongap_pred': 10,
                     'match_nongap_ref': 7,
                     'num_sites_nongap_ref': 10,
                 }}, ['regions']))

    files = [NoCloseStringIO() for i in range(9)]
    mocked_file = mocker.patch('analyze.filter_regions.open',
                               side_effect=files)

    mock_read = mocker.patch('analyze.filter_regions.Region_Reader')
    mock_read().__enter__().yield_fa.return_value = iter([
        ('r1', ['> seq', '> info'], ['atcg', 'x..']),
        ('r2', ['> seq', '> info'], ['attatt', 'aggcat', 'attatt',
                                     'ag-tat', np.array(list('CPPPPPPP'))]),
        ('r3', ['> seq', '> info'], ['actata', 'attatt', 'nng---',
                                     'ag-tat', np.array(list('CPPPPPPP'))])])
    mock_log = mocker.patch('analyze.filter_regions.log')

    filterer.run([0.99, 0.8, 0.1])

    assert mock_log.info.call_args_list == [
        mocker.call('pred'),
        mocker.call('pred2'),
    ]

    assert mock_read.call_count == 3  # called once during setup
    mock_read.assert_called_with('region_pred2.fa.gz', as_fa=True)

    assert mocked_file.call_args_list == [
        mocker.call('sweep.txt', 'w'),
        mocker.call('int_pred.txt', 'w'),
        mocker.call('int_int_pred.txt', 'w'),
        mocker.call('amb_pred.txt', 'w'),
        mocker.call('amb_int_pred.txt', 'w'),
        mocker.call('int_pred2.txt', 'w'),
        mocker.call('int_int_pred2.txt', 'w'),
        mocker.call('amb_pred2.txt', 'w'),
        mocker.call('amb_int_pred2.txt', 'w'),
    ]

    print(files[0].getvalue())
    assert files[0].getvalue() == (
        'threshold\tpredicted_state\talternative_states\tcount\n'
        '0.99\tpred\tpred,pred2\t1\n'
        '0.99\tpred\tpred\t1\n'
        '0.8\tpred\tpred,pred2\t1\n'
        '0.8\tpred\tpred\t1\n'
        '0.1\tpred\tpred,pred2\t1\n'
        '0.1\tpred\tpred\t1\n'
    )

    assert files[1].getvalue() == 'regions\nr2\t\nr3\t\n'
    assert files[2].getvalue() == (
        'regions\treason\n'
        'r1\tfraction gaps/masked in master = 0.6\n'
        'r2\t\n'
        'r3\t\n'
    )
    assert files[3].getvalue() == 'regions\nr3\t\n'  # pass filter 2
    assert files[4].getvalue() == (
        'regions\talternative_states\talternative_ids\talternative_P_counts\n'
        'r2\tpred,pred2\t1.0,1.0\t0,0\n'
        'r3\tpred\t1.0\t0\n'
    )


def test_filter_introgressed(filterer, mocker):
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

    assert filterer.filter_introgressed(region, '', 'ref') == \
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

    assert filterer.filter_introgressed(region, '', 'ref') == \
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

    assert filterer.filter_introgressed(region,
                                        np.array(list('CP')), 'ref') == \
        (False, 'count_P = 1')
    assert filterer.filter_introgressed(region,
                                        np.array(list('CCCCCCCCPPPPPPP')),
                                        'ref') == \
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

    assert filterer.filter_introgressed(region,
                                        np.array(list('CPPPPPPP')), 'ref') == \
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

    assert filterer.filter_introgressed(region,
                                        np.array(list('CPPPPPPP')), 'ref') == \
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

    assert filterer.filter_introgressed(region,
                                        np.array(list('CPPPPPPP')), 'ref') == \
        (True, '')


def test_filter_ambiguous(filterer, mocker):
    region = {'predicted_species': '1'}
    seqs = [list('attatt'),  # reference
            list('aggcat'),  # 4 / 5, p = 2
            list('a--tta'),  # 2 / 4, p = 1
            list('nng---'),  # no matches, '3' not in outputs
            list('attatt'),  # 2 / 5, p = 0
            list('ag-tat')]  # test sequence

    threshold = 0
    filt, states = filterer.filter_ambiguous(
        region, seqs, threshold, ['ref', '1', '2', '3', '4'])
    assert filt is False
    assert region['alternative_states'] == '1,2,4'
    assert region['alternative_ids'] == '0.8,0.5,0.4'
    assert region['alternative_P_counts'] == '2,1,0'
    assert states == ['1', '2', '4']

    threshold = 0.1
    filt, _ = filterer.filter_ambiguous(
        region, seqs, threshold, ['ref', '1', '2', '3', '4'])
    assert filt is False
    assert region['alternative_states'] == '1,2'
    assert region['alternative_ids'] == '0.8,0.5'
    assert region['alternative_P_counts'] == '2,1'

    threshold = 0.9
    filt, _ = filterer.filter_ambiguous(
        region, seqs, threshold, ['ref', '1', '2', '3', '4'])
    assert filt is True
    assert region['alternative_states'] == '1'
    assert region['alternative_ids'] == '0.8'
    assert region['alternative_P_counts'] == '2'


def test_filter_ambiguous_on_region_10817(filterer, mocker):

    fa = os.path.join(os.path.split(__file__)[0], 'r10817.fa')

    if os.path.exists(fa):
        headers, seqs = read_fasta.read_fasta(fa, gz=False)
        seqs = seqs[:-1]
        region = {'predicted_species': 'CBS432'}
        p, _ = filterer.filter_ambiguous(
            region, seqs, 0.98,
            ['S288c', 'CBS432', 'N_45', 'DBVPG6304', 'UWOPS91_917_1'])
        assert p is False
        assert region['alternative_states'] == (
            'CBS432,N_45')
        assert region['alternative_P_counts'] == '111,110'

    else:
        warnings.warn('Unable to test with datafile r10817.fa')


def test_filter_ambiguous_on_region_10805(filterer, mocker):

    fa = os.path.join(os.path.split(__file__)[0], 'r10805.fa')

    if os.path.exists(fa):
        headers, seqs = read_fasta.read_fasta(fa, gz=False)
        seqs = seqs[:-1]
        region = {'predicted_species': 'N_45'}
        p, _ = filterer.filter_ambiguous(
            region, seqs, 0.1,
            ['S288c', 'CBS432', 'N_45', 'DBVPG6304', 'UWOPS91_917_1'])
        assert p is False
        assert region['alternative_states'] == (
            'CBS432,N_45,UWOPS91_917_1,DBVPG6304')
        assert region['alternative_ids'] == (
            '0.9983805668016195,0.994331983805668,'
            '0.9642857142857143,0.9618506493506493')
        assert region['alternative_P_counts'] == '145,143,128,129'

        region = {'predicted_species': 'N_45'}
        p, _ = filterer.filter_ambiguous(
            region, seqs, 0.98,
            ['S288c', 'CBS432', 'N_45', 'DBVPG6304', 'UWOPS91_917_1'])
        assert p is False
        assert region['alternative_states'] == 'CBS432,N_45'
        assert region['alternative_ids'] == (
            '0.9983805668016195,0.994331983805668')
        assert region['alternative_P_counts'] == '145,143'

    else:
        warnings.warn('Unable to test with datafile r10805.fa')


@pytest.fixture
def filter_sweep():
    return filter_regions.Filter_Sweep(None, [])


def test_filter_sweep_context(mocker):
    # no file, no list
    mock_open = mocker.patch('analyze.filter_regions.open')
    fs = filter_regions.Filter_Sweep(None, [])
    mock_open.assert_not_called()

    fs.__enter__()
    assert fs.sweep_writer is None
    mock_open.assert_not_called()
    fs.__exit__(None, None, None)
    mock_open.return_value.close.assert_not_called()

    # file, no list
    mock_open = mocker.patch('analyze.filter_regions.open')
    fs = filter_regions.Filter_Sweep('sweep.txt', [])
    mock_open.assert_not_called()

    fs.__enter__()
    assert fs.sweep_writer is None
    mock_open.assert_not_called()
    assert not fs.__exit__(None, None, 'trace')
    mock_open.return_value.close.assert_not_called()

    # file, list
    mock_open = mocker.patch('analyze.filter_regions.open')
    fs = filter_regions.Filter_Sweep('sweep.txt', [.99])
    mock_open.assert_not_called()

    fs.__enter__()
    mock_open.assert_called_once_with('sweep.txt', 'w')
    assert fs.sweep_writer is not None
    assert fs.__exit__(None, None, None)
    mock_open.return_value.close.assert_called_once()


def test_sweep_write_header(filter_sweep):
    output = StringIO()
    filter_sweep.sweep_writer = output

    filter_sweep.write_header()
    assert output.getvalue() == \
        'threshold\tpredicted_state\talternative_states\tcount\n'


def test_sweep_record(filter_sweep, mocker):
    mock_lambda = mocker.MagicMock(
        side_effect=[
            (0, ['s1']),
            (0, ['s1', 's2']),
            (0, ['s2', 's3']),
            (0, ['s4', 's3']),
        ])
    filter_sweep.thresholds = [1, 0.9, 0.8, 0.7]

    filter_sweep.record('test', mock_lambda)
    mock_lambda.assert_not_called()

    filter_sweep.sweep_writer = ''
    filter_sweep.record('test', mock_lambda)
    assert mock_lambda.call_args_list == [
        mocker.call(1),
        mocker.call(0.9),
        mocker.call(0.8),
        mocker.call(0.7)]

    assert filter_sweep.data_table == {
        1: {'test': {'s1': 1}},
        0.9: {'test': {'s1,s2': 1}},
        0.8: {'test': {'s2,s3': 1}},
        0.7: {'test': {'s3,s4': 1}},
    }


def test_sweep_write_results(filter_sweep):
    filter_sweep.data_table == {
        1: {'test': {'s1': 1}},
        0.9: {'test': {'s1,s2': 1}},
        0.8: {'test': {'s2,s3': 1}},
        0.7: {'test': {'s3,s4': 1}},
    }
    filter_sweep.thresholds = [1, 0.9, 0.8, 0.7, 0]

    filter_sweep.write_results([])

    output = StringIO()
    filter_sweep.sweep_writer = output

    filter_sweep.write_results(['state'])
    assert output.getvalue() == ''

    filter_sweep.write_results(['test'])
    assert output.getvalue() == (
        ''
    )


def test_record_data_hit(filter_sweep):
    filter_sweep.record_data_hit(0.9, 's1', ['k1'])
    assert filter_sweep.data_table == {0.9: {'s1': {'k1': 1}}}
    filter_sweep.record_data_hit(0.9, 's1', ['k1'])
    filter_sweep.record_data_hit(0.9, 's1', ['k1'])
    assert filter_sweep.data_table == {0.9: {'s1': {'k1': 3}}}
    filter_sweep.record_data_hit(0.9, 's1', ['k2'])
    assert filter_sweep.data_table == {
        0.9: {
            's1': {'k1': 3, 'k2': 1}
        }
    }
    filter_sweep.record_data_hit(0.9, 's2', ['k2'])
    assert filter_sweep.data_table == {
        0.9: {
            's1': {'k1': 3, 'k2': 1},
            's2': {'k2': 1}
        }
    }
    filter_sweep.record_data_hit(0.8, 's2', ['k2'])
    assert filter_sweep.data_table == {
        0.9: {
            's1': {'k1': 3, 'k2': 1},
            's2': {'k2': 1}
        },
        0.8: {
            's2': {'k2': 1}
        }
    }
    filter_sweep.record_data_hit(0.9, 's2', ['k2', 'k3'])
    assert filter_sweep.data_table == {
        0.9: {
            's1': {'k1': 3, 'k2': 1},
            's2': {'k2': 1, 'k2,k3': 1}
        },
        0.8: {
            's2': {'k2': 1}
        }
    }


@pytest.fixture
def filter_writer():
    config = Configuration()
    config.set(introgressed='int_{state}.txt',
               introgressed_intermediate='int_int_{state}.txt',
               ambiguous='amb_{state}.txt',
               ambiguous_intermediate='amb_int_{state}.txt')
    return filter_regions.Filter_Writers(config)


def test_filter_writer_init(filter_writer):
    assert filter_writer.files == {
        'introgressed': 'int_{state}.txt',
        'introgressed_int': 'int_int_{state}.txt',
        'ambiguous': 'amb_{state}.txt',
        'ambiguous_int': 'amb_int_{state}.txt'
    }
    assert filter_writer.writers is None
    assert filter_writer.headers is None


def test_filter_writer_context(filter_writer, mocker):
    mock_open = mocker.patch('analyze.filter_regions.open')
    with filter_writer.open_state('s1', ['h1']) as filter_writer:
        assert mock_open.call_args_list == [
            mocker.call('int_s1.txt', 'w'),
            mocker.call('int_int_s1.txt', 'w'),
            mocker.call('amb_s1.txt', 'w'),
            mocker.call('amb_int_s1.txt', 'w')]
        mock_open.return_value.close.assert_not_called()
        assert filter_writer.headers == {
            'introgressed': ['h1'],
            'introgressed_int': ['h1', 'reason'],
            'ambiguous': ['h1'],
            'ambiguous_int': ['h1', 'alternative_states',
                              'alternative_ids', 'alternative_P_counts']
        }

    assert mock_open.return_value.close.call_count == 4
    assert filter_writer.writers is None
    assert filter_writer.headers is None
    mock_open.reset_mock()

    with filter_writer.open_state('s2', ['h2']) as filter_writer:
        assert mock_open.call_args_list == [
            mocker.call('int_s2.txt', 'w'),
            mocker.call('int_int_s2.txt', 'w'),
            mocker.call('amb_s2.txt', 'w'),
            mocker.call('amb_int_s2.txt', 'w')]

        mock_open.return_value.close.assert_not_called()

        assert filter_writer.headers == {
            'introgressed': ['h2'],
            'introgressed_int': ['h2', 'reason'],
            'ambiguous': ['h2'],
            'ambiguous_int': ['h2', 'alternative_states',
                              'alternative_ids', 'alternative_P_counts']
        }

    assert mock_open.return_value.close.call_count == 4


def test_filter_writers_write_headers(filter_writer):
    filter_writer.write_headers()  # nop

    filter_writer.writers = {
        'introgressed': StringIO(),
        'introgressed_int': StringIO(),
        'ambiguous': StringIO(),
        'ambiguous_int': StringIO()
    }

    filter_writer.write_headers()  # nop

    filter_writer.headers = {
        'introgressed': ['h1'],
        'introgressed_int': ['h2', 'h3'],
        'ambiguous': ['h4'],
        'ambiguous_int': ['h5']
    }

    filter_writer.write_headers()
    assert filter_writer.writers['introgressed'].getvalue() == 'h1\n'
    assert filter_writer.writers['introgressed_int'].getvalue() == 'h2\th3\n'
    assert filter_writer.writers['ambiguous'].getvalue() == 'h4\n'
    assert filter_writer.writers['ambiguous_int'].getvalue() == 'h5\n'


def test_write_filtered_line(filter_writer):
    # single value, first field is ignored
    output = StringIO()
    filter_writer.write_filtered_line(output, 'r1', {'chr': 'I'}, ['', 'chr'])

    assert output.getvalue() == 'r1\tI\n'

    # no value
    output = StringIO()
    filter_writer.write_filtered_line(output, 'r1', {}, [])

    assert output.getvalue() == 'r1\t\n'

    # two values
    output = StringIO()
    filter_writer.write_filtered_line(output, 'r1',
                                      {'a': 'b', 'c': 'd'},
                                      ['', 'c', 'a'])

    assert output.getvalue() == 'r1\td\tb\n'
