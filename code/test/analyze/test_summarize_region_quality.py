import analyze.summarize_region_quality as summarize
from io import StringIO
import pytest
from pytest import approx
import numpy as np
from numpy.testing import assert_array_equal as aae
from analyze.introgression_configuration import Configuration


@pytest.fixture
def summarizer():
    return summarize.Summarizer(Configuration())


def test_states_to_process(summarizer, mocker):
    summarizer.config.add_config({
        'analysis_params':
        {'reference': {'name': 'S288c'},
         'known_states': [
             {'name': 'CBS432'},
             {'name': 'N_45'},
         ],
         'unknown_states': [{'name': 'unknown'}]
         }
    })
    summarizer.config.set('states')

    assert summarizer.states_to_process() == \
        (0, 'S288c CBS432 N_45 unknown'.split())

    mock_warn = mocker.patch('analyze.summarize_region_quality.log.warning')
    assert summarizer.states_to_process('N_45 asdf S288c'.split()) == \
        (0, 'N_45 S288c'.split())
    mock_warn.assert_called_with("state 'asdf' was not found as a state")

    with pytest.raises(ValueError) as e:
        summarizer.states_to_process('asdf qwer'.split())
    assert 'No valid states were found to process' in str(e)

    summarizer.config.add_config({
        'analysis_params': {'reference': {'name': 'N_45'}}})

    assert summarizer.states_to_process() == \
        (2, 'S288c CBS432 N_45 unknown'.split())


def test_run(summarizer, mocker):
    summarizer.config.add_config({
        'analysis_params':
        {'reference': {'name': 'S288c'},
         'known_states': [
             {'name': 'CBS432'},
             {'name': 'N_45'},
             {'name': 'DBVPG6304'},
             {'name': 'UWOPS91_917_1'}
         ],
         'unknown_states': [{'name': 'unknown'}]
         }
    })
    summarizer.config.set('symbols', 'states',
                          positions='positions.txt.gz',
                          labeled_blocks='dir/tag/blocks_{state}_labeled.txt',
                          quality_blocks='dir/tag/blocks_{state}_quality.txt',
                          alignment='dir/tag/blocks_{chrom}_{strain}.txt',
                          regions='dir/tag/regions/{state}.fa.gz',
                          region_index='dir/tag/regions/{state}.pkl',
                          masks='dir/masked/{strain}_chr{chrom}.txt')
    summarizer.config.chromosomes = ['I', 'II']
    summarizer.validate_arguments()
    # for region database
    mock_table = mocker.patch(
        'misc.read_table.open',
        mocker.mock_open(
            read_data='region_id\tstrain\tchromosome\t'
            'predicted_species\tstart\tend\tnum_sites_hmm\n'
            'r4\tyjm1381\tI\tS288c\t2\t5\t60\n'
            'r5\tyjm689\tI\tS288c\t3\t6\t56\n'
            'r6\tyjm1381\tI\tS288c\t3\t7\t18\n'
            'r7\tyjm689\tI\tS288c\t3\t5\t13728\n'
            'r8\tyjm1208\tI\tS288c\t3\t4\t20\n'
            'r9\tyjm1304\tII\tS288c\t3\t7\t16\n'
        ))

    # sequence analyzer masked sites
    mock_masked = mocker.patch.object(
        summarize.Sequence_Analyzer,
        'read_masked_intervals',
        return_value=[
            (0, 2),
            (4, 5),
        ])
    position_in = StringIO(
        'yjm1381\tI\t2\t3\t5\t7\n'
        'yjm689\tI\t3\t4\t5\t6\n'
        'yjm1464\tI\t1\t2\t3\t3\n'
    )
    region_out = StringIO()

    def new_close():
        pass

    mocker.patch.object(region_out, 'close', new_close)

    mocked_gzip = mocker.patch(
        'analyze.summarize_region_quality.gzip.open',
        side_effect=[position_in, region_out])

    mocked_file = mocker.patch('analyze.summarize_region_quality.open',
                               mocker.mock_open())
    mock_log = mocker.patch('analyze.summarize_region_quality.log')

    mocker.patch('analyze.summarize_region_quality.read_fasta.read_fasta',
                 return_value=('',
                               np.asarray([
                                   list('--gatcctag--'),
                                   list('-agatgcaag-c'),
                                   list('-agatgcaag-c'),
                                   list('-agatgcaag-c'),
                                   list('-a-attacagt-'),
                                   list('-a-atttcagt-'),
                               ])))

    summarizer.run(['unknown'])

    mock_masked.assert_any_call('dir/masked/UWOPS91_917_1_chrII.txt')
    mock_masked.assert_any_call('dir/masked/yjm1381_chrI.txt')

    assert mocked_gzip.call_args_list == [
        mocker.call('positions.txt.gz', 'rt'),
        mocker.call('dir/tag/regions/unknown.fa.gz', 'wt'),
    ]
    assert mock_log.debug.call_args_list == [
          mocker.call('reference index: 0'),
          mocker.call("states to analyze: ['unknown']"),
          mocker.call("known_states ['S288c', 'CBS432', 'N_45', "
                      "'DBVPG6304', 'UWOPS91_917_1']"),
          mocker.call('Sequence_Analyzer init with:'),
          mocker.call('masks: dir/masked/{strain}_chr{chrom}.txt'),
          mocker.call('alignment: dir/tag/blocks_{chrom}_{strain}.txt'),
          mocker.call('yjm1381 I'),
          mocker.call('yjm689 I')
    ]
    assert mock_log.info.call_args_list == [
        mocker.call('Working on state unknown'),
        mocker.call('Working on chromosome I'),
        mocker.call('Working on chromosome II'),
    ]

    assert mocked_file.call_count == 2
    mocked_file.assert_any_call(
        'dir/tag/blocks_unknown_quality.txt', 'w')
    mocked_file.assert_any_call(
        'dir/tag/regions/unknown.pkl', 'wb')

    mock_table.assert_any_call(
        'dir/tag/blocks_unknown_labeled.txt', 'r')

    # just headers
    states = ['S288c', 'CBS432', 'N_45', 'DBVPG6304', 'UWOPS91_917_1']
    symbols = list('.-_npbcxNPBCX')
    assert mocked_file().write.call_args_list == [
        mocker.call(
            '\t'.join(
                ('region_id\tstrain\tchromosome\tpredicted_species\tstart'
                 '\tend\tnum_sites_hmm').split() +
                ['match_nongap_' + x for x in states] +
                ['num_sites_nongap_' + x for x in states] +
                ['match_hmm_' + x for x in states] +
                ['match_nonmask_' + x for x in states] +
                ['num_sites_nonmask_' + x for x in states] +
                ['count_' + x for x in symbols]
            ) + '\n'),

        mocker.call('r4\tyjm1381\tI\tS288c\t2\t5\t3\t1\t1\t1\t1\t3\t4\t4\t4\t4'
                    '\t4\t1\t1\t1\t1\t3\t0\t0\t0\t0\t1\t1\t0\t0\t0\t1\t0\t0'
                    '\t4\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n'),
        mocker.call('r5\tyjm689\tI\tS288c\t3\t6\t4\t1\t1\t1\t1\t3\t4\t4\t4\t4'
                    '\t4\t1\t1\t1\t1\t3\t1\t1\t1\t1\t2\t2\t1\t1\t1\t2\t1\t0'
                    '\t3\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n'),
        mocker.call('r6\tyjm1381\tI\tS288c\t3\t7\t3\t2\t2\t2\t2\t4\t5\t5\t5'
                    '\t5\t5\t1\t1\t1\t1\t3\t2\t2\t2\t2\t3\t3\t2\t2\t2\t3\t2'
                    '\t0\t3\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n'),
        mocker.call('r7\tyjm689\tI\tS288c\t3\t5\t3\t0\t0\t0\t0\t2\t3\t3\t3\t3'
                    '\t3\t0\t0\t0\t0\t2\t0\t0\t0\t0\t1\t1\t0\t0\t0\t1\t0\t0'
                    '\t3\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n'),
        mocker.call('r8\tyjm1208\tI\tS288c\t3\t4\t20\t0\t0\t0\t0\t0\t0\t0\t0'
                    '\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0'
                    '\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n'),
        mocker.call('r9\tyjm1304\tII\tS288c\t3\t7\t16\t0\t0\t0\t0\t0\t0\t0\t0'
                    '\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0'
                    '\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n'),
        mocker.ANY
    ]

    region_output = region_out.getvalue()
    region_out.close()

    assert region_output == (
        '#r4\n'
        '> S288c 2 5\ntcct\n'
        '> CBS432 3 6\ntgca\n'
        '> N_45 3 6\ntgca\n'
        '> DBVPG6304 3 6\ntgca\n'
        '> UWOPS91_917_1 2 5\nttac\n'
        '> yjm1381 2 5\ntttc\n'
        '> info\n____\n'
        '#r6\n'
        '> S288c 3 7\ncctag\n'
        '> CBS432 4 8\ngcaag\n'
        '> N_45 4 8\ngcaag\n'
        '> DBVPG6304 4 8\ngcaag\n'
        '> UWOPS91_917_1 3 7\ntacag\n'
        '> yjm1381 3 7\nttcag\n'
        '> info\n___..\n'
        '#r5\n'
        '> S288c 3 6\nccta\n'
        '> CBS432 4 7\ngcaa\n'
        '> N_45 4 7\ngcaa\n'
        '> DBVPG6304 4 7\ngcaa\n'
        '> UWOPS91_917_1 3 6\ntaca\n'
        '> yjm689 3 6\nttca\n'
        '> info\n___.\n'
        '#r7\n'
        '> S288c 3 5\ncct\n'
        '> CBS432 4 6\ngca\n'
        '> N_45 4 6\ngca\n'
        '> DBVPG6304 4 6\ngca\n'
        '> UWOPS91_917_1 3 5\ntac\n'
        '> yjm689 3 5\nttc\n'
        '> info\n___\n'
    )


def test_run_all_states(summarizer, mocker):
    summarizer.config.add_config({
        'analysis_params':
        {'reference': {'name': 'S288c'},
         'known_states': [
             {'name': 'CBS432'},
             {'name': 'N_45'},
             {'name': 'DBVPG6304'},
             {'name': 'UWOPS91_917_1'}
         ],
         'unknown_states': [{'name': 'unknown'}]
         }
    })
    summarizer.config.set('states',
                          'symbols',
                          positions='positions.txt.gz',
                          labeled_blocks='dir/tag/blocks_{state}_labeled.txt',
                          quality_blocks='dir/tag/blocks_{state}_quality.txt',
                          alignment='dir/tag/blocks_{chrom}_{strain}.txt',
                          regions='dir/tag/regions/{state}.fa.gz',
                          region_index='dir/tag/regions/{state}.pkl',
                          masks='dir/masked/{strain}_chr{chrom}.txt')
    summarizer.config.chromosomes = ['I', 'II']
    assert summarizer.validate_arguments()

    mock_log = mocker.patch('analyze.summarize_region_quality.log')

    with pytest.raises(FileNotFoundError) as e:
        summarizer.run()
    assert "No such file or directory: 'dir/masked/S288c_chrI.txt'" in str(e)

    assert mock_log.debug.call_args_list == [
        mocker.call('reference index: 0'),
        mocker.call("states to analyze: ['S288c', 'CBS432', 'N_45', "
                    "'DBVPG6304', 'UWOPS91_917_1', 'unknown']"),
        mocker.call("known_states ['S288c', 'CBS432', 'N_45', 'DBVPG6304', "
                    "'UWOPS91_917_1']"),
        mocker.call('Sequence_Analyzer init with:'),
        mocker.call('masks: dir/masked/{strain}_chr{chrom}.txt'),
        mocker.call('alignment: dir/tag/blocks_{chrom}_{strain}.txt')]


@pytest.fixture
def flag():
    return summarize.Flag_Info()


def test_flag_info_init(flag):
    assert flag.__dict__ == {
        'gap_any': None,
        'mask_any': None,
        'unseq_any': None,
        'hmm': None,
        'gap': None,
        'mask': None,
        'unseq': None,
        'match': None
    }


def test_intialize_flags(flag):
    flag.initialize_flags(3, 2)

    aae(flag.gap_any, np.zeros((3), bool))
    aae(flag.mask_any, np.zeros((3), bool))
    aae(flag.unseq_any, np.zeros((3), bool))
    aae(flag.gap, np.zeros((3, 2), bool))
    aae(flag.mask, np.zeros((3, 2), bool))
    aae(flag.unseq, np.zeros((3, 2), bool))
    aae(flag.match, np.zeros((3, 2), bool))

    flag.initialize_flags(5, 3)

    aae(flag.gap_any, np.zeros((5), bool))
    aae(flag.mask_any, np.zeros((5), bool))
    aae(flag.unseq_any, np.zeros((5), bool))
    aae(flag.gap, np.zeros((5, 3), bool))
    aae(flag.mask, np.zeros((5, 3), bool))
    aae(flag.unseq, np.zeros((5, 3), bool))
    aae(flag.match, np.zeros((5, 3), bool))


def test_add_sequence_flags(flag):
    flag.initialize_flags(3, 2)

    other = summarize.Flag_Info()
    other.gap = np.array([0, 0, 1], bool)
    other.unseq = np.array([1, 0, 0], bool)
    other.hmm = np.array([1, 0, 1], bool)
    other.match = np.array([0, 1, 1], bool)

    flag.add_sequence_flags(other, 0)
    aae(flag.hmm, np.array([1, 0, 1], bool))
    aae(flag.gap_any, np.array([0, 0, 1], bool))
    aae(flag.unseq_any, np.array([1, 0, 0], bool))

    aae(flag.gap, np.array(
        [
            [0, 0],
            [0, 0],
            [1, 0],
        ], bool))

    aae(flag.unseq, np.array(
        [
            [1, 0],
            [0, 0],
            [0, 0],
        ], bool))

    aae(flag.match, np.array(
        [
            [0, 0],
            [1, 0],
            [1, 0],
        ], bool))

    other = summarize.Flag_Info()
    other.gap = np.array([1, 0, 0], bool)
    other.unseq = np.array([0, 0, 1], bool)
    other.hmm = np.array([0, 1, 0], bool)
    other.match = np.array([0, 1, 1], bool)

    flag.add_sequence_flags(other, 1)
    aae(flag.hmm, np.array([1, 0, 1], bool))
    aae(flag.gap_any, np.array([1, 0, 1], bool))
    aae(flag.unseq_any, np.array([1, 0, 1], bool))

    aae(flag.gap, np.array(
        [
            [0, 1],
            [0, 0],
            [1, 0],
        ], bool))

    aae(flag.unseq, np.array(
        [
            [1, 0],
            [0, 0],
            [0, 1],
        ], bool))

    aae(flag.match, np.array(
        [
            [0, 0],
            [1, 1],
            [1, 1],
        ], bool))


def test_add_mask_flags(flag):
    flag.initialize_flags(3, 2)

    other = summarize.Flag_Info()
    other.mask = np.array([0, 0, 1], bool)

    flag.add_mask_flags(other, 0)
    aae(flag.mask_any, np.array([0, 0, 1], bool))

    aae(flag.mask, np.array(
        [
            [0, 0],
            [0, 0],
            [1, 0],
        ], bool))

    other = summarize.Flag_Info()
    other.mask = np.array([0, 1, 0], bool)

    flag.add_mask_flags(other, 1)
    aae(flag.mask_any, np.array([0, 1, 1], bool))

    aae(flag.mask, np.array(
        [
            [0, 0],
            [0, 1],
            [1, 0],
        ], bool))


def test_encode_info(flag):
    flag.initialize_flags(14, 3)
    flag.hmm = np.zeros((14), bool)

    flag.gap[0, 0] = True  # -
    flag.gap[11, 1] = True  # -
    flag.mask[1, 0] = True  # _
    flag.mask[12, 1] = True  # _
    flag.match[2, :] = True  # .
    flag.match[13, :] = True  # .
    flag.match[(3, 4, 5, 6), 0] = True  # b and c
    flag.match[(3, 4, 7, 8), 1] = True  # b and p
    # x is default
    flag.hmm[[4, 6, 8, 10]] = True  # capitalize

    s = flag.encode_info(master_ind=0, predict_ind=1)
    assert s == '-_.bBcCpPxX-_.'
    #            01234567890123

    flag.initialize_flags(0, 3)
    flag.hmm = np.zeros((0), bool)

    s = flag.encode_info(master_ind=0, predict_ind=1)
    assert s == ''


def test_encode_unknown_info(flag):
    flag.initialize_flags(5, 2)

    flag.gap_any[0] = True  # -
    flag.mask_any[1] = True  # _
    flag.match[2, :] = True  # .
    flag.match[3, 0] = True  # x
    flag.match[4, 1] = True  # X

    s = flag.encode_unknown_info(master_ind=0)
    assert s == '-_.xX'
    s = flag.encode_info(master_ind=0, predict_ind=3)
    assert s == '-_.xX'

    flag.initialize_flags(0, 2)
    s = flag.encode_unknown_info(master_ind=0)
    assert s == ''


@pytest.fixture
def region_db(mocker):
    mocker.patch(
        'misc.read_table.open',
        mocker.mock_open(
            read_data='region_id\tstrain\tchromosome\t'
            'predicted_species\tstart\tend\tnum_sites_hmm\n'
            'r4\tyjm1381\tI\tS288c\t24327\t26193\t60\n'
            'r5\tyjm689\tI\tS288c\t24327\t24444\t56\n'
            'r6\tyjm1381\tI\tS288c\t24612\t25439\t18\n'
            'r7\tyjm689\tI\tS288c\t24612\t138647\t13728\n'
            'r8\tyjm1208\tI\tS288c\t25395\t25448\t20\n'
        ))

    return summarize.Region_Database('labeled.txt',
                                     'I',
                                     ['s1', 's2'])


def test_region_init(mocker):
    mock_open = mocker.patch(
        'misc.read_table.open',
        mocker.mock_open(
            read_data='region_id\tstrain\tchromosome\t'
            'predicted_species\tstart\tend\tnum_sites_hmm\n'
            'r4\tyjm1381\tI\tS288c\t24327\t26193\t60\n'
            'r5\tyjm689\tI\tS288c\t24327\t24444\t56\n'
            'r6\tyjm1381\tI\tS288c\t24612\t25439\t18\n'
            'r7\tyjm689\tI\tS288c\t24612\t138647\t13728\n'
            'r8\tyjm1208\tI\tS288c\t25395\t25448\t20\n'
            'r9\tyjm1304\tII\tS288c\t25395\t25436\t16\n'
        ))

    db = summarize.Region_Database('labeled.txt',
                                   'I',
                                   ['s1', 's2'])
    mock_open.assert_called_with('labeled.txt', 'r')

    assert db.labels == ['region_id', 'strain', 'chromosome',
                         'predicted_species', 'start', 'end', 'num_sites_hmm',
                         'match_nongap_s1', 'match_nongap_s2',
                         'num_sites_nongap_s1', 'num_sites_nongap_s2',
                         'match_hmm_s1', 'match_hmm_s2', 'match_nonmask_s1',
                         'match_nonmask_s2', 'num_sites_nonmask_s1',
                         'num_sites_nonmask_s2', 'count_.', 'count_-',
                         'count__', 'count_n', 'count_p', 'count_b',
                         'count_c', 'count_x', 'count_N', 'count_P',
                         'count_B', 'count_C', 'count_X']

    assert db.data == {
        'yjm1381': {
            'region_id': ['r4', 'r6'],
            'strain': ['yjm1381', 'yjm1381'],
            'chromosome': ['I', 'I'],
            'predicted_species': ['S288c', 'S288c'],
            'start': ['24327', '24612'],
            'end': ['26193', '25439'],
            'num_sites_hmm': ['60', '18'],
            'match_nongap_s1': [0, 0],
            'num_sites_nongap_s1': [0, 0],
            'match_hmm_s1': [0, 0],
            'match_nonmask_s1': [0, 0],
            'num_sites_nonmask_s1': [0, 0],
            'match_nongap_s2': [0, 0],
            'num_sites_nongap_s2': [0, 0],
            'match_hmm_s2': [0, 0],
            'match_nonmask_s2': [0, 0],
            'num_sites_nonmask_s2': [0, 0],
            'count_.': [0, 0], 'count_-': [0, 0], 'count__': [0, 0],
            'count_n': [0, 0], 'count_p': [0, 0], 'count_b': [0, 0],
            'count_c': [0, 0], 'count_x': [0, 0], 'count_N': [0, 0],
            'count_P': [0, 0], 'count_B': [0, 0], 'count_C': [0, 0],
            'count_X': [0, 0]},
        'yjm689': {
            'region_id': ['r5', 'r7'],
            'strain': ['yjm689', 'yjm689'],
            'chromosome': ['I', 'I'],
            'predicted_species': ['S288c', 'S288c'],
            'start': ['24327', '24612'],
            'end': ['24444', '138647'],
            'num_sites_hmm': ['56', '13728'],
            'match_nongap_s1': [0, 0],
            'num_sites_nongap_s1': [0, 0],
            'match_hmm_s1': [0, 0],
            'match_nonmask_s1': [0, 0],
            'num_sites_nonmask_s1': [0, 0],
            'match_nongap_s2': [0, 0],
            'num_sites_nongap_s2': [0, 0],
            'match_hmm_s2': [0, 0],
            'match_nonmask_s2': [0, 0],
            'num_sites_nonmask_s2': [0, 0],
            'count_.': [0, 0], 'count_-': [0, 0], 'count__': [0, 0],
            'count_n': [0, 0], 'count_p': [0, 0], 'count_b': [0, 0],
            'count_c': [0, 0], 'count_x': [0, 0], 'count_N': [0, 0],
            'count_P': [0, 0], 'count_B': [0, 0], 'count_C': [0, 0],
            'count_X': [0, 0]},
        'yjm1208': {
            'region_id': ['r8'],
            'strain': ['yjm1208'],
            'chromosome': ['I'],
            'predicted_species': ['S288c'],
            'start': ['25395'],
            'end': ['25448'],
            'num_sites_hmm': ['20'],
            'match_nongap_s1': [0],
            'num_sites_nongap_s1': [0],
            'match_hmm_s1': [0],
            'match_nonmask_s1': [0],
            'num_sites_nonmask_s1': [0],
            'match_nongap_s2': [0],
            'num_sites_nongap_s2': [0],
            'match_hmm_s2': [0],
            'match_nonmask_s2': [0],
            'num_sites_nonmask_s2': [0],
            'count_.': [0], 'count_-': [0], 'count__': [0],
            'count_n': [0], 'count_p': [0], 'count_b': [0],
            'count_c': [0], 'count_x': [0], 'count_N': [0],
            'count_P': [0], 'count_B': [0], 'count_C': [0],
            'count_X': [0]}}

    assert db.info_string_symbols == list('.-_npbcxNPBCX')
    assert db.label_prefixes == [
        'match_nongap',
        'num_sites_nongap',
        'match_hmm',
        'match_nonmask',
        'num_sites_nonmask']


def test_has_strain(region_db):
    assert region_db.has_strain('yjm689')
    assert not region_db.has_strain('yjm688')


def test_get_entries(region_db):
    result = [('r4', 24327, 26193),
              ('r6', 24612, 25439)]
    for i, (r_id, start, end) in enumerate(region_db.get_entries('yjm1381')):
        entry = result[i]
        assert r_id == entry[0]
        assert start == entry[1]
        assert end == entry[2]

    result = [('r5', 24327, 24444),
              ('r7', 24612, 138647)]
    for i, (r_id, start, end) in enumerate(region_db.get_entries('yjm689')):
        entry = result[i]
        assert r_id == entry[0]
        assert start == entry[1]
        assert end == entry[2]

    result = [('r8', 25395, 25448)]
    for i, (r_id, start, end) in enumerate(region_db.get_entries('yjm1208')):
        entry = result[i]
        assert r_id == entry[0]
        assert start == entry[1]
        assert end == entry[2]

    with pytest.raises(ValueError) as e:
        list(region_db.get_entries('asdf'))
    assert 'Region Database does not contain strain asdf' in str(e)


def test_set_region(region_db):
    region_db.set_region('yjm1381',
                         0,
                         's1',
                         (10, 20),
                         (30, 40),
                         (50, 60))

    ds = region_db.data['yjm1381']

    assert ds['num_sites_hmm'][0] == 20
    assert ds[f'match_hmm_s1'][0] == 10
    assert ds[f'match_nongap_s1'][0] == 30
    assert ds[f'num_sites_nongap_s1'][0] == 40
    assert ds[f'match_nonmask_s1'][0] == 50
    assert ds[f'num_sites_nonmask_s1'][0] == 60

    region_db.set_region('yjm1381',
                         1,
                         's2',
                         (11, None),
                         (31, 41),
                         (51, 61))

    ds = region_db.data['yjm1381']

    # retained from initial value
    assert ds['num_sites_hmm'][1] == '18'
    assert ds[f'match_hmm_s2'][1] == 11
    assert ds[f'match_nongap_s2'][1] == 31
    assert ds[f'num_sites_nongap_s2'][1] == 41
    assert ds[f'match_nonmask_s2'][1] == 51
    assert ds[f'num_sites_nonmask_s2'][1] == 61


def test_generate_output(region_db):
    initial_lines = [
        'r4\tyjm1381\tI\tS288c\t24327\t26193\t60',
        'r5\tyjm689\tI\tS288c\t24327\t24444\t56',
        'r6\tyjm1381\tI\tS288c\t24612\t25439\t18',
        'r7\tyjm689\tI\tS288c\t24612\t138647\t13728',
        'r8\tyjm1208\tI\tS288c\t25395\t25448\t20',
    ]
    for i, line in enumerate(region_db.generate_output()):
        assert line == initial_lines[i] + '\t' + '\t'.join(['0']*23) + '\n'


def test_generate_header(region_db, mocker):
    assert region_db.generate_header() == (
        'region_id\tstrain\tchromosome\tpredicted_species\tstart\tend'
        '\tnum_sites_hmm\tmatch_nongap_s1\tmatch_nongap_s2'
        '\tnum_sites_nongap_s1\tnum_sites_nongap_s2\tmatch_hmm_s1'
        '\tmatch_hmm_s2\tmatch_nonmask_s1\tmatch_nonmask_s2'
        '\tnum_sites_nonmask_s1\tnum_sites_nonmask_s2\tcount_.\tcount_-'
        '\tcount__\tcount_n\tcount_p\tcount_b\tcount_c\tcount_x\tcount_N'
        '\tcount_P\tcount_B\tcount_C\tcount_X\n')

    mocker.patch(
        'misc.read_table.open',
        mocker.mock_open(
            read_data='region_id\tstrain\tchromosome\n'
            'r4\tyjm1381\tII\n'
        ))

    # fewer headers, states, same counts
    db = summarize.Region_Database('labeled.txt', 'II', ['st1'])
    assert db.generate_header() == (
        'region_id\tstrain\tchromosome\tmatch_nongap_st1'
        '\tnum_sites_nongap_st1\tmatch_hmm_st1\tmatch_nonmask_st1'
        '\tnum_sites_nonmask_st1\tcount_.\tcount_-\tcount__\tcount_n\tcount_p'
        '\tcount_b\tcount_c\tcount_x\tcount_N\tcount_P\tcount_B\tcount_C'
        '\tcount_X\n')


def test_update_counts(region_db):
    syms = list('.-_npbcxNPBCX')

    region_db.update_counts('yjm1381', 0, '.-_npbcxNPBCX')
    for s in syms:
        assert region_db.data['yjm1381'][f'count_{s}'][0] == 1

    region_db.update_counts('yjm1381', 1, '.-_npbcxNPBCX'*20)
    for s in syms:
        assert region_db.data['yjm1381'][f'count_{s}'][1] == 20

    region_db.update_counts('yjm689', 1, '._pcNBX'*40)
    for i, s in enumerate(syms):
        if i % 2 == 0:
            assert region_db.data['yjm689'][f'count_{s}'][1] == 40
        else:
            assert region_db.data['yjm689'][f'count_{s}'][1] == 0

    # same values, gets overwritten
    region_db.update_counts('yjm689', 1, '-nbxPC'*10)
    for i, s in enumerate(syms):
        if i % 2 == 0:
            assert region_db.data['yjm689'][f'count_{s}'][1] == 0
        else:
            assert region_db.data['yjm689'][f'count_{s}'][1] == 10


@pytest.fixture
def region_context(mocker):
    writer = summarize.Region_Writer(
        'test_region.gz', 'test_index.pkl', 's1 s2'.split())
    mock_gzip = mocker.patch('analyze.summarize_region_quality.gzip')
    mock_open = mocker.patch('analyze.summarize_region_quality.open')
    mock_pickle = mocker.patch('analyze.summarize_region_quality.pickle.dump')

    return (writer, mock_gzip, mock_open, mock_pickle)


def test_region_context(region_context, mocker):
    writer, mock_gzip, mock_open, mock_pickle = region_context

    assert writer.region_file == 'test_region.gz'
    assert writer.index_file == 'test_index.pkl'
    assert writer.index == {}
    assert writer.known_states == 's1 s2'.split()

    mock_gzip.assert_not_called()
    mock_open.assert_not_called()

    writer = writer.__enter__()
    assert writer.region_writer is not None
    mock_gzip.open.assert_called_once_with('test_region.gz', 'wt')
    mock_open.assert_not_called()

    writer.__exit__(None, None, None)
    mock_gzip.open.return_value.close.assert_called_once()
    mock_open.assert_called_once_with('test_index.pkl', 'wb')
    mock_pickle.assert_called_once_with({}, mocker.ANY)


def test_region_write_header(region_context):
    writer, _, _, _ = region_context
    output = StringIO()
    writer.region_writer = output
    writer.write_header('r123')
    writer.write_header('f13')
    writer.write_header('q3')
    writer.write_header('213')

    # note, 213 is changed to 13 when storing in index, overwriting f13
    assert output.getvalue() == '#r123\n#f13\n#q3\n#213\n'
    assert writer.index == {3: 11, 13: 15, 123: 0}


def test_region_write_sequences(region_context):
    writer, _, _, _ = region_context
    output = StringIO()
    writer.region_writer = output

    writer.write_sequences('strain',
                           [
                               list(range(10)),
                               list(range(2, 12)),
                               list(range(4, 14))
                           ],
                           np.array([
                               list('0123456789')*3,
                               list('123456789')*3,
                               list('23456789')*3,
                           ]),
                           (5, 9))
    assert output.getvalue() == (
        '> s1 5 9\n'
        '56789\n'
        '> s2 3 7\n'
        '67891\n'
        '> strain 1 5\n'
        '78923\n')


def test_region_write_info_string(region_context):
    writer, _, _, _ = region_context
    output = StringIO()
    writer.region_writer = output

    writer.write_info_string('this is my string')
    assert output.getvalue() == (
        '> info\n'
        'this is my string\n')


def test_position_reader_context(mocker):
    mock_gzip = mocker.patch('analyze.summarize_region_quality.gzip.open')
    reader = summarize.Position_Reader('test_file.gz')
    assert reader.position_file == 'test_file.gz'
    assert reader.last_position == 0

    reader = reader.__enter__()
    assert reader.reader is not None
    mock_gzip.assert_called_once_with('test_file.gz', 'rt')

    reader = reader.__exit__(None, None, None)
    mock_gzip.return_value.close.assert_called_once()


def test_position_reader_next_line():
    reader = summarize.Position_Reader('mock')
    positions = StringIO(
        'yjm1460\tI\t25957\t25958\t25961\t25963\n'
        'yjm1463\tI\t25665\t25668\t25670\t25676\n'
        'yjm1464\tI\t25665\t25668\t25670\t25676\n'
        'yjm1464\tII\t25665\t25668\t25670\t25676\n'
        'yjm1460\tIII\t25957\t25958\t25961\t25963\n'
    )
    reader.reader = positions

    assert reader.next_line() == \
        'yjm1460\tI\t25957\t25958\t25961\t25963\n'
    assert reader.last_position == 0

    assert reader.next_line() == \
        'yjm1463\tI\t25665\t25668\t25670\t25676\n'
    assert reader.last_position == 34

    assert reader.next_line() == \
        'yjm1464\tI\t25665\t25668\t25670\t25676\n'
    assert reader.last_position == 68

    assert reader.next_line() == \
        'yjm1464\tII\t25665\t25668\t25670\t25676\n'
    assert reader.last_position == 102

    assert reader.next_line() == \
        'yjm1460\tIII\t25957\t25958\t25961\t25963\n'
    assert reader.last_position == 137

    assert reader.next_line() == ''
    assert reader.last_position == 173


def test_position_reader_get_positions(region_db):
    # region_db contains yjm1381, 689, 1208
    reader = summarize.Position_Reader('mock')
    positions = StringIO(
        'yjm1381\tI\t25957\t25958\t25961\t25963\n'
        'yjm689\tI\t25665\t25668\t25670\t25676\n'
        'yjm1464\tI\t25665\t25668\t25670\t25676\n'
        'yjm1464\tII\t25665\t25668\t25670\t25676\n'
        'yjm1381\tIII\t25957\t25958\t25961\t25963\n'
    )
    reader.reader = positions

    results = [
        ('yjm1381', np.array([25957, 25958, 25961, 25963]), 0),
        ('yjm689', np.array([25665, 25668, 25670, 25676]), 34),
    ]

    for i, (strain, ps) in enumerate(reader.get_positions(region_db, 'I')):
        assert strain == results[i][0]
        aae(ps, results[i][1])
        assert reader.last_position == results[i][2]

    assert i == 1
    assert reader.last_position == 101

    i = None
    # won't run because chromosome is not in order (on II)
    for i, (strain, ps) in enumerate(reader.get_positions(region_db, 'III')):
        pass

    assert i is None
    assert reader.last_position == 101

    # won't return since strain not in regions, will change last position
    for i, (strain, ps) in enumerate(reader.get_positions(region_db, 'II')):
        pass

    # if loop has not run
    assert i is None
    assert reader.last_position == 136

    for i, (strain, ps) in enumerate(reader.get_positions(region_db, 'III')):
        assert strain == 'yjm1381'
        aae(ps, np.array([25957, 25958, 25961, 25963]))
        assert reader.last_position == 136

    assert i == 0
    assert reader.last_position == 172


def test_quality_writer_context(mocker):
    mock_open = mocker.patch('analyze.summarize_region_quality.open')
    writer = summarize.Quality_Writer('test_file.txt')
    assert writer.filename == 'test_file.txt'
    assert writer.first_write is True

    writer = writer.__enter__()
    assert writer.writer is not None
    mock_open.assert_called_once_with('test_file.txt', 'w')

    writer = writer.__exit__(None, None, None)
    mock_open.return_value.close.assert_called_once()


def test_quality_writer_write_quality(mocker):
    writer = summarize.Quality_Writer('test')
    output = StringIO()
    writer.writer = output

    assert writer.first_write is True
    mock_region = mocker.MagicMock()
    mock_region.generate_header.return_value = 'header1\n'
    mock_region.generate_output.return_value = [
        'region1\n',
        'region2\n',
    ]

    writer.write_quality(mock_region)
    assert writer.first_write is False

    mock_region = mocker.MagicMock()
    mock_region.generate_header.return_value = 'header2\n'
    mock_region.generate_output.return_value = [
        'a\tdifferent\tformat\n',
    ]

    writer.write_quality(mock_region)

    assert writer.first_write is False
    assert output.getvalue() == (
        'header1\n'
        'region1\n'
        'region2\n'
        'a\tdifferent\tformat\n'
    )


def test_sequence_analyzer_init():
    sa = summarize.Sequence_Analyzer('mask',
                                     'alignment',
                                     ['s1', 's2'],
                                     ['s1', 's2'],
                                     ['I', 'II'],
                                     {})

    assert sa.masks == 'mask'
    assert sa.alignments == 'alignment'
    assert sa.known_states == ['s1', 's2']
    assert sa.chromosomes == ['I', 'II']
    assert sa.symbols == {}


@pytest.fixture
def sa():
    symbols = {
        'match': '+',
        'mismatch': '-',
        'unknown': '?',
        'unsequenced': 'n',
        'gap': '-',
        'unaligned': '?',
        'masked': 'x'
    }
    return summarize.Sequence_Analyzer('', '', [], [], [], symbols)


def test_SA_build_masked_sites(sa, mocker):
    mock_open = mocker.patch(
        'analyze.summarize_region_quality.open',
        mocker.mock_open(
            read_data='>header\n'
            '0 - 2\n'
            '22 - 25\n'
            '32 - 33\n'
        ))

    sa.chromosomes = ['I', 'II']
    sa.known_states = ['s1', 's2']
    sa.interval_states = ['i1', 'i2']
    sa.masks = '{strain}_chr{chrom}_intervals.txt'
    sa.build_masked_sites()
    sites = sa.masked_sites

    assert mock_open.call_args_list == [
        mocker.call('i1_chrI_intervals.txt', 'r'),
        mocker.call('i2_chrI_intervals.txt', 'r'),
        mocker.call('i1_chrII_intervals.txt', 'r'),
        mocker.call('i2_chrII_intervals.txt', 'r'),
    ]

    expected = {
        'I': {
            's1': np.array([0, 1, 2, 22, 23, 24, 25, 32, 33]),
            's2': np.array([0, 1, 2, 22, 23, 24, 25, 32, 33]),
        },
        'II': {
            's1': np.array([0, 1, 2, 22, 23, 24, 25, 32, 33]),
            's2': np.array([0, 1, 2, 22, 23, 24, 25, 32, 33]),
        }
    }
    for chrom in sites:
        for state in sites[chrom]:
            aae(sites[chrom][state], expected[chrom][state])

    mock_open = mocker.patch(
        'analyze.summarize_region_quality.open',
        mocker.mock_open(
            read_data='>header\n'
        ))

    sa.build_masked_sites()
    sites = sa.masked_sites

    assert mock_open.call_args_list == [
        mocker.call('i1_chrI_intervals.txt', 'r'),
        mocker.call('i2_chrI_intervals.txt', 'r'),
        mocker.call('i1_chrII_intervals.txt', 'r'),
        mocker.call('i2_chrII_intervals.txt', 'r'),
    ]

    expected = {
        'I': {
            's1': np.array([]),
            's2': np.array([]),
        },
        'II': {
            's1': np.array([]),
            's2': np.array([]),
        }
    }
    for chrom in sites:
        for state in sites[chrom]:
            aae(sites[chrom][state], expected[chrom][state])


def test_read_masked_sites(sa, mocker):
    mock_open = mocker.patch(
        'analyze.summarize_region_quality.open',
        mocker.mock_open(
            read_data='>header\n'
            '0 - 2\n'
            '22 - 25\n'
            '32 - 33\n'
        ))

    sa.masks = '{chrom}_{strain}_mock'
    result = sa.read_masked_sites('I', 'str')
    assert mock_open.call_args_list == [
        mocker.call('I_str_mock', 'r')
    ]
    aae(result, np.array('0 1 2 22 23 24 25 32 33'.split(), dtype=int))


def test_convert_intervals_to_sites(sa):
    sites = sa.convert_intervals_to_sites([])
    assert sites == approx([])

    sites = sa.convert_intervals_to_sites([(1, 2)])
    assert sites == approx([1, 2])

    sites = sa.convert_intervals_to_sites([(1, 2), (4, 6)])
    assert sites == approx([1, 2, 4, 5, 6])


def test_read_masked_intervals(sa, mocker):
    lines = StringIO('')
    mocked_file = mocker.patch('analyze.summarize_region_quality.open',
                               return_value=lines)
    intervals = sa.read_masked_intervals('mocked')
    mocked_file.assert_called_with('mocked', 'r')
    assert intervals == []

    lines = StringIO('I am a header')
    mocked_file = mocker.patch('analyze.summarize_region_quality.open',
                               return_value=lines)
    intervals = sa.read_masked_intervals('mocked')
    mocked_file.assert_called_with('mocked', 'r')
    assert intervals == []

    lines = StringIO('I am a header\n'
                     'short and stout')
    mocked_file = mocker.patch('analyze.summarize_region_quality.open',
                               return_value=lines)
    with pytest.raises(ValueError):
        intervals = sa.read_masked_intervals('mocked')
    mocked_file.assert_called_with('mocked', 'r')

    lines = StringIO('I am a header\n'
                     '1 and 2')
    mocked_file = mocker.patch('analyze.summarize_region_quality.open',
                               return_value=lines)
    intervals = sa.read_masked_intervals('mocked')
    assert intervals == [(1, 2)]


def test_get_stats(sa):
    hmm, nongap, nonmask = sa.get_stats(
        np.array(list('abc')),
        np.array(list('abd')),
        0,
        [1, 2, 5],
        ([0], []))

    assert hmm[0] == 1
    assert hmm[1] == 2
    assert hmm[2].gap == approx([False] * 3)
    assert hmm[2].hmm == approx([False, True, True])
    assert hmm[2].match == approx([True, True, False])
    assert hmm[2].unseq == approx([False] * 3)

    assert nongap[0] == 1
    assert nongap[1] == 1

    assert nonmask[0] == 1
    assert nonmask[1] == 2
    assert nonmask[2].mask == approx([True, False, False])


def test_process_alignment(sa, region_context, mocker):
    sa.known_states = ['k1', 'k2']
    sa.masked_sites = {
        'I': {
            'k1': np.array([0, 1, 2]),
            'k2': np.array([2, 4, 6]),
        }}
    region_writer, _, _, _ = region_context
    output = StringIO()
    region_writer.region_writer = output
    region_writer.known_states = ['k1', 'k2']
    mocker.patch('analyze.summarize_region_quality.read_fasta.read_fasta',
                 return_value=('',
                               np.asarray([
                                   list('--gatcctag--'),
                                   list('-agatgcaag-c'),
                                   list('-a-att-cagt-'),
                               ])))
    mocker.patch.object(summarize.Sequence_Analyzer, 'read_masked_sites',
                        return_value=np.array(
                            '2 3 4'.split(), int))

    # for region database
    mocker.patch(
        'misc.read_table.open',
        mocker.mock_open(
            read_data='region_id\tstrain\tchromosome\t'
            'predicted_species\tstart\tend\tnum_sites_hmm\n'
            'r4\ts1\tI\tS288c\t1\t3\t60\n'
            'r5\tyjm689\tI\tS288c\t24327\t24444\t56\n'
            'r6\ts1\tI\tS288c\t5\t7\t18\n'
            'r7\tyjm689\tI\tS288c\t24612\t138647\t13728\n'
            'r8\tyjm1208\tI\tS288c\t25395\t25448\t20\n'
        ))
    rd = summarize.Region_Database('labeled_file.txt',
                                   'I',
                                   ['k1', 'k2'])

    sa.process_alignment(0,
                         1,
                         'I',
                         's1',
                         np.array([1, 3, 5, 7]),
                         rd,
                         region_writer)

    # these are inputs with modified hmm (last column)
    region_input = ['r4\ts1\tI\tS288c\t1\t3\t2\t',
                    'r5\tyjm689\tI\tS288c\t24327\t24444\t56\t',
                    'r6\ts1\tI\tS288c\t5\t7\t2\t',
                    'r7\tyjm689\tI\tS288c\t24612\t138647\t13728\t',
                    'r8\tyjm1208\tI\tS288c\t25395\t25448\t20\t']
    # these are the counts
    region_output = [
        '\t'.join('2 2 3 3 1 1 0 0 0 0 0 0 3 '
                  '0 0 0 0 0 0 0 0 0 0'.split()) + '\n',
        '\t'.join(['0' for _ in range(23)]) + '\n',
        '\t'.join('2 2 3 3 1 1 2 2 2 2 2 0 1 '
                  '0 0 0 0 0 0 0 0 0 0'.split()) + '\n',
        '\t'.join(['0' for _ in range(23)]) + '\n',
        '\t'.join(['0' for _ in range(23)]) + '\n',
    ]

    for i, line in enumerate(rd.generate_output()):
        assert line == region_input[i] + region_output[i]

    region_writer_output = (
        '#r4\n'
        '> k1 1 3\n'
        'atc\n'
        '> k2 2 4\n'
        'atg\n'
        '> s1 1 3\n'
        'att\n'
        '> info\n'
        '___\n'
        '#r6\n'
        '> k1 5 7\n'
        'tag\n'
        '> k2 6 8\n'
        'aag\n'
        '> s1 4 6\n'
        'cag\n'
        '> info\n'
        '_..\n'
    )
    assert output.getvalue() == region_writer_output


def test_get_indices(sa, mocker):
    sa.known_states = ['k1', 'k2']
    sa.masked_sites = {
        'I': {
            'k1': np.array([0, 1, 2]),
            'k2': np.array([2, 4, 6]),
        }}
    sa.alignments = 'align_{chrom}_{strain}'
    sa.masks = 'mask_{chrom}_{strain}'
    mock_fasta = mocker.patch(
        'analyze.summarize_region_quality.read_fasta.read_fasta',
        return_value=('',
                      np.asarray([
                          list('--gatcctag--'),
                          list('-agatgcaag-c'),
                          list('-a-att-cagt-'),
                      ])))
    mock_mask = mocker.patch.object(summarize.Sequence_Analyzer,
                                    'read_masked_sites',
                                    return_value=np.array(
                                        '2 3 4'.split(), int))

    seq, align, mask = sa.get_indices('I', 's1')

    mock_fasta.assert_called_once_with('align_I_s1')
    mock_mask.assert_called_once_with('I', 's1')

    aae(seq, np.asarray([
        list('--gatcctag--'),
        list('-agatgcaag-c'),
        list('-a-att-cagt-'),
    ]))

    result = [
        np.array(range(2, 10)),
        np.array(list(range(1, 10)) + [11]),
        np.array([1, 3, 4, 5, 7, 8, 9, 10])
    ]
    for i, a in enumerate(align):
        aae(a, result[i])

    result = [
        np.array([2, 3, 4]),
        np.array([3, 5, 7]),
        np.array([4, 5, 7]),
    ]
    for i, m in enumerate(mask):
        aae(m, result[i])


def test_get_slice(sa):
    alignment = np.array([1, 2, 3, 4, 5])
    ps_align = np.array([2, 4])

    start, end = sa.get_slice(1, 3, alignment, ps_align)
    assert start == 2
    assert end == 4

    with pytest.raises(ValueError) as e:
        sa.get_slice(0, 3, alignment, ps_align)
    assert 'Slice not found in position alignment' in str(e)

    with pytest.raises(ValueError) as e:
        sa.get_slice(1, 4, alignment, ps_align)
    assert 'Slice not found in position alignment' in str(e)


def test_index_alignment_by_reference(sa):
    output = sa.index_alignment_by_reference(np.array(list('abc')))
    assert output == approx([0, 1, 2])

    output = sa.index_alignment_by_reference(np.array(list('a-b-c')))
    assert output == approx([0, 2, 4])


def test_seq_id_unmasked(sa):
    match, sites, info = sa.seq_id_unmasked(np.array(list('abd')),
                                            np.array(list('abc')),
                                            0, [], [])
    assert match == 2
    assert sites == 3
    assert info.mask == approx([False, False, False])

    match, sites, info = sa.seq_id_unmasked(np.array(list('abd')),
                                            np.array(list('abc')),
                                            0, [0], [])
    assert match == 1
    assert sites == 2
    assert info.mask == approx([True, False, False])

    match, sites, info = sa.seq_id_unmasked(np.array(list('abd')),
                                            np.array(list('abc')),
                                            2, [0], [1])
    assert match == 2
    assert sites == 3
    assert info.mask == approx([False, False, False])

    match, sites, info = sa.seq_id_unmasked(np.array(list('abd')),
                                            np.array(list('abc')),
                                            0, [0], [1])
    assert match == 0
    assert sites == 1
    assert info.mask == approx([True, True, False])


def test_seq_id_hmm(sa):
    match, sites, info = sa.seq_id_hmm(np.array(list('abd')),
                                       np.array(list('abc')),
                                       0, [1, 2, 5])
    assert match == 1  # only count matches in included sites
    assert sites == 2  # included, not matching
    assert info.gap == approx([False] * 3)
    assert info.hmm == approx([False, True, True])
    assert info.match == approx([True, True, False])
    assert info.unseq == approx([False] * 3)

    match, sites, info = sa.seq_id_hmm(np.array(list('n-d')),
                                       np.array(list('--c')),
                                       1, [3, 5])
    assert match == 0
    assert sites == 1
    assert info.gap == approx([True, True, False])
    assert info.hmm == approx([False, False, True])
    assert info.match == approx([False, True, False])
    assert info.unseq == approx([True, False, False])

    with pytest.raises(ValueError) as e:
        match, sites, d = sa.seq_id_hmm(np.array(list('n-d')),
                                        np.array(list('--c')),
                                        1, [2, 5])
    assert ('Need to skip site specified as included '
            f'seq1: -, seq2: -, index: 1') in str(e)

    with pytest.raises(ValueError) as e:
        match, sites, d = sa.seq_id_hmm(np.array(list('n-d')),
                                        np.array(list('--c')),
                                        1, [1, 5])
    assert ('Need to skip site specified as included '
            f'seq1: n, seq2: -, index: 0') in str(e)
