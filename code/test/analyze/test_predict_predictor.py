from analyze import predict
from hmm import hmm_bw as hmm
import pytest
from io import StringIO
from collections import defaultdict
import random
import numpy as np
from analyze.introgression_configuration import Configuration


@pytest.fixture
def config():
    config = Configuration()
    config.add_config({
            'analysis_params':
            {'reference': {'name': 'S288c'},
             'known_states': [
                 {'name': 'CBS432'},
                 {'name': 'N_45'},
                 {'name': 'DBVPG6304'},
                 {'name': 'UWOPS91_917_1'},
             ],
             'unknown_states': [{'name': 'unknown'}]
             }
        })

    return config


@pytest.fixture
def predictor(config):
    result = predict.Predictor(config)
    config.set_states()
    return result


def test_predictor(predictor):
    assert predictor.config.known_states ==\
        'S288c CBS432 N_45 DBVPG6304 UWOPS91_917_1'.split()
    assert predictor.config.unknown_states == ['unknown']


def test_run_prediction_no_pos(predictor, config, mocker, capsys):
    config.chromosomes = ['I', 'II']
    config.blocks = 'blocks{state}.txt'
    config.prefix = 'prefix'
    config.strains = ['s1', 's2']
    config.hmm_initial = 'hmm_initial.txt'
    config.hmm_trained = 'hmm_trained.txt'
    config.probabilities = 'probs.txt'
    config.positions = None
    config.alignment = 'prefix_{strain}_chr{chrom}.maf'
    config.known_states = 'S288c CBS432 N_45 DBVP UWOP'.split()
    config.unknown_states = ['unknown']
    config.states = config.known_states + config.unknown_states
    config.threshold = 'viterbi'
    config.config = {
        'analysis_params':
        {'reference': {'name': 'S288c'},
         'known_states': [
             {'name': 'CBS432',
              'expected_length': 10000,
              'expected_fraction': 0.025},
             {'name': 'N_45',
              'expected_length': 10000,
              'expected_fraction': 0.025},
             {'name': 'DBVP',
              'expected_length': 10000,
              'expected_fraction': 0.025},
             {'name': 'UWOP',
              'expected_length': 10000,
              'expected_fraction': 0.025},
         ],
         'unknown_states': [{'name': 'unknown',
                             'expected_length': 1000,
                             'expected_fraction': 0.01},
                            ]
         }
    }
    mock_files = [mocker.MagicMock() for i in range(8)]
    mocker.patch('analyze.predict.open',
                 side_effect=mock_files)
    mock_gzip = mocker.patch('analyze.predict.gzip.open')
    mocker.patch('analyze.predict.log')
    mocker.patch('analyze.predict.os.path.exists', return_value=True)
    mocker.patch('analyze.predict.read_fasta',
                 return_value=(None,
                               [list('NNENNENNEN'),  # S288c
                                list('NNNENEENNN'),  # CBS432
                                list('NN-NNEENNN'),  # N_45
                                list('NEENN-ENEN'),  # DBVPG6304
                                list('ENENNEENEN'),  # UWOPS..
                                list('NNENNEENEN'),  # predicted
                                ]
                               ))

    mock_log_hmm = mocker.patch('hmm.hmm_bw.log.info')

    predictor.run_prediction(only_poly_sites=True)

    # check hmm output
    assert mock_log_hmm.call_args_list[-3:] == \
        [mocker.call('Iteration 8'),
         mocker.call('Iteration 9'),
         mocker.call('finished in 10 iterations')]

    assert mock_gzip.call_args_list == [mocker.call('probs.txt', 'wt')]

    # probs and pos interspersed
    print(mock_gzip.return_value.__enter__().write.call_args_list)
    assert mock_gzip.return_value.__enter__().write.call_args_list == \
        [
            mocker.call('s1\tI\t'),
            mocker.ANY,
            mocker.call('\n'),
            mocker.call('s2\tI\t'),
            mocker.ANY,
            mocker.call('\n'),
            mocker.call('s1\tII\t'),
            mocker.ANY,
            mocker.call('\n'),
            mocker.call('s2\tII\t'),
            mocker.ANY,
            mocker.call('\n'),

        ]


def test_run_prediction_full(predictor, config, mocker):
    config.chromosomes = ['I', 'II']
    config.blocks = 'blocks{state}.txt'
    config.prefix = 'prefix'
    config.strains = ['s1', 's2']
    config.hmm_initial = 'hmm_initial.txt'
    config.hmm_trained = 'hmm_trained.txt'
    config.probabilities = 'probs.txt'
    config.positions = 'pos.txt'
    config.alignment = 'prefix_{strain}_chr{chrom}.maf'
    config.known_states = 'S288c CBS432 N_45 DBVP UWOP'.split()
    config.unknown_states = ['unknown']
    config.states = config.known_states + config.unknown_states
    config.threshold = 'viterbi'
    config.config = {
        'analysis_params':
        {'reference': {'name': 'S288c'},
         'known_states': [
             {'name': 'CBS432',
              'expected_length': 10000,
              'expected_fraction': 0.025},
             {'name': 'N_45',
              'expected_length': 10000,
              'expected_fraction': 0.025},
             {'name': 'DBVP',
              'expected_length': 10000,
              'expected_fraction': 0.025},
             {'name': 'UWOP',
              'expected_length': 10000,
              'expected_fraction': 0.025},
         ],
         'unknown_states': [{'name': 'unknown',
                             'expected_length': 1000,
                             'expected_fraction': 0.01},
                            ]
         }
    }
    mock_files = [mocker.MagicMock() for i in range(8)]
    mock_open = mocker.patch('analyze.predict.open',
                             side_effect=mock_files)
    mock_gzip = mocker.patch('analyze.predict.gzip.open')
    mock_log = mocker.patch('analyze.predict.log')
    mocker.patch('analyze.predict.os.path.exists', return_value=True)
    mock_fasta = mocker.patch('analyze.predict.read_fasta',
                              return_value=(None,
                                            [list('NNENNENNEN'),  # S288c
                                             list('NNNENEENNN'),  # CBS432
                                             list('NN-NNEENNN'),  # N_45
                                             list('NEENN-ENEN'),  # DBVPG6304
                                             list('ENENNEENEN'),  # UWOPS..
                                             list('NNENNEENEN'),  # predicted
                                             ]
                                            ))
    mock_log_hmm = mocker.patch('hmm.hmm_bw.log.info')

    predictor.run_prediction(only_poly_sites=True)

    # check hmm output
    assert mock_log_hmm.call_args_list[-3:] == \
        [mocker.call('Iteration 8'),
         mocker.call('Iteration 9'),
         mocker.call('finished in 10 iterations')]

    mock_open.assert_has_calls([
        mocker.call('hmm_initial.txt', 'w'),
        mocker.call('hmm_trained.txt', 'w'),
        mocker.call('blocksS288c.txt', 'w'),
        mocker.call('blocksCBS432.txt', 'w'),
        mocker.call('blocksN_45.txt', 'w'),
        mocker.call('blocksDBVP.txt', 'w'),
        mocker.call('blocksUWOP.txt', 'w'),
        mocker.call('blocksunknown.txt', 'w')])

    # hmm_initial
    mock_files[0].__enter__().write.assert_has_calls(
        [
            mocker.call('strain\tchromosome\t'),
            mocker.ANY,
            mocker.call('\n'),
            mocker.call('s1\tI\t'),
            mocker.ANY,
            mocker.call('\n'),
            mocker.call('s2\tI\t'),
            mocker.ANY,
            mocker.call('\n'),
            mocker.call('s1\tII\t'),
            mocker.ANY,
            mocker.call('\n'),
            mocker.call('s2\tII\t'),
            mocker.ANY,
            mocker.call('\n')
        ])
    # trained
    mock_files[1].__enter__().write.assert_has_calls(
        [
            mocker.call('strain\tchromosome\t'),
            mocker.ANY,
            mocker.call('\n'),
            mocker.call('s1\tI\t'),
            mocker.ANY,
            mocker.call('\n'),
            mocker.call('s2\tI\t'),
            mocker.ANY,
            mocker.call('\n'),
            mocker.call('s1\tII\t'),
            mocker.ANY,
            mocker.call('\n'),
            mocker.call('s2\tII\t'),
            mocker.ANY,
            mocker.call('\n')
        ])
    # check initial probability (5th write, dereference to get string...)
    hmm_entry = mock_files[1].__enter__().write.\
        call_args_list[4][0][0].split('\t')
    assert hmm_entry[0] == '1.0'
    assert hmm_entry[1] == '0.0'
    assert hmm_entry[2] == '0.0'
    assert hmm_entry[3] == '0.0'
    assert hmm_entry[4] == '0.0'
    assert hmm_entry[5] == '0.0'
    assert hmm_entry[6] == '0.0'

    # blocks S288c
    mock_files[2].__enter__().write.assert_has_calls(
        [
            mocker.call('strain\tchromosome\tpredicted_species'
                        '\tstart\tend\tnum_sites_hmm\n'),
            mocker.call('s1\tI\tS288c\t0\t1\t2'),
            mocker.call('\n'),
            mocker.call('s2\tI\tS288c\t0\t1\t2'),
            mocker.call('\n'),
            mocker.call('s1\tII\tS288c\t0\t1\t2'),
            mocker.call('\n'),
            mocker.call('s2\tII\tS288c\t0\t1\t2'),
            mocker.call('\n')
        ])
    # blocks CBS432
    mock_files[3].__enter__().write.assert_has_calls(
        [
            mocker.call('strain\tchromosome\tpredicted_species'
                        '\tstart\tend\tnum_sites_hmm\n'),
            mocker.call(''),
            mocker.call(''),
            mocker.call(''),
            mocker.call('')])
    # blocks N_45
    mock_files[4].__enter__().write.assert_has_calls(
        [
            mocker.call('strain\tchromosome\tpredicted_species'
                        '\tstart\tend\tnum_sites_hmm\n'),
            mocker.call(''),
            mocker.call(''),
            mocker.call(''),
            mocker.call('')])
    # blocks DBVP
    mock_files[5].__enter__().write.assert_has_calls(
        [
            mocker.call('strain\tchromosome\tpredicted_species'
                        '\tstart\tend\tnum_sites_hmm\n'),
            mocker.call(''),
            mocker.call(''),
            mocker.call(''),
            mocker.call('')])
    # blocks UWOP
    mock_files[6].__enter__().write.assert_has_calls(
        [
            mocker.call('strain\tchromosome\tpredicted_species'
                        '\tstart\tend\tnum_sites_hmm\n'),
            mocker.call('s1\tI\tUWOP\t3\t8\t3'),
            mocker.call('\n'),
            mocker.call('s2\tI\tUWOP\t3\t8\t3'),
            mocker.call('\n'),
            mocker.call('s1\tII\tUWOP\t3\t8\t3'),
            mocker.call('\n'),
            mocker.call('s2\tII\tUWOP\t3\t8\t3'),
            mocker.call('\n')
        ])
    # blocks unknown
    mock_files[7].__enter__().write.assert_has_calls(
        [
            mocker.call('strain\tchromosome\tpredicted_species'
                        '\tstart\tend\tnum_sites_hmm\n'),
            mocker.call(''),
            mocker.call(''),
            mocker.call(''),
            mocker.call('')])

    mock_gzip.assert_any_call('probs.txt', 'wt')
    mock_gzip.assert_any_call('pos.txt', 'wt')

    # probs and pos interspersed
    mock_gzip.return_value.__enter__().write.assert_has_calls(
        [
            mocker.call('s1\tI\t0\t1\t3\t6\t8\n'),
            mocker.call('s1\tI\t'),
            mocker.ANY,
            mocker.call('\n'),
            mocker.call('s2\tI\t0\t1\t3\t6\t8\n'),
            mocker.call('s2\tI\t'),
            mocker.ANY,
            mocker.call('\n'),
            mocker.call('s1\tII\t0\t1\t3\t6\t8\n'),
            mocker.call('s1\tII\t'),
            mocker.ANY,
            mocker.call('\n'),
            mocker.call('s2\tII\t0\t1\t3\t6\t8\n'),
            mocker.call('s2\tII\t'),
            mocker.ANY,
            mocker.call('\n'),

        ])

    mock_fasta.assert_has_calls([
        mocker.call('prefix_s1_chrI.maf'),
        mocker.call('prefix_s2_chrI.maf'),
        mocker.call('prefix_s1_chrII.maf'),
        mocker.call('prefix_s2_chrII.maf')
    ])

    mock_log.info.assert_has_calls([
        mocker.call('working on: s1 I (1 of 4)'),
        mocker.call('working on: s2 I (2 of 4)'),
        mocker.call('working on: s1 II (3 of 4)'),
        mocker.call('working on: s2 II (4 of 4)')
    ])


def test_write_hmm_header(predictor, config):
    config.states = []
    predictor.emission_symbols = []
    writer = StringIO()
    predictor.write_hmm_header(writer)
    assert writer.getvalue() == 'strain\tchromosome\t\n'

    config.states = ['s1', 's2', 'u1']
    predictor.emission_symbols = ['-', '+']
    writer = StringIO()
    predictor.write_hmm_header(writer)

    header = 'strain\tchromosome\t'
    header += '\t'.join(
        ['init_{}'.format(s) for s in ['s1', 's2', 'u1']] +
        ['emis_{}_{}'.format(s, sym)
         for s in ['s1', 's2', 'u1']
         for sym in ['-', '+']] +
        ['trans_{}_{}'.format(s, s2)
         for s in ['s1', 's2', 'u1']
         for s2 in ['s1', 's2', 'u1']])

    assert writer.getvalue() == header + '\n'


def test_write_hmm(predictor):
    predictor.emission_symbols = list('abc')
    output = StringIO()

    hm = hmm.HMM()

    # empty hmm
    predictor.write_hmm(hm, output, 'strain', 'I')
    assert output.getvalue() == 'strain\tI\t\n'

    hm.set_hidden_states(list('abc'))
    hm.set_initial_p([0, 1, 0])
    hm.set_transitions([[0, 1, 0], [1, 0, 0], [0, 0, 1]])
    hm.set_emissions([{'a': 1, 'b': 0, 'c': 0},
                      {'a': 0, 'b': 0, 'c': 1},
                      {'a': 0, 'b': 1, 'c': 0},
                      ])

    output = StringIO()
    predictor.write_hmm(hm, output, 'strain', 'I')

    result = 'strain\tI\t'
    result += '\t'.join(list('010')) + '\t'  # init
    result += '\t'.join(list('100001010')) + '\t'  # emis
    result += '\t'.join(list('010100001')) + '\n'  # trans
    assert output.getvalue() == result


def test_write_blocks_header(predictor):
    writer = StringIO()
    predictor.write_blocks_header(writer)

    assert writer.getvalue() == '\t'.join(['strain',
                                           'chromosome',
                                           'predicted_species',
                                           'start',
                                           'end',
                                           'num_sites_hmm']) + '\n'


def test_write_blocks(predictor):
    output = StringIO()
    block = []
    pos = [i * 2 for i in range(20)]
    predictor.write_blocks(block,
                           pos,
                           output, 'test', 'I', 'pred')

    assert output.getvalue() == ''

    output = StringIO()
    block = [(0, 1), (4, 6), (10, 8)]
    pos = [i * 2 for i in range(20)]
    predictor.write_blocks(block,
                           pos,
                           output, 'test', 'I', 'pred')

    result = "\n".join(
        ["\t".join(['test', 'I', 'pred',
                    str(pos[s]), str(pos[e]), str(e - s + 1)])
         for s, e in block]) + "\n"

    assert output.getvalue() == result


def test_write_positions(predictor):
    output = StringIO()
    predictor.write_positions([0, 1, 3, 5, 7], output, 'test', 'I')
    assert output.getvalue() == "{}\t{}\t{}\n".format(
        "test",
        "I",
        "\t".join([str(i) for i in (0, 1, 3, 5, 7)]))


def test_write_state_probs(predictor):
    output = StringIO()
    predictor.config.states = []
    predictor.write_state_probs([{}], output, 'strain', 'I')

    assert output.getvalue() == 'strain\tI\t\n'

    output = StringIO()
    predictor.config.states = list('abc')
    predictor.write_state_probs([
        [0, 0, 1],
        [1, 0, 0],
        [0, 1, 1],
    ], output, 'strain', 'I')

    assert output.getvalue() == \
        ('strain\tI\t'
         'a:0.00000,1.00000,0.00000\t'
         'b:0.00000,0.00000,1.00000\t'
         'c:1.00000,0.00000,1.00000\n')


def test_process_path(predictor, config, hm):
    probs = hm.posterior_decoding()[0]
    config.set_threshold(0.8)
    config.states = 'N E'.split()
    config.known_states = 'N E'.split()
    path, probability = predictor.process_path(hm)
    assert (probability == probs).all()
    assert path == 'E E N E E N E E N N'.split()

    config.set_threshold('viterbi')
    path, probability = predictor.process_path(hm)

    assert (probability == probs).all()
    assert path == 'E E N E E N E E N E'.split()


def test_convert_to_blocks(predictor):
    random.seed(0)
    states = [str(i) for i in range(10)]
    help_test_convert_blocks(states, list('1'), predictor)
    help_test_convert_blocks(states, list('12'), predictor)
    help_test_convert_blocks(states, list('1111'), predictor)

    for test in range(10):
        seq = [str(random.randint(0, 9)) for i in range(100)]
        help_test_convert_blocks(states, seq, predictor)


def help_test_convert_blocks(states, seq, predictor):
    predictor.config.states = states
    blocks = predictor.convert_to_blocks(seq)

    nseq = np.array(seq, int)
    # add element to the end to catch repeats on last index
    nseq = np.append(nseq, nseq[-1]+1)
    diff = np.diff(nseq)
    locs = np.nonzero(diff)[0]
    lens = np.diff(locs)
    lens = np.append(locs[0]+1, lens)

    current = 0
    result = defaultdict(list)
    for i, l in enumerate(locs):
        result[seq[l]].append((current, current + lens[i] - 1))
        current += lens[i]

    for k in blocks:
        assert blocks[k] == result[k]


def test_read_blocks(mocker):
    block_in = StringIO('''
''')

    mocked_file = mocker.patch('analyze.predict.open',
                               return_value=block_in)
    output = predict.read_blocks('mocked')

    mocked_file.assert_called_with('mocked', 'r')
    assert list(output.keys()) == []

    block_in = StringIO('''header
test\tI\tpred\t100\t200\t10
''')

    mocked_file = mocker.patch('analyze.predict.open',
                               return_value=block_in)
    output = predict.read_blocks('mocked')

    assert len(output) == 1
    assert output['test']['I'] == [(100, 200, 10)]

    block_in = StringIO('''header
test\tI\tpred\t100\t200\t10
test\tI\tpred\t200\t200\t30
test\tI\tpred\t300\t400\t40
test\tII\tpred\t300\t400\t40
test2\tIII\tpred\t300\t400\t47
''')

    mocked_file = mocker.patch('analyze.predict.open',
                               return_value=block_in)
    output = predict.read_blocks('mocked')

    assert len(output) == 2
    assert len(output['test']) == 2
    assert len(output['test2']) == 1
    assert output['test']['I'] == [
        (100, 200, 10),
        (200, 200, 30),
        (300, 400, 40),
    ]
    assert output['test']['II'] == [(300, 400, 40)]
    assert output['test2']['III'] == [(300, 400, 47)]


def test_read_blocks_labeled(mocker):
    block_in = StringIO('''
''')

    mocked_file = mocker.patch('analyze.predict.open',
                               return_value=block_in)
    output = predict.read_blocks('mocked', labeled=True)

    mocked_file.assert_called_with('mocked', 'r')
    assert list(output.keys()) == []

    block_in = StringIO('''header
r1\ttest\tI\tpred\t100\t200\t10
''')

    mocked_file = mocker.patch('analyze.predict.open',
                               return_value=block_in)
    output = predict.read_blocks('mocked', labeled=True)

    assert len(output) == 1
    assert output['test']['I'] == [('r1', 100, 200, 10)]

    block_in = StringIO('''header
r1\ttest\tI\tpred\t100\t200\t10
r2\ttest\tI\tpred\t200\t200\t30
r3\ttest\tI\tpred\t300\t400\t40
r4\ttest\tII\tpred\t300\t400\t40
r5\ttest2\tIII\tpred\t300\t400\t47
''')

    mocked_file = mocker.patch('analyze.predict.open',
                               return_value=block_in)
    output = predict.read_blocks('mocked', labeled=True)

    assert len(output) == 2
    assert len(output['test']) == 2
    assert len(output['test2']) == 1
    assert output['test']['I'] == [
        ('r1', 100, 200, 10),
        ('r2', 200, 200, 30),
        ('r3', 300, 400, 40),
    ]
    assert output['test']['II'] == [('r4', 300, 400, 40)]
    assert output['test2']['III'] == [('r5', 300, 400, 47)]
