from analyze import predict
from hmm import hmm_bw as hmm
import pytest
from pytest import approx
from collections import defaultdict
import random
import numpy as np
from analyze.introgression_configuration import Configuration


@pytest.fixture
def config():
    return Configuration()


@pytest.fixture
def default_builder(config):
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
             {'name': 'DBVPG6304',
              'expected_length': 10000,
              'expected_fraction': 0.025},
             {'name': 'UWOPS91_917_1',
              'expected_length': 10000,
              'expected_fraction': 0.025},
         ],
         'unknown_states': [{'name': 'unknown',
                             'expected_length': 1000,
                             'expected_fraction': 0.01},
                            ]
         }
    }
    builder = predict.HMM_Builder(config)
    config.set('states')
    builder.set_expected_values()
    builder.update_expected_length(1e5)
    return builder


@pytest.fixture
def builder(config):
    return predict.HMM_Builder(config)


def test_builder(builder):
    assert builder.symbols == {
        'match': '+',
        'mismatch': '-',
        'unknown': '?',
        'unsequenced': 'n',
        'gap': '-',
        'unaligned': '?',
        'masked': 'x'
    }


def test_init(mocker, config):
    mock_log = mocker.patch('analyze.predict.log')
    predict.HMM_Builder(config)
    # no config, all warnings
    mock_log.warning.has_calls([
        mocker.call("Symbol for match unset in config, using default '+'"),
        mocker.call("Symbol for mismatch unset in config, using default '-'"),
        mocker.call("Symbol for unknown unset in config, using default '?'"),
        mocker.call("Symbol for unsequenced unset in config, "
                    "using default 'n'"),
        mocker.call("Symbol for gap unset in config, using default '-'"),
        mocker.call("Symbol for unaligned unset in config, using default '?'"),
        mocker.call("Symbol for masked unset in config, using default 'x'")
    ])

    # config, same warnings as above along with unused
    mock_log = mocker.patch('analyze.predict.log')
    config.config = {'HMM_symbols': {'unused': 'X'}}
    predict.HMM_Builder(config)
    mock_log.warning.has_calls([
        mocker.call("Unused symbol in configuration: unused -> 'X'"),
        mocker.call("Symbol for mismatch unset in config, using default '-'"),
        mocker.call("Symbol for unknown unset in config, using default '?'"),
        mocker.call("Symbol for unsequenced unset in config, "
                    "using default 'n'"),
        mocker.call("Symbol for gap unset in config, using default '-'"),
        mocker.call("Symbol for unaligned unset in config, using default '?'"),
        mocker.call("Symbol for masked unset in config, using default 'x'")
    ])

    # overwrite
    mock_log = mocker.patch('analyze.predict.log')
    config.config = {'HMM_symbols': {'masked': 'X'}}
    predict.HMM_Builder(config)
    mock_log.debug.has_calls([
        mocker.call("Overwriting default symbol for masked with 'X'")
    ])


def test_update_emission_symbols(builder):
    assert builder.update_emission_symbols(1) == ['+', '-']
    assert builder.update_emission_symbols(3) == ['+++',
                                                  '++-',
                                                  '+-+',
                                                  '+--',
                                                  '-++',
                                                  '-+-',
                                                  '--+',
                                                  '---',
                                                  ]


def test_get_symbol_freqs(builder):
    sequence = '-++ +-+ ++- ---'.split()
    symbol_test_helper(sequence, builder)
    symbol_test_helper(['+'], builder)
    # get all len 10 symbols
    syms = builder.update_emission_symbols(10)

    random.seed(0)
    for i in range(10):
        sequence = [random.choice(syms) for j in range(100)]
        symbol_test_helper(sequence, builder)


def symbol_test_helper(sequence, builder):
    symb, weigh = builder.get_symbol_freqs(np.array(sequence))

    num_states = len(sequence[0])
    num_sites = len(sequence)

    individual_symbol_freqs = []
    for s in range(num_states):
        d = defaultdict(int)
        for i in range(num_sites):
            d[sequence[i][s]] += 1
        for sym in d:
            d[sym] /= num_sites
        individual_symbol_freqs.append(d)

    symbol_freqs = defaultdict(int)
    for i in range(num_sites):
        symbol_freqs[sequence[i]] += 1
    for sym in symbol_freqs:
        symbol_freqs[sym] /= num_sites

    # for each state, how often seq matches that state relative to
    # others
    weighted_match_freqs = []
    for s in range(num_states):
        weighted_match_freqs.append(
            individual_symbol_freqs[s][builder.symbols['match']])

    weighted_match_freqs /= np.sum(weighted_match_freqs)

    assert symb == symbol_freqs
    assert weigh == approx(weighted_match_freqs)


def test_set_expected_values(builder, config):
    config.config = {
        'analysis_params':
        {'reference': {'name': 'S288c'},
         'known_states': [
             {'name': 'CBS432',
              'expected_length': 10,
              'expected_fraction': 0.01},
             {'name': 'N_45',
              'expected_length': 10,
              'expected_fraction': 0.01},
             {'name': 'DBVPG6304',
              'expected_length': 10,
              'expected_fraction': 0.01},
             {'name': 'UWOPS91_917_1',
              'expected_length': 10,
              'expected_fraction': 0.01},
         ],
         'unknown_states': [{'name': 'unknown',
                             'expected_length': 10,
                             'expected_fraction': 0.01},
                            ]
         }
    }
    config.set('states')
    builder.set_expected_values()
    assert builder.expected_lengths == {
        'CBS432': 10,
        'N_45': 10,
        'DBVPG6304': 10,
        'UWOPS91_917_1': 10,
        'unknown': 10}

    assert builder.expected_fractions == {
        'S288c': 0.95,
        'CBS432': 0.01,
        'N_45': 0.01,
        'DBVPG6304': 0.01,
        'UWOPS91_917_1': 0.01,
        'unknown': 0.01}
    assert builder.ref_fraction == 0.96
    assert builder.other_sum == 0.004
    assert builder.ref_state == 'S288c'


def test_update_expected_length(builder, config):
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
             {'name': 'DBVPG6304',
              'expected_length': 10000,
              'expected_fraction': 0.025},
             {'name': 'UWOPS91_917_1',
              'expected_length': 10000,
              'expected_fraction': 0.025},
         ],
         'unknown_states': [{'name': 'unknown',
                             'expected_length': 1000,
                             'expected_fraction': 0.01},
                            ]
         }
    }
    config.set('states')
    builder.set_expected_values()

    assert builder.expected_lengths == {
        'CBS432': 10000,
        'N_45': 10000,
        'DBVPG6304': 10000,
        'UWOPS91_917_1': 10000,
        'unknown': 1000}

    assert builder.expected_fractions == {
        'S288c': 0.89,
        'CBS432': 0.025,
        'N_45': 0.025,
        'DBVPG6304': 0.025,
        'UWOPS91_917_1': 0.025,
        'unknown': 0.01}

    assert builder.ref_fraction == 0.9
    assert builder.other_sum == 1e-5
    assert builder.ref_state == 'S288c'

    builder.update_expected_length(1e5)
    assert builder.expected_lengths['S288c'] == 45000


def test_initial_probabilities(default_builder):
    probs = default_builder.initial_probabilities(
        [0.1, 0.2, 0.3, 0.4, 0.5])

    assert default_builder.expected_fractions == (
        {'DBVPG6304': 0.025,
         'UWOPS91_917_1': 0.025,
         'unknown': 0.01,
         'CBS432': 0.025,
         'N_45': 0.025,
         'S288c': 0.89})

    p = [0.1 + (0.89 - 0.1) * 0.9,
         0.2 + (0.025 - 0.2) * 0.9,
         0.3 + (0.025 - 0.3) * 0.9,
         0.4 + (0.025 - 0.4) * 0.9,
         0.5 + (0.025 - 0.5) * 0.9,
         0.01]

    p = p / np.sum(p, dtype=np.float)

    assert probs == approx(p)


def test_emission_probabilities(default_builder):
    # normal mode, 5 known_states
    symbols = default_builder.update_emission_symbols(5)
    emis = default_builder.emission_probabilities(symbols)

    iter_emis = iter_emission(default_builder, symbols)
    is_approx_equal_list_dict(emis, iter_emis)

    # more unknowns
    default_builder.unknown_states.append('test')
    emis = default_builder.emission_probabilities(symbols)
    iter_emis = iter_emission(default_builder, symbols)
    is_approx_equal_list_dict(emis, iter_emis)

    # no unknowns
    default_builder.unknown_states = []
    emis = default_builder.emission_probabilities(symbols)
    iter_emis = iter_emission(default_builder, symbols)
    is_approx_equal_list_dict(emis, iter_emis)

    # too many symbols
    symbols = default_builder.update_emission_symbols(6)
    emis = default_builder.emission_probabilities(symbols)
    iter_emis = iter_emission(default_builder, symbols)
    is_approx_equal_list_dict(emis, iter_emis)


def iter_emission(builder, symbols):
    probs = {'-+': 0.9,
             '++': 0.09,
             '--': 0.009,
             '+-': 0.001}
    mismatch_bias = 0.99

    known_len = len(builder.known_states)
    for k in probs:
        probs[k] *= 2**(known_len - 2)

    emis = []
    # using older, iterative version
    for s in range(known_len):
        emis.append(defaultdict(float))
        for symbol in symbols:
            key = symbol[0] + symbol[s]
            emis[s][symbol] = probs[key]

        emis[s] = mynorm(emis[s])

    symbol_len = len(symbols[0])
    for s in range(len(builder.unknown_states)):
        emis.append(defaultdict(float))
        for symbol in symbols:
            match_count = symbol.count('+')
            mismatch_count = symbol_len - match_count
            emis[s + known_len][symbol] = \
                (match_count * (1 - mismatch_bias)
                 + mismatch_count * mismatch_bias)
        emis[s + known_len] = mynorm(emis[s + known_len])

    return emis


def is_approx_equal_list_dict(actual, expected):
    for i in range(len(actual)):
        for k in actual[i]:
            assert actual[i][k] == approx(expected[i][k]),\
                "failed at i={}, k={}".format(i, k)


def mynorm(d):
    total = float(sum(d.values()))
    return {k: v/total for k, v in d.items()}


def test_transition_probabilities(default_builder):
    trans = default_builder.transition_probabilities()

    iter_trans = iter_transition(default_builder)
    for i in range(len(trans)):
        assert trans[i] == approx(iter_trans[i])


def iter_transition(builder):
    states = builder.known_states + builder.unknown_states
    expected_frac = builder.expected_fractions
    expected_length = builder.expected_lengths
    trans = []
    for i in range(len(states)):
        state_from = states[i]
        trans.append([])
        scale_other = 1 / (1 - expected_frac[state_from])
        for j in range(len(states)):
            state_to = states[j]
            if state_from == state_to:
                trans[i].append(1 - 1./expected_length[state_from])
            else:
                trans[i].append(1./expected_length[state_from] *
                                expected_frac[state_to] * scale_other)

        trans[i] /= np.sum(trans[i])

    return trans


def test_build_initial_hmm(default_builder):
    symbols = default_builder.update_emission_symbols(5)
    hm = default_builder.build_initial_hmm(
        symbols)

    assert default_builder.expected_fractions == (
        {'DBVPG6304': 0.025,
         'UWOPS91_917_1': 0.025,
         'unknown': 0.01,
         'CBS432': 0.025,
         'N_45': 0.025,
         'S288c': 0.89})

    p = [0.2 + (0.89 - 0.2) * 0.9,
         0.2 + (0.025 - 0.2) * 0.9,
         0.2 + (0.025 - 0.2) * 0.9,
         0.2 + (0.025 - 0.2) * 0.9,
         0.2 + (0.025 - 0.2) * 0.9,
         0.01]

    p = p / np.sum(p, dtype=np.float)
    assert hm.initial_p == approx(p)

    iter_emis = iter_emission(default_builder, symbols)
    hm2 = hmm.HMM()
    hm2.set_emissions(iter_emis)
    assert hm.emissions == approx(hm2.emissions)

    iter_trans = iter_transition(default_builder)
    for i in range(len(hm.transitions)):
        assert hm.transitions[i] == approx(iter_trans[i])


def test_run_hmm(default_builder, capsys, mocker):
    seqs = [list('NNENNENNEN'),  # S2288c
            list('NNNENEENNN'),  # CBS432
            list('NN-NNEENNN'),  # N_45
            list('NEENN-ENEN'),  # DBVPG6304
            list('ENENNEENEN'),  # UWOPS..
            list('NNENNEENEN'),  # predicted
            ]
    mock_fasta = mocker.patch('analyze.predict.read_fasta',
                              return_value=(None, seqs))
    mock_log_hmm = mocker.patch('hmm.hmm_bw.log.info')

    hmm_init, hmm, positions = default_builder.run_hmm('MOCKED', True)

    mock_fasta.called_with('MOCKED')

    # check hmm output
    assert mock_log_hmm.call_args_list[-3:] == \
        [mocker.call('Iteration 8'),
         mocker.call('Iteration 9'),
         mocker.call('finished in 10 iterations')]

    # ps are locations of polymorphic sites, not counting missing '-'
    assert positions == approx([0, 1, 3, 6, 8])
    assert np.array_equal(hmm.initial_p, np.array([1, 0, 0, 0, 0, 0]))
    np.testing.assert_allclose(
        hmm_init.initial_p,
        np.array([0.8212314, 0.03825122, 0.04350912,
                  0.04350912, 0.04350912, 0.00999001]))


def test_encode_sequence(builder, mocker):
    mock_fasta = mocker.patch('analyze.predict.read_fasta',
                              return_value=(None,
                                            [
                                                list('abcd'),
                                                list('abed'),
                                                list('bbcf'),
                                            ]))

    seq_coded, positions, len_pred = builder.encode_sequence('test', True)
    assert (seq_coded == '-- +- --'.split()).all()
    assert (positions == [0, 2, 3]).all()
    assert len_pred == 4
    mock_fasta.called_with('test')

    seq_coded, positions, len_pred = builder.encode_sequence('test2', False)
    assert (seq_coded == '-- ++ +- --'.split()).all()
    assert (positions == [0, 1, 2, 3]).all()
    assert len_pred == 4
    mock_fasta.called_with('test2')


def test_ungap_and_code(builder):
    # nothing in prediction
    sequence, positions = builder.ungap_and_code(
        '---',  # predicted reference string
        ['abc', 'def', 'ghi'],  # several references
        0)  # reference index
    assert positions == approx([])
    assert sequence == approx([])

    # one match
    sequence, positions = builder.ungap_and_code(
        'a--',
        ['abc', 'def', 'ghi'],
        0)
    assert positions == approx([0])
    assert sequence == ['+--']

    # no match from refs
    sequence, positions = builder.ungap_and_code(
        'a--',
        ['abc', 'def', '-hi'],
        0)
    assert positions == approx([])
    assert sequence == approx([])

    # two matches
    sequence, positions = builder.ungap_and_code(
        'ae-',
        ['abc', 'def', 'gei'],
        0)
    assert positions == approx([0, 1])
    assert (sequence == ['+--', '-++']).all()

    # mess with ref index
    sequence, positions = builder.ungap_and_code(
        'a--e-',
        ['a--bc', 'deeef', 'geeei'],
        0)
    assert positions == approx([0, 1])
    assert (sequence == ['+--', '-++']).all()
    sequence, positions = builder.ungap_and_code(
        'a--e-',
        ['a--bc', 'deeef', 'geeei'],
        1)
    assert positions == approx([0, 3])
    assert (sequence == ['+--', '-++']).all()

    sequence, positions = builder.ungap_and_code(
        'a---ef--i',
        ['ab-dhfghi',
         'a-cceeg-i',
         'a-ceef-hh'],
        0)

    assert (sequence == '+++ -++ +-+ ++-'.split()).all()
    assert positions == approx([0, 3, 4, 7])


def test_poly_sites(builder):
    sequence, positions = builder.poly_sites(
        np.array('+++ -++ +-+ ++-'.split()),
        np.array([0, 3, 4, 7])
    )
    assert (sequence == '-++ +-+ ++-'.split()).all()
    assert positions == approx([3, 4, 7])
