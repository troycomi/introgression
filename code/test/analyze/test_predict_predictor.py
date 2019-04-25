from analyze import predict
from hmm import hmm_bw as hmm
import pytest
from pytest import approx
from io import StringIO
from collections import defaultdict
import random
import numpy as np


@pytest.fixture
def predictor():
    result = predict.Predictor(
        configuration={
            'analysis_params':
            {'reference': {'name': 'S228c'},
             'known_states': [
                 {'name': 'CBS432'},
                 {'name': 'N_45'},
                 {'name': 'DBVPG6304'},
                 {'name': 'UWOPS91_917_1'},
             ],
             'unknown_states': [{'name': 'unknown'}]
             }
        }
    )
    return result


def test_predictor(predictor):
    assert predictor.known_states ==\
        'S228c CBS432 N_45 DBVPG6304 UWOPS91_917_1'.split()
    assert predictor.unknown_states == ['unknown']


def test_set_chromosomes(predictor):
    with pytest.raises(ValueError) as e:
        predictor.set_chromosomes()
    assert 'No chromosomes specified in config file!' in str(e)

    predictor.config = {'chromosomes': ['I']}
    predictor.set_chromosomes()
    assert predictor.chromosomes == ['I']


def test_set_blocks_file(predictor):
    with pytest.raises(ValueError) as e:
        predictor.set_blocks_file('blocks_file')
    assert '{state} not found in blocks_file' in str(e)

    predictor.set_blocks_file('blocks_file{state}')
    assert predictor.blocks == 'blocks_file{state}'

    with pytest.raises(ValueError) as e:
        predictor.set_blocks_file()
    assert 'No block file provided' in str(e)

    predictor.config = {'paths': {'analysis': {'block_files': 'blocks_file'}}}
    with pytest.raises(ValueError) as e:
        predictor.set_blocks_file()
    assert '{state} not found in blocks_file' in str(e)

    predictor.config = {'paths': {'analysis': {'block_files':
                                               'blocks_file{state}'}}}
    predictor.set_blocks_file()
    assert predictor.blocks == 'blocks_file{state}'


def test_set_prefix(predictor):
    predictor.known_states = ['s1']
    predictor.set_prefix()
    assert predictor.prefix == 's1'

    predictor.known_states = 's1 s2'.split()
    predictor.set_prefix()
    assert predictor.prefix == 's1_s2'

    predictor.set_prefix('prefix')
    assert predictor.prefix == 'prefix'

    predictor.known_states = []
    with pytest.raises(ValueError) as e:
        predictor.set_prefix()
    assert 'Unable to build prefix, no known states provided' in str(e)


def test_set_threshold(predictor):
    with pytest.raises(ValueError) as e:
        predictor.set_threshold()
    assert 'No threshold provided' in str(e)

    predictor.config = {'analysis_params': {'threshold': 'asdf'}}
    with pytest.raises(ValueError) as e:
        predictor.set_threshold()
    assert 'Unsupported threshold value: asdf' in str(e)

    predictor.set_threshold(0.05)
    assert predictor.threshold == 0.05

    predictor.config = {'analysis_params':
                        {'threshold': 'viterbi'}}
    predictor.set_threshold()
    assert predictor.threshold == 'viterbi'


def test_set_strains(predictor, mocker):
    mock_find = mocker.patch.object(predict.Predictor, 'find_strains')

    predictor.set_strains()
    mock_find.called_with(None)

    with pytest.raises(ValueError) as e:
        predictor.config = {'paths': {'test_strains': ['test']}}
        predictor.set_strains()
    assert '{strain} not found in test' in str(e)

    with pytest.raises(ValueError) as e:
        predictor.config = {'paths': {'test_strains': ['test{strain}']}}
        predictor.set_strains()
    assert '{chrom} not found in test{strain}' in str(e)

    predictor.config = {'paths': {'test_strains':
                                  ['test{strain}{chrom}']}}
    predictor.set_strains()
    mock_find.called_with(['test{strain}{chrom}'])

    predictor.set_strains('test{strain}{chrom}')
    mock_find.called_with(['test{strain}{chrom}'])


def test_find_strains(predictor, mocker):
    with pytest.raises(ValueError) as e:
        predictor.find_strains()
    assert ('Unable to find strains in config and '
            'no test_strains provided') in str(e)

    predictor.config = {'strains': ['test2', 'test1']}
    predictor.find_strains()
    # sorted
    assert predictor.strains == 'test1 test2'.split()

    predictor.config = {}
    predictor.chromosomes = ['I']

    # too many chroms for s1
    mock_glob = mocker.patch('analyze.predict.glob.iglob',
                             side_effect=[[
                                 'test_prefix_s1_c1.fa',
                                 'test_prefix_s2_c1.fa',
                                 'test_prefix_s1_c2.fa',
                                 'test_prefix.fa',
                             ]])
    mock_log = mocker.patch('analyze.predict.log')
    with pytest.raises(ValueError) as e:
        predictor.find_strains(['test_prefix_{strain}_{chrom}.fa'])

    assert 'Strain s1 has incorrect number of chromosomes. Expected 1 found 2'\
        in str(e)
    mock_glob.assert_called_with('test_prefix_*_*.fa')
    mock_log.info.assert_called_with('searching for test_prefix_*_*.fa')
    assert mock_log.debug.call_args_list == \
        [mocker.call("matched with ('s1', 'c1')"),
         mocker.call("matched with ('s2', 'c1')"),
         mocker.call("matched with ('s1', 'c2')"),
         ]

    # no matches
    mock_glob = mocker.patch('analyze.predict.glob.iglob',
                             side_effect=[[
                                 'test_prefix.fa',
                             ]])
    mock_log = mocker.patch('analyze.predict.log')
    with pytest.raises(ValueError) as e:
        predictor.find_strains(['test_prefix_{strain}_{chrom}.fa'])
    assert ('Found no chromosome sequence files in '
            "['test_prefix_{strain}_{chrom}.fa']") in str(e)
    mock_glob.assert_called_with('test_prefix_*_*.fa')
    mock_log.info.assert_called_with('searching for test_prefix_*_*.fa')
    assert mock_log.debug.call_args_list == []

    # correct, with second test_strains
    mock_glob = mocker.patch('analyze.predict.glob.iglob',
                             side_effect=[
                                 [
                                     'test_prefix_s1_c1.fa',
                                     'test_prefix_s2_c1.fa',
                                     'test_prefix.fa',
                                 ],
                                 ['test_prefix_c2_s3.fa']
                             ])
    mock_log = mocker.patch('analyze.predict.log')
    predictor.find_strains(['test_prefix_{strain}_{chrom}.fa',
                            'test_prefix_{chrom}_{strain}.fa'])
    assert mock_glob.call_args_list == \
        [mocker.call('test_prefix_*_*.fa'),
         mocker.call('test_prefix_*_*.fa')]
    assert mock_log.info.call_args_list ==\
        [mocker.call('searching for test_prefix_*_*.fa'),
         mocker.call('searching for test_prefix_*_*.fa')]
    assert mock_log.debug.call_args_list == \
        [mocker.call("matched with ('s1', 'c1')"),
         mocker.call("matched with ('s2', 'c1')"),
         mocker.call("matched with ('s3', 'c2')"),
         ]
    assert predictor.strains == ['s1', 's2', 's3']


def test_set_output_files(predictor):
    with pytest.raises(ValueError) as e:
        predictor.set_output_files('', '', '', '', '')
    assert 'No initial hmm file provided' in str(e)

    with pytest.raises(ValueError) as e:
        predictor.set_output_files('init', '', '', '', '')
    assert 'No trained hmm file provided' in str(e)

    with pytest.raises(ValueError) as e:
        predictor.set_output_files('init', 'trained', 'pos', 'prob', '')
    assert 'No alignment file provided' in str(e)

    with pytest.raises(ValueError) as e:
        predictor.set_output_files('init', 'trained', 'pos', 'prob', 'align')
    assert '{prefix} not found in align' in str(e)

    with pytest.raises(ValueError) as e:
        predictor.set_output_files('init', 'trained', 'pos', 'prob',
                                   'align{prefix}')
    assert '{strain} not found in align{prefix}' in str(e)

    with pytest.raises(ValueError) as e:
        predictor.set_output_files('init', 'trained', 'pos', 'prob',
                                   'align{prefix}{strain}')
    assert '{chrom} not found in align{prefix}{strain}' in str(e)

    predictor.prefix = 'pre'
    predictor.set_output_files('init', 'trained', 'pos', 'prob',
                               'align{prefix}{strain}{chrom}')
    assert predictor.hmm_initial == 'init'
    assert predictor.hmm_trained == 'trained'
    assert predictor.positions == 'pos'
    assert predictor.probabilities == 'prob'
    assert predictor.alignment == 'alignpre{strain}{chrom}'

    predictor.set_output_files('init', 'trained', '', 'prob',
                               'align{prefix}{strain}{chrom}')
    assert predictor.hmm_initial == 'init'
    assert predictor.hmm_trained == 'trained'
    assert predictor.positions is None
    assert predictor.probabilities == 'prob'
    assert predictor.alignment == 'alignpre{strain}{chrom}'

    with pytest.raises(ValueError) as e:
        predictor.config = {'paths': {'analysis': {'hmm_initial': 'init'}}}
        predictor.set_output_files('', '', '', '', '')
    assert 'No trained hmm file provided' in str(e)

    with pytest.raises(ValueError) as e:
        predictor.config = {'paths': {'analysis': {'hmm_initial': 'init',
                                                   'hmm_trained': 'trained',
                                                   'positions': 'pos'
                                                   }}}
        predictor.set_output_files('', '', '', '', '')
    assert 'No probabilities file provided' in str(e)

    with pytest.raises(ValueError) as e:
        predictor.config = {'paths': {'analysis': {'hmm_initial': 'init',
                                                   'hmm_trained': 'trained',
                                                   'positions': 'pos',
                                                   'probabilities': 'prob'
                                                   }}}
        predictor.set_output_files('', '', '', '', '')
    assert 'No alignment file provided' in str(e)

    predictor.config = {'paths': {'analysis': {
        'hmm_initial': 'init',
        'hmm_trained': 'trained',
        'positions': 'pos',
        'probabilities': 'prob',
        'alignment': 'align{prefix}{strain}{chrom}'
    }}}
    predictor.set_output_files('', '', '', '', '')

    assert predictor.hmm_initial == 'init'
    assert predictor.hmm_trained == 'trained'
    assert predictor.positions == 'pos'
    assert predictor.probabilities == 'prob'
    assert predictor.alignment == 'alignpre{strain}{chrom}'

    predictor.config = {'paths': {'analysis': {
        'hmm_initial': 'init',
        'hmm_trained': 'trained',
        'probabilities': 'prob',
        'alignment': 'align{prefix}{strain}{chrom}'
    }}}
    predictor.set_output_files('', '', '', '', '')

    assert predictor.hmm_initial == 'init'
    assert predictor.hmm_trained == 'trained'
    assert predictor.positions is None
    assert predictor.probabilities == 'prob'
    assert predictor.alignment == 'alignpre{strain}{chrom}'


def test_validate_arguments(predictor):
    predictor.chromosomes = 1
    predictor.blocks = 1
    predictor.prefix = 1
    predictor.strains = 1
    predictor.hmm_initial = 1
    predictor.hmm_trained = 1
    predictor.probabilities = 1
    predictor.alignment = 1
    predictor.known_states = 1
    predictor.unknown_states = 1
    predictor.threshold = 1
    predictor.config = {
        'analysis_params':
        {'reference': {'name': 'S228c'},
         'known_states': [
             {'name': 'CBS432',
              'expected_length': 1,
              'expected_fraction': 0.01},
             {'name': 'N_45',
              'expected_length': 1,
              'expected_fraction': 0.01},
             {'name': 'DBVPG6304',
              'expected_length': 1,
              'expected_fraction': 0.01},
             {'name': 'UWOPS91_917_1',
              'expected_length': 1,
              'expected_fraction': 0.01},
         ],
         'unknown_states': [{'name': 'unknown',
                             'expected_length': 1,
                             'expected_fraction': 0.01},
                            ]
         }
    }

    assert predictor.validate_arguments()

    args = [
        'chromosomes',
        'blocks',
        'prefix',
        'strains',
        'hmm_initial',
        'hmm_trained',
        'probabilities',
        'alignment',
        'known_states',
        'unknown_states',
        'threshold'
    ]

    for arg in args:
        predictor.__dict__[arg] = None
        with pytest.raises(ValueError) as e:
            predictor.validate_arguments()
        assert ('Failed to validate Predictor, '
                f'required argument {arg} was unset') in str(e)
        predictor.__dict__[arg] = 1

    predictor.config = {
        'analysis_params':
        {'reference': {'name': 'S228c'},
         'unknown_states': [{'name': 'unknown',
                             'expected_length': 1,
                             'expected_fraction': 0.01},
                            ]
         }
    }
    with pytest.raises(ValueError) as e:
        predictor.validate_arguments()
    assert 'Configuration did not provide any known_states' in str(e)

    predictor.config = {
        'analysis_params':
        {'known_states': [
             {'name': 'CBS432',
              'expected_length': 1,
              'expected_fraction': 0.01},
             {'name': 'N_45',
              'expected_length': 1,
              'expected_fraction': 0.01},
         ],
         'unknown_states': [{'name': 'unknown',
                             'expected_length': 1,
                             'expected_fraction': 0.01},
                            ]
         }
    }
    with pytest.raises(ValueError) as e:
        predictor.validate_arguments()
    assert 'Configuration did not specify a reference strain' in str(e)

    predictor.config = {
        'analysis_params':
        {'reference': {'name': 'S228c'},
         'known_states': [
             {'name': 'CBS432',
              'expected_fraction': 0.01},
             {'name': 'N_45',
              'expected_length': 1,
              'expected_fraction': 0.01},
         ],
         'unknown_states': [{'name': 'unknown',
                             'expected_length': 1,
                             'expected_fraction': 0.01},
                            ]
         }
    }
    with pytest.raises(ValueError) as e:
        predictor.validate_arguments()
    assert 'CBS432 did not provide an expected_length' in str(e)

    predictor.config = {
        'analysis_params':
        {'reference': {'name': 'S228c'},
         'known_states': [
             {'name': 'CBS432',
              'expected_length': 1,
              'expected_fraction': 0.01},
             {'name': 'N_45',
              'expected_length': 1,
              },
         ],
         'unknown_states': [{'name': 'unknown',
                             'expected_length': 1,
                             'expected_fraction': 0.01},
                            ]
         }
    }
    with pytest.raises(ValueError) as e:
        predictor.validate_arguments()
    assert 'N_45 did not provide an expected_fraction' in str(e)

    predictor.config = {
        'analysis_params':
        {'reference': {'name': 'S228c'},
         'known_states': [
             {'name': 'CBS432',
              'expected_length': 1,
              'expected_fraction': 0.01},
             {'name': 'N_45',
              'expected_length': 1,
              'expected_fraction': 0.01},
         ],
         'unknown_states': [{'name': 'unknown',
                             'expected_fraction': 0.01},
                            ]
         }
    }
    with pytest.raises(ValueError) as e:
        predictor.validate_arguments()
    assert 'unknown did not provide an expected_length' in str(e)

    predictor.config = {
        'analysis_params':
        {'reference': {'name': 'S228c'},
         'known_states': [
             {'name': 'CBS432',
              'expected_length': 1,
              'expected_fraction': 0.01},
             {'name': 'N_45',
              'expected_length': 1,
              'expected_fraction': 0.01},
         ],
         'unknown_states': [{'name': 'unknown',
                             'expected_length': 1,
                             },
                            ]
         }
    }
    with pytest.raises(ValueError) as e:
        predictor.validate_arguments()
    assert 'unknown did not provide an expected_fraction' in str(e)


def test_run_prediction_no_pos(predictor, mocker, capsys):
    predictor.chromosomes = ['I', 'II']
    predictor.blocks = 'blocks{state}.txt'
    predictor.prefix = 'prefix'
    predictor.strains = ['s1', 's2']
    predictor.hmm_initial = 'hmm_initial.txt'
    predictor.hmm_trained = 'hmm_trained.txt'
    predictor.probabilities = 'probs.txt'
    predictor.alignment = 'prefix_{strain}_chr{chrom}.maf'
    predictor.known_states = 'S228c CBS432 N_45 DBVP UWOP'.split()
    predictor.unknown_states = ['unknown']
    predictor.states = predictor.known_states + predictor.unknown_states
    predictor.threshold = 'viterbi'
    predictor.config = {
        'analysis_params':
        {'reference': {'name': 'S228c'},
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
    mocker.patch('analyze.predict.read_fasta',
                 return_value=(None,
                               [list('NNENNENNEN'),  # S228c
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


def test_run_prediction_full(predictor, mocker):
    predictor.chromosomes = ['I', 'II']
    predictor.blocks = 'blocks{state}.txt'
    predictor.prefix = 'prefix'
    predictor.strains = ['s1', 's2']
    predictor.hmm_initial = 'hmm_initial.txt'
    predictor.hmm_trained = 'hmm_trained.txt'
    predictor.probabilities = 'probs.txt'
    predictor.positions = 'pos.txt'
    predictor.alignment = 'prefix_{strain}_chr{chrom}.maf'
    predictor.known_states = 'S228c CBS432 N_45 DBVP UWOP'.split()
    predictor.unknown_states = ['unknown']
    predictor.states = predictor.known_states + predictor.unknown_states
    predictor.threshold = 'viterbi'
    predictor.config = {
        'analysis_params':
        {'reference': {'name': 'S228c'},
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
    mock_fasta = mocker.patch('analyze.predict.read_fasta',
                              return_value=(None,
                                            [list('NNENNENNEN'),  # S228c
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
        mocker.call('blocksS228c.txt', 'w'),
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

    # blocks S228c
    mock_files[2].__enter__().write.assert_has_calls(
        [
            mocker.call('strain\tchromosome\tpredicted_species'
                        '\tstart\tend\tnum_sites_hmm\n'),
            mocker.call('s1\tI\tS228c\t0\t1\t2'),
            mocker.call('\n'),
            mocker.call('s2\tI\tS228c\t0\t1\t2'),
            mocker.call('\n'),
            mocker.call('s1\tII\tS228c\t0\t1\t2'),
            mocker.call('\n'),
            mocker.call('s2\tII\tS228c\t0\t1\t2'),
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
        mocker.call('working on: s1 I'),
        mocker.call('working on: s2 I'),
        mocker.call('working on: s1 II'),
        mocker.call('working on: s2 II')
    ])


def test_write_hmm_header(predictor):
    predictor.known_states = []
    predictor.unknown_states = []
    predictor.emission_symbols = []
    writer = StringIO()
    predictor.write_hmm_header(writer)
    assert writer.getvalue() == 'strain\tchromosome\t\n'

    predictor.known_states = ['s1', 's2']
    predictor.unknown_states = ['u1']
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
    predictor.states = []
    predictor.write_state_probs([{}], output, 'strain', 'I')

    assert output.getvalue() == 'strain\tI\t\n'

    output = StringIO()
    predictor.states = list('abc')
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


def test_process_path(predictor, hm):
    probs = hm.posterior_decoding()[0]
    predictor.set_threshold(0.8)
    predictor.states = 'N E'.split()
    predictor.known_states = 'N E'.split()
    path, probability = predictor.process_path(hm)
    assert (probability == probs).all()
    assert path == 'E E N E E N E E N N'.split()

    predictor.set_threshold('viterbi')
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
    predictor.states = states
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
