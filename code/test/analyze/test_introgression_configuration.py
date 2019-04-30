from analyze.introgression_configuration import Configuration
import pytest


@pytest.fixture()
def config():
    return Configuration()


def test_set_log_file(config):
    config.set_log_file()
    assert config.log_file is None

    config.set_log_file('test')
    assert config.log_file == 'test'

    config.config = {'paths': {'log_file': 'log'}}
    config.set_log_file()
    assert config.log_file == 'log'

    config.set_log_file('test')
    assert config.log_file == 'test'


def test_set_chromosomes(config):
    with pytest.raises(ValueError) as e:
        config.set_chromosomes()
    assert 'No chromosomes specified in config file!' in str(e)

    config.config = {'chromosomes': ['I']}
    config.set_chromosomes()
    assert config.chromosomes == ['I']


def test_get_states(config):
    assert config.get_states() == ([], [])

    config.config = {
            'analysis_params': {
                'known_states': [
                    {'name': 'k1'},
                    {'name': 'k2'},
                    {'name': 'k3'},
                ],
                'unknown_states': [
                    {'name': 'u1'},
                    {'name': 'u2'},
                ]
            }
        }
    assert config.get_states() == ('k1 k2 k3'.split(), 'u1 u2'.split())

    config.config = {
            'analysis_params': {
                'reference': {'name': 'ref'},
                'unknown_states': [
                    {'name': 'u1'},
                    {'name': 'u2'},
                ]
            }
        }
    assert config.get_states() == ('ref'.split(), 'u1 u2'.split())

    config.config = {
            'analysis_params': {
                'reference': {'name': 'ref'},
                'known_states': [
                    {'name': 'k1'},
                    {'name': 'k2'},
                    {'name': 'k3'},
                ],
                'unknown_states': [
                    {'name': 'u1'},
                    {'name': 'u2'},
                ]
            }
        }
    assert config.get_states() == ('ref k1 k2 k3'.split(), 'u1 u2'.split())


def test_set_states(config):
    config.config = {
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
        }

    config.set_states()
    assert config.known_states ==\
        'S288c CBS432 N_45 DBVPG6304 UWOPS91_917_1'.split()
    assert config.unknown_states ==\
        'unknown'.split()
    assert config.states ==\
        'S288c CBS432 N_45 DBVPG6304 UWOPS91_917_1 unknown'.split()

    config.set_states([])
    assert config.known_states ==\
        'S288c CBS432 N_45 DBVPG6304 UWOPS91_917_1'.split()
    assert config.unknown_states ==\
        'unknown'.split()
    assert config.states ==\
        'S288c CBS432 N_45 DBVPG6304 UWOPS91_917_1 unknown'.split()

    config.set_states('testing 123'.split())
    assert config.states == ['testing', '123']

    config.config = {}

    with pytest.raises(ValueError) as e:
        config.set_states()
    assert 'No states specified' in str(e)


def test_set_threshold(config):
    with pytest.raises(ValueError) as e:
        config.set_threshold()
    assert 'No threshold provided' in str(e)

    config.config = {'analysis_params': {'threshold': 'asdf'}}
    with pytest.raises(ValueError) as e:
        config.set_threshold()
    assert 'Unsupported threshold value: asdf' in str(e)

    config.set_threshold(0.05)
    assert config.threshold == 0.05

    config.config = {'analysis_params':
                     {'threshold': 'viterbi'}}
    config.set_threshold()
    assert config.threshold == 'viterbi'


def test_set_labeled_blocks_file(config):
    with pytest.raises(ValueError) as e:
        config.set_labeled_blocks_file('blocks_file')
    assert '{state} not found in blocks_file' in str(e)

    config.set_labeled_blocks_file('blocks_file{state}')
    assert config.labeled_blocks == 'blocks_file{state}'

    with pytest.raises(ValueError) as e:
        config.set_labeled_blocks_file()
    assert 'No labeled block file provided' in str(e)

    config.config = {'paths': {'analysis':
                               {'labeled_block_files': 'blocks_file'}}}
    with pytest.raises(ValueError) as e:
        config.set_labeled_blocks_file()
    assert '{state} not found in blocks_file' in str(e)

    config.config = {'paths': {'analysis': {'labeled_block_files':
                                            'blocks_file{state}'}}}
    config.set_labeled_blocks_file()
    assert config.labeled_blocks == 'blocks_file{state}'


def test_set_blocks_file(config):
    with pytest.raises(ValueError) as e:
        config.set_blocks_file('blocks_file')
    assert '{state} not found in blocks_file' in str(e)

    config.set_blocks_file('blocks_file{state}')
    assert config.blocks == 'blocks_file{state}'

    with pytest.raises(ValueError) as e:
        config.set_blocks_file()
    assert 'No block file provided' in str(e)

    config.config = {'paths': {'analysis': {'block_files': 'blocks_file'}}}
    with pytest.raises(ValueError) as e:
        config.set_blocks_file()
    assert '{state} not found in blocks_file' in str(e)

    config.config = {'paths': {'analysis': {'block_files':
                                            'blocks_file{state}'}}}
    config.set_blocks_file()
    assert config.blocks == 'blocks_file{state}'


def test_set_prefix(config):
    config.known_states = ['s1']
    config.set_prefix()
    assert config.prefix == 's1'

    config.known_states = 's1 s2'.split()
    config.set_prefix()
    assert config.prefix == 's1_s2'

    config.set_prefix('prefix')
    assert config.prefix == 'prefix'

    config.known_states = []
    with pytest.raises(ValueError) as e:
        config.set_prefix()
    assert 'Unable to build prefix, no known states provided' in str(e)


def test_set_strains(config, mocker):
    mock_find = mocker.patch.object(Configuration, 'find_strains')

    config.set_strains()
    mock_find.called_with(None)

    with pytest.raises(ValueError) as e:
        config.config = {'paths': {'test_strains': ['test']}}
        config.set_strains()
    assert '{strain} not found in test' in str(e)

    with pytest.raises(ValueError) as e:
        config.config = {'paths': {'test_strains': ['test{strain}']}}
        config.set_strains()
    assert '{chrom} not found in test{strain}' in str(e)

    config.config = {'paths': {'test_strains':
                               ['test{strain}{chrom}']}}
    config.set_strains()
    mock_find.called_with(['test{strain}{chrom}'])

    config.set_strains('test{strain}{chrom}')
    mock_find.called_with(['test{strain}{chrom}'])


def test_find_strains(config, mocker):
    with pytest.raises(ValueError) as e:
        config.find_strains()
    assert ('Unable to find strains in config and '
            'no test_strains provided') in str(e)

    config.config = {'strains': ['test2', 'test1']}
    config.find_strains()
    # sorted
    assert config.strains == 'test1 test2'.split()

    config.config = {}
    config.chromosomes = ['I']

    # too many chroms for s1
    mock_glob = mocker.patch('analyze.introgression_configuration.glob.iglob',
                             side_effect=[[
                                 'test_prefix_s1_cII.fa',
                                 'test_prefix_s2_cII.fa',
                                 'test_prefix_s1_cIII.fa',
                                 'test_prefix.fa',
                             ]])
    mock_log = mocker.patch('analyze.introgression_configuration.log')
    with pytest.raises(ValueError) as e:
        config.find_strains(['test_prefix_{strain}_c{chrom}.fa'])

    assert "Strain s1 is missing chromosomes. Unable to find chromosome 'I'"\
        in str(e)
    mock_glob.assert_called_with('test_prefix_*_c*.fa')
    mock_log.info.assert_called_with('searching for test_prefix_*_c*.fa')
    assert mock_log.debug.call_args_list == \
        [mocker.call("matched with ('s1', 'II')"),
         mocker.call("matched with ('s2', 'II')"),
         mocker.call("matched with ('s1', 'III')"),
         ]

    # no matches
    mock_glob = mocker.patch('analyze.introgression_configuration.glob.iglob',
                             side_effect=[[
                                 'test_prefix.fa',
                             ]])
    mock_log = mocker.patch('analyze.introgression_configuration.log')
    with pytest.raises(ValueError) as e:
        config.find_strains(['test_prefix_{strain}_{chrom}.fa'])
    assert ('Found no chromosome sequence files in '
            "['test_prefix_{strain}_{chrom}.fa']") in str(e)
    mock_glob.assert_called_with('test_prefix_*_*.fa')
    mock_log.info.assert_called_with('searching for test_prefix_*_*.fa')
    assert mock_log.debug.call_args_list == []

    # correct, with second test_strains, extra chromosomes
    mock_glob = mocker.patch('analyze.introgression_configuration.glob.iglob',
                             side_effect=[
                                 [
                                     'test_prefix_s1_cI.fa',
                                     'test_prefix_s2_cI.fa',
                                     'test_prefix_s2_cII.fa',
                                     'test_prefix.fa',
                                 ],
                                 ['test_prefix_cI_s3.fa']
                             ])
    mock_log = mocker.patch('analyze.introgression_configuration.log')
    config.find_strains(['test_prefix_{strain}_c{chrom}.fa',
                         'test_prefix_c{chrom}_{strain}.fa'])
    assert mock_glob.call_args_list == \
        [mocker.call('test_prefix_*_c*.fa'),
         mocker.call('test_prefix_c*_*.fa')]
    assert mock_log.info.call_args_list ==\
        [mocker.call('searching for test_prefix_*_c*.fa'),
         mocker.call('searching for test_prefix_c*_*.fa')]
    assert mock_log.debug.call_args_list == \
        [mocker.call("matched with ('s1', 'I')"),
         mocker.call("matched with ('s2', 'I')"),
         mocker.call("matched with ('s2', 'II')"),
         mocker.call("matched with ('s3', 'I')"),
         ]
    assert config.strains == ['s1', 's2', 's3']


def test_set_predict_files(config):
    with pytest.raises(ValueError) as e:
        config.set_predict_files('', '', '', '', '')
    assert 'No initial hmm file provided' in str(e)

    with pytest.raises(ValueError) as e:
        config.set_predict_files('init', '', '', '', '')
    assert 'No trained hmm file provided' in str(e)

    with pytest.raises(ValueError) as e:
        config.set_predict_files('init', 'trained', 'pos', 'prob', '')
    assert 'No alignment file provided' in str(e)

    with pytest.raises(ValueError) as e:
        config.set_predict_files('init', 'trained', 'pos', 'prob', 'align')
    assert '{prefix} not found in align' in str(e)

    with pytest.raises(ValueError) as e:
        config.set_predict_files('init', 'trained', 'pos', 'prob',
                                 'align{prefix}')
    assert '{strain} not found in align{prefix}' in str(e)

    with pytest.raises(ValueError) as e:
        config.set_predict_files('init', 'trained', 'pos', 'prob',
                                 'align{prefix}{strain}')
    assert '{chrom} not found in align{prefix}{strain}' in str(e)

    config.prefix = 'pre'
    config.set_predict_files('init', 'trained', 'pos', 'prob',
                             'align{prefix}{strain}{chrom}')
    assert config.hmm_initial == 'init'
    assert config.hmm_trained == 'trained'
    assert config.positions == 'pos'
    assert config.probabilities == 'prob'
    assert config.alignment == 'alignpre{strain}{chrom}'

    with pytest.raises(ValueError) as e:
        config.config = {'paths': {'analysis': {'hmm_initial': 'init'}}}
        config.set_predict_files('', '', '', '', '')
    assert 'No trained hmm file provided' in str(e)

    with pytest.raises(ValueError) as e:
        config.config = {'paths': {'analysis': {'hmm_initial': 'init',
                                                'hmm_trained': 'trained',
                                                'positions': 'pos'
                                                }}}
        config.set_predict_files('', '', '', '', '')
    assert 'No probabilities file provided' in str(e)

    with pytest.raises(ValueError) as e:
        config.config = {'paths': {'analysis': {'hmm_initial': 'init',
                                                'hmm_trained': 'trained',
                                                'positions': 'pos',
                                                'probabilities': 'prob'
                                                }}}
        config.set_predict_files('', '', '', '', '')
    assert 'No alignment file provided' in str(e)

    config.config = {'paths': {'analysis': {
        'hmm_initial': 'init',
        'hmm_trained': 'trained',
        'positions': 'pos',
        'probabilities': 'prob',
        'alignment': 'align{prefix}{strain}{chrom}'
    }}}
    config.set_predict_files('', '', '', '', '')

    assert config.hmm_initial == 'init'
    assert config.hmm_trained == 'trained'
    assert config.positions == 'pos'
    assert config.probabilities == 'prob'
    assert config.alignment == 'alignpre{strain}{chrom}'
