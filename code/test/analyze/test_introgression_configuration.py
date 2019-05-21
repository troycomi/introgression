from analyze.introgression_configuration import (
    Configuration, Variable)
import pytest


@pytest.fixture()
def config():
    return Configuration()


def test_set(config):
    # unknown key
    with pytest.raises(ValueError) as e:
        config.set(asdf=None)
    assert 'Unknown variable to set: asdf' in str(e)

    # chromosomes
    with pytest.raises(ValueError) as e:
        config.set('chromosomes')
    assert 'No chromosomes provided' in str(e)

    config.config = {'chromosomes': ['I']}
    config.set('chromosomes')
    assert config.chromosomes == ['I']

    # log file
    config.set(log_file='')
    assert config.log_file is None

    config.set(log_file='test')
    assert config.log_file == 'test'

    config.config = {'paths': {'log_file': 'log'}}
    config.set(log_file='')
    assert config.log_file == 'log'

    config.set(log_file='test')
    assert config.log_file == 'test'


def test_set_state_files(config):
    state_files = [
        'blocks',
        'labeled_blocks',
        'quality_blocks',
        'introgressed',
        'introgressed_intermediate',
        'ambiguous',
        'ambiguous_intermediate',
    ]
    for sf in state_files:
        with pytest.raises(ValueError) as e:
            config.set(**{sf: None})
        assert f'No {sf} provided' in str(e)

        with pytest.raises(ValueError) as e:
            config.set(**{sf: 'test'})
        assert '{state} not found in test' in str(e)

        config.set(**{sf: 'test{state}'})
        assert config.__dict__[sf] == 'test{state}'

        config.config = {'paths': {'analysis': {sf: 'test2{state}'}}}
        config.set(**{sf: None})
        assert config.__dict__[sf] == 'test2{state}'


def test_set_nonwild_files(config):
    nonwild_files = [
        'hmm_initial',
        'hmm_trained',
        'positions'
    ]
    for nwf in nonwild_files:
        with pytest.raises(ValueError) as e:
            config.set(**{nwf: None})
        assert f'No {nwf} provided' in str(e)

        config.set(**{nwf: 'test'})
        assert config.__dict__[nwf] == 'test'

        config.config = {'paths': {'analysis': {nwf: 'test2'}}}
        config.set(**{nwf: None})
        assert config.__dict__[nwf] == 'test2'


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


def test_get_interval_states(config):
    assert config.get_interval_states() == []

    config.config = {
            'analysis_params': {
                'reference': {'name': 'ref'},
                'known_states': [
                    {'name': 'k1'},
                    {'name': 'k2'},
                    {'name': 'k3'},
                ],
            }
        }
    assert config.get_interval_states() == 'ref k1 k2 k3'.split()

    config.config = {
            'analysis_params': {
                'known_states': [
                    {'name': 'k1'},
                    {'name': 'k2'},
                    {'name': 'k3'},
                ],
            }
        }
    assert config.get_interval_states() == 'k1 k2 k3'.split()

    config.config = {
            'analysis_params': {
                'reference': {'name': 'ref'},
            }
        }
    assert config.get_interval_states() == 'ref'.split()

    config.config = {
            'analysis_params': {
                'reference': {'name': 'ref',
                              'interval_name': 'int_ref'},
                'known_states': [
                    {'name': 'k1',
                     'interval_name': 'i1'},
                    {'name': 'k2'},
                    {'name': 'k3',
                     'interval_name': 'i3'},
                ],
            }
        }
    assert config.get_interval_states() == 'int_ref i1 k2 i3'.split()


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

    config.set('states')
    assert config.known_states ==\
        'S288c CBS432 N_45 DBVPG6304 UWOPS91_917_1'.split()
    assert config.unknown_states ==\
        'unknown'.split()
    assert config.states ==\
        'S288c CBS432 N_45 DBVPG6304 UWOPS91_917_1 unknown'.split()

    config.set(states=[])
    assert config.known_states ==\
        'S288c CBS432 N_45 DBVPG6304 UWOPS91_917_1'.split()
    assert config.unknown_states ==\
        'unknown'.split()
    assert config.states ==\
        'S288c CBS432 N_45 DBVPG6304 UWOPS91_917_1 unknown'.split()

    config.set(states='testing 123'.split())
    assert config.states == ['testing', '123']

    config.config = {}

    with pytest.raises(ValueError) as e:
        config.set('states')
    assert 'No states specified' in str(e)


def test_set_threshold(config):
    with pytest.raises(ValueError) as e:
        config.set('threshold')
    assert 'No threshold provided' in str(e)

    config.config = {'analysis_params': {'threshold': 'asdf'}}
    with pytest.raises(ValueError) as e:
        config.set('threshold')
    assert 'Unsupported threshold value: asdf' in str(e)

    config.set(threshold=0.05)
    assert config.threshold == 0.05

    config.config = {'analysis_params':
                     {'threshold': 'viterbi'}}
    config.set('threshold')
    assert config.threshold == 'viterbi'


def test_set_prefix(config):
    config.known_states = ['s1']
    config.set('prefix')
    assert config.prefix == 's1'

    config.known_states = 's1 s2'.split()
    config.set('prefix')
    assert config.prefix == 's1_s2'

    config.set(prefix='prefix')
    assert config.prefix == 'prefix'

    config.known_states = []
    with pytest.raises(ValueError) as e:
        config.set('prefix')
    assert 'Unable to build prefix, no known states provided' in str(e)


def test_set_strains(config, mocker):
    mock_find = mocker.patch.object(Configuration, 'find_strains')

    config.set('strains')
    mock_find.called_with(None)

    with pytest.raises(ValueError) as e:
        config.config = {'paths': {'test_strains': ['test']}}
        config.set('strains')
    assert '{strain} not found in test' in str(e)

    with pytest.raises(ValueError) as e:
        config.config = {'paths': {'test_strains': ['test{strain}']}}
        config.set('strains')
    assert '{chrom} not found in test{strain}' in str(e)

    config.config = {'paths': {'test_strains':
                               ['test{strain}{chrom}']}}
    config.set('strains')
    mock_find.called_with(['test{strain}{chrom}'])

    config.set(strains='test{strain}{chrom}')
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


def test_set_alignment(config):
    config.set(alignment='align{strain}{chrom}')
    assert config.alignment == 'align{strain}{chrom}'

    with pytest.raises(AttributeError) as e:
        config.set(alignment='align{prefix}{strain}{chrom}')
    assert "'Configuration' object has no attribute 'prefix'" in str(e)

    config.prefix = 'prefix'
    config.set(alignment='align{prefix}{strain}{chrom}')
    assert config.alignment == 'alignprefix{strain}{chrom}'


def test_set_masked_file(config):
    with pytest.raises(ValueError) as e:
        config.set('masks')
    assert 'No masks provided' in str(e)

    with pytest.raises(ValueError) as e:
        config.set(masks='mask')
    assert '{strain} not found in mask' in str(e)

    with pytest.raises(ValueError) as e:
        config.set(masks='mask{strain}')
    assert '{chrom} not found in mask{strain}' in str(e)

    config.set(masks='mask{strain}{chrom}')
    assert config.masks == 'mask{strain}{chrom}'

    config.config = {'paths': {'analysis':
                               {'masked_intervals': 'msk{strain}{chrom}'}}}
    config.set('masks')
    assert config.masks == 'msk{strain}{chrom}'


def test_set_filter_threshold(config):
    with pytest.raises(ValueError) as e:
        config.set('filter_threshold')
    assert 'No filter_threshold provided' in str(e)

    config.set(filter_threshold=0.9)
    assert config.filter_threshold == 0.9

    config.config = {'analysis_params': {'filter_threshold': 0.8}}
    config.set('filter_threshold')
    assert config.filter_threshold == 0.8

    with pytest.raises(ValueError) as e:
        config.set(filter_threshold='test')
    assert 'Filter threshold is not a valid number' in str(e)


@pytest.fixture
def variable():
    return Variable('test')


def test_variable_init(variable):
    assert variable.name == 'test'
    assert variable.config_path == 'test'
    assert variable.nullable is False
    assert variable.wildcards is None

    var2 = Variable('test2', 'test.path', True, 'wild')
    assert var2.name == 'test2'
    assert var2.config_path == 'test.path'
    assert var2.nullable is True
    assert var2.wildcards == 'wild'


def test_variable_parse(variable):
    with pytest.raises(ValueError) as e:
        variable.parse(None)
    assert 'No test provided' in str(e)

    assert variable.parse('test', {}) == 'test'
    assert variable.parse(None, {'test': 'test'}) == 'test'

    variable.config_path = 'test.path'
    assert variable.parse(None, {'test': {'path': 'test'}}) == 'test'

    variable.nullable = True
    assert variable.parse(None) is None
    assert variable.parse('test') == 'test'

    variable.wildcards = 'state'
    with pytest.raises(ValueError) as e:
        variable.parse('test')
    assert '{state} not found in test' in str(e)

    assert variable.parse('test{state}') == 'test{state}'
