import pytest
from click.testing import CliRunner
import analyze.main as main
import yaml
from analyze.id_regions import ID_producer


'''
Unit tests for the id_regions command of main.py when all parameters are
provided by the config file
'''


@pytest.fixture
def runner():
    return CliRunner()


def test_empty(runner):
    result = runner.invoke(
        main.cli,
        'id-regions')
    assert result.exit_code != 0
    assert str(result.exception) == 'No chromosomes specified in config file!'


def test_chroms(runner, mocker):
    mock_log = mocker.patch('analyze.main.log.info')
    with runner.isolated_filesystem():
        with open('config.yaml', 'w') as f:
            yaml.dump(
                {
                    'chromosomes': 'I II III'.split()
                }, f)

        result = runner.invoke(
            main.cli,
            '--config config.yaml id-regions')

        assert result.exit_code != 0
        assert str(result.exception) == 'No states specified'
        assert mock_log.call_args_list == [
            mocker.call('Verbosity set to WARNING'),
            mocker.call('Read in 1 config file'),
            mocker.call('Found 3 chromosomes in config')
        ]


def test_states(runner, mocker):
    mock_log = mocker.patch('analyze.main.log.info')
    with runner.isolated_filesystem():
        with open('config.yaml', 'w') as f:
            yaml.dump(
                {
                    'chromosomes': 'I II III'.split(),
                    'analysis_params': {
                        'known_states': [
                            {'name': 's1'},
                            {'name': 's2'}],
                    },
                }, f)

        result = runner.invoke(
            main.cli,
            '--config config.yaml id-regions')

        assert result.exit_code != 0
        assert str(result.exception) == 'No block file provided'
        assert mock_log.call_args_list == [
            mocker.call('Verbosity set to WARNING'),
            mocker.call('Read in 1 config file'),
            mocker.call('Found 3 chromosomes in config'),
            mocker.call('Found 2 states to process'),
        ]


def test_block_file(runner, mocker):
    mock_log = mocker.patch('analyze.main.log.info')
    with runner.isolated_filesystem():
        with open('config.yaml', 'w') as f:
            yaml.dump(
                {
                    'chromosomes': 'I II III'.split(),
                    'analysis_params': {
                        'known_states': [
                            {'name': 's1'},
                            {'name': 's2'}],
                    },
                    'paths': {'analysis': {
                        'block_files': 'block_{state}.txt',
                    }}
                }, f)

        result = runner.invoke(
            main.cli,
            '--config config.yaml id-regions')

        assert result.exit_code != 0
        assert str(result.exception) == 'No labeled block file provided'
        assert mock_log.call_args_list == [
            mocker.call('Verbosity set to WARNING'),
            mocker.call('Read in 1 config file'),
            mocker.call('Found 3 chromosomes in config'),
            mocker.call('Found 2 states to process'),
            mocker.call("Input blocks file is 'block_{state}.txt'"),
        ]


def test_labeled_block_file(runner, mocker):
    mock_log = mocker.patch('analyze.main.log.info')
    with runner.isolated_filesystem():
        with open('config.yaml', 'w') as f:
            yaml.dump(
                {
                    'chromosomes': 'I II III'.split(),
                    'analysis_params': {
                        'known_states': [
                            {'name': 's1'},
                            {'name': 's2'}],
                    },
                    'paths': {'analysis': {
                        'block_files': 'block_{state}.txt',
                        'labeled_block_files': 'labeled_block_{state}.txt',
                    }}
                }, f)

        mock_id = mocker.patch.object(ID_producer, 'add_ids')

        result = runner.invoke(
            main.cli,
            '--config config.yaml id-regions')

        assert result.exit_code == 0
        assert mock_log.call_args_list == [
            mocker.call('Verbosity set to WARNING'),
            mocker.call('Read in 1 config file'),
            mocker.call('Found 3 chromosomes in config'),
            mocker.call('Found 2 states to process'),
            mocker.call("Input blocks file is 'block_{state}.txt'"),
            mocker.call("Output blocks file is 'labeled_block_{state}.txt'"),
        ]

        mock_id.called_once()
