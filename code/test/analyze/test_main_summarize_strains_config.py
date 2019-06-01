import pytest
from click.testing import CliRunner
import analyze.main as main
import yaml
from analyze.summarize_strain_states import Strain_Summarizer


'''
Unit tests for the summarize_strains command of main.py when parameters are
provided by args
'''


@pytest.fixture
def runner():
    return CliRunner()


def test_empty(runner, mocker):
    mock_log = mocker.patch('analyze.main.log.info')
    with runner.isolated_filesystem():
        result = runner.invoke(
            main.cli,
            'summarize-strains')

        assert result.exit_code != 0
        assert str(result.exception) == 'No states specified'
        assert mock_log.call_args_list == [
            mocker.call('Verbosity set to WARNING'),
            mocker.call('Read in 0 config files'),
        ]


def test_states(runner, mocker):
    mock_log = mocker.patch('analyze.main.log.info')
    with runner.isolated_filesystem():
        with open('config.yaml', 'w') as f:
            yaml.dump(
                {
                    'analysis_params': {
                        'known_states': [
                            {'name': 's1'},
                            {'name': 's2'}],
                    },
                }, f)

        result = runner.invoke(
            main.cli,
            '--config config.yaml summarize-strains')

        assert result.exit_code != 0
        assert str(result.exception) == 'No introgressed_intermediate provided'
        assert mock_log.call_args_list == [
            mocker.call('Verbosity set to WARNING'),
            mocker.call('Read in 1 config file'),
        ]


def test_files(runner, mocker):
    mock_log = mocker.patch('analyze.main.log.info')
    mocker.patch.object(Strain_Summarizer, 'run')
    with runner.isolated_filesystem():
        with open('config.yaml', 'w') as f:
            yaml.dump(
                {
                    'analysis_params': {
                        'known_states': [
                            {'name': 's1'},
                            {'name': 's2'}],
                        'filter_threshold': 0.9,
                    },
                    'paths': {'analysis': {
                        'introgressed_intermediate': 'int_int_{state}.txt',
                        'ambiguous_intermediate': 'amb_int_{state}.txt',
                        'strain_info': 'strain_info.txt',
                        'state_counts': 'state_counts.txt'
                    }}
                }, f)
        result = runner.invoke(
            main.cli,
            '--config config.yaml summarize-strains '
        )

        assert result.exit_code == 0
        log = [
            mocker.call('Verbosity set to WARNING'),
            mocker.call('Read in 1 config file'),
            mocker.call("Introgressed intermediate file is "
                        "'int_int_{state}.txt'"),
            mocker.call("Ambiguous intermediate file is "
                        "'amb_int_{state}.txt'"),
            mocker.call("Strain information from 'strain_info.txt'"),
            mocker.call("State counts saved to 'state_counts.txt'"),
        ]
        assert mock_log.call_args_list == log
