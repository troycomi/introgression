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


def test_introgressed(runner, mocker):
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
            '--config config.yaml summarize-strains '
            '--introgress-inter int_int_{state}.txt '
        )

        assert result.exit_code != 0
        assert str(result.exception) == 'No ambiguous_intermediate provided'
        assert mock_log.call_args_list == [
            mocker.call('Verbosity set to WARNING'),
            mocker.call('Read in 1 config file'),
        ]


def test_ambiguous(runner, mocker):
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
            '--config config.yaml summarize-strains '
            '--introgress-inter int_int_{state}.txt '
            '--ambiguous-inter amb_int_{state}.txt '
        )

        assert result.exit_code != 0
        assert str(result.exception) == 'No strain_info provided'
        assert mock_log.call_args_list == [
            mocker.call('Verbosity set to WARNING'),
            mocker.call('Read in 1 config file'),
        ]


def test_info(runner, mocker):
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
            '--config config.yaml summarize-strains '
            '--introgress-inter int_int_{state}.txt '
            '--ambiguous-inter amb_int_{state}.txt '
            '--strain-info strain_info.txt '
        )

        assert result.exit_code != 0
        assert str(result.exception) == 'No state_counts provided'
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
                    },
                }, f)

        result = runner.invoke(
            main.cli,
            '--config config.yaml summarize-strains '
            '--introgress-inter int_int_{state}.txt '
            '--ambiguous-inter amb_int_{state}.txt '
            '--strain-info strain_info.txt '
            '--state-counts state_counts.txt '
        )

        assert result.exit_code == 0
        assert mock_log.call_args_list == [
            mocker.call('Verbosity set to WARNING'),
            mocker.call('Read in 1 config file'),
            mocker.call("Introgressed intermediate file is "
                        "'int_int_{state}.txt'"),
            mocker.call("Ambiguous intermediate file is "
                        "'amb_int_{state}.txt'"),
            mocker.call("Strain information from 'strain_info.txt'"),
            mocker.call("State counts saved to 'state_counts.txt'"),
        ]
