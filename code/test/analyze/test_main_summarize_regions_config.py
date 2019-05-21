import pytest
from click.testing import CliRunner
import analyze.main as main
import yaml
from analyze.summarize_region_quality import Summarizer


'''
Unit tests for the summarize_regions command of main.py when parameters are
provided by config
'''


@pytest.fixture
def runner():
    return CliRunner()


def test_empty(runner, mocker):
    mock_log = mocker.patch('analyze.main.log.info')
    with runner.isolated_filesystem():
        with open('config.yaml', 'w') as f:
            yaml.dump(
                {
                    'analysis_params': {
                        'reference':
                            {'name': 'r1'},
                        'known_states': [
                            {'name': 's1'},
                            {'name': 's2'}],
                        'unknown_states': [
                            {'name': 'u1'}]
                    },
                    'chromosomes': 'I II III'.split(),
                }, f)

        result = runner.invoke(
            main.cli,
            '--config config.yaml summarize-regions')

        assert result.exit_code != 0
        assert str(result.exception) == 'No labeled_blocks provided'
        assert mock_log.call_args_list == [
            mocker.call('Verbosity set to WARNING'),
            mocker.call('Read in 1 config file'),
            mocker.call('Found 3 chromosomes in config'),
        ]


def test_labeled(runner, mocker):
    mock_log = mocker.patch('analyze.main.log.info')
    with runner.isolated_filesystem():
        with open('config.yaml', 'w') as f:
            yaml.dump(
                {
                    'analysis_params': {
                        'reference':
                            {'name': 'r1'},
                        'known_states': [
                            {'name': 's1'},
                            {'name': 's2'}],
                        'unknown_states': [
                            {'name': 'u1'}]
                    },
                    'chromosomes': 'I II III'.split(),
                    'paths': {'analysis': {
                        'labeled_blocks': '{state}lbl.txt',
                    }}
                }, f)

        result = runner.invoke(
            main.cli,
            '--config config.yaml summarize-regions')

        assert result.exit_code != 0
        assert str(result.exception) == 'No quality_blocks provided'
        assert mock_log.call_args_list == [
            mocker.call('Verbosity set to WARNING'),
            mocker.call('Read in 1 config file'),
            mocker.call('Found 3 chromosomes in config'),
            mocker.call("Labeled blocks file is '{state}lbl.txt'"),
        ]


def test_quality(runner, mocker):
    mock_log = mocker.patch('analyze.main.log.info')
    with runner.isolated_filesystem():
        with open('config.yaml', 'w') as f:
            yaml.dump(
                {
                    'analysis_params': {
                        'reference':
                            {'name': 'r1'},
                        'known_states': [
                            {'name': 's1'},
                            {'name': 's2'}],
                        'unknown_states': [
                            {'name': 'u1'}]
                    },
                    'chromosomes': 'I II III'.split(),
                    'paths': {'analysis': {
                        'labeled_blocks': '{state}lbl.txt',
                        'quality_blocks': '{state}qual.txt',
                    }}
                }, f)

        result = runner.invoke(
            main.cli,
            '--config config.yaml summarize-regions')

        assert result.exit_code != 0
        assert str(result.exception) == 'No masks provided'
        assert mock_log.call_args_list == [
            mocker.call('Verbosity set to WARNING'),
            mocker.call('Read in 1 config file'),
            mocker.call('Found 3 chromosomes in config'),
            mocker.call("Labeled blocks file is '{state}lbl.txt'"),
            mocker.call("Quality file is '{state}qual.txt'"),
        ]


def test_masked(runner, mocker):
    mock_log = mocker.patch('analyze.main.log.info')
    with runner.isolated_filesystem():
        with open('config.yaml', 'w') as f:
            yaml.dump(
                {
                    'analysis_params': {
                        'reference':
                            {'name': 'r1'},
                        'known_states': [
                            {'name': 's1'},
                            {'name': 's2'}],
                        'unknown_states': [
                            {'name': 'u1'}]
                    },
                    'chromosomes': 'I II III'.split(),
                    'paths': {'analysis': {
                        'labeled_blocks': '{state}lbl.txt',
                        'quality_blocks': '{state}qual.txt',
                        'masked_intervals': '{strain}_{chrom}mask.txt',
                    }}
                }, f)

        result = runner.invoke(
            main.cli,
            '--config config.yaml summarize-regions')

        assert result.exit_code != 0
        assert str(result.exception) == 'No alignment provided'
        assert mock_log.call_args_list == [
            mocker.call('Verbosity set to WARNING'),
            mocker.call('Read in 1 config file'),
            mocker.call('Found 3 chromosomes in config'),
            mocker.call("Labeled blocks file is '{state}lbl.txt'"),
            mocker.call("Quality file is '{state}qual.txt'"),
            mocker.call("Mask file is '{strain}_{chrom}mask.txt'"),
        ]


def test_alignment(runner, mocker):
    mock_log = mocker.patch('analyze.main.log.info')
    with runner.isolated_filesystem():
        with open('config.yaml', 'w') as f:
            yaml.dump(
                {
                    'analysis_params': {
                        'reference':
                            {'name': 'r1'},
                        'known_states': [
                            {'name': 's1'},
                            {'name': 's2'}],
                        'unknown_states': [
                            {'name': 'u1'}]
                    },
                    'chromosomes': 'I II III'.split(),
                    'paths': {'analysis': {
                        'labeled_blocks': '{state}lbl.txt',
                        'quality_blocks': '{state}qual.txt',
                        'masked_intervals': '{strain}_{chrom}mask.txt',
                        'alignment': '{strain}_{chrom}_align.txt',
                    }}
                }, f)

        # no prefix
        result = runner.invoke(
            main.cli,
            '--config config.yaml summarize-regions')

        assert result.exit_code != 0
        assert str(result.exception) == 'No positions provided'
        assert mock_log.call_args_list == [
            mocker.call('Verbosity set to WARNING'),
            mocker.call('Read in 1 config file'),
            mocker.call('Found 3 chromosomes in config'),
            mocker.call("Labeled blocks file is '{state}lbl.txt'"),
            mocker.call("Quality file is '{state}qual.txt'"),
            mocker.call("Mask file is '{strain}_{chrom}mask.txt'"),
            mocker.call("Alignment file is '{strain}_{chrom}_align.txt'"),
        ]


def test_positions(runner, mocker):
    mock_log = mocker.patch('analyze.main.log.info')
    with runner.isolated_filesystem():
        with open('config.yaml', 'w') as f:
            yaml.dump(
                {
                    'analysis_params': {
                        'reference':
                            {'name': 'r1'},
                        'known_states': [
                            {'name': 's1'},
                            {'name': 's2'}],
                        'unknown_states': [
                            {'name': 'u1'}]
                    },
                    'chromosomes': 'I II III'.split(),
                    'paths': {'analysis': {
                        'labeled_blocks': '{state}lbl.txt',
                        'quality_blocks': '{state}qual.txt',
                        'masked_intervals': '{strain}_{chrom}mask.txt',
                        'alignment': '{strain}_{chrom}_align.txt',
                        'positions': 'pos.txt',
                    }}
                }, f)

        result = runner.invoke(
            main.cli,
            '--config config.yaml summarize-regions')

        assert result.exit_code != 0
        assert str(result.exception) == 'No regions provided'
        assert mock_log.call_args_list == [
            mocker.call('Verbosity set to WARNING'),
            mocker.call('Read in 1 config file'),
            mocker.call('Found 3 chromosomes in config'),
            mocker.call("Labeled blocks file is '{state}lbl.txt'"),
            mocker.call("Quality file is '{state}qual.txt'"),
            mocker.call("Mask file is '{strain}_{chrom}mask.txt'"),
            mocker.call("Alignment file is '{strain}_{chrom}_align.txt'"),
            mocker.call("Positions file is 'pos.txt'"),
        ]


def test_region(runner, mocker):
    mock_log = mocker.patch('analyze.main.log.info')
    with runner.isolated_filesystem():
        with open('config.yaml', 'w') as f:
            yaml.dump(
                {
                    'analysis_params': {
                        'reference':
                            {'name': 'r1'},
                        'known_states': [
                            {'name': 's1'},
                            {'name': 's2'}],
                        'unknown_states': [
                            {'name': 'u1'}]
                    },
                    'chromosomes': 'I II III'.split(),
                    'paths': {'analysis': {
                        'labeled_blocks': '{state}lbl.txt',
                        'quality_blocks': '{state}qual.txt',
                        'masked_intervals': '{strain}_{chrom}mask.txt',
                        'alignment': '{strain}_{chrom}_align.txt',
                        'positions': 'pos.txt',
                        'regions': 'region{state}.gz',
                    }}
                }, f)

        result = runner.invoke(
            main.cli,
            '--config config.yaml summarize-regions')

        assert result.exit_code != 0
        assert str(result.exception) == 'No region_index provided'
        assert mock_log.call_args_list == [
            mocker.call('Verbosity set to WARNING'),
            mocker.call('Read in 1 config file'),
            mocker.call('Found 3 chromosomes in config'),
            mocker.call("Labeled blocks file is '{state}lbl.txt'"),
            mocker.call("Quality file is '{state}qual.txt'"),
            mocker.call("Mask file is '{strain}_{chrom}mask.txt'"),
            mocker.call("Alignment file is '{strain}_{chrom}_align.txt'"),
            mocker.call("Positions file is 'pos.txt'"),
        ]


def test_run(runner, mocker):
    mock_log = mocker.patch('analyze.main.log.info')
    mock_summarize = mocker.patch.object(Summarizer, 'run')
    with runner.isolated_filesystem():
        with open('config.yaml', 'w') as f:
            yaml.dump(
                {
                    'analysis_params': {
                        'reference':
                            {'name': 'r1'},
                        'known_states': [
                            {'name': 's1'},
                            {'name': 's2'}],
                        'unknown_states': [
                            {'name': 'u1'}]
                    },
                    'chromosomes': 'I II III'.split(),
                    'paths': {'analysis': {
                        'labeled_blocks': '{state}lbl.txt',
                        'quality_blocks': '{state}qual.txt',
                        'masked_intervals': '{strain}_{chrom}mask.txt',
                        'alignment': '{strain}_{chrom}_align.txt',
                        'positions': 'pos.txt',
                        'regions': 'region{state}.gz',
                        'region_index': 'ind{state}.pkl',
                    }}
                }, f)

        result = runner.invoke(
            main.cli,
            '--config config.yaml summarize-regions')

        assert result.exit_code == 0
        assert mock_log.call_args_list == [
            mocker.call('Verbosity set to WARNING'),
            mocker.call('Read in 1 config file'),
            mocker.call('Found 3 chromosomes in config'),
            mocker.call("Labeled blocks file is '{state}lbl.txt'"),
            mocker.call("Quality file is '{state}qual.txt'"),
            mocker.call("Mask file is '{strain}_{chrom}mask.txt'"),
            mocker.call("Alignment file is '{strain}_{chrom}_align.txt'"),
            mocker.call("Positions file is 'pos.txt'"),
            mocker.call("Region file is 'region{state}.gz'"),
            mocker.call("Region index file is 'ind{state}.pkl'"),
        ]
        mock_summarize.assert_called_once_with([])
