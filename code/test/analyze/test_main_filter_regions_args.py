import pytest
from click.testing import CliRunner
import analyze.main as main
import yaml
from analyze.filter_regions import Filterer


'''
Unit tests for the filter_regions command of main.py when parameters are
provided by args
'''


@pytest.fixture
def runner():
    return CliRunner()


def test_threshold(runner, mocker):
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
            '--config config.yaml filter-regions '
            '--thresh 0.9'
        )

        assert result.exit_code != 0
        assert str(result.exception) == 'No introgressed provided'
        assert mock_log.call_args_list == [
            mocker.call('Verbosity set to WARNING'),
            mocker.call('Read in 1 config file'),
            mocker.call('Filter threshold set to \'0.9\''),
        ]


def test_filter_files(runner, mocker):
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

        files = ['--introgress-filter int_{state}.txt',
                 '--introgress-inter int_int_{state}.txt',
                 '--ambiguous-filter amb_{state}.txt',
                 '--ambiguous-inter amb_int_{state}.txt',
                 '--filter-sweep filter.txt',
                 ]
        results = [
            'No introgressed_intermediate provided',
            'No ambiguous provided',
            'No ambiguous_intermediate provided',
            'No regions provided',  # sweep is not required
            'No regions provided',
        ]
        for i, expected in enumerate(results):
            result = runner.invoke(
                main.cli,
                '--config config.yaml filter-regions '
                '--thresh 0.9 ' +
                ' '.join(files[0:i+1])
            )

            assert result.exit_code != 0
            assert str(result.exception) == expected
            log = [
                mocker.call('Verbosity set to WARNING'),
                mocker.call('Read in 1 config file'),
                mocker.call('Filter threshold set to \'0.9\''),
            ]
            if i >= 3:
                log += [
                    mocker.call('Introgressed filtered file '
                                'is \'int_{state}.txt\''),
                    mocker.call('Introgressed intermediate file '
                                'is \'int_int_{state}.txt\''),
                    mocker.call('Ambiguous filtered file '
                                'is \'amb_{state}.txt\''),
                    mocker.call('Ambiguous intermediate file '
                                'is \'amb_int_{state}.txt\''),
                ]
            # filter sweep line is not printed if it is unset
            if i == 4:
                log += [mocker.call('Filter sweep file is \'filter.txt\'')]
            assert mock_log.call_args_list == log
            mock_log.reset_mock()


def test_region_files(runner, mocker):
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
                '--config config.yaml filter-regions '
                '--thresh 0.9 '
                '--introgress-filter int_{state}.txt '
                '--introgress-inter int_int_{state}.txt '
                '--ambiguous-filter amb_{state}.txt '
                '--ambiguous-inter amb_int_{state}.txt '
                '--filter-sweep filter.txt '
                '--region region_{state}.gz '
            )

            assert result.exit_code != 0
            assert str(result.exception) == 'No region_index provided'
            assert mock_log.call_args_list == [
                mocker.call('Verbosity set to WARNING'),
                mocker.call('Read in 1 config file'),
                mocker.call('Filter threshold set to \'0.9\''),
                mocker.call('Introgressed filtered file '
                            'is \'int_{state}.txt\''),
                mocker.call('Introgressed intermediate file '
                            'is \'int_int_{state}.txt\''),
                mocker.call('Ambiguous filtered file '
                            'is \'amb_{state}.txt\''),
                mocker.call('Ambiguous intermediate file '
                            'is \'amb_int_{state}.txt\''),
                mocker.call('Filter sweep file is \'filter.txt\'')
            ]
            mock_log.reset_mock()

            result = runner.invoke(
                main.cli,
                '--config config.yaml filter-regions '
                '--thresh 0.9 '
                '--introgress-filter int_{state}.txt '
                '--introgress-inter int_int_{state}.txt '
                '--ambiguous-filter amb_{state}.txt '
                '--ambiguous-inter amb_int_{state}.txt '
                '--filter-sweep filter.txt '
                '--region region_{state}.gz '
                '--region-index region_{state}.pkl '
            )

            assert result.exit_code != 0
            assert str(result.exception) == 'No quality_blocks provided'
            assert mock_log.call_args_list == [
                mocker.call('Verbosity set to WARNING'),
                mocker.call('Read in 1 config file'),
                mocker.call('Filter threshold set to \'0.9\''),
                mocker.call('Introgressed filtered file '
                            'is \'int_{state}.txt\''),
                mocker.call('Introgressed intermediate file '
                            'is \'int_int_{state}.txt\''),
                mocker.call('Ambiguous filtered file '
                            'is \'amb_{state}.txt\''),
                mocker.call('Ambiguous intermediate file '
                            'is \'amb_int_{state}.txt\''),
                mocker.call('Filter sweep file is \'filter.txt\''),
                mocker.call('Region file is \'region_{state}.gz\''),
                mocker.call('Region index file is \'region_{state}.pkl\'')
            ]


def test_quality_files(runner, mocker):
    mock_log = mocker.patch('analyze.main.log.info')
    mock_run = mocker.patch.object(Filterer, 'run')
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
                '--config config.yaml filter-regions '
                '--thresh 0.9 '
                '--introgress-filter int_{state}.txt '
                '--introgress-inter int_int_{state}.txt '
                '--ambiguous-filter amb_{state}.txt '
                '--ambiguous-inter amb_int_{state}.txt '
                '--filter-sweep filter.txt '
                '--region region_{state}.gz '
                '--region-index region_{state}.pkl '
                '--quality quality_{state}.txt '
            )

            assert result.exit_code == 0
            assert mock_log.call_args_list == [
                mocker.call('Verbosity set to WARNING'),
                mocker.call('Read in 1 config file'),
                mocker.call('Filter threshold set to \'0.9\''),
                mocker.call('Introgressed filtered file '
                            'is \'int_{state}.txt\''),
                mocker.call('Introgressed intermediate file '
                            'is \'int_int_{state}.txt\''),
                mocker.call('Ambiguous filtered file '
                            'is \'amb_{state}.txt\''),
                mocker.call('Ambiguous intermediate file '
                            'is \'amb_int_{state}.txt\''),
                mocker.call('Filter sweep file is \'filter.txt\''),
                mocker.call('Region file is \'region_{state}.gz\''),
                mocker.call('Region index file is \'region_{state}.pkl\''),
                mocker.call('Quality file is \'quality_{state}.txt\''),
                mocker.call('Threshold sweep with: []'),
            ]

            mock_run.assert_called_once_with([])


def test_thresholds_files(runner, mocker):
    mock_log = mocker.patch('analyze.main.log.info')
    mock_run = mocker.patch.object(Filterer, 'run')
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
                '--config config.yaml filter-regions '
                '--thresh 0.9 '
                '--introgress-filter int_{state}.txt '
                '--introgress-inter int_int_{state}.txt '
                '--ambiguous-filter amb_{state}.txt '
                '--ambiguous-inter amb_int_{state}.txt '
                '--filter-sweep filter.txt '
                '--region region_{state}.gz '
                '--region-index region_{state}.pkl '
                '--quality quality_{state}.txt '
                '1.0 .99 .98 .1 .01'
            )

            assert result.exit_code == 0
            assert mock_log.call_args_list == [
                mocker.call('Verbosity set to WARNING'),
                mocker.call('Read in 1 config file'),
                mocker.call('Filter threshold set to \'0.9\''),
                mocker.call('Introgressed filtered file '
                            'is \'int_{state}.txt\''),
                mocker.call('Introgressed intermediate file '
                            'is \'int_int_{state}.txt\''),
                mocker.call('Ambiguous filtered file '
                            'is \'amb_{state}.txt\''),
                mocker.call('Ambiguous intermediate file '
                            'is \'amb_int_{state}.txt\''),
                mocker.call('Filter sweep file is \'filter.txt\''),
                mocker.call('Region file is \'region_{state}.gz\''),
                mocker.call('Region index file is \'region_{state}.pkl\''),
                mocker.call('Quality file is \'quality_{state}.txt\''),
                mocker.call('Threshold sweep with: '
                            '[1.0, 0.99, 0.98, 0.1, 0.01]'),
            ]

            mock_run.assert_called_once_with([1.0, 0.99, 0.98, 0.1, 0.01])
