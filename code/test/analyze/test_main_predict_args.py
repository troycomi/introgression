import pytest
from click.testing import CliRunner
import analyze.main as main
import yaml
from analyze import predict
from pathlib import Path


'''
Unit tests for the predict command of main.py when parameters are
provided by the arguments primarily
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
                    'chromosomes': 'I II III'.split(),
                }, f)

        result = runner.invoke(
            main.cli,
            '--config config.yaml predict --threshold 0.05')

        assert result.exit_code != 0
        assert str(result.exception) == 'No block file provided'
        assert mock_log.call_args_list == [
            mocker.call('Verbosity set to WARNING'),
            mocker.call('Read in 1 config file'),
            mocker.call('Found 3 chromosomes in config'),
            mocker.call("Threshold value is '0.05'")
        ]


def test_block(runner, mocker):
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
            '--config config.yaml predict --threshold viterbi '
            '--blocks blocks_{state}.txt'
        )

        assert result.exit_code != 0
        assert str(result.exception) == \
            'Unable to find strains in config and no test_strains provided'
        assert mock_log.call_args_list == [
            mocker.call('Verbosity set to WARNING'),
            mocker.call('Read in 1 config file'),
            mocker.call('Found 3 chromosomes in config'),
            mocker.call("Threshold value is 'viterbi'"),
            mocker.call("Output blocks file is 'blocks_{state}.txt'"),
            mocker.call("Prefix is 's1_s2'"),
        ]


def test_prefix(runner, mocker):
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
            '--config config.yaml predict --threshold viterbi '
            '--blocks blocks_{state}.txt --prefix s1_s2'
        )

        assert result.exit_code != 0
        assert str(result.exception) == \
            'Unable to find strains in config and no test_strains provided'
        assert mock_log.call_args_list == [
            mocker.call('Verbosity set to WARNING'),
            mocker.call('Read in 1 config file'),
            mocker.call('Found 3 chromosomes in config'),
            mocker.call("Threshold value is 'viterbi'"),
            mocker.call("Output blocks file is 'blocks_{state}.txt'"),
            mocker.call("Prefix is 's1_s2'"),
        ]


def test_test_strains(runner, mocker):
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

        Path('s1_chrI.fa').touch()
        Path('s1_chrII.fa').touch()
        Path('s1_chrIII.fa').touch()
        Path('s2_chrI.fa').touch()
        Path('s2_chrII.fa').touch()
        Path('s2_chrIII.fa').touch()

        result = runner.invoke(
            main.cli,
            '--config config.yaml predict --threshold viterbi '
            '--blocks blocks_{state}.txt --prefix s1_s2 '
            '--test-strains {strain}_chr{chrom}.fa'
        )

        assert result.exit_code != 0
        assert str(result.exception) == \
            'No initial hmm file provided'

        print(mock_log.call_args_list)
        assert mock_log.call_args_list == [
            mocker.call('Verbosity set to WARNING'),
            mocker.call('Read in 1 config file'),
            mocker.call('Found 3 chromosomes in config'),
            mocker.call("Threshold value is 'viterbi'"),
            mocker.call("Output blocks file is 'blocks_{state}.txt'"),
            mocker.call("Prefix is 's1_s2'"),
            mocker.call('searching for *_chr*.fa'),
            mocker.call('Found 1 test strain'),
            mocker.call('Found 2 unique strains'),
        ]


def test_outputs(runner, mocker):
    mock_log = mocker.patch('analyze.main.log.info')
    mock_calls = [
        mocker.call('Verbosity set to WARNING'),
        mocker.call('Read in 1 config file'),
        mocker.call('Found 3 chromosomes in config'),
        mocker.call("Threshold value is 'viterbi'"),
        mocker.call("Output blocks file is 'blocks_{state}.txt'"),
        mocker.call("Prefix is 's1_s2'"),
        mocker.call('No test_strains provided'),
        mocker.call('Found 2 unique strains'),
    ]

    with runner.isolated_filesystem():
        with open('config.yaml', 'w') as f:
            yaml.dump(
                {
                    'chromosomes': 'I II III'.split(),
                    'strains': 'str1 str2 str1'.split(),
                    'analysis_params': {
                        'known_states': [
                            {'name': 's1'},
                            {'name': 's2'}],
                    },
                }, f)

        mock_log.reset_mock()
        result = runner.invoke(
            main.cli,
            '--config config.yaml predict --threshold viterbi '
            '--blocks blocks_{state}.txt --prefix s1_s2 '
            '--hmm-initial hmm_init.txt '
        )

        assert result.exit_code != 0
        assert str(result.exception) == \
            'No trained hmm file provided'
        assert mock_log.call_args_list == mock_calls

    with runner.isolated_filesystem():
        with open('config.yaml', 'w') as f:
            yaml.dump(
                {
                    'chromosomes': 'I II III'.split(),
                    'strains': 'str1 str2 str1'.split(),
                    'analysis_params': {
                        'known_states': [
                            {'name': 's1'},
                            {'name': 's2'}],
                    },
                }, f)

        mock_log.reset_mock()
        result = runner.invoke(
            main.cli,
            '--config config.yaml predict --threshold viterbi '
            '--blocks blocks_{state}.txt --prefix s1_s2 '
            '--hmm-initial hmm_init.txt '
            '--hmm-trained hmm_trained.txt '
        )

        assert result.exit_code != 0
        assert str(result.exception) == \
            'No probabilities file provided'
        assert mock_log.call_args_list == mock_calls

    with runner.isolated_filesystem():
        with open('config.yaml', 'w') as f:
            yaml.dump(
                {
                    'chromosomes': 'I II III'.split(),
                    'strains': 'str1 str2 str1'.split(),
                    'analysis_params': {
                        'known_states': [
                            {'name': 's1'},
                            {'name': 's2'}],
                    },
                }, f)

        mock_log.reset_mock()
        result = runner.invoke(
            main.cli,
            '--config config.yaml predict --threshold viterbi '
            '--blocks blocks_{state}.txt --prefix s1_s2 '
            '--hmm-initial hmm_init.txt '
            '--hmm-trained hmm_trained.txt '
            '--probabilities probs.txt.gz '
        )

        assert result.exit_code != 0
        assert str(result.exception) == \
            'No alignment file provided'
        assert mock_log.call_args_list == mock_calls

    with runner.isolated_filesystem():
        with open('config.yaml', 'w') as f:
            yaml.dump(
                {
                    'chromosomes': 'I II III'.split(),
                    'strains': 'str1 str2 str1'.split(),
                    'analysis_params': {
                        'known_states': [
                            {'name': 's1'},
                            {'name': 's2'}],
                    },
                }, f)

        mock_log.reset_mock()
        result = runner.invoke(
            main.cli,
            '--config config.yaml predict --threshold viterbi '
            '--blocks blocks_{state}.txt --prefix s1_s2 '
            '--hmm-initial hmm_init.txt '
            '--hmm-trained hmm_trained.txt '
            '--probabilities probs.txt.gz '
            '--alignment {prefix}_{strain}_chr{chrom}.maf '
        )

        assert result.exit_code != 0
        assert str(result.exception) == \
            's1 did not provide an expected_length'
        assert mock_log.call_args_list == mock_calls + [
            mocker.call("Hmm_initial file is 'hmm_init.txt'"),
            mocker.call("Hmm_trained file is 'hmm_trained.txt'"),
            mocker.call("Positions file is 'None'"),
            mocker.call("Probabilities file is 'probs.txt.gz'"),
            mocker.call("Alignment file is 's1_s2_{strain}_chr{chrom}.maf'"),
            mocker.call("Only considering polymorphic sites"),
        ]

    mock_predict = mocker.patch.object(predict.Predictor, 'run_prediction')
    with runner.isolated_filesystem():
        with open('config.yaml', 'w') as f:
            yaml.dump(
                {
                    'chromosomes': 'I II III'.split(),
                    'strains': 'str1 str2 str1'.split(),
                    'analysis_params': {
                        'known_states': [
                            {'name': 's1'},
                            {'name': 's2'}],
                    },
                }, f)

        mock_log.reset_mock()
        result = runner.invoke(
            main.cli,
            '--config config.yaml predict --threshold viterbi '
            '--blocks blocks_{state}.txt --prefix s1_s2 '
            '--hmm-initial hmm_init.txt '
            '--hmm-trained hmm_trained.txt '
            '--probabilities probs.txt.gz '
            '--positions pos.txt.gz '
            '--alignment {prefix}_{strain}_chr{chrom}.maf '
        )

        assert result.exit_code == 0
        assert mock_log.call_args_list == mock_calls + [
            mocker.call("Hmm_initial file is 'hmm_init.txt'"),
            mocker.call("Hmm_trained file is 'hmm_trained.txt'"),
            mocker.call("Positions file is 'pos.txt.gz'"),
            mocker.call("Probabilities file is 'probs.txt.gz'"),
            mocker.call("Alignment file is 's1_s2_{strain}_chr{chrom}.maf'"),
            mocker.call("Only considering polymorphic sites"),
        ]
        mock_predict.called_once_with(True)

        mock_predict.reset_mock()
        result = runner.invoke(
            main.cli,
            '--config config.yaml predict --threshold viterbi '
            '--blocks blocks_{state}.txt --prefix s1_s2 '
            '--hmm-initial hmm_init.txt '
            '--hmm-trained hmm_trained.txt '
            '--probabilities probs.txt.gz '
            '--positions pos.txt.gz '
            '--alignment {prefix}_{strain}_chr{chrom}.maf '
            '--all-sites'
        )

        mock_predict.called_once_with(False)
