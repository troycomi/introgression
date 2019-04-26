import pytest
from click.testing import CliRunner
import analyze.main as main
import yaml
import logging as log


@pytest.fixture
def runner():
    return CliRunner()


def test_main_cli_configs(runner, mocker):
    result = runner.invoke(main.cli)
    assert result.exit_code == 0

    with runner.isolated_filesystem():
        mock_clean = mocker.patch('analyze.main.config_utils.clean_config',
                                  side_effect=lambda x: x)
        mock_echo = mocker.patch('analyze.main.click.echo_via_pager')
        mock_log_info = mocker.patch('analyze.main.log.info')
        mock_log_debug = mocker.patch('analyze.main.log.debug')
        mock_log_lvl = mocker.patch('analyze.main.log.basicConfig')

        with open('config1.yaml', 'w') as f:
            yaml.dump({'test': '123'}, f)
        with open('config2.yaml', 'w') as f:
            yaml.dump({'test': '23', 'test2': '34'}, f)

        result = runner.invoke(
            main.cli,
            '--config config1.yaml --config config2.yaml'.split())
        assert result.exit_code == 0
        mock_clean.assert_called_with(
            {'test': '23',
             'test2': '34'})

        # since no subcommand was called
        mock_echo.assert_called_once()

        mock_log_lvl.assert_called_once_with(level=log.WARNING)
        assert mock_log_info.call_args_list == [
            mocker.call('Verbosity set to WARNING'),
            mocker.call('Read in 2 config files')
        ]
        assert mock_log_debug.call_args_list == [
            mocker.call('Cleaned config:\ntest - 23\ntest2 - 34\n')
        ]


def test_main_cli_verbosity(runner, mocker):
    mock_log_info = mocker.patch('analyze.main.log.info')
    mock_log_lvl = mocker.patch('analyze.main.log.basicConfig')

    result = runner.invoke(
        main.cli,
        '-v')
    assert result.exit_code == 0
    mock_log_lvl.assert_called_once_with(level=log.CRITICAL)
    assert mock_log_info.call_args_list == [
        mocker.call('Verbosity set to CRITICAL'),
        mocker.call('Read in 0 config files')
    ]

    mock_log_info.reset_mock()
    mock_log_lvl.reset_mock()
    result = runner.invoke(
        main.cli,
        '-vv')
    assert result.exit_code == 0
    mock_log_lvl.assert_called_once_with(level=log.ERROR)
    assert mock_log_info.call_args_list == [
        mocker.call('Verbosity set to ERROR'),
        mocker.call('Read in 0 config files')
    ]

    mock_log_info.reset_mock()
    mock_log_lvl.reset_mock()
    result = runner.invoke(
        main.cli,
        '-vvv')
    assert result.exit_code == 0
    mock_log_lvl.assert_called_once_with(level=log.WARNING)
    assert mock_log_info.call_args_list == [
        mocker.call('Verbosity set to WARNING'),
        mocker.call('Read in 0 config files')
    ]

    mock_log_info.reset_mock()
    mock_log_lvl.reset_mock()
    result = runner.invoke(
        main.cli,
        '-vvvv')
    assert result.exit_code == 0
    mock_log_lvl.assert_called_once_with(level=log.INFO)
    assert mock_log_info.call_args_list == [
        mocker.call('Verbosity set to INFO'),
        mocker.call('Read in 0 config files')
    ]

    mock_log_info.reset_mock()
    mock_log_lvl.reset_mock()
    result = runner.invoke(
        main.cli,
        '-vvvvv')
    assert result.exit_code == 0
    mock_log_lvl.assert_called_once_with(level=log.DEBUG)
    assert mock_log_info.call_args_list == [
        mocker.call('Verbosity set to DEBUG'),
        mocker.call('Read in 0 config files')
    ]

    mock_log_info.reset_mock()
    mock_log_lvl.reset_mock()
    result = runner.invoke(
        main.cli,
        '-vvvvvv')
    assert result.exit_code == 0
    mock_log_lvl.assert_called_once_with(level=log.DEBUG)
    assert mock_log_info.call_args_list == [
        mocker.call('Verbosity set to DEBUG'),
        mocker.call('Read in 0 config files')
    ]
