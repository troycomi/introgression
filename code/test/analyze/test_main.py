import pytest
from click.testing import CliRunner
import analyze.main as main
import yaml


@pytest.fixture
def runner():
    return CliRunner()


def test_main_cli(runner, mocker):
    result = runner.invoke(main.cli)
    assert result.exit_code == 0

    with runner.isolated_filesystem():
        clean = mocker.patch('analyze.main.config_utils.clean_config',
                             return_value=dict())
        with open('config1.yaml', 'w') as f:
            yaml.dump({'test': '123'}, f)
        with open('config2.yaml', 'w') as f:
            yaml.dump({'test': '23', 'test2': '34'}, f)

        result = runner.invoke(
            main.cli,
            '--config config1.yaml --config config2.yaml'.split())
        assert result.exit_code == 0
        clean.assert_called_with(
            {'test': '23',
             'test2': '34'})
