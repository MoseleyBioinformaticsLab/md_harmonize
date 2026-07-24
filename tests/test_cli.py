import pytest
import md_harmonize.__main__ as cli
from md_harmonize import __version__


def test_version(mocker):
    mocker.patch('sys.argv', ['md_harmonize', '--version'])
    print_mock: mocker.MagicMock = mocker.patch('builtins.print')
    with pytest.raises(SystemExit):
        cli.main()
    print_mock.assert_called_with(__version__)
