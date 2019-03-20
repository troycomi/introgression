from misc.region_reader import Region_Reader
import pytest
from pytest import approx
from io import StringIO
import numpy as np


def test_init(mocker):
    # fail on filename existing
    mocker.patch('os.path.exists', return_value=False)
    with pytest.raises(ValueError) as e:
        Region_Reader('test')
    assert 'test not found' in str(e)

    # fail on filename format
    mocker.patch('os.path.exists', return_value=True)
    with pytest.raises(ValueError) as e:
        Region_Reader('test')

    # fail on pickle
    mocker.patch('os.path.exists', side_effect=[True, False])
    with pytest.raises(ValueError) as e:
        Region_Reader('test.fa.gz')

    # success, with defaults
    mocker.patch('os.path.exists', side_effect=[True, True])
    r = Region_Reader('test.fa.gz')
    assert r.region_file == 'test.fa.gz'
    assert r.pickle == 'test.pkl'
    assert r.as_fa is False
    assert r.suppress_header is True
    assert r.num_lines == 14, 'Suppress header did not change num_lines'

    # non defaults
    mocker.patch('os.path.exists', side_effect=[True, True])
    r = Region_Reader('test1.fa.gz',
                      as_fa=True,
                      suppress_header=False,
                      num_lines=4)
    assert r.region_file == 'test1.fa.gz'
    assert r.pickle == 'test1.pkl'
    assert r.as_fa is True
    assert r.suppress_header is False
    assert r.num_lines == 4


@pytest.fixture
def r(mocker):
    mocker.patch('os.path.exists', side_effect=[True, True])
    return Region_Reader('test.fa.gz')


def test_read_region(r, capsys):
    # get fa, don't suppress header
    r.region_reader = StringIO('header 1\n'
                               'line 1\n'
                               '#header\n'
                               'header 2\n'
                               'line 2\n')
    r.num_lines = 2
    r.as_fa = True
    r.suppress_header = False
    r.index = {1: 16}
    header, seqs = r.read_region('r1')
    assert header == ['header 2']
    assert seqs == approx(np.asarray(['line 2']))
    soe = capsys.readouterr()
    assert soe.out == '#header\n'
    assert soe.err == ''

    # print, suppress header
    r.region_reader = StringIO('header 1\n'
                               'line 1\n'
                               '#header\n'
                               'header 2\n'
                               'line 2\n')
    r.num_lines = 2
    r.as_fa = False
    r.suppress_header = True
    r.read_region('1')
    soe = capsys.readouterr()
    assert soe.out == 'header 2\nline 2\n'
    assert soe.err == ''


def test_read_location(r, capsys):
    # get fa, don't suppress header
    r.region_reader = StringIO('header 1\n'
                               'line 1\n'
                               '#header\n'
                               'header 2\n'
                               'line 2\n')
    r.num_lines = 2
    r.as_fa = True
    r.suppress_header = False
    header, seqs = r.read_location(16)
    assert header == ['header 2']
    assert seqs == approx(np.asarray(['line 2']))
    soe = capsys.readouterr()
    assert soe.out == '#header\n'
    assert soe.err == ''

    # print, suppress header
    r.region_reader = StringIO('header 1\n'
                               'line 1\n'
                               '#header\n'
                               'header 2\n'
                               'line 2\n')
    r.num_lines = 2
    r.as_fa = False
    r.suppress_header = True
    r.read_location(16)
    soe = capsys.readouterr()
    assert soe.out == 'header 2\nline 2\n'
    assert soe.err == ''


def test_convert_region(r):
    with pytest.raises(ValueError) as e:
        r.convert_region('z123')
    assert 'z123 could not be parsed' in str(e)

    with pytest.raises(ValueError) as e:
        r.convert_region('zr123')
    assert 'zr123 could not be parsed' in str(e)

    assert r.convert_region('123') == 123
    assert r.convert_region('r123') == 123


def test_decode_region(r):
    index = {1: 2, 10: 3, 100: 4}
    r.index = index

    # raise key error
    with pytest.raises(KeyError) as e:
        r.decode_region(3)
    assert 'r3 not found in index' in str(e)

    assert r.decode_region(1) == 2
    assert r.decode_region(10) == 3
    assert r.decode_region(100) == 4


def test_encode_fa(r):
    # outside of file
    r.region_reader = StringIO('')
    with pytest.raises(ValueError) as e:
        r.encode_fa(100)
    assert '100 outside of file' in str(e)

    r.region_reader = StringIO('header 1\n'
                               'line 1\n'
                               'header 2\n'
                               'line 2\n')
    r.num_lines = 4
    header, seqs = r.encode_fa(0)
    assert header == ['header 1', 'header 2']
    assert seqs == approx(np.asarray(['line 1', 'line 2']))

    r.region_reader = StringIO('header 1\n'
                               'line 1\n'
                               'header 2\n'
                               'line 2\n')
    r.num_lines = 3
    header, seqs = r.encode_fa(0)
    assert header == ['header 1', 'header 2']
    assert seqs == approx(np.asarray(['line 1']))


def test_print_region(r, capsys):
    # outside of file
    r.region_reader = StringIO('')
    r.print_region(100)
    soe = capsys.readouterr()
    assert soe.out == ''
    assert soe.err == '100 outside of file\n'

    # outside of file on second position
    r.region_reader = StringIO('a test\n')
    r.num_lines = 2
    r.print_region(0)
    soe = capsys.readouterr()
    assert soe.err == '0 outside of file\n'
    assert soe.out == 'a test\n'

    # normal
    r.region_reader = StringIO('header\n'
                               'line 1\n'
                               'line 2\n'
                               'header\n'
                               'line 3\n')
    r.num_lines = 1
    r.print_region(0)
    r.region_reader = StringIO('header\n'
                               'line 1\n')
    soe = capsys.readouterr()
    assert soe.err == ''
    assert soe.out == 'header\n'

    # normal
    r.region_reader = StringIO('head 1\n'
                               'line 1\n'
                               'line 2\n'
                               'head 2\n'
                               'line 3\n')
    r.num_lines = 2
    r.print_region(0)
    soe = capsys.readouterr()
    assert soe.err == ''
    assert soe.out == 'head 1\nline 1\n'
