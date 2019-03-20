from analyze import extract_region as ex
import pytest
from io import StringIO
from misc.region_reader import Region_Reader


def compare_args(args, non_defaults):
    defaults = {'regions': [],
                'filename': '',
                'list_sort': False,
                'suppress_header': False}
    for k, v in non_defaults.items():
        defaults[k] = v

    assert len(defaults) == len(args)
    for k in defaults:
        assert defaults[k] == args[k]


def test_parse_args():
    # simple case
    args = ex.parse_args('--filename test.fa.gz r123'.split())
    compare_args(args, {'filename': 'test.fa.gz', 'regions': ['r123']})

    # add flags
    args = ex.parse_args('--filename test.fa.gz r123 \
                         --list_sort --suppress_header'.split())
    compare_args(args, {'filename': 'test.fa.gz',
                        'regions': ['r123'],
                        'list_sort': True,
                        'suppress_header': True})

    # more regions
    args = ex.parse_args('--filename test.fa.gz r123 r4 \
                         r246 r1 r2 r3'.split())
    compare_args(args, {'filename': 'test.fa.gz',
                        'regions': 'r123 r4 r246 r1 r2 r3'.split()})


def test_validate_args(mocker):
    # fail on filename existing
    mocker.patch('os.path.exists', return_value=False)
    with pytest.raises(ValueError) as e:
        ex.validate_args({'filename': 'test', 'suppress_header': False})
    assert 'test not found' in str(e)

    # fail on filename format
    mocker.patch('os.path.exists', return_value=True)
    with pytest.raises(ValueError) as e:
        ex.validate_args({'filename': 'test',
                          'suppress_header': False})

    # fail on pickle
    mocker.patch('os.path.exists', side_effect=[True, False])
    with pytest.raises(ValueError) as e:
        ex.validate_args({'filename': 'test.fa.gz',
                          'suppress_header': False})

    # fail on regions
    mocker.patch('os.path.exists', side_effect=[True, True])
    with pytest.raises(ValueError) as e:
        ex.validate_args({'filename': 'test.fa.gz',
                          'regions': ['z123'],
                          'suppress_header': False})

    # fail on regions
    mocker.patch('os.path.exists', side_effect=[True, True])
    with pytest.raises(ValueError) as e:
        ex.validate_args({'filename': 'test.fa.gz',
                          'regions': ['rz123'],
                          'suppress_header': False})

    # fail on regions
    mocker.patch('os.path.exists', side_effect=[True, True])
    with pytest.raises(ValueError) as e:
        ex.validate_args({'filename': 'test.fa.gz',
                          'regions': 'r123 12 z2'.split(),
                          'suppress_header': False})
    assert 'z2 could not be parsed' in str(e)

    # success!
    mocker.patch('os.path.exists', side_effect=[True, True])
    args, reader = ex.validate_args({'filename': 'test.fa.gz',
                                     'regions': 'r123 12 42'.split(),
                                     'suppress_header': False})
    assert reader.pickle == 'test.pkl'
    assert args['regions'] == [123, 12, 42]


@pytest.fixture
def r(mocker):
    mocker.patch('os.path.exists', side_effect=[True, True])
    return Region_Reader('test.fa.gz')


def test_decode_regions(r):
    r.index = {1: 2, 10: 3, 100: 4}

    # raise key error
    with pytest.raises(KeyError) as e:
        ex.decode_regions([1, 3], r, True)
    assert 'r3 not found in index' in str(e)

    result = ex.decode_regions([1, 1, 100, 10], r, True)
    assert result == [2, 2, 4, 3]

    result = ex.decode_regions([1, 1, 100, 10], r, False)
    assert result == [2, 2, 3, 4]


def test_write_regions(r, capsys):
    # empty regions
    r.region_reader = StringIO('')
    r.suppress_header = True
    r.num_lines = 1
    ex.write_regions(r, [])
    assert capsys.readouterr().out == ''

    # outside of file
    r.region_reader = StringIO('')
    r.suppress_header = False
    r.num_lines = 2
    ex.write_regions(r, [100])
    soe = capsys.readouterr()
    assert soe.out == ''
    assert soe.err == '100 outside of file\n'

    # outside of file on second
    r.region_reader = StringIO('a test\n')
    r.suppress_header = False
    r.num_lines = 2
    ex.write_regions(r, [0])
    soe = capsys.readouterr()
    assert soe.err == '0 outside of file\n'
    assert soe.out == 'a test\n'

    # normal, no header
    r.region_reader = StringIO('header\n'
                               'line 1\n'
                               'line 2\n'
                               'header\n'
                               'line 3\n')
    r.suppress_header = True
    r.num_lines = 1
    ex.write_regions(r, [0, 21, 0])
    soe = capsys.readouterr()
    assert soe.err == ''
    assert soe.out == 'line 1\nline 3\nline 1\n'

    # normal, with header
    r.region_reader = StringIO('head 1\n'
                               'line 1\n'
                               'line 2\n'
                               'head 2\n'
                               'line 3\n')
    r.suppress_header = False
    r.num_lines = 1
    ex.write_regions(r, [0, 21, 0])
    soe = capsys.readouterr()
    assert soe.out == 'head 1\nline 1\nhead 2\nline 3\nhead 1\nline 1\n'
    assert soe.err == ''
