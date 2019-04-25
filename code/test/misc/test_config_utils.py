import pytest
from misc.config_utils import (clean_config, clean_list,
                               merge_lists, merge_dicts,
                               get_nested, check_wildcards,
                               get_states, validate)


def test_simple():
    config = {'base_name': 'base',
              'test': '__BASE_NAME__/test.txt',
              'test2': '__BASE_NAME__/__BASE_NAME__/test.txt',
              'test3': '__BASE_NAME____BASE_NAME__/test.txt',
              'test4': '__BASE_NAME__/__TEST__/test.txt',
              'test5': '__BASE_NAME__/__TEST__/test__DIGIT__.txt',
              'test6': '__BASE_NAME__//test__DIGIT__.txt',
              'test7': '//test///test////test/',
              'digit': 10,
              }
    config = clean_config(config)
    assert config == {'base_name': 'base',
                      'test': 'base/test.txt',
                      'test2': 'base/base/test.txt',
                      'test3': 'basebase/test.txt',
                      'test4': 'base/base/test.txt/test.txt',
                      'test5': 'base/base/test.txt/test10.txt',
                      'test6': 'base/test10.txt',
                      'test7': '/test/test/test/',
                      'digit': 10,
                      }


def test_circular():
    config = {'base_name': '__TEST2__',
              'test': '__BASE_NAME__/test.txt',
              'test2': '__BASE_NAME__/test.txt',
              }
    with pytest.raises(Exception) as e:
        clean_config(config)
    assert 'Failed to dereference all keys' in str(e)


def test_nest_dict():
    config = {
        'base_name': 'base',
        'test': '__BASE_NAME__/test.txt',
        'dict2': {
            'test2': '__BASE_NAME__/test2.txt',
            'base_name': 'base2',
            'dict3': {
                'base_name': 'base3',
                'test3': '__BASE_NAME__/test3.txt'
            },
            'dict4': {
                'test4': '__BASE_NAME__/test3.txt'
            }
        },
        'test5': '__BASE_NAME__/test5.txt',
        'dict5': {
            'base_name': '__BASE_NAME__',
            'test6': '__BASE_NAME__/test_5.txt'
        }
    }
    config = clean_config(config)
    assert config == {
        'base_name': 'base',
        'test': 'base/test.txt',
        'dict2': {
            'test2': 'base2/test2.txt',
            'base_name': 'base2',
            'dict3': {
                'base_name': 'base3',
                'test3': 'base3/test3.txt'
            },
            'dict4': {
                'test4': 'base2/test3.txt'
            }},
        'test5': 'base/test5.txt',
        'dict5': {
            'base_name': 'base',
            'test6': 'base/test_5.txt'
        }
    }


def test_clean_list():
    config = [
        'test',
        {'base': 'base',
         'test2': '__BASE__/test2'},
        {'base': 'base2',
         'test3': '__BASE__/test3'},
        [1, 2, 3],
        17
    ]
    assert clean_list(config) == [
        'test',
        {'base': 'base',
         'test2': 'base/test2'},
        {'base': 'base2',
         'test3': 'base2/test3'},
        [1, 2, 3],
        17
    ]

    with pytest.raises(Exception) as e:
        clean_list(['__NOT_FOUND__'])
    assert 'Failed to dereference list entry: "__NOT_FOUND__"' in str(e)

    assert clean_list(['__BASE__/test'], {'base': 'base'}) == ['base/test']


def test_merge_lists():
    assert merge_lists([], list('abc')) == list('abc')
    assert merge_lists(list('abc'), list('abc')) == list('abc')
    assert merge_lists(list('abc'), list('d')) == list('abcd')
    assert merge_lists(list('abc'), []) == list('abc')
    assert merge_lists(list('abc'), [1, 2, 3]) == ['a', 'b', 'c', 1, 2, 3]
    assert merge_lists([{'a': 1, 'b': 2}, 1], [{'a': 1, 'b': 2}, 3]) ==\
        [{'a': 1, 'b': 2}, 1, 3]
    assert merge_lists([{'a': 1, 'b': 2}, 1], [{'a': 1, 'b': 3}, 3]) ==\
        [{'a': 1, 'b': 2}, 1, {'a': 1, 'b': 3}, 3]


def test_merge_dicts():
    assert merge_dicts({1: 1, 2: 2}, {}) == {1: 1, 2: 2}
    assert merge_dicts({}, {1: 1, 2: 2}) == {1: 1, 2: 2}
    assert merge_dicts({1: 3}, {1: 1, 2: 2}) == {1: 1, 2: 2}
    # only new value type matters
    assert merge_dicts({1: 3, 2: {}}, {1: 1, 2: 2}) == {1: 1, 2: 2}
    # nested dict
    assert merge_dicts({1: 3, 2: {3: 4}}, {1: 1, 2: {3: 3}}) == \
        {1: 1, 2: {3: 3}}
    # nested list, just overwrite
    assert merge_dicts({1: 3, 2: {3: [1, 2]}}, {1: 1, 2: {3: [3, 4]}}) == \
        {1: 1, 2: {3: [3, 4]}}


def test_get_nested():
    assert get_nested({'a': 1}, 'a') == 1
    assert get_nested({'a': 1}, 'b') is None
    assert get_nested({'a': {'b': 2}}, 'a.b') == 2
    assert get_nested({'a': {'b': 2}}, 'a.c') is None
    assert get_nested({'a': {'b': {'c': 3}}}, 'a.b.c') == 3
    assert get_nested(None, 'key') is None


def test_check_wildcards(mocker):
    assert check_wildcards('{test}.txt', 'test')
    assert check_wildcards('{test}{string}.txt', 'test,string')

    mock_log = mocker.patch('misc.config_utils.log.exception')
    with pytest.raises(ValueError) as e:
        check_wildcards('test.txt', 'test')

    mock_log.assert_called_with('{test} not found in test.txt')
    assert '{test} not found in test.txt' in str(e)


def test_get_states():
    assert get_states({}) == ([], [])
    assert get_states(
        {
            'analysis_params': {
                'known_states': [
                    {'name': 'k1'},
                    {'name': 'k2'},
                    {'name': 'k3'},
                ],
                'unknown_states': [
                    {'name': 'u1'},
                    {'name': 'u2'},
                ]
            }
        }) == ('k1 k2 k3'.split(), 'u1 u2'.split())
    assert get_states(
        {
            'analysis_params': {
                'reference': {'name': 'ref'},
                'unknown_states': [
                    {'name': 'u1'},
                    {'name': 'u2'},
                ]
            }
        }) == ('ref'.split(), 'u1 u2'.split())
    assert get_states(
        {
            'analysis_params': {
                'reference': {'name': 'ref'},
                'known_states': [
                    {'name': 'k1'},
                    {'name': 'k2'},
                    {'name': 'k3'},
                ],
                'unknown_states': [
                    {'name': 'u1'},
                    {'name': 'u2'},
                ]
            }
        }) == ('ref k1 k2 k3'.split(), 'u1 u2'.split())


def test_validate(mocker):
    assert validate({}, '', '', 'test') == 'test'
    assert validate({'path': 'test'}, 'path', '') == 'test'
    assert validate({'path': 'test'}, 'path', '', '') == 'test'
    mock_log = mocker.patch('misc.config_utils.log.exception')
    with pytest.raises(ValueError) as e:
        validate({'path': 'test'}, 'path2', 'except', '')
    assert 'except' in str(e)
    mock_log.assert_called_with('except')
