import pytest
from misc.clean_config import clean_config, print_dict, clean_list
from yaml import load


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
