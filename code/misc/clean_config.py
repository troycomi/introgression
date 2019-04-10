import re
from copy import copy
from typing import Dict, List


'''
clean_config.py

Helper functions for performing replacements on yaml config files
'''


def clean_config(config: Dict,
                 valid_replacements: Dict[str, str] = None) -> Dict:
    '''
    Performs subsitution of variables in string recursively replacing
    strings of the form __.+__ with the matching key.  Nested variables
    with the same name replace parent values.
    config is the possibly nested dict with values to replace
    valid_replacements are the valid entries for performing replacements
    '''
    result = {}
    if valid_replacements is None:
        valid_replacements = dict()
    len_values = len(config)
    while config:
        # want to look at valid replacements first,
        # to possibly replace their values
        keys = config.keys()
        keys = list([k for k in keys if k in valid_replacements] +
                    [k for k in keys if k not in valid_replacements])

        for key in keys:
            value = config[key]
            if isinstance(value, str):
                value = replace_entry(value, valid_replacements)
                if value is None:
                    continue  # don't remove
                result[key] = value
                valid_replacements[key] = value

            elif isinstance(value, dict):
                result[key] = clean_config(value,
                                           copy(valid_replacements))

            elif isinstance(value, list):
                result[key] = clean_list(value,
                                         valid_replacements)

            else:
                result[key] = value
                valid_replacements[key] = str(value)

            config.pop(key)

        if len_values == len(config):
            raise Exception('Failed to dereference all keys, remaining '
                            f'values are:\n {print_dict(config)}')

        len_values = len(config)

    return result


def clean_list(config: List,
               valid_replacements: Dict[str, str] = None) -> List:
    '''
    Performs substitution on list of config objects
    '''
    result = []
    for value in config:
        if isinstance(value, str):
            output = replace_entry(value, valid_replacements)
            if output is None:
                raise Exception(f'Failed to dereference list entry: "{value}"')
            result.append(output)

        elif isinstance(value, list):
            result.append(clean_list(value, valid_replacements))

        elif isinstance(value, dict):
            result.append(clean_config(value, copy(valid_replacements)))

        else:
            result.append(value)

    return result


def replace_entry(value: str, valid_replacements: Dict[str, str]) -> str:
    '''
    Replace instances of __.+__ with the key in valid_replacements
    If valid replacements is none or the key is not found, return None
    Else return the (possibly) substituted string with all instances of /+
    replaced with / (common in path replacements)
    '''
    replacements = re.findall('__(.+?)__', value)
    for replacement in set(replacements):
        replace = replacement.lower()
        if valid_replacements is None or replace not in valid_replacements:
            return None
        value = re.sub(f'__{replacement}__',
                       valid_replacements[replace],
                       value)
    return re.sub('/+', '/', value)


def print_dict(d: Dict, lvl: int = 0) -> str:
    '''
    Return pretty representation of the dictionary d.
    lvl is the starting amount to indent the line
    '''
    result = ''
    for k, v in d.items():
        if isinstance(v, dict):
            result += '  ' * lvl + f'{k} -\n'
            result += print_dict(v, lvl+1)
        elif isinstance(v, list):
            result += '  ' * lvl + f'{k} -\n'
            result += print_list(v, lvl+1)
        else:
            result += '  ' * lvl + f'{k} - {v}\n'
    return result


def print_list(l: List, lvl: int = 0) -> str:
    '''
    Return pretty representation of the list l.
    lvl is the startin amount to indent the line
    '''
    result = ''
    for i, v in enumerate(l):
        result += '  ' * lvl + f'{i}:\n'
        if isinstance(v, dict):
            result += print_dict(v, lvl+1)
        elif isinstance(v, list):
            result += print_list(v, lvl+1)
        else:
            result += '  ' * lvl + f'{v},\n'
    return result
