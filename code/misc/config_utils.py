import re
from copy import copy
from typing import Dict, List
import logging as log


'''
config_utils.py

Helper functions for working with yaml config files
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
        if isinstance(v, dict):
            result += '  ' * lvl + f'{i}:\n' + print_dict(v, lvl+1)
        elif isinstance(v, list):
            result += '  ' * lvl + f'{i}:\n' + print_list(v, lvl+1)
        else:
            result += '  ' * lvl + f'{i}:\t{v},\n'
    return result


def merge_dicts(parent: Dict, new: Dict) -> Dict:
    '''
    Merge the new dict into parent.  Existing items are overwritten,
    dicts are merged recursively, lists are combined as sets.
    '''

    for k, v in new.items():
        if k in parent:
            if isinstance(v, dict):
                parent[k] = merge_dicts(parent[k], v)

            else:
                parent[k] = v
        else:
            parent[k] = v

    return parent


def merge_lists(parent: List, new: List) -> List:
    '''
    Merge new list into parent.  If new item isn't in list, add it.
    Overwriting and nesting is not supported as it seems ill-defined.
    '''
    for i, v in enumerate(new):
        if v not in parent:
            parent.append(v)

    return parent


def get_nested(config: Dict, keys: str):
    '''
    Return the value of the nested keys, or none if the key is invalid
    keys is a period separated list of keys as a string
    '''
    if config is None:
        return None
    keys = keys.split('.')
    value = config
    try:
        for k in keys:
            value = value[k]
    except KeyError:
        return None
    return value


def check_wildcards(path: str, wildcards: str) -> bool:
    '''
    Check if the supplied path contains all required wildcards
    wildcards are provided as a comma separated list string
    returns true if all wildcards are present in path, e.g. {wildcard} in path
    else raises a ValueError with the unfound wildcard
    '''
    for wildcard in wildcards.split(','):
        if f'{{{wildcard}}}' not in path:
            err = f'{{{wildcard}}} not found in {path}'
            log.exception(err)
            raise ValueError(err)

    return True


def validate(config: Dict,
             path: str,
             exception: str,
             value: str = None):
    '''
    validate the supplied value, raising exception if no value is found
    config: the config dictionary to lookup
    path: the path in nested config dict
    exception: string to display if no value is found
    value: starting value. values of None or '' will cause lookup into config
    '''

    if value is None or value == '':
        value = get_nested(config, path)

    if value is None:
        log.exception(exception)
        raise ValueError(exception)

    return value
