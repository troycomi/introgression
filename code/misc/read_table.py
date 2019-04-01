import gzip
from typing import List, Dict, Tuple


def read_table_rows(fn: str,
                    sep: str,
                    header: bool = True,
                    key_ind: int = 0) -> Tuple[
                        Dict[str, Dict[str, List[str]]],
                        List[str]]:
    '''
    Read the text file of tabular data by rows
    fn: filename to read
    sep: the column delimiter
    header: flag to indicate a header is present
     If a header is provided, labels are returned from the first row.
     Return value becomes a dictionary of dictionaries, keyed first
     by the key_ind, then the column label
    key_ind: the column index to use as keys in output
    returns dictionary of rows keyed by key_ind and labels
    '''

    reader = None
    if fn.endswith('.gz'):
        reader = gzip.open(fn, 'rt')
    else:
        reader = open(fn, 'r')

    labels = None
    if header:
        labels = reader.readline()[:-1].split(sep)
    else:
        labels = None

    table = {}
    for line in reader:
        region = line[:-1].split(sep)
        if header:
            table[region[key_ind]] = \
                dict(zip(labels[:key_ind] + labels[key_ind + 1:],
                         region[:key_ind] + region[key_ind + 1:]))
        else:
            table[region[key_ind]] = region[:key_ind] + region[key_ind + 1:]

    reader.close()
    return table, labels


def read_table_columns(fn: str,
                       sep: str,
                       group_by: str = None,
                       **filter_output) -> Tuple[Dict, List[str]]:
    '''
    Reads sep delimited file to generate dictionary of columns, keyed by labels
    Optionally, a column to group by can be specified, changing the return
    value to a dictionary, keyed by column values, of dictionaries of columns
    Further, filters specified as column=value can limit output to just
    rows with matching values
    '''

    reader = None
    if fn.endswith('.gz'):
        reader = gzip.open(fn, 'rt')
    else:
        reader = open(fn, 'r')

    labels = reader.readline()[:-1].split(sep)

    filters = [None] * len(labels)
    if filter_output:
        for i, l in enumerate(labels):
            if l in filter_output:
                filters[i] = filter_output[l]
    else:
        filters = None

    if group_by not in labels:
        group_by = None
    else:
        group_by = labels.index(group_by)

    if group_by is None:
        table = dict(zip(labels, [[] for l in labels]))
        for line in reader:
            line = line[:-1].split(sep)
            if filters is None or \
                    all([v is None or v == l
                         for v, l in zip(filters, line)]):
                for label, l in zip(labels, line):
                    table[label].append(l)

    else:
        table = {}
        for line in reader:
            line = line[:-1].split(sep)
            if filters is None or \
                    all([v is None or v == l
                         for v, l in zip(filters, line)]):
                key = line[group_by]
                if key not in table:
                    table[key] = dict(zip(labels, [[] for l in labels]))
                for label, l in zip(labels, line):
                    table[key][label].append(l)

    reader.close()
    return table, labels
