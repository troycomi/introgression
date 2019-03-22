#!/usr/bin/env python3
import argparse
from misc.region_reader import Region_Reader
from typing import List, Tuple


def main():
    '''
    Main method to read in arguments from stdin and perform lookup with
    Region_Reader
    '''
    args = parse_args()
    args, reader = validate_args(args)
    with reader:
        locations = decode_regions(args['regions'],
                                   reader, args['list_sort'])
        write_regions(reader, locations)


def parse_args(args: List[str] = None) -> argparse.Namespace:
    '''
    Read in input arguments or the supplied list of strings
    Returns a dictionary of options
    '''
    parser = argparse.ArgumentParser(
        description='retrieve regions from indexed file.')

    parser.add_argument('regions',
                        nargs='+',
                        help='one or more region ids to retrieve')
    parser.add_argument('--filename',
                        required=True,
                        help='fa.gz file to look for regions')
    parser.add_argument('--list_sort',
                        action='store_true',
                        help='sort regions by the input order. Defualt sort by'
                        ' disk location')
    parser.add_argument('--suppress_header',
                        action='store_true',
                        help='suppress printing of header line in stdout')

    return vars(parser.parse_args(args))


def validate_args(args: argparse.Namespace) -> Tuple[argparse.Namespace,
                                                     Region_Reader]:
    '''
    Performs checks and conversions of input, raises ValueErrors if invalid
    '''
    reader = Region_Reader(args['filename'],
                           as_fa=False,
                           suppress_header=args['suppress_header'],
                           num_lines=15)

    args['regions'] = [reader.convert_region(r) for r in args['regions']]

    return args, reader


def decode_regions(regions: List[int],
                   reader: Region_Reader,
                   retain_sort: bool) -> List[int]:
    '''
    Converts list of regions to file locations based on index dictionary
    Retain_sort controls if the output list order is determined by the
    region order or the disk location (i.e. values of index dict)
    '''

    result = [reader.decode_region(r) for r in regions]

    if retain_sort:
        return result
    else:
        return sorted(result)


def write_regions(reader: Region_Reader, locations: List[int]) -> None:
    '''
    Writes the regions specified by index to stdout
    If print_header is false, ignore first line after location
    '''
    for location in locations:
        reader.read_location(location)


if __name__ == '__main__':
    main()
