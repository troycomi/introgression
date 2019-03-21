import pickle
import gzip
import os
import sys
import numpy as np


class Region_Reader():
    def __init__(self, region_file,
                 as_fa=False,
                 suppress_header=True,
                 num_lines=14):
        '''
        Checks for valid filename and existance of corresponding pickle
        as_fa: if true will return headers and sequences as read_fasta does
        suppress_header: if true will not print the #region_id line
        num_lines: number of lines to print once seek to index. Does not
        include region header line.
        '''
        if not os.path.exists(region_file):
            raise ValueError(f'{region_file} not found')

        if region_file[-6:] != '.fa.gz':
            raise ValueError(f'{region_file} expected to be .fa.gz')

        pickle = region_file[:-6] + '.pkl'
        if not os.path.exists(pickle):
            raise ValueError(f'{pickle} not found with region file')

        self.region_file = region_file
        self.pickle = pickle
        self.as_fa = as_fa
        self.suppress_header = suppress_header
        self.num_lines = num_lines

    def __enter__(self):
        self.region_reader = gzip.open(self.region_file, 'rt')
        self.index = pickle.load(open(self.pickle, 'rb'))
        return self

    def __exit__(self, type, value, traceback):
        self.region_reader.close()

    def __repr__(self):
        return (
            f'region_file = {self.region_file}\n'
            f'pickle = {self.pickle}\n'
            f'as_fa = {self.as_fa}\n'
            f'suppress_header = {self.suppress_header}\n'
            f'num_lines = {self.num_lines}\n'
        )

    def read_region(self, region_name):
        '''
        read the supplied region name, either printing to stdout or returning
        (headers, seqs) tuple depending on as_fa value
        '''
        region = self.convert_region(region_name)
        location = self.decode_region(region)
        return self.read_location(location)

    def read_location(self, location):
        '''
        helper method used in extract_region for directly handling locations
        '''
        self.region_reader.seek(location)

        if self.suppress_header is True:
            self.region_reader.readline()
        else:
            print(self.region_reader.readline(), end='')

        if self.as_fa:
            return self.encode_fa(location)
        else:
            self.print_region(location)

    def convert_region(self, region_name):
        '''
        Checks that region is a digit that starts with r
        If so, returns the integer value of the region for decoding
        '''
        r = region_name
        if r[0] == 'r':
            r = r[1:]
        if not r.isdigit():
            raise ValueError(f'{region_name} could not be parsed')
        return int(r)

    def decode_region(self, region_number):
        '''
        Convert region to disk location.
        Raises key error if region doesn't exist
        '''
        try:
            result = self.index[region_number]
        except KeyError as e:
            raise KeyError(f'r{e} not found in index')

        return result

    def yield_fa(self):
        '''
        repeatedly yield tuples of region, headers, sequences from fa file
        assumes file position starts at header for region
        suppress_header is taken as true (will not print)
        '''
        while True:
            region = self.region_reader.readline()[1:-1]
            try:
                header, seq = self.encode_fa(region)
                yield (region, header, seq)
            except ValueError:
                break

    def encode_fa(self, location):
        '''
        Reads the region file entry and returns headers, seqs
        Assumes even numbered lines are headers, odd are sequences
        '''
        headers = []
        seqs = []
        for i in range(self.num_lines):
            line = self.region_reader.readline()
            if line == '':
                raise ValueError(f'{location} outside of file')
            if i % 2 == 0:  # header
                headers.append(line[:-1])
            else:
                seqs.append(line[:-1])

        return headers, np.asarray(seqs)

    def print_region(self, location):
        '''
        reads the region file entry, printing to stdout
        '''
        for i in range(self.num_lines):
            line = self.region_reader.readline()
            if line == '':
                print(f'{location} outside of file', file=sys.stderr)
                break
            else:
                print(line, end='')
