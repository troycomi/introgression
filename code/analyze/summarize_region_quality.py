from __future__ import annotations
import bisect
import gzip
import numpy as np
import pickle
from contextlib import ExitStack
import click
import logging as log
from collections import Counter
from misc import read_fasta
from misc import read_table
from misc import seq_functions
from typing import List, Tuple, Dict
from analyze.introgression_configuration import Configuration


cen_starts = [151465, 238207, 114385, 449711, 151987, 148510,
              496920, 105586, 355629, 436307, 440129, 150828,
              268031, 628758, 326584, 555957]
cen_starts = [x-1 for x in cen_starts]

cen_ends = [151582, 238323, 114501, 449821, 152104, 148627,
            497038, 105703, 355745, 436425, 440246, 150947,
            268149, 628875, 326702, 556073]
cen_ends = [x-1 for x in cen_ends]

tel_coords = [1, 801, 229411, 230218,
              1, 6608, 812379, 813184,
              1, 1098, 315783, 316620,
              1, 904, 1524625, 1531933,
              1, 6473, 569599, 576874,
              1, 5530, 269731, 270161,
              1, 781, 1083635, 1090940,
              1, 5505, 556105, 562643,
              1, 7784, 439068, 439888,
              1, 7767, 744902, 745751,
              1, 807, 665904, 666816,
              1, 12085, 1064281, 1078177,
              1, 6344, 923541, 924431,
              1, 7428, 783278, 784333,
              1, 847, 1083922, 1091291,
              1, 7223, 942396, 948010]
tel_coords = [x-1 for x in tel_coords]

tel_left_starts = [tel_coords[i] for i in range(0, len(tel_coords), 4)]
tel_left_ends = [tel_coords[i] for i in range(1, len(tel_coords), 4)]
tel_right_starts = [tel_coords[i] for i in range(2, len(tel_coords), 4)]
tel_right_ends = [tel_coords[i] for i in range(3, len(tel_coords), 4)]

chromosomes = ('I II III IV V '
               'VI VII VIII IX X '
               'XI XII XIII XIV XV XVI').split()


def distance_from_telomere(start, end, chrm):

    assert start <= end, str(start) + ' ' + str(end)

    i = chromosomes.index(chrm)
    # region entirely on left arm
    if end <= cen_starts[i]:
        return start - tel_left_ends[i]
    # region entirely on right arm
    if start >= cen_ends[i]:
        return tel_right_starts[i] - end
    # region overlaps centromere: return minimum distance from either telomere
    return min(start - tel_left_ends[i], tel_right_starts[i] - end)


def distance_from_centromere(start, end, chrm):

    assert start <= end, str(start) + ' ' + str(end)

    i = chromosomes.index(chrm)
    # region entirely on left arm
    if end <= cen_starts[i]:
        return cen_starts[i] - end
    # region entirely on right arm
    if start >= cen_ends[i]:
        return start - cen_ends[i]
    # region overlaps centromere: return 0
    return 0


def write_region_summary_plus(fn, regions, fields):
    f = open(fn, 'w')
    f.write('region_id\t' + '\t'.join(fields) + '\n')
    keys = sorted(regions.keys(), key=lambda x: int(x[1:]))
    for region_id in keys:
        f.write(region_id + '\t')
        f.write('\t'.join([str(regions[region_id][field])
                           for field in fields]))
        f.write('\n')
    f.close()


def gap_columns(seqs):
    g = 0
    for i in range(len(seqs[0])):
        for seq in seqs:
            if seq[i] == '-':  # gp.gap_symbol:
                g += 1
                break
    return g


def longest_consecutive(s, c):
    max_consecutive = 0
    current_consecutive = 0
    in_segment = False
    for i in range(len(s)):
        if s[i] == c:
            current_consecutive += 1
            in_segment = True
        else:
            if in_segment:
                max_consecutive = max(max_consecutive, current_consecutive)
                current_consecutive = 0
            in_segment = False
    return max_consecutive


def masked_columns(seqs):
    # return two things:
    # - number of columns that are masked in any sequence
    # - above, but excluding columns with gaps
    num_seqs = len(seqs)
    num_sites = len(seqs[0])
    mask_total = 0
    mask_non_gap_total = 0
    for ps in range(num_sites):
        mask = False
        gap = False
        for s in range(num_seqs):
            if seqs[s][ps] == '-':  # gp.gap_symbol:
                gap = True
            elif seqs[s][ps] == 'x':  # gp.masked_symbol:
                mask = True
        if mask:
            mask_total += 1
            if not gap:
                mask_non_gap_total += 1
    return mask_total, mask_non_gap_total


def index_by_reference(ref_seq, seq):
    # return dictionary keyed by reference index, with value the
    # corresponding index in non-reference sequence

    d = {}
    ri = 0
    si = 0
    for i in range(len(ref_seq)):
        if ref_seq[i] != '-':  # gp.gap_symbol:
            d[ri] = si
            ri += 1
        if seq[i] != '-':  # gp.gap_symbol:
            si += 1
    return d


def num_sites_between(sites, start, end):
    # sites are sorted
    i = bisect.bisect_left(sites, start)
    j = bisect.bisect_right(sites, end)
    return j - i, sites[i:j]


class Summarizer():
    '''
    Summarize region quality of each region
    '''
    def __init__(self, configuration: Configuration):
        self.config = configuration

    def validate_arguments(self):
        '''
        Check that all required instance variables are set to perform a
        summarize run. Returns true if valid, raises value error otherwise
        '''
        args = [
            'chromosomes',
            'labeled_blocks',
            'quality_blocks',
            'masks',
            'alignment',
            'positions',
            'regions',
            'region_index',
            'known_states',
            'unknown_states',
            'states',
            'symbols'
        ]
        variables = self.config.__dict__
        for arg in args:
            if arg not in variables or variables[arg] is None:
                err = ('Failed to validate Summarizer, required argument '
                       f"'{arg}' was unset")
                log.exception(err)
                raise ValueError(err)

        reference = self.config.get('analysis_params.reference')
        if reference is None:
            err = f'Configuration did not specify a reference strain'
            log.exception(err)
            raise ValueError(err)

        return True

    def run(self, states: List[str] = None):
        '''
        Summarize region quality of each region for the states specified
        '''
        ref_ind, states = self.states_to_process(states)

        log.debug(f'reference index: {ref_ind}')
        log.debug(f'states to analyze: {states}')

        known_states = self.config.known_states
        log.debug(f'known_states {known_states}')

        analyzer = Sequence_Analyzer(
            self.config.masks,
            self.config.alignment,
            self.config.known_states,
            self.config.interval_states,
            self.config.chromosomes,
            self.config.symbols)

        log.debug(f'Sequence_Analyzer init with:')
        log.debug(f'masks: {self.config.masks}')
        log.debug(f'alignment: {self.config.alignment}')

        analyzer.build_masked_sites()

        for ind, state in enumerate(states):
            log.info(f'Working on state {state}')
            state_ind = self.config.states.index(state)

            with Position_Reader(
                    self.config.positions
                                 ) as positions,\
                    Region_Writer(
                        self.config.regions.format(state=state),
                        self.config.region_index.format(state=state),
                        known_states
                    ) as region_writer,\
                    Quality_Writer(
                        self.config.quality_blocks.format(state=state)
                    ) as quality_writer,\
                    ExitStack() as stack:

                progress_bar = None
                if self.config.log_file:
                    progress_bar = stack.enter_context(
                        click.progressbar(
                            length=len(self.config.chromosomes),
                            label=f'State {ind+1} of {len(states)}'))

                for chrm in self.config.chromosomes:
                    log.info(f'Working on chromosome {chrm}')
                    region = Region_Database(
                        self.config.labeled_blocks.format(state=state),
                        chrm,
                        known_states)

                    for strain, ps in positions.get_positions(region, chrm):
                        log.debug(f'{strain} {chrm}')

                        analyzer.process_alignment(ref_ind,
                                                   state_ind,
                                                   chrm,
                                                   strain,
                                                   ps,
                                                   region,
                                                   region_writer)

                    quality_writer.write_quality(region)

                    if progress_bar:
                        progress_bar.update(1)

    def states_to_process(self,
                          states: List[str] = None) -> Tuple[int,
                                                             List[str]]:
        '''
        Set the states to summarize to the values passed in.
        If no values are specified, run all states in config
        Checks if states are in config, warning if a state is not
        found and raising an error if none of the states are in config.
        '''
        reference = self.config.get('analysis_params.reference.name')
        ref_ind = self.config.states.index(reference)

        if states is None or states == []:
            to_process = self.config.states

        else:
            to_process = []
            for s in states:
                if s in self.config.states:
                    to_process.append(s)
                else:
                    log.warning(f"state '{s}' was not found as a state")

            if to_process == []:
                err = 'No valid states were found to process'
                log.exception(err)
                raise ValueError(err)

        return ref_ind, to_process


class Flag_Info():
    '''
    Collection of boolean flags for sequence summary
    '''
    def __init__(self):
        self.gap_any = None
        self.mask_any = None
        self.unseq_any = None
        self.hmm = None
        self.gap = None
        self.mask = None
        self.unseq = None
        self.match = None

    def initialize_flags(self, number_sequences: int, number_states: int):
        '''
        Initialize internal flags to np arrays of false
        '''
        self.gap_any = np.zeros((number_sequences), bool)
        self.mask_any = np.zeros((number_sequences), bool)
        self.unseq_any = np.zeros((number_sequences), bool)
        self.gap = np.zeros((number_sequences, number_states), bool)
        self.mask = np.zeros((number_sequences, number_states), bool)
        self.unseq = np.zeros((number_sequences, number_states), bool)
        self.match = np.zeros((number_sequences, number_states), bool)

    def add_sequence_flags(self, other: Flag_Info, state: int):
        '''
        Join the other flag info with this info by replacing values
        in the gap, unseq, and match arrays and performing OR with anys
        '''
        # only write the first time
        if state == 0:
            self.hmm = other.hmm

        self.gap_any = np.logical_or(self.gap_any, other.gap)
        self.unseq_any = np.logical_or(self.unseq_any, other.unseq)

        self.gap[:, state] = other.gap
        self.unseq[:, state] = other.unseq
        self.match[:, state] = other.match

    def add_mask_flags(self, other: Flag_Info, state: int):
        '''
        Join the other flag info with this by replacing values in mask and
        performing an OR with mask_any
        '''
        self.mask_any = np.logical_or(self.mask_any, other.mask)
        self.mask[:, state] = other.mask

    def encode_info(self,
                    master_ind: int,
                    predict_ind: int) -> str:
        '''
        Summarize info flags into a string. master_ind is the index of
        the master reference state. predict_ind is the index of the predicted
        state.  The return string is encoded for each position as:
         '-': if either master or predict has a gap
         '_': if either master or predict is masked
         '.': if any state has a match
         'b': both predict and master match
         'c': master matches but not predict
         'p': predict matches but not master
         'x': no other condition applies
         if the position is in the hmm_flag
          it will be capitalized for x, p, c, or b
        in order of precidence, e.g. if a position satisfies both '-' and '.',
        it will be '-'.
        '''

        if predict_ind >= self.match.shape[1]:
            return self.encode_unknown_info(master_ind)

        decoder = np.array(list('xXpPcCbB._-'))
        indices = np.zeros(self.match.shape[0], int)

        indices[self.match[:, predict_ind]] += 2  # x to p if true
        indices[self.match[:, master_ind]] += 4  # x to c, p to b
        indices[self.hmm] += 1  # to upper

        matches = np.all(self.match, axis=1)
        indices[matches] = 8  # .
        indices[np.any(
            self.mask[:, [master_ind, predict_ind]],
            axis=1)] = 9  # _
        indices[np.any(
            self.gap[:, [master_ind, predict_ind]],
            axis=1)] = 10  # -

        return ''.join(decoder[indices])

    def encode_unknown_info(self,
                            master_ind: int) -> str:
        '''
        Summarize info dictionary into a string for unknown state.
        master_ind is the index of the master reference state.
        The return string is encoded as each position as:
         '-': if any state has a gap
         '_': if any state has a mask
         '.': all states match
         'x': master matches
         'X': no other condition applies
        in order of precidence, e.g. if a position satisfies both '-' and '.',
        it will be '-'.
        '''

        # used with indices to decode result
        decoder = np.array(list('Xx._-'))
        indices = np.zeros(self.gap_any.shape, int)

        indices[self.match[:, master_ind]] = 1  # x
        matches = np.all(self.match, axis=1)
        indices[matches] = 2  # .
        indices[self.mask_any] = 3  # _
        indices[self.gap_any] = 4  # -

        return ''.join(decoder[indices])


class Sequence_Analyzer():
    '''
    Performs handling of masking, reading, and analyzing sequence data for
    summarizing the sequences
    '''
    def __init__(self,
                 mask_file: str,
                 alignment_file: str,
                 known_states: List,
                 interval_states: List,
                 chromosomes: List,
                 symbols: Dict):
        self.masks = mask_file
        self.alignments = alignment_file
        self.known_states = known_states
        self.interval_states = interval_states
        self.chromosomes = chromosomes
        self.symbols = symbols

    def build_masked_sites(self):
        '''
        Read in all intervals files and return dictionary of intervals,
        keyed first by chromosome, then state
        '''
        result = {}
        for chrom in self.chromosomes:
            result[chrom] = {}
            for state, name in zip(self.known_states, self.interval_states):
                result[chrom][state] = self.read_masked_sites(chrom, name)

        self.masked_sites = result

    def read_masked_sites(self, chrom: str, strain: str) -> np.array:
        filename = self.masks.format(chrom=chrom, strain=strain)
        intervals = self.read_masked_intervals(filename)
        sites = self.convert_intervals_to_sites(intervals)
        return sites

    def convert_intervals_to_sites(self,
                                   intervals: List[Tuple]) -> np.array:
        '''
        Given a list of start, end positions, returns a 1D np.array of sites
        contained in the intervals List
        convert_intervals_to_sites([(1, 2), (4, 6)]) -> [1, 2, 4, 5, 6]
        '''
        sites = []
        for start, end in intervals:
            sites += range(start, end + 1)
        return np.array(sites, dtype=int)

    def read_masked_intervals(self,
                              filename: str) -> List[Tuple[int, int]]:
        '''
        Read the interval file provided and return start and end sequences
        as a list of tuples of 2 ints
        '''
        with open(filename, 'r') as reader:
            reader.readline()  # header
            intervals = []
            for line in reader:
                line = line.split()
                intervals.append((int(line[0]), int(line[2])))

        return intervals

    def get_stats(self,
                  current_sequence,
                  other_sequence,
                  slice_start,
                  aligned_index_positions,
                  masked_site):
        '''
        Helper function to perform analyses on the sequences returning
        the results of seq_id_hmm, seq_id, and seq_id_unmasked
        '''

        # only alignment columns used by HMM (polymorphic, no
        # gaps in any strain)
        hmm_stats = self.seq_id_hmm(other_sequence,
                                    current_sequence,
                                    slice_start,
                                    aligned_index_positions)

        # all alignment columns, excluding ones with gaps in
        # these two sequences
        nongap_stats = seq_functions.seq_id(other_sequence,
                                            current_sequence)

        # all alignment columns, excluding ones with gaps or
        # masked bases or unsequenced in *these two sequences*
        nonmask_stats = self.seq_id_unmasked(other_sequence,
                                             current_sequence,
                                             slice_start,
                                             masked_site[0],
                                             masked_site[1])

        return hmm_stats, nongap_stats, nonmask_stats

    def seq_id_hmm(self,
                   seq1: np.array,
                   seq2: np.array,
                   offset: int,
                   include_sites: List[int]) -> Tuple[
                       int, int, Flag_Info]:
        '''
        Compare two sequences and provide statistics of their overlap
        considering only the included sites.
        Takes the two sequences to consider, an offset of the included sites,
        and a list of the included sites.
        Returns:
        -the total number of matching sites, where seq1[i] == seq2[i] and
         i is an element in included_sites - offset
        -the total number of sites considered in the included sites, e.g. where
         included_sites - offset >= 0 and < len(seq)
        -a Flag_Info object with:
         -gap: true where seq1 or seq1 == gap_symbol
         -unseq: true where seq1 or seq1 == unsequenced_symbol
         -hmm: true where hmm[i] is in included_sites - offset
         -match: true where seq1 == seq2, regardless of symbol
        '''
        sites = np.array(include_sites) - offset

        info = Flag_Info()
        info.gap = np.logical_or(seq1 == self.symbols['gap'],
                                 seq2 == self.symbols['gap'])
        info.unseq = np.logical_or(seq1 == self.symbols['unsequenced'],
                                   seq2 == self.symbols['unsequenced'])
        info.match = seq1 == seq2
        info.hmm = np.zeros(info.match.shape, bool)
        sites = sites[np.logical_and(sites < len(info.match), sites >= 0)]
        info.hmm[sites] = True

        total_sites = np.sum(info.hmm)
        total_match = np.sum(np.logical_and(info.hmm, info.match))

        # check all included are not gapped or skipped
        include_in_skip = np.logical_and(
            info.hmm, np.logical_or(
                info.unseq, info.gap))
        if np.any(include_in_skip):
            ind = np.where(include_in_skip)[0][0]
            err = ('Need to skip site specified as included '
                   f'seq1: {seq1[ind]}, seq2: {seq2[ind]}, index: {ind}')
            log.exception(err)
            raise ValueError(err)

        return total_match, total_sites, info

    def seq_id_unmasked(self,
                        seq1: np.array,
                        seq2: np.array,
                        offset: int,
                        exclude_sites1: List[int],
                        exclude_sites2: List[int]) -> Tuple[
                            int, int, Flag_Info]:
        '''
        Compare two sequences and provide statistics of their overlap
        considering only the included sites.
        Takes two sequences, an offset applied to each excluded sites list
        Returns:
         -total number of matching sites in non-excluded sites. A position is
          excluded if it is an element of either excluded site list - offset,
          or it is a gap or unsequenced symbol in either sequence.
         -total number of non-excluded sites
         A Flag_Info object with:
          -mask_flag: a boolean array that is true if the position is in
           either excluded list - offset
        '''
        info = Flag_Info()
        info.gap = np.logical_or(seq1 == self.symbols['gap'],
                                 seq2 == self.symbols['gap'])
        info.unseq = np.logical_or(seq1 == self.symbols['unsequenced'],
                                   seq2 == self.symbols['unsequenced'])
        exclude_sites1 = np.array(exclude_sites1)
        exclude_sites2 = np.array(exclude_sites2)

        # convert offset excluded sites to boolean array
        info.mask = np.zeros(seq1.shape, bool)
        if exclude_sites1.size != 0:
            sites1 = exclude_sites1 - offset
            sites1 = sites1[np.logical_and(sites1 < len(info.gap),
                                           sites1 >= 0)]
            info.mask[sites1] = True

        if exclude_sites2.size != 0:
            sites2 = exclude_sites2 - offset
            sites2 = sites2[np.logical_and(sites2 < len(info.gap),
                                           sites2 >= 0)]
            info.mask[sites2] = True

        # find sites that are not masked, gapped, or unsequenced
        sites = np.logical_not(
            np.logical_or(
                info.mask,
                np.logical_or(
                    info.gap, info.unseq)))

        # determine totals
        total_sites = np.sum(sites)
        total_match = np.sum(
            np.logical_and(
                seq1 == seq2,
                sites))

        return total_match, total_sites, info

    def process_alignment(self,
                          reference_index: int,
                          state_index: int,
                          chromosome: str,
                          strain: str,
                          positions: np.array,
                          region: Region_Database,
                          region_writer: Region_Writer):
        '''
        Analyze the alignment of a given strain, chromosome, and position.
        Result is stored in the provided region database
        '''
        sequences, alignments, masked_sites = self.get_indices(chromosome,
                                                               strain)

        # convert position indices from indices in master reference to
        # indices in alignment
        ps_align = alignments[reference_index][positions]

        for i, (r_id, start, end) in enumerate(region.get_entries(strain)):
            start, end = self.get_slice(start, end,
                                        alignments[reference_index],
                                        ps_align)

            info = Flag_Info()
            info.initialize_flags(
                end - start + 1,
                len(self.known_states))

            for ind, state in enumerate(self.known_states):
                hmm, nongap, nonmask = self.get_stats(
                    sequences[-1][start:end + 1],
                    sequences[ind][start:end + 1],
                    start,
                    ps_align,
                    (masked_sites[ind],
                     masked_sites[-1]))

                region.set_region(strain, i, state,
                                  hmm,
                                  nongap,
                                  nonmask)

                info.add_sequence_flags(hmm[2], ind)
                info.add_mask_flags(nonmask[2], ind)

            info_string = info.encode_info(reference_index, state_index)

            region_writer.write_header(r_id)
            region_writer.write_sequences(
                strain,
                alignments,
                sequences,
                (start, end))
            region_writer.write_info_string(info_string)

            # and keep track of each symbol count
            region.update_counts(strain, i, info_string)

    def get_indices(self, chromosome: str, strain: str) -> Tuple:
        '''
        Get the sequences and different indices for the provided
        chromosome and strain
        Returned tuple contains:
        -sequences as np.array
        -index alignment list of indices for each sequence
        -masked_sites, index aligned for each sequence
        '''
        _, sequences = read_fasta.read_fasta(
            self.alignments.format(chrom=chromosome, strain=strain))

        # to go from index in reference seq to index in alignment
        alignments = [
            self.index_alignment_by_reference(seq)
            for seq in sequences
        ]

        masked = self.read_masked_sites(chromosome, strain)

        masked_sites = [
            alignments[ind][self.masked_sites[chromosome][state]]
            for ind, state in enumerate(self.known_states)
        ] + [alignments[-1][masked]]  # for strain

        return sequences, alignments, masked_sites

    def index_alignment_by_reference(self, sequence: np.array) -> np.array:
        '''
        Find locations of non-gapped sites in sequence
        want a way to go from reference sequence coordinate to index in
        alignment
        '''
        return np.where(sequence != self.symbols['gap'])[0]

    def get_slice(self,
                  start: int,
                  end: int,
                  alignment: np.array,
                  ps_align: np.array) -> Tuple[int, int]:
        '''
        Get start and end positions of index aligned sequence.
        Checks that positions are valid (in ps_align), and raises
        value errors otherwise
        '''
        # index of start and end of region in aligned sequences
        slice_start, slice_end = alignment[[start, end]]

        if not np.in1d([slice_start, slice_end], ps_align).all():
            err = 'Slice not found in position alignment'
            log.exception(err)
            raise ValueError(err)

        return slice_start, slice_end


class Region_Database():
    '''
    Contains data and logic for regions data during summarizing
    '''
    def __init__(self,
                 labeled_file: str,
                 chromosome: str,
                 known_states: List[str]):
        '''
        Read in labeled file and store resulting table and labels
        '''
        self.info_string_symbols = list('.-_npbcxNPBCX')

        self.label_prefixes = ['match_nongap',
                               'num_sites_nongap',
                               'match_hmm',
                               'match_nonmask',
                               'num_sites_nonmask']

        self.data, self.labels = read_table.read_table_columns(
            labeled_file,
            sep='\t',
            group_by='strain',
            chromosome=chromosome)

        if self.labels[0] != 'region_id':
            err = 'Unexpected labeled format'
            log.exception(err)
            raise ValueError(err)

        for strain, data in self.data.items():
            n = len(data['region_id'])

            for s in known_states:
                for lbl in self.label_prefixes:
                    data[f'{lbl}_{s}'] = [0] * n

            for s in self.info_string_symbols:
                data['count_' + s] = [0] * n

        self.labels += [f'{lbl}_{st}' for lbl in self.label_prefixes
                        for st in known_states]
        self.labels += ['count_' + x for x in self.info_string_symbols]

    def has_strain(self, strain: str) -> bool:
        '''
        Checks if the strain is in this database
        '''
        return strain in self.data

    def get_entries(self, strain: str) -> Tuple[str, int, int]:
        '''
        returns an iterator for the region entries of the strain
        with region id (string), start (int) and end (int) positions
        '''
        if not self.has_strain(strain):
            err = f'Region Database does not contain strain {strain}'
            log.exception(err)
            raise ValueError(err)

        r_ids = self.data[strain]['region_id']
        starts = self.data[strain]['start']
        ends = self.data[strain]['end']
        for i in range(len(r_ids)):
            yield (r_ids[i], int(starts[i]), int(ends[i]))

    def set_region(self,
                   strain: str,
                   index: int,
                   state: str,
                   hmm, nongap, nonmask):
        '''
        Set the region state with the provided values.
        hmm, nongap and nonmask are tuples of the (match, total) values
        '''
        ds = self.data[strain]
        MATCH, TOTAL = 0, 1
        if hmm[TOTAL] is not None:
            ds['num_sites_hmm'][index] = hmm[TOTAL]

        ds[f'match_hmm_{state}'][index] = hmm[MATCH]

        ds[f'match_nongap_{state}'][index] = nongap[MATCH]
        ds[f'num_sites_nongap_{state}'][index] = nongap[TOTAL]

        ds[f'match_nonmask_{state}'][index] = nonmask[MATCH]
        ds[f'num_sites_nonmask_{state}'][index] = nonmask[TOTAL]

    def update_counts(self,
                      strain: str,
                      index: int,
                      info_string: str):
        '''
        Update the counts variables based on the provided info string
        '''
        counts = Counter(info_string)
        for sym in self.info_string_symbols:
            self.data[strain]['count_' + sym][index] = counts[sym]

    def generate_output(self):
        '''
        Yield lines for writing to the quality output file.
        To save memory, this effectively deletes the data structure!
        Outputs are tab delimited, sorted by region_id
        '''
        # reorganize output as list of tuples ordered by label
        output = []
        # have to store this as dict changes during iterations
        strains = list(self.data.keys())
        for strain in strains:
            # pop to limit memory usage
            d = self.data.pop(strain)
            output += list(zip(*[d[l] for l in self.labels]))

        # sort by region id (index 0, remove r #[1:])
        for entry in sorted(output, key=lambda e: int(e[0][1:])):
            yield '\t'.join([str(e) for e in entry]) + '\n'

    def generate_header(self):
        '''
        Generate a header line for the region database
        '''
        return '\t'.join(self.labels) + '\n'


class Region_Writer():
    '''
    Controls the writing of region files and indices
    '''
    def __init__(self,
                 region_file: str,
                 index_file: str,
                 known_states: List[str]):
        self.region_file = region_file
        self.index_file = index_file
        self.index = {}
        self.known_states = known_states

    def __enter__(self):
        self.region_writer = gzip.open(self.region_file, 'wt')

        return self

    def __exit__(self, type, value, traceback):
        self.region_writer.close()

        if traceback is None:
            # write index
            with open(self.index_file, 'wb') as index_writer:
                pickle.dump(self.index, index_writer)
            return True

        else:
            return False

    def write_header(self, region_id: str):
        '''
        Add a header line with the region id
        '''
        self.index[int(region_id[1:])] = self.region_writer.tell()
        self.region_writer.write(f'#{region_id}\n')

    def write_sequences(self,
                        strain: str,
                        alignments: List,
                        sequences: np.array,
                        indices: Tuple):
        '''
        Write sequences to region file
        '''
        start, end = indices
        names = self.known_states + [strain]
        for sj, name in enumerate(names):
            startj = bisect.bisect_left(alignments[sj], start)
            endj = bisect.bisect_left(alignments[sj], end)

            self.region_writer.write(f'> {name} {startj} {endj}\n')

            self.region_writer.write(''.join(
                sequences[sj][start:end+1]) + '\n')

    def write_info_string(self, info_string: str):
        '''
        Write info string with header to region file
        '''
        # write info string
        self.region_writer.write('> info\n')
        self.region_writer.write(info_string + '\n')


class Position_Reader():
    '''
    Read in position file, yielding positions until no longer on current
    chromosome
    '''

    def __init__(self, position_file):
        self.position_file = position_file
        self.last_position = 0

    def __enter__(self):
        self.reader = gzip.open(self.position_file, 'rt')
        return self

    def __exit__(self, type, value, traceback):
        self.reader.close()
        return traceback is None

    def get_positions(self,
                      region: Region_Database,
                      chromosome: str) -> Tuple[str, np.array]:
        self.reader.seek(self.last_position)
        line = self.next_line()
        while line != '':
            line = line.split('\t')

            chrm = line[1]
            if chrm != chromosome:
                break

            strain = line[0]
            if not region.has_strain(strain):
                line = self.next_line()
                continue

            yield strain, np.array(line[2:], dtype=int)

            line = self.next_line()

    def next_line(self) -> str:
        self.last_position = self.reader.tell()
        line = self.reader.readline()
        return line


class Quality_Writer():
    '''
    Control writing of quality file from region database
    '''
    def __init__(self, quality_filename):
        self.filename = quality_filename
        self.first_write = True

    def __enter__(self):
        self.writer = open(self.filename, 'w')
        return self

    def __exit__(self, type, value, traceback):
        self.writer.close()
        return traceback is None

    def write_quality(self, region: Region_Database):
        '''
        Writes header if needed and region database values
        '''
        if self.first_write is True:
            self.writer.write(region.generate_header())
            self.first_write = False

        for line in region.generate_output():
            self.writer.write(line)
