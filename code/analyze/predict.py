import copy
import gzip
import os
import itertools
import click
from collections import defaultdict, Counter
from hmm import hmm_bw
from sim import sim_predict
from sim import sim_process
import numpy as np
from typing import List, Dict, Tuple, TextIO
from contextlib import ExitStack
import logging as log
from misc.read_fasta import read_fasta
from analyze.introgression_configuration import Configuration


# TODO remove gp references for symbols. pass args or fold into object?
def process_predict_args(arg_list: List[str]) -> Dict:
    '''
    Parses arguments from argv, producing dictionary of parsed values
    '''

    import global_params as gp
    d = {}
    i = 0

    d['tag'] = arg_list[i]
    i += 1

    d['improvement_frac'] = float(arg_list[i])
    i += 1

    d['threshold'] = 'viterbi'
    try:
        d['threshold'] = float(arg_list[i])
    except ValueError:
        pass
    i += 1

    # expected length of introgressed tracts and fraction of sequence
    # introgressed
    expected_tract_lengths = {}
    expected_frac = {}

    d['known_states'] = gp.alignment_ref_order
    for ref in gp.alignment_ref_order[1:]:
        expected_tract_lengths[ref] = float(arg_list[i])
        i += 1
        expected_frac[ref] = float(arg_list[i])
        i += 1

    d['unknown_states'] = []
    while i < len(arg_list):
        state = arg_list[i]
        d['unknown_states'].append(state)
        i += 1
        expected_tract_lengths[state] = float(arg_list[i])
        i += 1
        expected_frac[state] = float(arg_list[i])
        i += 1

    d['states'] = d['known_states'] + d['unknown_states']

    expected_frac[d['states'][0]] = 0
    expected_frac[d['states'][0]] = 1 - sum(expected_frac.values())
    d['expected_frac'] = expected_frac

    # calculate these based on remaining bases, but after we know
    # which chromosome we're looking at
    expected_tract_lengths[d['states'][0]] = 0
    d['expected_tract_lengths'] = expected_tract_lengths
    d['expected_num_tracts'] = {}
    d['expected_bases'] = {}

    return d


class Predictor():
    '''
    Predictor class
    Stores all variables needed to run an HMM prediction
    '''
    def __init__(self, configuration: Configuration):
        self.config = configuration

    def run_prediction(self, only_poly_sites=True):
        '''
        Run prediction with this predictor object
        '''
        self.config.validate_predict_arguments()

        hmm_builder = HMM_Builder(self.config)
        hmm_builder.set_expected_values()
        self.emission_symbols = \
            hmm_builder.update_emission_symbols(len(self.config.known_states))

        with open(self.config.hmm_initial, 'w') as initial, \
                open(self.config.hmm_trained, 'w') as trained, \
                gzip.open(self.config.probabilities, 'wt') as probabilities, \
                ExitStack() as stack:

            self.write_hmm_header(initial)
            self.write_hmm_header(trained)

            if self.config.positions is not None:
                positions = stack.enter_context(
                    gzip.open(self.config.positions, 'wt'))
            else:
                positions = None

            block_writers = {state:
                             stack.enter_context(
                                 open(self.config.blocks.format(
                                     state=state), 'w'))
                             for state in
                             self.config.states}
            for writer in block_writers.values():
                self.write_blocks_header(writer)

            counter = 0
            total = len(self.config.chromosomes) * len(self.config.strains)
            progress_bar = None
            if self.config.log_file:  # logging to file
                progress_bar = stack.enter_context(
                    click.progressbar(
                        length=total,
                        label='Running prediction'))

            for chrom in self.config.chromosomes:
                for strain in self.config.strains:
                    counter += 1
                    log.info(f'working on: {strain} {chrom} '
                             f'({counter} of {total})')

                    # get sequences and encode
                    alignment_file = self.config.alignment.format(
                        strain=strain, chrom=chrom)

                    if not os.path.exists(alignment_file):
                        log.info(f'skipping, file {alignment_file} not found')
                    else:
                        hmm_initial, hmm_trained, pos = hmm_builder.run_hmm(
                            alignment_file, only_poly_sites)

                        self.write_hmm(hmm_initial, initial, strain, chrom)
                        self.write_hmm(hmm_trained, trained, strain, chrom)

                        # process and threshold hmm result
                        predicted_states, probs = self.process_path(
                            hmm_trained)
                        state_blocks = self.convert_to_blocks(predicted_states)

                        if positions is not None:
                            self.write_positions(pos, positions, strain, chrom)

                        for state, block in state_blocks.items():
                            self.write_blocks(block,
                                              pos,
                                              block_writers[state],
                                              strain,
                                              chrom,
                                              state)

                        self.write_state_probs(probs, probabilities,
                                               strain, chrom)

                    if progress_bar:
                        progress_bar.update(1)

    def write_hmm_header(self, writer: TextIO) -> None:
        '''
        Write the header line for an hmm file to the provided textIO object
        Output is tab delimited with:
        strain chromosome initial_probs emissions transitions
        '''

        writer.write('strain\tchromosome\t')

        states = self.config.states

        writer.write('\t'.join(
            [f'init_{s}' for s in states] +  # initial
            [f'emis_{s}_{symbol}'
             for s in states
             for symbol in self.emission_symbols] +  # emissions
            [f'trans_{s1}_{s2}'
             for s1 in states
             for s2 in states]))  # transitions

        writer.write('\n')

    def write_hmm(self,
                  hmm: hmm_bw.HMM,
                  writer: TextIO,
                  strain: str,
                  chrm: str):
        '''
        Write information on the provided hmm as a line to the supplied textIO
        object.
        Output is tab delimited with:
        strain chromosome initial_probs emissions transitions
        '''
        writer.write(f'{strain}\t{chrm}\t')

        states = len(hmm.hidden_states)
        writer.write('\t'.join(
            [f'{p}' for p in hmm.initial_p] +  # initial
            [f'{hmm.emissions[i, hmm.symbol_to_ind[symbol]]}'
             if symbol in hmm.symbol_to_ind else '0.0'
             for i in range(states)
             for symbol in self.emission_symbols] +  # emission
            [f'{hmm.transitions[i, j]}'
             for i in range(states)
             for j in range(states)]  # transition
        ))
        writer.write('\n')

    def write_blocks_header(self, writer: TextIO) -> None:
        '''
        Write header line to tab delimited block file:
        strain chromosome predicted_species start end num_sites_hmm
        '''
        # NOTE: num_sites_hmm represents the sites considered by the HMM,
        # so it might exclude non-polymorphic sites in addition to gaps
        writer.write('\t'.join(['strain',
                                'chromosome',
                                'predicted_species',
                                'start',
                                'end',
                                'num_sites_hmm'])
                     + '\n')

    def write_blocks(self,
                     state_seq_blocks: List[Tuple[int, int]],
                     positions: np.array,
                     writer: TextIO,
                     strain: str,
                     chrm: str,
                     species_pred: str) -> None:
        '''
        Write entry into tab delimited block file, with columns:
        strain chromosome predicted_species start end num_sites_hmm
        '''
        writer.write('\n'.join(
            ['\t'.join([strain,
                        chrm,
                        species_pred,
                        str(positions[start]),
                        str(positions[end]),
                        str(end - start + 1)])
             for start, end in state_seq_blocks]))
        if state_seq_blocks:  # ensure ends with \n
            writer.write('\n')

    def write_positions(self,
                        positions: np.array,
                        writer: TextIO,
                        strain: str,
                        chrm: str) -> None:
        '''
        Write the positions of the specific strain, chromosome as a line to the
        provided textIO object
        '''
        writer.write(f'{strain}\t{chrm}\t' +
                     '\t'.join([str(x) for x in positions]) + '\n')

    def write_state_probs(self,
                          probs: Dict[str, List[float]],
                          writer: TextIO,
                          strain: str,
                          chrm: str) -> None:
        '''
        Write the probability of each state to the supplied textIO object
        Output is tab delimited with:
        strain chrom state1:prob1,prob2,...,probn state2...
        '''
        writer.write(f'{strain}\t{chrm}\t')

        writer.write('\t'.join(
            [f'{state}:' +
             ','.join([f'{site[i]:.5f}' for site in probs])
             for i, state in enumerate(self.config.states)]))

        writer.write('\n')

    def process_path(self, hmm: hmm_bw.HMM) -> Tuple[List[str], np.array]:
        '''
        Process the hmm path based the the predictor threshold value
        Return the predicted states and the probabilities of the master
        reference sequence
        '''
        probabilities = hmm.posterior_decoding()[0]

        # posterior
        if type(self.config.threshold) is float:
            path, path_probs = sim_process.get_max_path(probabilities,
                                                        hmm.hidden_states)
            path_t = sim_process.threshold_predicted(
                path,
                path_probs,
                self.config.threshold,
                self.config.known_states[0])

            return path_t, probabilities

        else:
            predicted = sim_predict.convert_predictions(hmm.viterbi(),
                                                        self.config.states)
            return predicted, probabilities

    def convert_to_blocks(self,
                          state_seq: List[str]) -> Dict[
                              str, List[Tuple[int, int]]]:
        '''
        Convert a list of sequences into a structure of start and end positions
        Return structure is a dict keyed on species with values of Lists of
        each block, which is a tuple with start and end positions
        '''
        # single individual state sequence
        blocks = {}
        for state in self.config.states:
            blocks[state] = []
        prev_species = state_seq[0]
        block_start = 0
        block_end = 0
        for i in range(len(state_seq)):
            if state_seq[i] == prev_species:
                block_end = i
            else:
                blocks[prev_species].append((block_start, block_end))
                block_start = i
                block_end = i
                prev_species = state_seq[i]
        # add last block
        if prev_species not in blocks:
            blocks[prev_species] = []
        blocks[prev_species].append((block_start, block_end))

        return blocks


class HMM_Builder():
    def __init__(self, configuration: Configuration):
        self.config = configuration
        self.config.set_HMM_symbols()
        self.symbols = self.config.symbols
        self.config.set_convergence()

    def update_emission_symbols(self, repeats: int):
        '''
        Generate all permutations of match and mismatch symbols with
        repeats number of characters, in lexigraphical order.
        Sets internal state and returns the emission symbols
        '''
        syms = [self.symbols['match'], self.symbols['mismatch']]
        emis_symbols = [''.join(x) for x in
                        itertools.product(syms,
                                          repeat=repeats)]
        emis_symbols.sort()
        self.emission_symbols = emis_symbols
        return emis_symbols

    def get_symbol_freqs(self, sequence: np.array) -> Tuple[Dict, List]:
        '''
        Calculate metrics from the provided, coded sequence
        Returns:
        the fraction of each matching pattern (e.g. +--++)
        the weighted fraction of matches for each species
        '''

        weighted = []

        symbols = defaultdict(int, Counter(sequence))
        total = len(sequence)
        for k in symbols:
            symbols[k] /= total

        sequence = np.array([list(s) for s in sequence])

        # look along species
        for s in np.transpose(sequence):
            s = ''.join(s)
            counts = Counter(s)
            weighted.append(counts[self.symbols['match']])

        total = sum(weighted)
        weighted = [w / total for w in weighted]
        return symbols, weighted

    def set_expected_values(self):
        '''
        Get expected lengths and fractions for each state.
        Assumes config has been validated by Predictor prior to running
        '''
        self.expected_lengths = {}
        self.expected_fractions = {}
        known_states = self.config.get('analysis_params.known_states')
        for state in known_states:
            self.expected_lengths[state['name']] = state['expected_length']
            self.expected_fractions[state['name']] = state['expected_fraction']

        unknown_states = self.config.get('analysis_params.unknown_states')
        for state in unknown_states:
            self.expected_lengths[state['name']] = state['expected_length']
            self.expected_fractions[state['name']] = state['expected_fraction']

        reference = self.config.get('analysis_params.reference')
        # expected fraction of reference is the remainder after other states
        # are specified
        self.expected_fractions[reference['name']] =\
            1 - sum(self.expected_fractions.values())

        self.ref_state = self.config.get('analysis_params.reference.name')
        self.known_states = self.config.known_states
        self.unknown_states = self.config.unknown_states

        # have to remove effect of unknown of these values for later
        self.ref_fraction = self.expected_fractions[self.ref_state] + \
            sum([self.expected_fractions[s] for s in self.unknown_states])
        # sum of fraction / length, or 1 / tract length
        self.other_sum = sum([self.expected_fractions[s['name']] /
                              self.expected_lengths[s['name']]
                              for s in known_states])

    def update_expected_length(self, total_length: int):
        '''
        Updates the expected length for the reference state
        based on the provided total_length of the sequence.
        This is the expected length of a single tract, determined as the sum
        of the total length (sequence length * fraction) divided by the number
        of tracts (sequence length * 1 / other's tracts). The + 1 assumes that
        the sequence will start and end with the reference.
        '''
        self.expected_lengths[self.ref_state] = (
            total_length * self.ref_fraction /
            (total_length * self.other_sum + 1))

    def initial_probabilities(self,
                              weighted_match_freqs: List[float]) -> np.array:
        '''
        Estimate the initial probability of being in each state
        based on the number of states and their expected fractions
        Returns the initial probability of each state
        '''

        init = []
        expectation_weight = .9
        for s, state in enumerate(self.known_states):
            expected = self.expected_fractions[state]
            estimated = weighted_match_freqs[s]
            init.append(expected * expectation_weight +
                        estimated * (1 - expectation_weight))

        for state in self.unknown_states:
            expected_frac = self.expected_fractions[state]
            init.append(expected_frac)

        return init / np.sum(init)

    def emission_probabilities(self,
                               symbols: List[str]) -> List[Dict]:
        '''
        Estimate initial emission probabilities
        Return estimates as list of default dict of probabilities
        '''

        match = self.symbols['match']
        mismatch = self.symbols['mismatch']
        probabilities = {
            mismatch + match: 0.9,
            match + match: 0.09,
            mismatch + mismatch: 0.009,
            match + mismatch: 0.001,
        }

        mismatch_bias = .99

        num_per_category = 2 ** (len(self.known_states) - 2)
        for key in probabilities:
            probabilities[key] *= num_per_category

        # for known states
        symbol_array = np.array([list(s) for s in symbols], dtype='<U1')
        # for unknown states
        symbol_length = symbol_array.shape[1]
        number_matches = (symbol_array == match).sum(axis=1)
        # combine first column with rest to generate probabilities
        first_column = np.tile(symbol_array[:, 0:1],
                               (1, len(self.known_states)))
        symbol_array = np.core.defchararray.add(
            first_column, symbol_array[:, 0:len(self.known_states)])
        # index into probabilities and normalize
        emissions = np.vectorize(probabilities.__getitem__)(symbol_array)
        emissions /= sum(emissions)

        # convert to match * (1-bias) + mismatch * bias, simplified
        number_matches = (number_matches + mismatch_bias *
                          (symbol_length - 2 * number_matches))
        number_matches /= sum(number_matches)
        # repeat for each unknown state
        number_matches = np.transpose(
            np.tile(number_matches, (len(self.unknown_states), 1)))

        # convert result into default dict
        result = [defaultdict(float,
                              {k: v for k, v in
                               zip(symbols, emissions[:, i])})
                  for i in range(emissions.shape[1])]
        result.extend([defaultdict(float,
                                   {k: v for k, v in
                                    zip(symbols, number_matches[:, i])})
                       for i in range(number_matches.shape[1])])

        return result

    def transition_probabilities(self) -> np.array:
        '''
        Estimate initial transition probabilities
        '''

        # doesn't depend on sequence observations but maybe it should?

        # also should we care about number of tracts rather than fraction
        # of genome? maybe theoretically, but that number is a lot more
        # suspect

        states = self.config.states

        fractions = np.array([self.expected_fractions[s] for s in states])
        lengths = 1/np.array([self.expected_lengths[s] for s in states])

        # general case,
        # trans[i,j] = 1/ length[i] * expected[j] * 1 /(1 - fraction[i])
        transitions = np.outer(
            np.multiply(lengths, 1/(1-fractions)),
            fractions)
        # when i == j, trans[i,j] = 1 - 1/length[i]
        np.fill_diagonal(transitions, 1-lengths)

        # normalize
        return transitions / transitions.sum(axis=1)[:, None]

    def build_initial_hmm(self, seq: np.array) -> hmm_bw.HMM:
        '''
        Build a HMM object initialized based on expected values and sequence
        '''

        # get frequencies of individual symbols (e.g. '+') and all full
        # combinations of symbols (e.g. '+++-')
        (symbol_freqs,
         weighted_match_freqs) = self.get_symbol_freqs(seq)

        # new Hidden Markov Model
        hmm = hmm_bw.HMM()

        hmm.set_initial_p(self.initial_probabilities(weighted_match_freqs))
        hmm.set_emissions(self.emission_probabilities(symbol_freqs.keys()))
        hmm.set_transitions(self.transition_probabilities())
        return hmm

    def run_hmm(self,
                alignment_file: str,
                only_poly_sites: bool = True) -> Tuple[hmm_bw.HMM,
                                                       hmm_bw.HMM,
                                                       np.array]:
        '''
        Runs the hmm training, returning the initial and trained HMM along
        with the positions of hmm importance
        '''
        coded_sequence, positions, len_seq = \
            self.encode_sequence(alignment_file, only_poly_sites)

        self.update_expected_length(len_seq)
        # set initial hmm parameters based on combination of (1) initial
        # expectations (length of introgressed tract and fraction of
        # genome/total number tracts and bases) and (2) number of sites at
        # which predict seq matches each reference
        hmm = self.build_initial_hmm(coded_sequence)

        # set states and initial probabilties
        hmm.set_hidden_states(self.known_states + self.unknown_states)

        # copy before setting observations to save memory
        hmm_init = copy.deepcopy(hmm)

        # set obs
        hmm.set_observations([coded_sequence])

        # Baum-Welch parameter estimation
        hmm.train(self.config.convergence)

        return hmm_init, hmm, positions

    def encode_sequence(self,
                        alignment_file: str,
                        only_poly_sites: bool = True) -> Tuple[
                            np.array,
                            np.array,
                            int]:
        '''
        open the supplied alignment file, encode, and return the coded
        sequence along with the positions.  If only_poly_sites is True,
        also filter out non-polymorphic sites.
        Returns the encoded sequence, positions, and length of original seq
        '''
        _, sequences = read_fasta(alignment_file)

        references = sequences[:-1]
        predicted = sequences[-1]

        seq_coded, positions = self.ungap_and_code(predicted, references)
        if only_poly_sites:
            seq_coded, positions = self.poly_sites(seq_coded, positions)

        return seq_coded, positions, len(predicted)

    def ungap_and_code(self,
                       predict_seq: str,
                       ref_seqs: List[str],
                       index_ref: int = 0) -> Tuple[np.array, np.array]:
        '''
        Remove any sequence locations where a gap is present and code
        into matching or mismatching sequence
        Returns the coded sequences, by default an array of + where matching, -
        where mismatching.  Also return the positions where the sequences are
        not gapped.
        '''
        # index_ref is index of reference strain to index relative to
        # build character array
        sequences = np.array([list(predict_seq)] +
                             [list(r) for r in ref_seqs])

        isbase = sequences != self.symbols['gap']

        # make boolean for valid characters
        isvalid = np.logical_and(sequences != self.symbols['gap'],
                                 sequences != self.symbols['unsequenced'])

        # positions are where everything is valid, index where the reference is
        # valid.  The +1 removes the predict sequence at index 0
        positions = np.where(
            np.all(isvalid[:, isbase[index_ref+1, :]], axis=0))[0]

        matches = np.where(sequences[0] == sequences[1:],
                           self.symbols['match'],
                           self.symbols['mismatch'])

        matches = np.fromiter((''.join(row)
                               for row in np.transpose(
                                   matches[:, np.all(isvalid, axis=0)])),
                              dtype=f'U{len(sequences) - 1}')

        return matches, positions

    def poly_sites(self,
                   sequences: np.array,
                   positions: np.array) -> Tuple[np.array, np.array]:
        '''
        Remove all sequences where the sequence is all match_symbol
        Returns the filtered sequence and position
        '''
        seq_len = len(sequences[0])
        # check if seq only contains match_symbol
        retain = np.vectorize(
            lambda x: x.count(self.symbols['match']) != seq_len)(sequences)
        indices = np.where(retain)[0]

        ps_poly = positions[indices]
        seq_poly = sequences[indices]

        return seq_poly, ps_poly


def read_positions(filename: str) -> Dict[str, Dict[str, List[int]]]:
    '''
    Read in positions from the provided filename, returning a dictionary
    keyed first by the strain, then chromosome.  Returned positions are
    lists of ints
    '''
    with gzip.open(filename, 'rt') as reader:
        result = defaultdict({})
        for line in reader:
            line = line.split()
            strain, chrm = line[0:2]
            positions = [int(x) for x in line[2:]]
            result[strain][chrm] = positions
    return result


def read_blocks(filename: str,
                labeled: bool = False) -> Dict[
                    str, Dict[str, Tuple[int, int, int, str]]]:
    '''
    Read in the supplied block file, returning a dict keyed on strain,
    then chromosome.  Values are tuples of start, end, and number of postions
    for the block.
    If labeled is true, values contain the region_id as last element
    '''
    with open(filename, 'r') as reader:
        reader.readline()  # header
        result = defaultdict(lambda: defaultdict(list))
        for line in reader:
            tokens = line.split()
            if labeled:
                (region_id, strain, chrm, species,
                 start, end, number_non_gap) = tokens
                item = (region_id, int(start), int(end), int(number_non_gap))
            else:
                (strain, chrm, species,
                 start, end, number_non_gap) = tokens
                item = (int(start), int(end), int(number_non_gap))
            result[strain][chrm].append(item)
    return result
