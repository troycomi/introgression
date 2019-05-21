from misc import seq_functions
import numpy as np
from typing import List, Dict, TextIO, Tuple
import click
import logging as log
from contextlib import ExitStack, contextmanager
from misc import read_table
from misc.region_reader import Region_Reader
from analyze.introgression_configuration import Configuration


class Filterer():
    def __init__(self, configuration: Configuration):
        self.config = configuration

    def filter_introgressed(self,
                            region: Dict,
                            info: str,
                            reference_species: str) -> Tuple[bool, str]:
        '''
        filtering out things that we can't call introgressed in general
        with confidence (i.e. doesn't seem like a strong case against
        being S288c)
        Return true if the region passes the filter, or false with a string
        specifying which filter failed
        Tests:
        -fraction of gaps masked in reference > 0.5
        -fraction of gaps masked in predicted species > 0.5
        -number of matches to predicted > 7
        -number of matches to predicted > number matches to reference
        -divergence with predicted species
        '''

        predicted_species = region['predicted_species']

        aligned_length = (int(region['end']) - int(region['start']) + 1)

        # FILTER: fraction gaps + masked
        fraction_gaps_masked_threshold = .5
        # num_sites_nonmask_x is number of sites at which neither
        # reference x nor the test sequence is masked or has a gap or
        # unsequenced character
        fraction_gaps_masked_r = \
            1 - int(region['num_sites_nonmask_' +
                           reference_species])/aligned_length
        fraction_gaps_masked_s = \
            1 - int(region['num_sites_nonmask_' +
                           predicted_species])/aligned_length

        if fraction_gaps_masked_r > fraction_gaps_masked_threshold:
            return False, f'fraction gaps/masked in master = '\
                f'{fraction_gaps_masked_r}'
        if fraction_gaps_masked_s > fraction_gaps_masked_threshold:
            return False, f'fraction gaps/masked in predicted = '\
                f'{fraction_gaps_masked_s}'

        # FILTER: number sites analyzed by HMM that match predicted (P)
        # reference (C)
        count_P = np.sum(info == 'P')
        count_C = np.sum(info == 'C')
        number_match_only_threshold = 7
        if count_P < number_match_only_threshold:
            return False, f'count_P = {count_P}'
        if count_P <= count_C:
            return False, f'count_P = {count_P} and count_C = {count_C}'

        # FILTER: divergence with predicted reference and master reference
        # (S288c)
        id_predicted = float(region['match_nongap_' + predicted_species]) / \
            float(region['num_sites_nongap_' + predicted_species])
        id_master = float(region['match_nongap_' + reference_species]) / \
            float(region['num_sites_nongap_' + reference_species])

        if id_master >= id_predicted:
            return False, f'id with master = {id_master} '\
                f'and id with predicted = {id_predicted}'
        if id_master < .7:
            return False, f'id with master = {id_master}'

        return True, ''

    def filter_ambiguous(self,
                         region: Dict,
                         seqs: np.array,
                         threshold: float,
                         refs: List[str]) -> Tuple[bool,
                                                   List[str],
                                                   List[float],
                                                   List[int]]:
        '''
        filter out things we can't assign to one species specifically;
        return the other reasonable alternatives if we're filtering
        it out
        Returns:
        True if the region passes the filter
        Fails the filter if number of matches and fraction matching are >= more
        than one state for the region
        Region is updated with:
        A list of likely species for the region
        A list of fraction of matching sequence for each species
        A list of total matching sites
        '''

        s = region['predicted_species']

        ids = {}
        P_counts = {}

        seqs = np.asarray(seqs)
        # skip any gap or unsequenced in ref or test
        # also skip if ref and test equal (later test ri == test but not ref)
        symbols = self.config.symbols
        skip = np.any(
            (seqs[0] == symbols['gap'],
             seqs[0] == symbols['unsequenced'],
             seqs[-1] == symbols['gap'],
             seqs[-1] == symbols['unsequenced'],
             seqs[0] == seqs[-1]),
            axis=0)

        for ri, ref in enumerate(refs):
            if ri == 0:
                continue
            r_match, r_total = seq_functions.seq_id(seqs[-1], seqs[ri])
            if r_total != 0:
                ids[ref] = r_match / r_total
                P_counts[ref] = np.sum(
                    np.logical_and(
                        np.logical_not(skip),
                        seqs[ri] == seqs[-1]))

        alts = {}
        for r in ids.keys():
            # TODO should threshold be the same for both?
            if ids[r] >= threshold * ids[s] and \
               P_counts[r] >= threshold * P_counts[s]:
                alts[r] = (ids[r], P_counts[r])

        alt_states = sorted(alts.keys(),
                            key=lambda x: alts[x][0],
                            reverse=True)
        region['alternative_states'] = ','.join(alt_states)

        alt_ids = [alts[state][0] for state in alt_states]
        region['alternative_ids'] = ','.join(
            [str(x) for x in alt_ids])

        alt_P_counts = [alts[state][1] for state in alt_states]
        region['alternative_P_counts'] = ','.join(
            [str(x) for x in alt_P_counts])

        return len(alts) <= 1, alt_states

    def validate_arguments(self):
        args = [
            'introgressed',
            'introgressed_intermediate',
            'ambiguous',
            'ambiguous_intermediate',
            'filter_threshold',
            'known_states',
            'regions',
            'region_index',
            'symbols',
            'quality_blocks'
        ]
        variables = self.config.__dict__
        for arg in args:
            if arg not in variables or variables[arg] is None:
                err = ('Failed to validate Filterer, required argument '
                       f"'{arg}' was unset")
                log.exception(err)
                raise ValueError(err)

        if 'filter_sweep' not in variables or \
                variables['filter_sweep'] is None:
            log.warning(f"'filter_sweep' was unset and will not be run")

    def run(self, thresholds=[]):
        '''
        Filter region files based on thresold in config and sweep
        with the supplied threshold list
        '''
        self.validate_arguments()
        known_states = self.config.known_states
        log.debug(f'Known states: {known_states}')

        with Filter_Sweep(self.config.filter_sweep, thresholds) as sweeper,\
                ExitStack() as stack:

            progress_bar = None
            if self.config.log_file:
                progress_bar = stack.enter_context(
                    click.progressbar(
                        length=len(known_states[1:]),
                        label='Filtering'))

            sweeper.write_header()
            writers = Filter_Writers(self.config)

            for species_from in known_states[1:]:

                log.info(species_from)

                region_summary, fields = read_table.read_table_rows(
                    self.config.quality_blocks.format(state=species_from),
                    '\t')

                with writers.open_state(species_from, fields) as writers,\
                        Region_Reader(self.config.regions.format(
                            state=species_from), as_fa=True) as region_reader:

                    writers.write_headers()

                    for region_id, _, seqs in region_reader.yield_fa():
                        region = region_summary[region_id]
                        seqs, info_string = seqs[:-1], seqs[-1]

                        # filtering stage 1: things that we're confident in
                        # calling not S288c
                        passes, reason = self.filter_introgressed(
                            region,
                            info_string,
                            known_states[0])
                        region['reason'] = reason

                        writers.write_introgressed(region_id, region, passes)

                        if passes:
                            sweeper.record(
                                species_from,
                                lambda thresh: self.filter_ambiguous(
                                    region, seqs, thresh, known_states))

                            passes, _ = self.filter_ambiguous(
                                 region, seqs,
                                 self.config.filter_threshold, known_states)
                            writers.write_ambiguous(region_id, region, passes)

                if progress_bar:
                    progress_bar.update(1)

            sweeper.write_results(known_states[1:])


class Filter_Sweep():
    def __init__(self,
                 sweep_file: str,
                 thresholds: List[float]):
        self.sweep_file = sweep_file
        self.sweep_writer = None
        self.thresholds = thresholds
        self.data_table = {}

    def __enter__(self):
        if self.sweep_file is not None and self.thresholds != []:
            self.sweep_writer = open(self.sweep_file, 'w')

        return self

    def __exit__(self, type, value, traceback):
        if self.sweep_writer:
            self.sweep_writer.close()

        return traceback is None

    def write_header(self):
        '''
        Write the header for the sweep filter file
        '''
        if self.sweep_writer:
            self.sweep_writer.write(
                'threshold\tpredicted_state\talternative_states\tcount\n')

    def record(self, species_from, thresh_lambda):
        '''
        Record the thresholds for this filter sweep object.
        The thresh lambda is an anonymous function that takes a threshold
        and returns a tuple with the value at index 1 being the alternative
        states. Filter_ambiguous is what this is meant for.
        '''
        if self.sweep_writer is None:
            return

        for thresh in self.thresholds:
            _, states = thresh_lambda(thresh)
            self.record_data_hit(thresh, species_from, states)

    def record_data_hit(self, threshold: float, species: str, states: List):
        '''
        adds an entry to the data table or increments if exists
        '''
        key = ','.join(sorted(states))
        if threshold not in self.data_table:
            self.data_table[threshold] = {}

        if species not in self.data_table[threshold]:
            self.data_table[threshold][species] = {}

        if key not in self.data_table[threshold][species]:
            self.data_table[threshold][species][key] = 0

        self.data_table[threshold][species][key] += 1

    def write_results(self, states):
        if self.sweep_writer is None:
            return

        for thresh in self.thresholds:
            for species in states:
                if thresh in self.data_table and \
                        species in self.data_table[thresh]:
                    d = self.data_table[thresh][species]
                    for key, value in d.items():
                        self.sweep_writer.write(
                            f'{thresh}\t{species}\t{key}\t{value}\n')


class Filter_Writers():
    '''
    Writes the filter and intermediate files
    '''
    def __init__(self, config):
        self.files = {
            'introgressed': config.introgressed,
            'introgressed_int': config.introgressed_intermediate,
            'ambiguous': config.ambiguous,
            'ambiguous_int': config.ambiguous_intermediate
        }
        self.headers = None
        self.writers = None

    @contextmanager
    def open_state(self, state: str, fields: List):
        '''
        Open output files for the particular state
        '''
        self.headers = {
            'introgressed': fields,
            'introgressed_int': fields + ['reason'],
            'ambiguous': fields,
            'ambiguous_int': fields + ['alternative_states',
                                       'alternative_ids',
                                       'alternative_P_counts']
        }

        self.writers = {k: open(v.format(state=state), 'w')
                        for k, v in self.files.items()}

        yield self

        for writer in self.writers.values():
            writer.close()

        self.headers = None
        self.writers = None

    def write_headers(self):
        if self.headers is None or self.writers is None:
            return

        for key, writer in self.writers.items():
            writer.write('\t'.join(self.headers[key]) + '\n')

    def write_filtered_line(self,
                            writer: TextIO,
                            region_id: str,
                            region: Dict,
                            fields: List) -> None:
        '''
        Write the region id and values in "region" dict to open file writer
        '''
        writer.write(f'{region_id}\t')
        writer.write('\t'.join([str(region[field]) for field in fields[1:]]))
        writer.write('\n')

    def write_introgressed(self,
                           region_id: str,
                           region: Dict,
                           passes: bool):
        self.write_filtered_line(
            self.writers['introgressed_int'],
            region_id,
            region,
            self.headers['introgressed_int'])

        if passes:
            self.write_filtered_line(
                self.writers['introgressed'],
                region_id,
                region,
                self.headers['introgressed'])

    def write_ambiguous(self,
                        region_id: str,
                        region: Dict,
                        passes: bool):
        self.write_filtered_line(
            self.writers['ambiguous_int'],
            region_id,
            region,
            self.headers['ambiguous_int'])

        if passes:
            self.write_filtered_line(
                self.writers['ambiguous'],
                region_id,
                region,
                self.headers['ambiguous'])
