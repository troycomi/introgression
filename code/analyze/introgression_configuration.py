import glob
import re
from typing import Tuple, Dict, List
import logging as log
from misc.config_utils import (get_nested, clean_config, merge_dicts,
                               print_dict, validate, check_wildcards)


class Configuration():
    def __init__(self):
        self.config = {}
        self.log_file = None

    def add_config(self, configuration: Dict):
        '''
        merge the provided configuration dictionary with this object.
        Cleans configuration
        '''
        self.config = clean_config(
            merge_dicts(self.config, configuration))

    def get_states(self) -> Tuple[List, List]:
        '''
        Build lists of known and unknown states from the analysis params
        '''

        ref = get_nested(self.config, 'analysis_params.reference.name')
        if ref is None:
            ref = []
        else:
            ref = [ref]

        known = get_nested(self.config, 'analysis_params.known_states')
        if known is None:
            known = []

        known_states = ref + [s['name'] for s in known]

        unknown = get_nested(self.config, 'analysis_params.unknown_states')
        if unknown is None:
            unknown = []

        unknown_states = [s['name'] for s in unknown]

        return known_states, unknown_states

    def get_interval_states(self) -> List:
        '''
        Build list of interval states, typically just known names
        but if the state has an interval name, use that
        '''
        ref = get_nested(self.config, 'analysis_params.reference')

        # set with name or empty list
        if ref is None:
            ref = []
        else:
            ref = [ref]

        known = get_nested(self.config, 'analysis_params.known_states')
        if known is None:
            known = []

        return [s['interval_name']
                if 'interval_name' in s
                else s['name']
                for s in ref + known]

    def set_states(self, states: List[str] = None):
        '''
        Set the states for which to perform region naming
        '''
        if states is None or states == []:
            self.known_states, self.unknown_states = self.get_states()
            self.states = self.known_states + self.unknown_states
        else:
            self.states = states

        self.interval_states = self.get_interval_states()

        if self.states == []:
            err = 'No states specified'
            log.exception(err)
            raise ValueError(err)

    def set_log_file(self, log_file: str = ''):
        '''
        sets log file based on provided value or config
        '''
        if log_file == '':
            self.log_file = get_nested(self.config, 'paths.log_file')
        else:
            self.log_file = log_file

    def set_chromosomes(self):
        '''
        Gets the chromosome list from config, raising a ValueError
        if undefined.
        '''
        self.chromosomes = validate(
            self.config,
            'chromosomes',
            'No chromosomes specified in config file!')

    def set_threshold(self, threshold: str = None):
        '''
        Set the threshold. Checks if set and converts to float if possible.
        Failing float casting, will store a string if it is 'viterbi',
        otherwise throws a ValueError
        '''
        self.threshold = validate(
            self.config,
            'analysis_params.threshold',
            'No threshold provided',
            threshold)
        try:
            self.threshold = float(self.threshold)
        except ValueError:
            if self.threshold != 'viterbi':
                err = f'Unsupported threshold value: {self.threshold}'
                log.exception(err)
                raise ValueError(err)

    def set_blocks_file(self, blocks: str = None):
        '''
        Set the block wildcard filename.  Checks for appropriate wildcards
        '''
        self.blocks = validate(
            self.config,
            'paths.analysis.blocks',
            'No block file provided',
            blocks)

        check_wildcards(self.blocks, 'state')

    def set_labeled_blocks_file(self, blocks: str = None):
        '''
        Set the labeled block wildcard filename.
        Checks for appropriate wildcards
        '''
        self.labeled_blocks = validate(
            self.config,
            'paths.analysis.labeled_blocks',
            'No labeled block file provided',
            blocks)

        check_wildcards(self.labeled_blocks, 'state')

    def set_quality_file(self, quality: str = None):
        '''
        Set the quality block wildcard filename.
        Checks for appropriate wildcards
        '''
        self.quality_blocks = validate(
            self.config,
            'paths.analysis.quality',
            'No quality block file provided',
            quality)

        check_wildcards(self.quality_blocks, 'state')

    def set_masked_file(self, masks: str = None):
        '''
        Set the masked interval block wildcard filename.
        Checks for appropriate wildcards
        '''
        self.masks = validate(
            self.config,
            'paths.analysis.masked_intervals',
            'No masked interval file provided',
            masks)

        check_wildcards(self.masks, 'strain,chrom')

    def set_prefix(self, prefix: str = ''):
        '''
        Set prefix string of the predictor to the supplied value or
        build it from the known states
        '''
        if prefix == '':
            if self.known_states == []:
                err = 'Unable to build prefix, no known states provided'
                log.exception(err)
                raise ValueError(err)

            self.prefix = '_'.join(self.known_states)
        else:
            self.prefix = prefix

    def set_strains(self, test_strains: str = ''):
        '''
        build the strains to perform prediction on
        '''
        if test_strains == '':
            test_strains = get_nested(self.config, 'paths.test_strains')
        else:
            # need to support list for test strains
            test_strains = [test_strains]

        if test_strains is not None:
            for test_strain in test_strains:
                check_wildcards(test_strain, 'strain,chrom')

        self.find_strains(test_strains)

    def find_strains(self, test_strains: List[str] = None):
        '''
        Helper method to get strains supplied in config, or from test_strains
        '''
        strains = get_nested(self.config, 'strains')
        self.test_strains = test_strains

        if strains is None:
            if test_strains is None:
                err = ('Unable to find strains in config and '
                       'no test_strains provided')
                log.exception(err)
                raise ValueError(err)

            # try to build strains from wildcards in test_strains
            strains = {}
            for test_strain in test_strains:
                # find matching files
                strain_glob = test_strain.format(
                    strain='*',
                    chrom='*')
                log.info(f'searching for {strain_glob}')
                for fname in glob.iglob(strain_glob):
                    # extract wildcard matches
                    match = re.match(
                        test_strain.format(
                            strain='(?P<strain>.*?)',
                            chrom='(?P<chrom>[^_]*?)'
                        ),
                        fname)
                    if match:
                        log.debug(
                            f'matched with {match.group("strain", "chrom")}')
                        strain, chrom = match.group('strain', 'chrom')
                        if strain not in strains:
                            strains[strain] = set()
                        strains[strain].add(chrom)

            if len(strains) == 0:
                err = ('Found no chromosome sequence files '
                       f'in {test_strains}')
                log.exception(err)
                raise ValueError(err)

            # check if requested chromosomes are within the list of chroms
            chrom_set = set(self.chromosomes)
            for strain, chroms in strains.items():
                if not chrom_set.issubset(chroms):
                    not_found = chrom_set.difference(chroms).pop()
                    err = (f'Strain {strain} is missing chromosomes. '
                           f'Unable to find chromosome \'{not_found}\'')
                    log.exception(err)
                    raise ValueError(err)

            self.strains = list(sorted(strains.keys()))

        else:  # strains set in config
            self.strains = list(sorted(set(strains)))

    def set_predict_files(self,
                          hmm_initial: str,
                          hmm_trained: str,
                          positions: str,
                          probabilities: str,
                          alignment: str):
        '''
        Set output files from provided values or config.
        Raises value errors if a file is not provided.
        Checks alignment for all wildcards and replaces prefix.
        '''
        self.hmm_initial = validate(self.config,
                                    'paths.analysis.hmm_initial',
                                    'No initial hmm file provided',
                                    hmm_initial)

        self.hmm_trained = validate(self.config,
                                    'paths.analysis.hmm_trained',
                                    'No trained hmm file provided',
                                    hmm_trained)

        self.set_positions(positions)

        self.probabilities = validate(self.config,
                                      'paths.analysis.probabilities',
                                      'No probabilities file provided',
                                      probabilities)

        self.set_alignment(alignment)

    def set_alignment(self, alignment: str):
        '''
        Set the alignment file, checking wildcards prefix, strain and chrom.
        If prefix is present, it is substituted, otherwise checks just
        strain and chrom
        '''
        alignment = validate(self.config,
                             'paths.analysis.alignment',
                             'No alignment file provided',
                             alignment)
        if '{prefix}' in alignment:
            check_wildcards(alignment, 'prefix,strain,chrom')
            self.alignment = alignment.replace('{prefix}', self.prefix)
        else:
            check_wildcards(alignment, 'strain,chrom')
            self.alignment = alignment

    def set_positions(self, positions: str):
        '''
        Sets the position file
        '''
        self.positions = validate(self.config,
                                  'paths.analysis.positions',
                                  'No positions file provided',
                                  positions)

    def set_regions_files(self,
                          regions: str = None,
                          region_index: str = None):
        '''
        Set the region and pickle wildcard filename. Checks for state wildcards
        '''
        self.regions = validate(
            self.config,
            'paths.analysis.regions',
            'No region file provided',
            regions)
        check_wildcards(self.regions, 'state')

        self.region_index = validate(
            self.config,
            'paths.analysis.region_index',
            'No region index file provided',
            region_index)
        check_wildcards(self.region_index, 'state')

    def set_HMM_symbols(self):
        '''
        Set symbols based on config values, using defaults if unset
        '''
        self.symbols = {
            'match': '+',
            'mismatch': '-',
            'unknown': '?',
            'unsequenced': 'n',
            'gap': '-',
            'unaligned': '?',
            'masked': 'x'
        }
        config_symbols = get_nested(self.config, 'HMM_symbols')
        if config_symbols is not None:
            for k, v in config_symbols.items():
                if k not in self.symbols:
                    log.warning("Unused symbol in configuration: "
                                f"{k} -> '{v}'")
                else:
                    self.symbols[k] = v
                    log.debug(f"Overwriting default symbol for {k} with '{v}'")

            for k, v in self.symbols.items():
                if k not in config_symbols:
                    log.warning(f'Symbol for {k} unset in config, '
                                f"using default '{v}'")

        else:
            for k, v in self.symbols.items():
                log.warning(f'Symbol for {k} unset in config, '
                            f"using default '{v}'")

    def set_convergence(self):
        '''
        Set convergence for HMM training, using default if unset
        '''
        self.convergence = get_nested(self.config,
                                      'analysis_params.convergence_threshold')
        if self.convergence is None:
            log.warning('No value set for convergence_threshold, using '
                        'default of 0.001')
            self.convergence = 0.001

    def get(self, key: str):
        '''
        Get nested key from underlying dictionary. Returning none if any
        key is not in dict
        '''
        return get_nested(self.config, key)

    def __repr__(self):
        return ('Config file:\n' +
                print_dict(self.config) +
                '\nSettings:\n' +
                print_dict({k: v for k, v in self.__dict__.items()
                            if k != 'config'})
                )
