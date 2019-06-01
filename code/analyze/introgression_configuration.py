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

        # these are very regular variables with state as a wildcard
        state_files = [
            'blocks',
            'labeled_blocks',
            'quality_blocks',
            'introgressed',
            'introgressed_intermediate',
            'ambiguous',
            'ambiguous_intermediate',
            'regions',
            'region_index',
        ]
        # no wildcards, non nullable
        nonwild_files = [
            'hmm_initial',
            'hmm_trained',
            'positions',
            'probabilities',
            'strain_info',
            'state_counts',
        ]
        var_list = [
            Variable('chromosomes'),
            Threshold_Variable(),
            Convergence_Variable(),
            Symbols_Variable(),
            Filter_Threshold_Variable(),
            Variable('log_file', 'paths.log_file', nullable=True),
            Variable('filter_sweep', 'paths.analysis.filter_sweep',
                     nullable=True),
            Variable('masks', 'paths.analysis.masked_intervals',
                     wildcards='strain,chrom'),
        ] + [
            Variable(n, f'paths.analysis.{n}', wildcards='state')
            for n in state_files
        ] + [
            Variable(n, f'paths.analysis.{n}')
            for n in nonwild_files
        ]

        self.variables = {v.name: v for v in var_list}
        # these require too much state from configuration to split out
        self.other_parsers = {
            'states': self._set_states,
            'prefix': self._set_prefix,
            'strains': self._set_strains,
            'alignment': self._set_alignment
        }

    def add_config(self, configuration: Dict):
        '''
        merge the provided configuration dictionary with this object.
        Cleans configuration
        '''
        self.config = clean_config(
            merge_dicts(self.config, configuration))

    def set(self, *args, **kwargs):
        '''
        Set the supplied variable to the value provided.
        If just a name is provided, set the value with a value of None
        '''
        kwargs.update({a: None for a in args})
        for key, value in kwargs.items():
            if key in self.variables:
                variable = self.variables[key]
                self.__dict__[key] = variable.parse(value, self.config)

            elif key in self.other_parsers:
                self.other_parsers[key](value)

            else:
                err = f'Unknown variable to set: {key}'
                log.exception(err)
                raise ValueError(err)

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

    def _set_states(self, states: List[str] = None):
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

    def _set_prefix(self, prefix: str = ''):
        '''
        Set prefix string of the predictor to the supplied value or
        build it from the known states
        '''
        if not prefix:
            if self.known_states == []:
                err = 'Unable to build prefix, no known states provided'
                log.exception(err)
                raise ValueError(err)

            self.prefix = '_'.join(self.known_states)
        else:
            self.prefix = prefix

    def _set_strains(self, test_strains: str = ''):
        '''
        build the strains to perform prediction on
        '''
        if not test_strains:
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

    def _set_alignment(self, alignment: str):
        '''
        Set the alignment file, checking wildcards prefix, strain and chrom.
        If prefix is present, it is substituted, otherwise checks just
        strain and chrom
        '''
        alignment = validate(self.config,
                             'paths.analysis.alignment',
                             'No alignment provided',
                             alignment)

        check_wildcards(alignment, 'strain,chrom')
        if '{prefix}' in alignment:
            self.alignment = alignment.replace('{prefix}', self.prefix)
        else:
            self.alignment = alignment

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
                            if k != 'config' and k != 'variables'
                            and k != 'other_parsers'})
                )


class Variable():
    def __init__(self, name, config_path=None, nullable=False, wildcards=None):
        self.name = name
        if config_path:
            self.config_path = config_path
        else:
            self.config_path = name

        self.nullable = nullable
        self.wildcards = wildcards

    def parse(self, value, config={}):
        if self.nullable:
            if not value:
                value = get_nested(config, self.config_path)

        else:
            value = validate(config, self.config_path,
                             f'No {self.name} provided', value)

        if self.wildcards:
            check_wildcards(value, self.wildcards)

        return value


class Threshold_Variable(Variable):
    def __init__(self):
        super().__init__('threshold', 'analysis_params.threshold')

    def parse(self, value, config={}):
        value = super().parse(value, config)

        try:
            value = float(value)

        except ValueError:
            if value != 'viterbi':
                err = f'Unsupported threshold value: {value}'
                log.exception(err)
                raise ValueError(err)

        return value


class Filter_Threshold_Variable(Variable):
    def __init__(self):
        super().__init__('filter_threshold',
                         'analysis_params.filter_threshold')

    def parse(self, value, config={}):
        value = super().parse(value, config)

        try:
            value = float(value)

        except (ValueError, TypeError):
            err = 'Filter threshold is not a valid number'
            log.exception(err)
            raise ValueError(err)

        return value


class Convergence_Variable(Variable):
    def __init__(self):
        super().__init__('convergence',
                         'analysis_params.convergence_threshold',
                         nullable=True)

    def parse(self, value, config={}):
        value = super().parse(value, config)

        try:
            value = float(value)

        except (ValueError, TypeError):
            log.warning('No value set for convergence_threshold, using '
                        'default of 0.001')
            value = 0.001

        return value


class Symbols_Variable(Variable):
    def __init__(self):
        super().__init__('symbols', '')

    def parse(self, value, config):
        '''
        Set symbols based on config values, using defaults if unset
        '''
        symbols = {
            'match': '+',
            'mismatch': '-',
            'unknown': '?',
            'unsequenced': 'n',
            'gap': '-',
            'unaligned': '?',
            'masked': 'x'
        }
        config_symbols = get_nested(config, 'HMM_symbols')
        if config_symbols is not None:
            for k, v in config_symbols.items():
                if k not in symbols:
                    log.warning("Unused symbol in configuration: "
                                f"{k} -> '{v}'")
                else:
                    symbols[k] = v
                    log.debug(f"Overwriting default symbol for {k} with '{v}'")

            for k, v in symbols.items():
                if k not in config_symbols:
                    log.warning(f'Symbol for {k} unset in config, '
                                f"using default '{v}'")

        else:
            for k, v in symbols.items():
                log.warning(f'Symbol for {k} unset in config, '
                            f"using default '{v}'")

        return symbols
