from analyze.introgression_configuration import Configuration
import logging as log
import itertools
from misc import read_table
from typing import List
from contextlib import ExitStack
import click


class Strain_Summarizer():
    def __init__(self, configuration: Configuration):
        self.config = configuration

    def validate_arguments(self):
        '''
        Check that all required instance variables are set to perform a
        strain summary run. Returns true if valid, raises value error otherwise
        '''
        args = [
            'known_states',
            'introgressed_intermediate',
            'ambiguous_intermediate',
            'strain_info',
            'state_counts',
        ]
        variables = self.config.__dict__
        for arg in args:
            if arg not in variables or variables[arg] is None:
                err = ('Failed to validate strain summarizer,'
                       f" required argument '{arg}' was unset")
                log.exception(err)
                raise ValueError(err)

        return True

    def run(self):
        '''
        Generate summary information for the state of
        each position in the sequence
        '''
        self.validate_arguments()

        summary = Summary_Table()

        states = self.config.known_states[1:]
        with ExitStack() as stack:
            progress_bar = None
            if self.config.log_file:
                progress_bar = stack.enter_context(
                    click.progressbar(
                        length=len(states),
                        label='State'))
            for species_from in states:

                log.info(species_from)

                regions1, _ = read_table.read_table_rows(
                    self.config.introgressed_intermediate.format(
                        state=species_from), '\t')
                regions2, _ = read_table.read_table_rows(
                    self.config.ambiguous_intermediate.format(
                        state=species_from), '\t')

                for region_id in regions1:
                    region1 = regions1[region_id]

                    strain = region1['strain']
                    length = int(region1['end']) - int(region1['start']) + 1

                    summary.set_region(strain, species_from, length)
                    summary.region_found()

                    if region1['reason'] != '':  # failed filter
                        continue

                    summary.region_passes_filter1()

                    region2 = regions2[region_id]
                    summary.record_alt_species(
                        region2['alternative_states'].split(','))

                if progress_bar:
                    progress_bar.update(1)

            with open(self.config.strain_info, 'r') as reader:
                summary.add_strain_info(reader)

            with open(self.config.state_counts, 'w') as writer:
                summary.write_summary(states, writer)


class Summary_Table():
    def __init__(self):
        self.table = {}

    def set_region(self, strain, species, length):
        self.strain = strain
        self.species = species
        self.length = length

    def record_element(self,
                       strain: str,
                       key: str,
                       count: int = 1):
        '''
        Increment the count of table[strain][key], adding new values as needed
        '''

        if strain not in self.table:
            self.table[strain] = {}

        t = self.table[strain]
        if key not in t:
            t[key] = 0

        t[key] += count

    def record_region(self,
                      strain: str,
                      species: str,
                      length: int,
                      suffix: str = "",
                      update_total: bool = True):
        '''
        Record a region of provided length.
        '''
        if suffix and suffix[0] != '_':
            suffix = '_' + suffix

        self.record_element(strain, f'num_regions_{species}{suffix}', 1)
        self.record_element(strain, f'num_bases_{species}{suffix}', length)
        if update_total:
            self.record_element(strain, f'num_bases_total{suffix}', length)
            self.record_element(strain, f'num_regions_total{suffix}', 1)

    def record_alt_species(self, alt_states: List):
        for species in alt_states:
            self.record_alt(species)

        if len(alt_states) == 1:
            self.record_region(self.strain, self.species,
                               self.length, '_filtered2')
        else:
            self.record_element(self.strain,
                                ('num_bases_' +
                                 '_or_'.join(sorted(alt_states)) +
                                 '_filtered2i'),
                                self.length)

        self.record_element(self.strain,
                            f'num_bases_{len(alt_states)}_filtered2i',
                            self.length)

    def region_found(self):
        self.record_region(self.strain, self.species, self.length)

    def region_passes_filter1(self):
        self.record_region(self.strain, self.species,
                           self.length, '_filtered1')

    def record_alt(self, alt_species):
        self.record_region(self.strain, alt_species,
                           self.length, '_filtered2_inclusive',
                           self.species == alt_species)

    def add_strain_info(self, reader):
        for line in reader:
            strain, _, _, geo, env, pop = line[:-1].split('\t')
            strain = strain.lower()
            if strain in self.table:
                d = self.table[strain]
                d['population'] = pop
                d['geographic_origin'] = geo
                d['environmental_origin'] = env

    def write_summary(self, states, writer):
        fields = self.get_fields(states)

        # write header
        writer.write('strain\t' + '\t'.join(fields) + '\n')

        for strain in sorted(self.table.keys()):
            row = self.table[strain]
            entries = [row[field]
                       if field in row
                       else 0
                       for field in fields]

            entries = [str(s) for s in [strain] + entries]

            writer.write('\t'.join(entries) + '\n')

    def get_fields(self, states):
        fields = ['population', 'geographic_origin', 'environmental_origin'] +\
            [f'num_{thing}_{state}{value}'
             for thing in ('regions', 'bases')
             for value in ('', '_filtered1',
                           '_filtered2', '_filtered2_inclusive')
             for state in states + ['total']
             ]

        r = sorted(states)
        for n in range(2, len(r)+1):
            fields += [f'num_bases_{"_or_".join(combo)}_filtered2i'
                       for combo in itertools.combinations(r, n)]
            fields += [f'num_bases_{n}_filtered2i']

        return fields
