from contextlib import ExitStack
from operator import itemgetter
from analyze.introgression_configuration import Configuration
from analyze.predict import read_blocks
import click


class ID_producer():
    '''
    ID_producer
    Adds unique region id to block files
    '''
    def __init__(self, configuration: Configuration):
        self.config = configuration

    def add_ids(self):
        '''
        Adds a unique region id to block files, producing labeled text files
        '''
        self.config.validate_id_regions_arguments()
        regions = dict(zip(self.config.chromosomes,
                           [[] for _ in self.config.chromosomes]))
        with ExitStack() as stack:
            writers = {}

            # Progress bars don't seem to show since these complete too fast
            progress_bar = None
            if self.config.log_file:
                progress_bar = stack.enter_context(
                    click.progressbar(
                        length=len(self.config.states),
                        label='Reading in states'))

            for state in self.config.states:
                # read in region as dict keyed by strain, chromosome:
                # (start, end, number non gapped)
                region = read_blocks(self.config.blocks.format(state=state))
                for strain, d_strain in region.items():
                    for chrm, d_chrm in d_strain.items():
                        for start, end, num in d_chrm:
                            regions[chrm].append(
                                (start, end, num, strain, state))

                # open writer
                writers[state] = stack.enter_context(
                    open(self.config.labeled_blocks.format(state=state), 'w'))
                writers[state].write(
                    'region_id\tstrain\tchromosome\tpredicted_species\t'
                    'start\tend\tnum_sites_hmm\n')

                if progress_bar:
                    progress_bar.update(1)
            id_counter = 1

            if progress_bar:
                progress_bar = stack.enter_context(
                    click.progressbar(
                        length=len(regions.keys()),
                        label='Adding regions'))

            for chrm, entries in regions.items():
                # sort by start, then strain
                for start, end, num, strain, state in \
                        sorted(entries, key=itemgetter(0, 3)):
                    writers[state].write(
                        f'r{id_counter}\t{strain}\t{chrm}\t{state}\t'
                        f'{start}\t{end}\t{num}\n')
                    id_counter += 1
                if progress_bar:
                    progress_bar.update(1)
