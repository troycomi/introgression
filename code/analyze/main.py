import click
import yaml
import glob
import re
import logging as log
from misc import config_utils
from misc.config_utils import (get_nested, check_wildcards, get_states,
                               validate)
from typing import List, Dict


# TODO also check for snakemake object?
@click.group(invoke_without_command=True)
@click.option('--config', '-c',
              multiple=True,
              type=click.File('r'),
              help='Base configuration yaml.')
@click.option('-v', '--verbosity', count=True, default=2)
@click.pass_context
def cli(ctx, config, verbosity):
    '''
    Main entry script to run analyze methods
    '''

    verbosity = 4 if verbosity > 4 else verbosity
    levelstr = ['CRITICAL', 'ERROR',
                'WARNING', 'INFO',
                'DEBUG'][verbosity]
    level = [log.CRITICAL, log.ERROR,
             log.WARNING, log.INFO,
             log.DEBUG][verbosity]

    log.basicConfig(level=level)
    log.info(f'Verbosity set to {levelstr}')

    ctx.ensure_object(dict)

    log.info(f'Reading in {len(config)} config files')
    for path in config:
        conf = yaml.safe_load(path)
        ctx.obj = config_utils.merge_dicts(ctx.obj, conf)

    ctx.obj = config_utils.clean_config(ctx.obj)
    log.debug('Cleaned config:\n' + config_utils.print_dict(ctx.obj))

    if ctx.invoked_subcommand is None:
        click.echo_via_pager(
            click.style(
                'No command supplied. Read in the following config:\n',
                fg='yellow') +
            config_utils.print_dict(ctx.obj))


@cli.command()
@click.pass_context
@click.option('--blocks', default='', help='Block file location with {state}')
@click.option('--prefix', default='', help='Prefix of test-strain files '
              'default to list of states joined with _.')
@click.option('--test-strains', default='',
              help='Test files location with {strain} and {chrom}')
@click.option('--hmm-initial', default='',
              help='Initial hmm parameter text file')
@click.option('--hmm-trained', default='',
              help='Trained hmm parameter text file')
@click.option('--positions', default='',
              help='Positions file, gzipped')
@click.option('--probabilities', default='',
              help='Probabilities file, gzipped')
@click.option('--alignment', default='',
              help='Alignment file location with '
              '{prefix}, {strain}, and {chrom}')
def predict(ctx,
            blocks,
            prefix,
            test_strains,
            hmm_initial,
            hmm_trained,
            positions,
            probabilities,
            alignment):
    config = ctx.obj

    chromosomes = validate(config,
                           'chromosomes',
                           'No chromosomes specified in config file!')

    blocks = validate(config,
                      'paths.analysis.block_files',
                      'No block file provided',
                      blocks)

    check_wildcards(blocks, 'state')
    log.info(f'output blocks file for predict is {blocks}')

    known, unknown = get_states(config)
    if prefix == '':
        prefix = '_'.join(known)

    log.info(f'prefix is {prefix}')

    if test_strains == '':
        test_strains = get_nested(config, 'paths.test_strains')
    else:
        # need to support list for test strains
        test_strains = [test_strains]
    for test_strain in test_strains:
        check_wildcards(test_strain, 'strain,chrom')

    log.info(f'found {len(test_strains)} test strains')

    strains = get_strains(config, test_strains, prefix, chromosomes)
    log.info(f'found {len(strains)} unique strains')

    hmm_initial = validate(config,
                           'paths.analysis.hmm_initial',
                           'No initial hmm file provided',
                           hmm_initial)
    log.info(f'hmm_initial is {hmm_initial}')

    hmm_trained = validate(config,
                           'paths.analysis.hmm_trained',
                           'No trained hmm file provided',
                           hmm_trained)
    log.info(f'hmm_trained is {hmm_trained}')

    positions = validate(config,
                         'paths.analysis.positions',
                         'No positions file provided',
                         positions)
    log.info(f'positions is {positions}')

    probabilities = validate(config,
                             'paths.analysis.probabilities',
                             'No probabilities file provided',
                             probabilities)
    log.info(f'probabilities is {probabilities}')

    alignment = validate(config,
                         'paths.analysis.alignment',
                         'No alignment file provided',
                         alignment)
    check_wildcards(alignment, 'prefix,strain,chrom')
    alignment = alignment.replace('{prefix}', prefix)
    log.info(f'alignment is {alignment}')


def get_strains(config: Dict,
                test_strains: List,
                prefix: str,
                chromosomes: List):
    '''
    Helper method to get strains supplied in config, or from test_strains
    '''
    strains = get_nested(config, 'strains')

    if strains is None:
        # try to build strains from wildcards in test_strains
        strains = {}
        for test_strain in test_strains:
            strain_glob = test_strain.format(
                prefix=prefix,
                strain='*',
                chrom='*')
            log.info(f'searching for {strain_glob}')
            for fname in glob.iglob(strain_glob):
                match = re.match(
                    test_strain.format(
                        prefix=prefix,
                        strain='(?P<strain>.*?)',
                        chrom='(?P<chrom>[^_]*?)'
                    ),
                    fname)
                if match:
                    log.debug(f'matched with {match.group("strain", "chrom")}')
                    strain, chrom = match.group('strain', 'chrom')
                    if strain not in strains:
                        strains[strain] = []
                    strains[strain].append(chrom)

        if len(strains) == 0:
            err = f'Found no chromosome sequence files in {test_strains}'
            log.exception(err)
            raise ValueError(err)

        for strain, chroms in strains.items():
            if len(chromosomes) != len(chroms):
                err = (f'Strain {strain} has incorrect number of chromosomes. '
                       f'Expected {len(chromosomes)} found {len(chroms)}')
                log.exception(err)
                raise ValueError(err)
    return list(sorted(strains.keys()))
