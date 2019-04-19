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

    predictor = predict.Predictor(config)
    predictor.set_chromosomes()

    predictor.set_blocks_file(blocks)
    log.info(f'output blocks file for predict is {predictor.blocks}')

    predictor.set_prefix(prefix)
    log.info(f'prefix is {predictor.prefix}')

    predictor.set_strains(test_strains)
    log.info(f'found {len(predictor.test_strains)} test strains')
    log.info(f'found {len(predictor.strains)} unique strains')

    predictor.set_output_files(hmm_initial,
                               hmm_trained,
                               positions,
                               probabilities,
                               alignment)
    log.info(f'hmm_initial is {predictor.hmm_initial}')
    log.info(f'hmm_trained is {predictor.hmm_trained}')
    log.info(f'positions is {positions}')
    log.info(f'probabilities is {predictor.probabilities}')
    log.info(f'alignment is {predictor.alignment}')

    predictor.validate_arguments()
    predictor.run_prediction()
