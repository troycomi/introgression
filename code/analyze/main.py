import click
import yaml
import logging as log
import analyze.predict
from analyze.introgression_configuration import Configuration
from analyze.id_regions import ID_producer


# TODO also check for snakemake object?
@click.group(invoke_without_command=True)
@click.option('--config', '-c',
              multiple=True,
              type=click.File('r'),
              help='Base configuration yaml.')
@click.option('-v', '--verbosity', count=True, default=3)
@click.option('--log-file',
              default='',
              help='Optional log file. If unset print to stdout.')
@click.pass_context
def cli(ctx, config, verbosity, log_file):
    '''
    Main entry script to run analyze methods
    '''

    verbosity -= 1
    verbosity = 4 if verbosity > 4 else verbosity
    levelstr, level = [
        ('CRITICAL', log.CRITICAL),
        ('ERROR', log.ERROR),
        ('WARNING', log.WARNING),
        ('INFO', log.INFO),
        ('DEBUG', log.DEBUG),
    ][verbosity]

    ctx.ensure_object(Configuration)

    confs = len(config)
    for path in config:
        conf = yaml.safe_load(path)
        ctx.obj.add_config(conf)

    ctx.obj.set_log_file(log_file)
    if ctx.obj.log_file is not None:
        log.basicConfig(level=level, filename=ctx.obj.log_file, filemode='w')
    else:
        log.basicConfig(level=level)
    log.info(f'Verbosity set to {levelstr}')

    log.info(f'Read in {confs} config file{"" if confs == 1 else "s"}')
    log.debug('Cleaned config:\n' + repr(ctx.obj))

    if ctx.invoked_subcommand is None:
        click.echo_via_pager(
            click.style(
                'No command supplied. Read in the following config:\n',
                fg='yellow') + repr(ctx.obj))


@cli.command()
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
@click.option('--threshold', default='',
              help='Threshold to apply to estimated path. Valid values are '
              'floats or `viterbi\'')
@click.option('--alignment', default='',
              help='Alignment file location with '
              '{prefix}, {strain}, and {chrom}')
@click.option('--only-poly-sites/--all-sites', default=True,
              help='Consider only polymorphic sites or all sites. '
              'Default is only polymorphic.')
@click.pass_context
def predict(ctx,
            blocks,
            prefix,
            test_strains,
            hmm_initial,
            hmm_trained,
            positions,
            probabilities,
            threshold,
            alignment,
            only_poly_sites):
    config = ctx.obj

    config.set_chromosomes()
    log.info(f'Found {len(config.chromosomes)} chromosomes in config')

    config.set_threshold(threshold)
    log.info(f'Threshold value is \'{config.threshold}\'')

    config.set_blocks_file(blocks)
    log.info(f'Output blocks file is \'{config.blocks}\'')

    config.set_states()
    config.set_prefix(prefix)
    log.info(f'Prefix is \'{config.prefix}\'')

    config.set_strains(test_strains)
    if config.test_strains is None:
        log.info(f'No test_strains provided')
    else:
        str_len = len(config.test_strains)
        log.info(f'Found {str_len} test strain'
                 f'{"" if str_len == 1 else "s"}')
    str_len = len(config.strains)
    log.info(f'Found {str_len} unique strain'
             f'{"" if str_len == 1 else "s"}')

    config.set_predict_files(hmm_initial,
                             hmm_trained,
                             positions,
                             probabilities,
                             alignment)
    log.info(f'Hmm_initial file is \'{config.hmm_initial}\'')
    log.info(f'Hmm_trained file is \'{config.hmm_trained}\'')
    log.info(f'Positions file is \'{config.positions}\'')
    log.info(f'Probabilities file is \'{config.probabilities}\'')
    log.info(f'Alignment file is \'{config.alignment}\'')

    predictor = analyze.predict.Predictor(config)
    if only_poly_sites:
        log.info('Only considering polymorphic sites')
    else:
        log.info('Considering all sites')
    predictor.run_prediction(only_poly_sites)


# accept multiple states and pass as list
@cli.command()
@click.option('--blocks', default='', help='Block file location with {state}')
@click.option('--labeled', default='', help='Block file location with {state}')
@click.option('--state', multiple=True, help='States to add ids to')
@click.pass_context
def id_regions(ctx, blocks, labeled, state):
    config = ctx.obj
    config.set_chromosomes()
    log.info(f'Found {len(config.chromosomes)} chromosomes in config')

    state = list(state)
    config.set_states(state)
    log.info(f'Found {len(config.states)} states to process')

    config.set_blocks_file(blocks)
    log.info(f'Input blocks file is \'{config.blocks}\'')

    config.set_labeled_blocks_file(labeled)
    log.info(f'Output blocks file is \'{config.labeled_blocks}\'')

    id_producer = ID_producer(config)
    id_producer.add_ids()


if __name__ == '__main__':
    cli()
