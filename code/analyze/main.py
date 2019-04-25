import click
import yaml
import logging as log
from misc import config_utils
import analyze.predict


# TODO also check for snakemake object?
@click.group(invoke_without_command=True)
@click.option('--config', '-c',
              multiple=True,
              type=click.File('r'),
              help='Base configuration yaml.')
@click.option('-v', '--verbosity', count=True, default=3)
@click.pass_context
def cli(ctx, config, verbosity):
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

    log.basicConfig(level=level)
    log.info(f'Verbosity set to {levelstr}')

    ctx.ensure_object(dict)

    confs = len(config)
    log.info(f'Reading in {confs} config file{"" if confs == 1 else "s"}')
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
@click.option('--threshold', default='',
              help='Threshold to apply to estimated path. Valid values are '
              'floats or `viterbi\'')
@click.option('--alignment', default='',
              help='Alignment file location with '
              '{prefix}, {strain}, and {chrom}')
@click.option('--only-poly-sites/--all-sites', default=True,
              help='Consider only polymorphic sites or all sites')
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

    predictor = analyze.predict.Predictor(config)
    predictor.set_chromosomes()
    log.info(f'Found {len(predictor.chromosomes)} chromosomes in config')

    predictor.set_threshold(threshold)
    log.info(f'Threshold value is \'{predictor.threshold}\'')

    predictor.set_blocks_file(blocks)
    log.info(f'Output blocks file is \'{predictor.blocks}\'')

    predictor.set_prefix(prefix)
    log.info(f'Prefix is \'{predictor.prefix}\'')

    predictor.set_strains(test_strains)
    if predictor.test_strains is None:
        log.info(f'No test_strains provided')
    else:
        str_len = len(predictor.test_strains)
        log.info(f'Found {str_len} test strain'
                 f'{"" if str_len == 1 else "s"}')
    log.info(f'Found {len(predictor.strains)} unique strains')

    predictor.set_output_files(hmm_initial,
                               hmm_trained,
                               positions,
                               probabilities,
                               alignment)
    log.info(f'Hmm_initial file is \'{predictor.hmm_initial}\'')
    log.info(f'Hmm_trained file is \'{predictor.hmm_trained}\'')
    log.info(f'Positions file is \'{predictor.positions}\'')
    log.info(f'Probabilities file is \'{predictor.probabilities}\'')
    log.info(f'Alignment file is \'{predictor.alignment}\'')

    predictor.run_prediction(only_poly_sites)


if __name__ == '__main__':
    cli()
