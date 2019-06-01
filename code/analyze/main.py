import click
import yaml
import logging as log
import analyze.predict
from analyze.introgression_configuration import Configuration
from analyze.id_regions import ID_producer
from analyze.summarize_region_quality import Summarizer
from analyze.filter_regions import Filterer
from analyze.summarize_strain_states import Strain_Summarizer


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

    ctx.obj.set(log_file=log_file)
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

    config.set('chromosomes')
    log.info(f'Found {len(config.chromosomes)} chromosomes in config')

    config.set(threshold=threshold)
    log.info(f'Threshold value is \'{config.threshold}\'')

    config.set(blocks=blocks)
    log.info(f'Output blocks file is \'{config.blocks}\'')

    config.set('states')
    config.set(prefix=prefix)
    log.info(f'Prefix is \'{config.prefix}\'')

    config.set(strains=test_strains)
    if config.test_strains is None:
        log.info(f'No test_strains provided')
    else:
        str_len = len(config.test_strains)
        log.info(f'Found {str_len} test strain'
                 f'{"" if str_len == 1 else "s"}')
    str_len = len(config.strains)
    log.info(f'Found {str_len} unique strain'
             f'{"" if str_len == 1 else "s"}')

    config.set(hmm_initial=hmm_initial,
               hmm_trained=hmm_trained,
               positions=positions,
               probabilities=probabilities,
               alignment=alignment)
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


@cli.command()
@click.option('--blocks', default='', help='Block file location with {state}')
@click.option('--labeled', default='', help='Block file location with {state}')
@click.option('--state', multiple=True, help='States to add ids to')
@click.pass_context
def id_regions(ctx, blocks, labeled, state):
    config = ctx.obj
    config.set('chromosomes')
    log.info(f'Found {len(config.chromosomes)} chromosomes in config')

    state = list(state)
    config.set(states=state)
    log.info(f'Found {len(config.states)} states to process')

    config.set(blocks=blocks)
    log.info(f'Input blocks file is \'{config.blocks}\'')

    config.set(labeled_blocks=labeled)
    log.info(f'Output blocks file is \'{config.labeled_blocks}\'')

    id_producer = ID_producer(config)
    id_producer.add_ids()


@cli.command()
@click.option('--state', multiple=True, help='States to summarize')
@click.option('--labeled', default='',
              help='Labeled block file with {state} '
              'Created during id_regions')
@click.option('--masks', default='',
              help='Mask file with {strain} and {chrom}')
@click.option('--alignment', default='',
              help='Alignment file with {prefix} [optional], '
              '{strain} and {chrom}')
@click.option('--positions', default='',
              help='Position file created during prediction')
@click.option('--quality', default='',
              help='Output quality file with {state}')
@click.option('--region', default='',
              help='Output region file with {state}, gzipped')
@click.option('--region-index', default='',
              help='Output region index file with {state}, pickled')
@click.pass_context
def summarize_regions(ctx,
                      state,
                      labeled,
                      quality,
                      masks,
                      alignment,
                      positions,
                      region,
                      region_index):
    config = ctx.obj

    config.set('states',
               'chromosomes')
    log.info(f'Found {len(config.chromosomes)} chromosomes in config')

    config.set(labeled_blocks=labeled)
    log.info(f'Labeled blocks file is \'{config.labeled_blocks}\'')

    config.set(quality_blocks=quality)
    log.info(f'Quality file is \'{config.quality_blocks}\'')

    config.set(masks=masks)
    log.info(f'Mask file is \'{config.masks}\'')

    config.set('prefix')
    config.set(alignment=alignment)
    log.info(f'Alignment file is \'{config.alignment}\'')

    config.set(positions=positions)
    log.info(f'Positions file is \'{config.positions}\'')

    config.set(regions=region, region_index=region_index)
    log.info(f'Region file is \'{config.regions}\'')
    log.info(f'Region index file is \'{config.region_index}\'')

    config.set('symbols')

    summarizer = Summarizer(config)
    summarizer.run(list(state))


@cli.command()
@click.option('--thresh', help='Threshold to apply to ambiguous filter',
              default=None, type=float)
@click.option('--introgress-filter', default='',
              help='Filtered block file location with {state}.'
              ' Contains only regions passing introgression filter')
@click.option('--introgress-inter', default='',
              help='Filtered block file location with {state}.'
              ' Contains all regions with reasons they failed filtering')
@click.option('--ambiguous-filter', default='',
              help='Filtered block file location with {state}.'
              ' Contains only regions passing ambiguous filter')
@click.option('--ambiguous-inter', default='',
              help='Filtered block file location with {state}.'
              ' Contains all regions passing introgressing filtering, '
              'with reasons they failed ambiguous filtering')
@click.option('--filter-sweep', default='',
              help='Contains summary results for applying ambiguous filter '
              'with various threshold values supplied as arguments.')
@click.option('--region', default='',
              help='Region file with {state}, gzipped')
@click.option('--region-index', default='',
              help='Region index file with {state}, pickled')
@click.option('--quality', default='',
              help='Quality file with {state}')
@click.argument('thresholds', nargs=-1, type=float)
@click.pass_context
def filter_regions(ctx,
                   thresh,
                   introgress_filter,
                   introgress_inter,
                   ambiguous_filter,
                   ambiguous_inter,
                   filter_sweep,
                   region,
                   region_index,
                   quality,
                   thresholds):
    config = ctx.obj  # type: Configuration
    config.set('states')

    config.set(filter_threshold=thresh)
    log.info(f"Filter threshold set to '{config.filter_threshold}'")

    config.set(introgressed=introgress_filter,
               introgressed_intermediate=introgress_inter,
               ambiguous=ambiguous_filter,
               ambiguous_intermediate=ambiguous_inter,
               filter_sweep=filter_sweep)
    log.info(f"Introgressed filtered file is '{config.introgressed}'")
    log.info('Introgressed intermediate file is '
             f"'{config.introgressed_intermediate}'")
    log.info(f"Ambiguous filtered file is '{config.ambiguous}'")
    log.info('Ambiguous intermediate file is '
             f"'{config.ambiguous_intermediate}'")
    if config.filter_sweep is not None:
        log.info(f"Filter sweep file is '{config.filter_sweep}'")

    config.set(regions=region,
               region_index=region_index)
    log.info(f'Region file is \'{config.regions}\'')
    log.info(f'Region index file is \'{config.region_index}\'')

    config.set(quality_blocks=quality)
    log.info(f'Quality file is \'{config.quality_blocks}\'')

    config.set('symbols')

    thresholds = list(thresholds)
    log.info(f'Threshold sweep with: {thresholds}')

    filterer = Filterer(config)
    filterer.run(thresholds)


@cli.command()
@click.option('--introgress-inter', default='',
              help='Filtered block file location with {state}.'
              ' Contains all regions with reasons they failed filtering')
@click.option('--ambiguous-inter', default='',
              help='Filtered block file location with {state}.'
              ' Contains all regions passing introgressing filtering, '
              'with reasons they failed ambiguous filtering')
@click.option('--strain-info', default='',
              help='Tab separated table with strain name, alternate name, '
              'location, envionment, and population')
@click.option('--state-counts', default='',
              help='Output state summary file')
@click.pass_context
def summarize_strains(ctx,
                      introgress_inter,
                      ambiguous_inter,
                      strain_info,
                      state_counts):
    config = ctx.obj  # type: Configuration
    config.set('states')
    config.set(introgressed_intermediate=introgress_inter,
               ambiguous_intermediate=ambiguous_inter,
               strain_info=strain_info,
               state_counts=state_counts)
    log.info('Introgressed intermediate file is '
             f"'{config.introgressed_intermediate}'")
    log.info('Ambiguous intermediate file is '
             f"'{config.ambiguous_intermediate}'")
    log.info(f"Strain information from '{config.strain_info}'")
    log.info(f"State counts saved to '{config.state_counts}'")
    strain_summarizer = Strain_Summarizer(config)
    strain_summarizer.run()


if __name__ == '__main__':
    cli()
