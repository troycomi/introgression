from analyze import summarize_strain_states
from analyze.summarize_strain_states import Summary_Table
from analyze.introgression_configuration import Configuration
import pytest
from io import StringIO


@pytest.fixture
def summarizer():
    config = Configuration()
    config.known_states = 'state1 state2 state3'.split()
    config.set(
        introgressed_intermediate='/dirtag/'
        'blocks_{state}_tag_filtered1intermediate.txt',
        ambiguous_intermediate='/dirtag/'
        'blocks_{state}_tag_filtered2intermediate.txt',
        strain_info='100_genomes_info.txt',
        state_counts='/dirtag/state_counts_by_strain.txt'
    )

    return summarize_strain_states.Strain_Summarizer(config)


def test_run(summarizer, mocker, capsys):
    mock_log = mocker.patch('analyze.summarize_strain_states.log.info')
    mock_read = mocker.patch(
        'analyze.summarize_strain_states.read_table.read_table_rows',
        side_effect=[
            ({'r1': {'strain': 's1', 'start': 10, 'end': 20,
                     'reason': ''},
              'r2': {'strain': 's2', 'start': 25, 'end': 40,
                     'reason': 'test'}},
             ['regions']),
            ({'r1': {'alternative_states': 'state2,state3'},
              'r2': {'not_called': ''}},
             ['regions']),
            ({'r1': {'strain': 's1', 'start': 35, 'end': 40,
                     'reason': ''},
              'r2': {'strain': 's2', 'start': 4, 'end': 8,
                     'reason': ''}},
             ['regions']),
            ({'r1': {'alternative_states': 'state3,state2'},
              'r2': {'alternative_states': 'state3'}},
             ['regions'])
        ])

    handle1 = mocker.MagicMock()
    handle1.__enter__.return_value.__iter__.return_value = \
        ("s1\tnothing\tnothing\tgeo1\tenv1\tpop1\n",
         "s2\tnothing\tnothing\tgeo2\tenv2\tpop2\n")
    handle2 = mocker.MagicMock()

    mocker.patch(
        'analyze.summarize_strain_states.open',
        side_effect=(handle1, handle2)
    )

    summarizer.run()

    assert mock_log.call_args_list == [
        mocker.call('state2'),
        mocker.call('state3'),
    ]

    assert handle1.called_with(mocker.ANY, 'r')

    assert handle2.called_with('/dirtag/state_counts_by_strain.txt', 'w')

    calls = handle2.__enter__().write.call_args_list
    assert handle2.__enter__().write.call_count == 3
    assert calls[0][0] == \
        ('strain\tpopulation\tgeographic_origin\tenvironmental_origin\t'
         'num_regions_state2\tnum_regions_state3\tnum_regions_total\t'
         'num_regions_state2_filtered1\tnum_regions_state3_filtered1\t'
         'num_regions_total_filtered1\tnum_regions_state2_filtered2\t'
         'num_regions_state3_filtered2\tnum_regions_total_filtered2\t'
         'num_regions_state2_filtered2_inclusive\t'
         'num_regions_state3_filtered2_inclusive\t'
         'num_regions_total_filtered2_inclusive\tnum_bases_state2\t'
         'num_bases_state3\tnum_bases_total\tnum_bases_state2_filtered1\t'
         'num_bases_state3_filtered1\tnum_bases_total_filtered1\t'
         'num_bases_state2_filtered2\tnum_bases_state3_filtered2\t'
         'num_bases_total_filtered2\tnum_bases_state2_filtered2_inclusive\t'
         'num_bases_state3_filtered2_inclusive\t'
         'num_bases_total_filtered2_inclusive\t'
         'num_bases_state2_or_state3_filtered2i\tnum_bases_2_filtered2i\n',)
    assert calls[1][0] == \
        ('s1\tpop1\tgeo1\tenv1\t1\t1\t2\t1\t1\t2\t0\t0\t0\t2\t2\t2\t11\t6\t'
         '17\t11\t6\t17\t0\t0\t0\t17\t17\t17\t17\t17\n',)
    assert calls[2][0] == \
        ('s2\tpop2\tgeo2\tenv2\t1\t1\t2\t0\t1\t1\t0\t1\t1\t0\t1\t1\t16\t5\t21'
         '\t0\t5\t5\t0\t5\t5\t0\t5\t5\t0\t0\n',)

    assert mock_read.call_count == 4
    fname = '/dirtag/blocks_state{s}_tag_filtered{f}intermediate.txt'
    mock_read.assert_has_calls([
        mocker.call(fname.format(s=2, f=1), '\t'),
        mocker.call(fname.format(s=2, f=2), '\t'),
        mocker.call(fname.format(s=3, f=1), '\t'),
        mocker.call(fname.format(s=3, f=2), '\t')])


@pytest.fixture
def table():
    return summarize_strain_states.Summary_Table()


def test_table_init(table):
    assert table.table == {}


def test_table_set_region(table):
    assert table.__dict__ == {'table': {}}
    table.set_region('strain', 'species', 100)
    assert table.__dict__ == {'table': {},
                              'strain': 'strain',
                              'species': 'species',
                              'length': 100}


def test_table_record_element(table):
    assert table.table == {}
    table.record_element('s1', 'k1')
    assert table.table == {'s1': {'k1': 1}}

    table.record_element('s1', 'k1', 10)
    assert table.table == {'s1': {'k1': 11}}

    table.record_element('s1', 'k2', -1)
    assert table.table == {'s1': {'k1': 11, 'k2': -1}}

    table.record_element('s2', 'k2', 0)
    assert table.table == {'s1': {'k1': 11, 'k2': -1},
                           's2': {'k2': 0}}


def test_table_record_region(table):
    table.record_region('s1', 'CBS', 100)
    assert table.table == {'s1': {
        'num_regions_CBS': 1,
        'num_regions_total': 1,
        'num_bases_CBS': 100,
        'num_bases_total': 100,
    }}

    table.record_region('s1', 'N', 50)
    assert table.table == {'s1': {
        'num_regions_CBS': 1,
        'num_regions_N': 1,
        'num_regions_total': 2,
        'num_bases_CBS': 100,
        'num_bases_N': 50,
        'num_bases_total': 150,
    }}

    table.record_region('s1', 'N', 50, '_filt')
    assert table.table == {'s1': {
        'num_regions_CBS': 1,
        'num_regions_N': 1,
        'num_regions_total': 2,
        'num_bases_CBS': 100,
        'num_bases_N': 50,
        'num_bases_total': 150,
        'num_regions_N_filt': 1,
        'num_regions_total_filt': 1,
        'num_bases_N_filt': 50,
        'num_bases_total_filt': 50,
    }}

    # add underscore when concatenating
    table.record_region('s1', 'N', 50, 'filt')
    assert table.table == {'s1': {
        'num_regions_CBS': 1,
        'num_regions_N': 1,
        'num_regions_total': 2,
        'num_bases_CBS': 100,
        'num_bases_N': 50,
        'num_bases_total': 150,
        'num_regions_N_filt': 2,
        'num_regions_total_filt': 2,
        'num_bases_N_filt': 100,
        'num_bases_total_filt': 100,
    }}

    # don't update total
    table.record_region('s1', 'N', 50, 'filt', False)
    assert table.table == {'s1': {
        'num_regions_CBS': 1,
        'num_regions_N': 1,
        'num_regions_total': 2,
        'num_bases_CBS': 100,
        'num_bases_N': 50,
        'num_bases_total': 150,
        'num_regions_N_filt': 3,
        'num_regions_total_filt': 2,
        'num_bases_N_filt': 150,
        'num_bases_total_filt': 100,
    }}

    # don't update total
    table.record_region('s2', 'N', 50, 'filt', False)
    assert table.table == {
        's1': {
            'num_regions_CBS': 1,
            'num_regions_N': 1,
            'num_regions_total': 2,
            'num_bases_CBS': 100,
            'num_bases_N': 50,
            'num_bases_total': 150,
            'num_regions_N_filt': 3,
            'num_regions_total_filt': 2,
            'num_bases_N_filt': 150,
            'num_bases_total_filt': 100,
        },
        's2': {
            'num_regions_N_filt': 1,
            'num_bases_N_filt': 50,
        }}


def test_table_record_alt_species(table):
    table.set_region('strain', 'species', 100)
    assert table.table == {}
    table.record_alt_species(['a1', 'a2'])
    assert table.table == {
        'strain': {
            'num_bases_2_filtered2i': 100,
            'num_bases_a1_filtered2_inclusive': 100,
            'num_bases_a1_or_a2_filtered2i': 100,
            'num_bases_a2_filtered2_inclusive': 100,
            'num_regions_a1_filtered2_inclusive': 1,
            'num_regions_a2_filtered2_inclusive': 1,
        }
    }

    table.set_region('strain', 'species', 50)
    table.record_alt_species(['a1'])
    assert table.table == {
        'strain': {
            'num_bases_2_filtered2i': 100,
            'num_bases_a1_filtered2_inclusive': 150,
            'num_bases_a1_or_a2_filtered2i': 100,
            'num_bases_a2_filtered2_inclusive': 100,
            'num_regions_a1_filtered2_inclusive': 2,
            'num_regions_a2_filtered2_inclusive': 1,
            'num_bases_species_filtered2': 50,
            'num_regions_species_filtered2': 1,
            'num_bases_total_filtered2': 50,
            'num_regions_total_filtered2': 1,
            'num_bases_1_filtered2i': 50,
        }
    }


def test_table_region_found(table, mocker):
    table.set_region('strain', 'species', 100)
    mock_record = mocker.patch.object(Summary_Table, 'record_region')
    table.region_found()
    mock_record.assert_called_once_with('strain', 'species', 100)


def test_table_region_passes_filter1(table, mocker):
    table.set_region('strain', 'species', 100)
    mock_record = mocker.patch.object(Summary_Table, 'record_region')
    table.region_passes_filter1()
    mock_record.assert_called_once_with('strain', 'species', 100, '_filtered1')


def test_table_record_alt(table, mocker):
    table.set_region('strain', 'species', 100)
    mock_record = mocker.patch.object(Summary_Table, 'record_region')
    table.record_alt('species2')
    mock_record.assert_called_once_with('strain', 'species2',
                                        100, '_filtered2_inclusive',
                                        False)
    mock_record.reset_mock()
    table.record_alt('species')
    mock_record.assert_called_once_with('strain', 'species',
                                        100, '_filtered2_inclusive',
                                        True)


def test_table_add_strain_info(table):
    table.table = {'s1': {}, 's2': {}}
    strain_info = StringIO(
        's1\t_\t_\tgeo1\tenv1\tpop1\n'  # overwrite below
        'S1\t_\t_\tgeo2\tenv2\tpop2\n'  # upper should cast lower
        's2\t_\t_\tgeo3\tenv3\tpop3\n'
        's3\t_\t_\tgeo4\tenv4\tpop4\n'  # not recorded
    )
    table.add_strain_info(strain_info)
    assert table.table == {
        's1': {
            'population': 'pop2',
            'geographic_origin': 'geo2',
            'environmental_origin': 'env2'
        },
        's2': {
            'population': 'pop3',
            'geographic_origin': 'geo3',
            'environmental_origin': 'env3'
        }
    }


def test_table_add_strain_info_empty(table):
    table.table = {'s1': {}, 's2': {}}
    strain_info = StringIO(
        's1\t\t_\tgeo1\tenv1\tpop1\n'  # should not fail
        's1\t_\t_\tgeo2\tenv2\tpop2\n'
        's2\t_\t_\tgeo3\tenv3\tpop3\n'
        's3\t_\t_\tgeo4\tenv4\tpop4\n'  # not recorded
    )
    table.add_strain_info(strain_info)
    assert table.table == {
        's1': {
            'population': 'pop2',
            'geographic_origin': 'geo2',
            'environmental_origin': 'env2'
        },
        's2': {
            'population': 'pop3',
            'geographic_origin': 'geo3',
            'environmental_origin': 'env3'
        }
    }


def test_table_write_summary(table):
    output = StringIO()
    table.write_summary([], output)
    assert output.getvalue() == (
        'strain\tpopulation\tgeographic_origin\tenvironmental_origin\t'
        'num_regions_total\tnum_regions_total_filtered1\t'
        'num_regions_total_filtered2\tnum_regions_total_filtered2_inclusive\t'
        'num_bases_total\tnum_bases_total_filtered1\tnum_bases_total_filtered2'
        '\tnum_bases_total_filtered2_inclusive\n'
    )

    table.table = {
        's1': {
            'num_regions_total': 15,
            'num_bases_total_filtered2_inclusive': 500,
            'NOTHING': 'NOT READ'
        },
        's2': {
            'population': 'pop2'
        }
    }
    output = StringIO()
    table.write_summary(['state'], output)
    assert output.getvalue() == (
        'strain\tpopulation\tgeographic_origin\tenvironmental_origin\t'
        'num_regions_state\tnum_regions_total\tnum_regions_state_filtered1\t'
        'num_regions_total_filtered1\tnum_regions_state_filtered2\t'
        'num_regions_total_filtered2\tnum_regions_state_filtered2_inclusive\t'
        'num_regions_total_filtered2_inclusive\tnum_bases_state\t'
        'num_bases_total\tnum_bases_state_filtered1\t'
        'num_bases_total_filtered1\t'
        'num_bases_state_filtered2\tnum_bases_total_filtered2\t'
        'num_bases_state_filtered2_inclusive\t'
        'num_bases_total_filtered2_inclusive\n'
        's1\t0\t0\t0\t0\t15\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t500\n'
        's2\tpop2\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n'
    )


def test_table_get_fields(table):
    fields = table.get_fields([])
    assert fields == [
        'population', 'geographic_origin', 'environmental_origin',
        'num_regions_total',
        'num_regions_total_filtered1',
        'num_regions_total_filtered2',
        'num_regions_total_filtered2_inclusive',
        'num_bases_total',
        'num_bases_total_filtered1',
        'num_bases_total_filtered2',
        'num_bases_total_filtered2_inclusive']

    fields = table.get_fields(['s1'])
    assert fields == [
        'population', 'geographic_origin', 'environmental_origin',
        'num_regions_s1', 'num_regions_total', 'num_regions_s1_filtered1',
        'num_regions_total_filtered1', 'num_regions_s1_filtered2',
        'num_regions_total_filtered2', 'num_regions_s1_filtered2_inclusive',
        'num_regions_total_filtered2_inclusive', 'num_bases_s1',
        'num_bases_total', 'num_bases_s1_filtered1',
        'num_bases_total_filtered1', 'num_bases_s1_filtered2',
        'num_bases_total_filtered2', 'num_bases_s1_filtered2_inclusive',
        'num_bases_total_filtered2_inclusive']

    fields = table.get_fields(['s1', 's2', 's3'])
    assert fields == [
        'population', 'geographic_origin', 'environmental_origin',
        'num_regions_s1', 'num_regions_s2', 'num_regions_s3',
        'num_regions_total', 'num_regions_s1_filtered1',
        'num_regions_s2_filtered1', 'num_regions_s3_filtered1',
        'num_regions_total_filtered1', 'num_regions_s1_filtered2',
        'num_regions_s2_filtered2', 'num_regions_s3_filtered2',
        'num_regions_total_filtered2', 'num_regions_s1_filtered2_inclusive',
        'num_regions_s2_filtered2_inclusive',
        'num_regions_s3_filtered2_inclusive',
        'num_regions_total_filtered2_inclusive', 'num_bases_s1',
        'num_bases_s2', 'num_bases_s3', 'num_bases_total',
        'num_bases_s1_filtered1', 'num_bases_s2_filtered1',
        'num_bases_s3_filtered1', 'num_bases_total_filtered1',
        'num_bases_s1_filtered2', 'num_bases_s2_filtered2',
        'num_bases_s3_filtered2', 'num_bases_total_filtered2',
        'num_bases_s1_filtered2_inclusive', 'num_bases_s2_filtered2_inclusive',
        'num_bases_s3_filtered2_inclusive',
        'num_bases_total_filtered2_inclusive', 'num_bases_s1_or_s2_filtered2i',
        'num_bases_s1_or_s3_filtered2i', 'num_bases_s2_or_s3_filtered2i',
        'num_bases_2_filtered2i', 'num_bases_s1_or_s2_or_s3_filtered2i',
        'num_bases_3_filtered2i']
