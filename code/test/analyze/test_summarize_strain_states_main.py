import analyze.summarize_strain_states_main as main


def test_main(mocker, capsys):
    mocker.patch(
        'analyze.summarize_strain_states_main.predict.process_predict_args',
        return_value={
            'known_states': ['state1', 'state2', 'state3'],
            'tag': 'tag'
        })
    mocker.patch(
        'analyze.summarize_strain_states_main.gp.analysis_out_dir_absolute',
        '/dir')
    mocker.patch(
        'analyze.summarize_strain_states_main.gp.alignment_ref_order',
        ['state1', 'state2', 'state3'])

    mock_read = mocker.patch(
        'analyze.summarize_strain_states_main.read_table.read_table_rows',
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
        'analyze.summarize_strain_states_main.open',
        side_effect=(handle1, handle2)
    )

    main.main()

    captured = capsys.readouterr().out
    assert captured == 'state2\nstate3\n'

    assert handle1.called_with(mocker.ANY, 'r')
    assert handle2.called_with('/dirtag/state_counts_by_strain.txt', 'w')
    calls = handle2.write.call_args_list
    assert handle2.write.call_count == 7
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
    assert calls[1][0] == ('s1\t',)
    assert calls[2][0] == \
        ('pop1\tgeo1\tenv1\t1\t1\t2\t1\t1\t2\t0\t0\t0\t2\t2\t2\t11\t6\t17\t'
         '11\t6\t17\t0\t0\t0\t17\t17\t17\t17\t17',)
    assert calls[3][0] == ('\n',)
    assert calls[4][0] == ('s2\t',)
    assert calls[5][0] == \
        ('pop2\tgeo2\tenv2\t1\t1\t2\t0\t1\t1\t0\t1\t1\t0\t1\t1\t16\t5\t21\t0\t'
         '5\t5\t0\t5\t5\t0\t5\t5\t0\t0',)
    assert calls[6][0] == ('\n',)

    assert mock_read.call_count == 4
    fname = '/dirtag/blocks_state{s}_tag_filtered{f}intermediate.txt'
    mock_read.assert_has_calls([
        mocker.call(fname.format(s=2, f=1), '\t'),
        mocker.call(fname.format(s=2, f=2), '\t'),
        mocker.call(fname.format(s=3, f=1), '\t'),
        mocker.call(fname.format(s=3, f=2), '\t')])
