from analyze import filter_2_main as main


def test_main(mocker, capsys):
    mocker.patch('sys.argv', ['', '0.1'])
    mocker.patch('analyze.filter_2_main.predict.process_predict_args',
                 return_value={
                     'known_states': ['state1', 'state2'],
                     'tag': 'tag'
                 })
    mocker.patch('analyze.filter_2_main.gp.analysis_out_dir_absolute',
                 '/dir')
    mocker.patch('analyze.filter_2_main.read_table.read_table_rows',
                 return_value=({'r1': {}, 'r2': {'a': 1}}, ['regions']))

    mocked_file = mocker.patch('analyze.filter_2_main.open')

    mock_read = mocker.patch('analyze.filter_2_main.Region_Reader')
    mock_read().__enter__().yield_fa.return_value = iter([
        ('r1', ['> seq', '> info'], ['atcg', 'x..']),
        ('r2', ['> seq', '> info'], ['atcg', 'x..'])])

    mock_filter = mocker.patch('analyze.filter_2_main.passes_filters2',
                               side_effect=[
                                   (False, ['1', '2'], [0.8, 0.5], [2, 1, 0]),
                                   (True, ['1'], [0.8], [2])
                               ])
    mock_write = mocker.patch('analyze.filter_2_main.write_filtered_line')

    main.main()

    captured = capsys.readouterr().out
    assert captured == 'state2\n'

    assert mock_read.call_count == 2
    mock_read.assert_called_with('/dirtag/regions/state2.fa.gz', as_fa=True)

    assert mocked_file.call_count == 2
    mocked_file.assert_any_call(
        '/dirtag/blocks_state2_tag_filtered2intermediate.txt', 'w')
    mocked_file.assert_any_call(
        '/dirtag/blocks_state2_tag_filtered2.txt', 'w')

    # just headers, capture others
    mocked_file().__enter__().write.assert_has_calls([
        mocker.call('regions\talternative_states\t'
                    'alternative_ids\talternative_P_counts\n'),
        mocker.call('regions\n')])

    assert mock_filter.call_count == 2
    # seems like this references the object, which changes after call
    mock_filter.assert_has_calls([
        mocker.call(
            {'alternative_states': '1,2',
             'alternative_ids': '0.8,0.5',
             'alternative_P_counts': '2,1,0'},
            ['atcg'], 0.1, ['state1', 'state2']),
        mocker.call(
            {'a': 1,
             'alternative_states': '1',
             'alternative_ids': '0.8',
             'alternative_P_counts': '2'},
            ['atcg'], 0.1, ['state1', 'state2'])])

    assert mock_write.call_count == 3
    mock_write.assert_has_calls([
        mocker.call(mocker.ANY, 'r1',
                    {'alternative_states': '1,2',
                     'alternative_ids': '0.8,0.5',
                     'alternative_P_counts': '2,1,0'},
                    ['regions', 'alternative_states',
                     'alternative_ids', 'alternative_P_counts']
                    ),
        mocker.call(mocker.ANY, 'r2',
                    {'a': 1,
                     'alternative_states': '1',
                     'alternative_ids': '0.8',
                     'alternative_P_counts': '2'},
                    ['regions', 'alternative_states',
                     'alternative_ids', 'alternative_P_counts']
                    ),
        mocker.call(mocker.ANY, 'r2',
                    {'a': 1,
                     'alternative_states': '1',
                     'alternative_ids': '0.8',
                     'alternative_P_counts': '2'},
                    ['regions']
                    )
    ])
