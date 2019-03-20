from analyze import filter_1_main as main


def test_main(mocker, capsys):
    mocker.patch('analyze.filter_1_main.predict.process_predict_args',
                 return_value={
                     'known_states': ['state1', 'state2'],
                     'tag': 'tag'
                 })
    mocker.patch('analyze.filter_1_main.gp.analysis_out_dir_absolute',
                '/dir')
    mocker.patch('analyze.filter_1_main.read_table.read_table_rows',
                 return_value=({'r1': {}, 'r2': {'a': 1}}, ['regions']))
    mocked_file = mocker.patch('analyze.filter_1_main.open')
    mock_fasta = mocker.patch('analyze.filter_1_main.read_fasta.read_fasta',
                              return_value=(['> seq', '> info'],
                                            ['atcg', 'x..']))
    mock_filter = mocker.patch('analyze.filter_1_main.passes_filters1',
                               side_effect=[(False, 'test'),  # r1
                                            (True, '')])  # r2
    mock_write = mocker.patch('analyze.filter_1_main.write_filtered_line')

    main.main()

    captured = capsys.readouterr().out
    assert captured == 'state2\n'

    assert mock_fasta.call_count == 2
    mock_fasta.assert_any_call('/dirtag/regions/r2.fa.gz', gz=True)
    mock_fasta.assert_any_call('/dirtag/regions/r1.fa.gz', gz=True)

    assert mocked_file.call_count == 2
    mocked_file.assert_any_call(
        '/dirtag/blocks_state2_tag_filtered1intermediate.txt', 'w')
    mocked_file.assert_any_call(
        '/dirtag/blocks_state2_tag_filtered1.txt', 'w')

    # just headers, capture others
    mocked_file().write.assert_has_calls([
        mocker.call('regions\treason\n'),
        mocker.call('regions\n')])

    assert mock_filter.call_count == 2
    # seems like this references the object, which changes after call
    mock_filter.assert_has_calls([
        mocker.call({'reason': 'test'}, 'x..'),
        mocker.call({'reason': '', 'a': 1}, 'x..')])

    assert mock_write.call_count == 3
    mock_write.assert_has_calls([
        mocker.call(mocker.ANY, 'r1', {'reason': 'test'},
                    ['regions', 'reason']),
        mocker.call(mocker.ANY, 'r2', {'a': 1, 'reason': ''},
                    ['regions', 'reason']),
        mocker.call(mocker.ANY, 'r2', {'a': 1, 'reason': ''},
                    ['regions']),
    ])
