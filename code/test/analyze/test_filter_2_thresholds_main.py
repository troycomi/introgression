from analyze import filter_2_thresholds_main as main


def test_main(mocker, capsys):
    mocker.patch('sys.argv', ['', '0.1'])
    mocker.patch(
        'analyze.filter_2_thresholds_main.predict.process_predict_args',
        return_value={
            'known_states': ['state1', 'state2'],
            'tag': 'tag'
        })
    mocker.patch(
        'analyze.filter_2_thresholds_main.thresholds',
        [0.99, 0.95])
    mocker.patch(
        'analyze.filter_2_thresholds_main.gp.analysis_out_dir_absolute',
        '/dir')
    mocker.patch('analyze.filter_2_thresholds_main.read_table.read_table_rows',
                 return_value=({'r1': {}, 'r2': {'a': 1}}, ['regions']))

    mocked_file = mocker.patch('analyze.filter_2_thresholds_main.open')
    mock_fasta = mocker.patch(
        'analyze.filter_2_thresholds_main.read_fasta.read_fasta',
        return_value=(['> seq', '> info'],
                      ['atcg', 'x..']))
    mock_filter = mocker.patch(
        'analyze.filter_2_thresholds_main.passes_filters2',
        side_effect=[
            (False, ['1', '2'], [0.8, 0.5], [2, 1, 0]),
            (True, ['1'], [0.8], [2]),
            (True, ['1'], [0.8], [2]),
            (False, ['1', '2'], [0.8, 0.5], [2, 1, 0])
        ])

    main.main()

    captured = capsys.readouterr().out
    assert captured == '0.99\n* state2\n0.95\n* state2\n'

    assert mock_fasta.call_count == 4
    mock_fasta.assert_has_calls([
        mocker.call('/dirtag/regions/r1.fa.gz', gz=True),
        mocker.call('/dirtag/regions/r2.fa.gz', gz=True),
        mocker.call('/dirtag/regions/r1.fa.gz', gz=True),
        mocker.call('/dirtag/regions/r2.fa.gz', gz=True)])

    assert mocked_file.call_count == 1
    mocked_file.assert_any_call(
        '/dirtag/filter_2_thresholds_tag.txt', 'w')

    mocked_file().write.assert_has_calls([
        mocker.call('threshold\tpredicted_state\talternative_states\tcount\n'),
        mocker.call('0.99\tstate2\t1\t1\n'),
        mocker.call('0.99\tstate2\t1,2\t1\n'),
        mocker.call('0.95\tstate2\t1\t1\n'),
        mocker.call('0.95\tstate2\t1,2\t1\n'),
        ])

    assert mock_filter.call_count == 4
    # seems like this references the object, which changes after call
    print(mock_filter.call_args_list)
    mock_filter.assert_has_calls([
        mocker.call({}, ['atcg'], 0.99),
        mocker.call({'a': 1}, ['atcg'], 0.99),
        mocker.call({}, ['atcg'], 0.95),
        mocker.call({'a': 1}, ['atcg'], 0.95),
        ])
