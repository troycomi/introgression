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
    mock_read = mocker.patch('analyze.filter_2_thresholds_main.Region_Reader')
    mock_read().__enter__().yield_fa.return_value = iter([
        ('r1', ['> seq', '> info'], ['atcg', 'x..']),
        ('r2', ['> seq', '> info'], ['atcg', 'x..'])])
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
    assert captured == '* state2\n'

    assert mock_read.call_count == 2
    mock_read.assert_called_with('/dirtag/regions/state2.fa.gz', as_fa=True)

    assert mocked_file.call_count == 1
    mocked_file.assert_any_call(
        '/dirtag/filter_2_thresholds_tag.txt', 'w')

    mocked_file().__enter__().write.assert_has_calls([
        mocker.call('threshold\tpredicted_state\talternative_states\tcount\n'),
        mocker.call('0.99\tstate2\t1,2\t1\n'),
        mocker.call('0.99\tstate2\t1\t1\n'),
        mocker.call('0.95\tstate2\t1\t1\n'),
        mocker.call('0.95\tstate2\t1,2\t1\n'),
        ])

    assert mock_filter.call_count == 4
    print(mock_filter.call_args_list)
    mock_filter.assert_has_calls([
        mocker.call({}, ['atcg'], 0.99),
        mocker.call({}, ['atcg'], 0.95),
        mocker.call({'a': 1}, ['atcg'], 0.99),
        mocker.call({'a': 1}, ['atcg'], 0.95),
        ])


def test_record_data_hit():
    dt = {}
    main.record_data_hit(dt, 0.9, 's1', 'k1')
    assert dt == {0.9: {'s1': {'k1': 1}}}
    main.record_data_hit(dt, 0.9, 's1', 'k1')
    main.record_data_hit(dt, 0.9, 's1', 'k1')
    assert dt == {0.9: {'s1': {'k1': 3}}}
    main.record_data_hit(dt, 0.9, 's1', 'k2')
    assert dt == {
        0.9: {
            's1': {'k1': 3, 'k2': 1}
        }
    }
    main.record_data_hit(dt, 0.9, 's2', 'k2')
    assert dt == {
        0.9: {
            's1': {'k1': 3, 'k2': 1},
            's2': {'k2': 1}
        }
    }
    main.record_data_hit(dt, 0.8, 's2', 'k2')
    assert dt == {
        0.9: {
            's1': {'k1': 3, 'k2': 1},
            's2': {'k2': 1}
        },
        0.8: {
            's2': {'k2': 1}
        }
    }
    main.record_data_hit(dt, 0.9, 's2', 'k2')
    assert dt == {
        0.9: {
            's1': {'k1': 3, 'k2': 1},
            's2': {'k2': 2}
        },
        0.8: {
            's2': {'k2': 1}
        }
    }
