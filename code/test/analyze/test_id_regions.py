from analyze import id_regions
import pytest
from analyze.introgression_configuration import Configuration


@pytest.fixture
def id_producer():
    config = Configuration()
    config.add_config({
        'analysis_params':
        {'reference': {'name': 'S288c'},
         'known_states': [
             {'name': 'CBS432'},
             {'name': 'N_45'},
             {'name': 'DBVPG6304'},
             {'name': 'UWOPS91_917_1'},
         ],
         'unknown_states': [{'name': 'unknown'}]
         }
    })
    config.set_states()
    result = id_regions.ID_producer(config)
    return result


def test_producer(id_producer):
    assert id_producer.config.known_states == \
        'S288c CBS432 N_45 DBVPG6304 UWOPS91_917_1'.split()
    assert id_producer.config.unknown_states == \
        'unknown'.split()
    assert id_producer.config.states == \
        'S288c CBS432 N_45 DBVPG6304 UWOPS91_917_1 unknown'.split()


def test_add_ids_empty(id_producer, mocker):
    id_producer.config.add_config({
        'chromosomes': ['I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'IX',
                        'X', 'XI', 'XII', 'XIII', 'XIV', 'XV', 'XVI'],
        'paths': {'analysis': {'blocks': 'dir/blocks_{state}.txt',
                               'labeled_blocks':
                               'dir/blocks_{state}_labeled.txt',
                               }}})

    id_producer.config.states = 'ref state1 unknown'.split()
    id_producer.config.set_blocks_file()
    id_producer.config.set_labeled_blocks_file()
    id_producer.config.set_chromosomes()

    mocker.patch('analyze.id_regions.read_blocks',
                 return_value={})

    mocked_file = mocker.patch('analyze.id_regions.open',
                               mocker.mock_open())

    id_producer.add_ids()

    assert mocked_file.call_count == 3
    mocked_file.assert_any_call('dir/blocks_ref_labeled.txt', 'w')
    mocked_file.assert_any_call('dir/blocks_state1_labeled.txt', 'w')
    mocked_file.assert_any_call('dir/blocks_unknown_labeled.txt', 'w')

    # just headers
    mocked_file().write.assert_has_calls([
        mocker.call('region_id\tstrain\tchromosome\tpredicted_species'
                    '\tstart\tend\tnum_sites_hmm\n')
    ]*3)


def test_add_ids(id_producer, mocker):
    id_producer.config.add_config({
        'chromosomes': ['I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'IX',
                        'X', 'XI', 'XII', 'XIII', 'XIV', 'XV', 'XVI'],
        'paths': {'analysis': {'blocks': 'dir/blocks_{state}.txt',
                               'labeled_blocks':
                               'dir/blocks_{state}_labeled.txt',
                               }}})

    id_producer.config.states = 'ref state1 unknown'.split()
    id_producer.config.set_blocks_file()
    id_producer.config.set_labeled_blocks_file()
    id_producer.config.set_chromosomes()

    regions = [
        {
            'strain1': {
                'I': [(10, 100, 10), (10, 100, 1)],
                'VI': [(10, 100, 10), (10, 100, 1)],
            },
            'strain2': {
                'V': [(10, 100, 10), (10, 100, 1)],
            },
            'strain3': {
                'III': [(10, 100, 10), (10, 100, 1)],
            }
        },
        {
            'strain1': {
                'IX': [(10, 100, 10), (10, 100, 1)],
            },
            'strain2': {
                'II': [(10, 100, 10), (10, 100, 1)],
            },
            'strain3': {
                'X': [(10, 100, 10), (10, 100, 1)],
            }
        },
        {}
    ]
    mocker.patch('analyze.id_regions.read_blocks',
                 side_effect=regions)

    mocked_file = mocker.patch('analyze.id_regions.open',
                               mocker.mock_open())

    id_producer.add_ids()

    assert mocked_file.call_count == 3
    mocked_file.assert_any_call('dir/blocks_ref_labeled.txt', 'w')
    mocked_file.assert_any_call('dir/blocks_state1_labeled.txt', 'w')
    mocked_file.assert_any_call('dir/blocks_unknown_labeled.txt', 'w')

    # headers
    calls = [
        mocker.call('region_id\tstrain\tchromosome\tpredicted_species'
                    '\tstart\tend\tnum_sites_hmm\n')
    ]*3 + [
        mocker.call('r1\tstrain1\tI\tref\t10\t100\t10\n'),
        mocker.call('r2\tstrain1\tI\tref\t10\t100\t1\n'),
        mocker.call('r3\tstrain2\tII\tstate1\t10\t100\t10\n'),
        mocker.call('r4\tstrain2\tII\tstate1\t10\t100\t1\n'),
        mocker.call('r5\tstrain3\tIII\tref\t10\t100\t10\n'),
        mocker.call('r6\tstrain3\tIII\tref\t10\t100\t1\n'),
        mocker.call('r7\tstrain2\tV\tref\t10\t100\t10\n'),
        mocker.call('r8\tstrain2\tV\tref\t10\t100\t1\n'),
        mocker.call('r9\tstrain1\tVI\tref\t10\t100\t10\n'),
        mocker.call('r10\tstrain1\tVI\tref\t10\t100\t1\n'),
        mocker.call('r11\tstrain1\tIX\tstate1\t10\t100\t10\n'),
        mocker.call('r12\tstrain1\tIX\tstate1\t10\t100\t1\n'),
        mocker.call('r13\tstrain3\tX\tstate1\t10\t100\t10\n'),
        mocker.call('r14\tstrain3\tX\tstate1\t10\t100\t1\n'),
    ]
    mocked_file().write.assert_has_calls(calls)


def test_validate_arguments(id_producer):
    with pytest.raises(ValueError) as e:
        id_producer.validate_arguments()
    assert ('Failed to validate ID Producer, '
            "required argument 'chromosomes' was unset") in str(e)

    config = id_producer.config
    config.chromosomes = 1
    config.blocks = 1
    config.labeled_blocks = 1
    config.states = 1

    assert id_producer.validate_arguments()
