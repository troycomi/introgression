# combine all chromosomal alignments into one master
# indexed relative to cerevisiae reference

from misc import read_maf
import global_params as gp

complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
              'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
              'N': 'N', 'n': 'n', '-': '-'}

flip = {'-': '+', '+': '-'}


def reverse_start(start, length, total_length):
    return total_length - start - length


def reverse_complement(s):
    r = []
    for b in s[::-1]:
        r.append(complement[b])
    return r


def forward_index(blocks):
    # go through all blocks and add a field for start relative to
    # forward strand, and sequence in forward direction
    for label in blocks.keys():
        for strain in blocks[label]['strains'].keys():

            start = blocks[label]['strains'][strain]['start']
            seq = blocks[label]['strains'][strain]['sequence']

            blocks[label]['strains'][strain]['forward_start'] = start
            blocks[label]['strains'][strain]['forward_sequence'] = seq

            if blocks[label]['strains'][strain]['strand'] == '-':
                blocks[label]['strains'][strain]['forward_sequence'] = \
                    seq[::-1]
                blocks[label]['strains'][strain]['forward_start'] = \
                    reverse_start(
                        start,
                        blocks[label]['strains'][strain]['length'],
                        blocks[label]['strains'][strain]['aligned_length'])

    return blocks


def master_forward(blocks, master):
    # make all master sequences go in forward direction (+) and flip
    # others as necessary
    for label in blocks.keys():
        if master in blocks[label]['strains']:
            if blocks[label]['strains'][master]['strand'] == '-':
                for strain in blocks[label]['strains'].keys():
                    aligned_length = \
                        blocks[label]['strains'][strain]['aligned_length']
                    seq = blocks[label]['strains'][strain]['sequence']
                    start = blocks[label]['strains'][strain]['start']
                    length = blocks[label]['strains'][strain]['length']

                    rc_seq = reverse_complement(seq)
                    rc_start = reverse_start(start, length, aligned_length)

                    blocks[label]['strains'][strain]['sequence'] = rc_seq
                    blocks[label]['strains'][strain]['strand'] = \
                        flip[blocks[label]['strains'][strain]['strand']]
                    blocks[label]['strains'][strain]['start'] = rc_start

    return blocks


def make_master(fn, master):

    # keyed by block label; most of info in each keyed by ['strains'][strain]
    blocks = read_maf.read_mugsy(fn)

    # flip all blocks so that master sequence is on + strand
    blocks = master_forward(blocks, master)
    # add fields giving index and sequence relative to + strand
    blocks = forward_index(blocks)

    # make sequences with alignment columns present in master
    n = blocks['1']['strains'][master]['aligned_length']
    all_strains = blocks['1']['strains'].keys()
    a = dict(zip(all_strains,
                 [[gp.unaligned_symbol] * n for s in all_strains]))

    # loop through all blocks
    for label in blocks.keys():
        # only care about aligned blocks that include master sequence
        if master in blocks[label]['strains']:
            absolute_ind = blocks[label]['strains'][master]['start']
            master_seq = blocks[label]['strains'][master]['sequence']
            block_length = len(master_seq)
            strains = blocks[label]['strains'].keys()

            # loop through all positions in the alignment block
            for relative_ind in range(block_length):
                # only care about positions that aren't gaps in master sequence
                if master_seq[relative_ind] != gp.gap_symbol:
                    # make sure we haven't already dealt with this position

                    # apparently mugsy sometimes aligns the same part
                    # of one genome to multiple parts of another
                    # genome. this is a problem.
                    assert a[master][absolute_ind] == gp.unaligned_symbol,\
                        absolute_ind
                    # loop through all the strains in this block
                    for strain in strains:
                        a[strain][absolute_ind] = \
                            blocks[label]['strains'][strain][
                                'forward_sequence'][relative_ind]
                    absolute_ind += 1

    for strain in all_strains:
        a[strain] = ''.join(a[strain])
        print(strain, a[strain].count(gp.unaligned_symbol))

    return a


def write_master(fn, a):
    f = open(fn, 'w')
    for strain in a.keys():
        f.write('> ' + strain + '\n')
        f.write(a[strain] + '\n')
    f.close()
