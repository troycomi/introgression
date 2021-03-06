

def write_coding_table(seqs_coded, files, rep, header):

    for i in files.keys():
        f = files[i]
        if header:
            f.write('site\tcode\trep\n')
        for j in range(len(seqs_coded[i])):
            f.write(str(j) + '\t' + seqs_coded[i][j] + '\t' + str(rep) + '\n')


def write_one_block_set(blocks_dic, ind, f, rep, suffix=''):

    for block_type in blocks_dic[ind].keys():
        for species in blocks_dic[ind][block_type].keys():
            for block in blocks_dic[ind][block_type][species]:
                f.write(block_type + suffix + '\t' + species + '\t' +
                        str(block[0]) + '\t' + str(block[1]) + '\t' +
                        str(rep) + '\n')


def write_blocks_table(blocks_dic, files, rep, ref_ind, header):

    for ind in blocks_dic.keys():
        f = files[int(ind)]
        if header:
            f.write('type\tspecies\tstart\tend\trep\n')
        # write blocks for the current individual
        write_one_block_set(blocks_dic, ind, f, rep)
        # and also the reference individual
        write_one_block_set(blocks_dic, ref_ind, f, rep, suffix='_ref')
