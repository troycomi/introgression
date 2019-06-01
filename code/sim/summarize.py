import global_params as gp


def read_list(l):
    return [x.strip() for x in l[1:-1].split(',')]


tag = None  # was not defined
gp_dir = '../'
outfilename = gp.sim_out_prefix + tag + '.txt'
results_filename = gp.sim_out_prefix + tag + '_summary.txt'
