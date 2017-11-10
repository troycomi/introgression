import sys
import os
import process_args
import sim_process
from sim_predict import *
sys.path.append('..')
import global_params as gp

##======
# read in simulation parameters
##======

sim_tag = sys.argv[2]
sim_args = process_args.process_args_by_tag(sys.argv[1], sim_tag)
predict_args, last_read = process_args(sys.argv, sim_args, i=2)

# need to read in all sim args so that we can find the one with the
# correct tag
all_sim_args = process_args.process_all_args(sys.argv[1])
# then read prediction-specific args and combine those with sim args
# for the correct tag
args, last_read = process_args(sys.argv, all_sim_args, i=2)

##======
# loop through all simulations predict introgression
##======

gp_dir = '../'
# for reading output from ms
ms_f = open(gp_dir + gp.sim_out_dir + '/ms/' + gp.sim_out_prefix + \
            args['tag'] + '.txt', 'r')
# summary output
out_f = open(gp_dir + gp.sim_out_dir + gp.sim_out_prefix + \
             args['tag'] + '_hmm_' + args['predict_tag'] + '.txt', 'w')
# summary output
out_init_f = open(gp_dir + gp.sim_out_dir + gp.sim_out_prefix + \
                  args['tag'] + '_hmm_init_' + args['predict_tag'] + '.txt', 'w')
# introgression output
introgression_f = open(gp_dir + gp.sim_out_dir + gp.sim_out_prefix + \
                       args['tag'] + '_introgressed_predicted_' + \
                       args['predict_tag'] + '.txt', 'w')

for i in range(args['num_reps']):
    
    print i

    ##======
    # read in simulated sequences
    ##======
    
    # trees, recomb_sites, seg_sites, positions, seqs
    sim = sim_process.read_one_sim(ms_f, args['num_sites'], args['num_samples'])

    ##======
    # predict introgressed/non-introgressed tracts
    ##======
    
    state_seq, probs, hmm, hmm_init = predict_introgressed(sim, args, train=True)
    state_seq_blocks = sim_process.convert_to_blocks(state_seq, \
                                                     args['states'])

    ##======
    # output
    ##======

    # summary info about HMM (before training)
    write_hmm_line(hmm_init, out_init_f, i==0) 

    # summary info about HMM (after training)
    write_hmm_line(hmm, out_f, i==0) 

    # specific locations of introgression (for comparing predictions
    # to)
    sim_process.write_introgression_blocks(state_seq_blocks, introgression_f, \
                                           i, args['states'])

    # TODO file format that gives probabilities
    
ms_f.close()
out_f.close()
introgression_f.close()

