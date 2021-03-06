import sys
import process_args
import sim_process
from compare_introgressed import (write_compare_header,
                                  count_bases,
                                  write_compare_line)
import global_params as gp

# takes in two sets of introgression calls on same set of data and
# compares them (for example, actual and predicted regions, or calls
# from two different prediction methods)

sim_tag = sys.argv[1]
suffix1 = sys.argv[2]  # actual, or first prediction method
suffix2 = sys.argv[3]  # prediction method
args = process_args.process_args_by_tag(sys.argv[4], sim_tag)

gp_dir = '../'
fn1 = gp.sim_out_dir_absolute + gp.sim_out_prefix + sim_tag + \
      '_introgressed_' + suffix1 + '.txt'
fn2 = gp.sim_out_dir_absolute + gp.sim_out_prefix + sim_tag + \
      '_introgressed_' + suffix2 + '.txt'
f_out = gp.sim_out_dir_absolute + gp.sim_out_prefix + sim_tag + \
      '_introgressed_compare_' + suffix1 + '_' + suffix2 + '.txt'

f1 = open(fn1, 'r')
f2 = open(fn2, 'r')
f_out = open(f_out, 'w')

line1 = f1.readline()
line2 = f2.readline()

write_compare_header(f_out, args['species'], suffix1, suffix2)
while line1 != '' and line2 != '':

    d1, rep1, line1 = sim_process.read_introgression_blocks(f1, line1,
                                                            args['species'])
    d2, rep2, line2 = sim_process.read_introgression_blocks(f2, line2,
                                                            args['species'])
    assert rep1 == rep2, str(rep1) + ' ' + str(rep2)
    print('rep', rep1)
    print(d1)
    print(d2)

    base_counts, avg_base_counts = count_bases(d1, d2, args,
                                               suffix1, suffix2, 'par')

    write_compare_line(avg_base_counts, f_out, args['species'],
                       suffix1, suffix2)

if line1 != '' or line2 != '':
    print('one of these files is incomplete and for some reason '
          'I\'m not bothering to tell you which!')

f_out.close()
