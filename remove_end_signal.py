##########################################################################################################
# Deal with end signal as softclip
# Author: Zhe Jiang
# Created on：2022/06/07
##########################################################################################################

##############################################################
# import packages
##############################################################

import pysam
import realignment_forward
import sys
from multiprocessing import Pool, Manager
import random
import string
import shutil
import os


##############################################################
# define functions
##############################################################

# define a function to solve end signals
def solve_end_sig(input_bam_ses, output_bam_ses, max_soft_clip_ses, end_length_ses):

    save_ses = pysam.set_verbosity(0)  # help suppressingHTTSlib's warning message
    bam_input_ses = pysam.AlignmentFile(input_bam_ses, "rb")
    pysam.set_verbosity(save_ses)

    # load all reads in allreads
    allreads_ses = bam_input_ses.fetch(until_eof=True)

    # create a bam file to store reads after realignment
    bam_header_ses = bam_input_ses.header  # record the header

    output_reads_ls_ses = []

    for read_ses in allreads_ses:

        cigar_change_result_ses = signal_cigar_change(read_ses.cigarstring, max_soft_clip_ses, end_length_ses)
            # cigar_change_result_ses[0]: True False, True is keep, False is drop
            # cigar_change_result_ses[1]: start site move
            # cigar_change_result_ses[2]: new cigar

#        print(cigar_change_result_ses)

        if cigar_change_result_ses[0]:

            out_read_ses = pysam.AlignedSegment()
            out_read_ses.query_name = read_ses.query_name  # name
            out_read_ses.flag = read_ses.flag  # flag
            out_read_ses.reference_id = read_ses.reference_id  # reference name
            out_read_ses.reference_start = read_ses.reference_start + cigar_change_result_ses[1]  # reference position
            out_read_ses.mapping_quality = read_ses.mapping_quality  # mapping quality
            out_read_ses.cigar = realignment_forward.cigar_change(cigar_change_result_ses[2])  # cigar
            out_read_ses.next_reference_id = read_ses.next_reference_id  # reference name of paired read (mate)
            out_read_ses.next_reference_start = read_ses.next_reference_start  # reference position of paired read (mate)
            out_read_ses.template_length = read_ses.template_length  # length between read and mate
            out_read_ses.query_sequence = read_ses.query_sequence  # squence
            out_read_ses.query_qualities = read_ses.query_qualities  # quality
            out_read_ses.tags = (("NH", read_ses.get_tag("NH")),  # set "NH" tag to store mapping times
                             ("SC", read_ses.get_tag("SC")))  # create a tag "SC" to store best result

            output_reads_ls_ses.append(out_read_ses)

    bam_output_ses = pysam.AlignmentFile(output_bam_ses, "wb", header = bam_header_ses)

    for out_read_ses in output_reads_ls_ses:
        bam_output_ses.write(out_read_ses)

    bam_output_ses.close()


# define a funtion to change cigar and start sites
def signal_cigar_change(input_cigar_scc, max_soft_clip_scc, end_length_scc):

    start_sites_dif_scc = 0                 # difference of start site
    soft_clip_scc = 0                       # soft clip statistcs
    length_before_gap_scc = 0               # length before deletion
    del_scc = 0                             # length of deletion signal
    cig_num_scc = ''
    len_cig_scc = len(input_cigar_scc)

    output_cigar_scc = ""
    temp_cigar_scc = ""

    keep_bool = True                        # True: keep the reads, False: drop the reads
    #soft_clip_bool = False                  # True: there is softclip prior

    # start part
    for i in range(len_cig_scc):
        if input_cigar_scc[i].isdigit():
            cig_num_scc += input_cigar_scc[i]
        else:
            if input_cigar_scc[i] == 'S':
                soft_clip_scc += int(cig_num_scc)
                #soft_clip_bool = True
                cig_num_scc = ''
            elif input_cigar_scc[i] == 'M':
                #if soft_clip_bool:
                soft_clip_scc += int(cig_num_scc)
                length_before_gap_scc += int(cig_num_scc)
                cig_num_scc = ''
            elif input_cigar_scc[i] == 'D':
                del_scc = int(cig_num_scc)
                cig_num_scc = ''
                # if the length before gap fit the end signal condition
                if length_before_gap_scc <= end_length_scc:
                    # start sites, move forward
                    start_sites_dif_scc += length_before_gap_scc
                    start_sites_dif_scc += del_scc
                    # change cigar
                    temp_cigar_scc = str(soft_clip_scc) + "S"
                    for j in range(i + 1, len_cig_scc):
                        temp_cigar_scc += input_cigar_scc[j]

                    # softclip fit the condition, keep
                    if soft_clip_scc <= max_soft_clip_scc:
                        keep_bool = True
                    else:
                        keep_bool = False
                    # not end signal
                else:
                    temp_cigar_scc = input_cigar_scc
                # softclip do not fit the condition, drop
                break
            elif input_cigar_scc[i] == 'I':
                soft_clip_scc += int(cig_num_scc)
                length_before_gap_scc += int(cig_num_scc)
                cig_num_scc = ''
            else:
                cig_num_scc = ''

    if temp_cigar_scc == "":
        temp_cigar_scc = input_cigar_scc

    # end part
    if keep_bool:
        len_cig_scc = len(temp_cigar_scc)
        current_ess = ""                    # store current S, M, D, I

        # initialize
        soft_clip_scc = 0
        length_before_gap_scc = 0
        cig_num_scc = ''
        #soft_clip_bool = False

        for i in range(len_cig_scc):
            if temp_cigar_scc[len_cig_scc - i - 1].isdigit():
                cig_num_scc += temp_cigar_scc[len_cig_scc - i - 1]
            else:
                if current_ess == "":
                    current_ess = temp_cigar_scc[len_cig_scc - i - 1]
                else:
                    if current_ess == "S":
                        soft_clip_scc += int(cig_num_scc[::-1])
                        #soft_clip_bool = True
                        cig_num_scc = ''
                    elif current_ess == "M":
                        #if soft_clip_bool:
                        soft_clip_scc += int(cig_num_scc[::-1])
                        length_before_gap_scc += int(cig_num_scc[::-1])
                        cig_num_scc = ''
                    elif current_ess == "I":
                        soft_clip_scc += int(cig_num_scc[::-1])
                        length_before_gap_scc += int(cig_num_scc[::-1])
                        cig_num_scc = ''
                    elif current_ess == "D":

                        # if the length before gap fit the end signal condition
                        if length_before_gap_scc <= end_length_scc:
                            # change cigar
                            output_cigar_scc = "S" + str(soft_clip_scc)[::-1]
                            for j in range(i, len_cig_scc):
                                output_cigar_scc += temp_cigar_scc[len_cig_scc - j - 1]
                            output_cigar_scc = output_cigar_scc[::-1]

                            # softclip fit the condition, keep
                            if soft_clip_scc <= max_soft_clip_scc:
                                keep_bool = True
                            else:
                                keep_bool = False
                        # not end signal
                        else:
                            output_cigar_scc = temp_cigar_scc
                        break

                    current_ess = temp_cigar_scc[len_cig_scc - i - 1]

    if output_cigar_scc == "":
        output_cigar_scc = temp_cigar_scc

    return keep_bool, start_sites_dif_scc, output_cigar_scc


# define a funtion to do multiprocessing
def res_multi_proc(input_bam_rmp, output_bam_rmp, max_soft_clip_rmp, end_length_rmp, temp_dic_rmp, core_num_rmp):

    mk_temp_result = realignment_forward.mk_temp_file(input_bam_rmp, temp_dic_rmp, 50000)

    temp_file_rmp = mk_temp_result[0]

    temp_out_file_rmp = []

    # create name list of temp output bam
    for ii in range(len(temp_file_rmp)):
        temp_out_file_rmp.append(temp_file_rmp[ii][:-4] + '.remove_end_sig.bam')

    po = Pool(core_num_rmp)
    for ii in range(len(temp_file_rmp)):
        po.apply_async(solve_end_sig, (temp_file_rmp[ii], temp_out_file_rmp[ii], max_soft_clip_rmp, end_length_rmp))

    po.close()
    po.join()

    # output file
    bam_input_all_outer_rmp = pysam.AlignmentFile(input_bam_rmp, "rb")
    bam_header_rmp = bam_input_all_outer_rmp.header

    # to create a temp file to store unsorted bam
    ran_str_rmp = ''.join(random.sample(string.ascii_letters + string.digits, 8))
    output_bam_file_temp_all_rmp = temp_dic_rmp + '/' + 'output_temp_' + ran_str_rmp + '.bam'

    bam_output_all_rmp = pysam.AlignmentFile(output_bam_file_temp_all_rmp, "wb", header = bam_header_rmp)

    bam_input_all_outer_rmp.close()

    for ii in range(len(temp_out_file_rmp)):
        bam_output_all_rmp = realignment_forward.bam_merge(temp_out_file_rmp[ii], bam_output_all_rmp)

    bam_output_all_rmp.close()

    realignment_forward.sort_bam_index(output_bam_file_temp_all_rmp, output_bam_rmp, core_num_rmp)

    shutil.rmtree(temp_dic_rmp)


##############################################################
# Main Program
##############################################################

if __name__ == "__main__":

    # input parameters
    core_num = 1                                # number of cores, default = 1
    max_soft_clip = 3                           # max length of output reads' softclip, default = 3
    end_length = 1                              # length of deletion to the terminal, default = 1

    # files
    input_bam = ""  # input bam
    output_bam = ""  # output bam

    temp_file_dic = 'tmp'                       # path of temp file directory

    for i in range(len(sys.argv)):
        input_para = sys.argv[i]
        if input_para == '-i':
            input_bam = sys.argv[i + 1]
        elif input_para == '-o':
            output_bam = sys.argv[i + 1]
        elif input_para == '-sc':
            max_soft_clip = int(sys.argv[i + 1])
        elif input_para == '-el':
            end_length = int(sys.argv[i + 1])
        elif input_para == '-t':
            core_num = int(sys.argv[i + 1])
        elif input_para == '-T':
            temp_file_dic = sys.argv[i + 1]
        elif input_para == '-h' or input_para == '--help':
            print('-----------------------------------------------------Welcome to use PRAISE realignment remove end signals-----------------------------------------------------')
            print('Usage:')
            print('    python remove_end_signal.py [Options]* {-i <bam_in>} [-o <bam_out>]')
            print('')
            print('    <bam_in>        Input bam file, need to be sorted by name')
            print('    <bam_out>       Output sorted bam file with its index (.bai file)')
            print('')
            print('[Options]')
            print('')
            print('  Input:')
            print('    -t [INT]        cores number, default =', core_num)
            print('    -T [Directory]  directory of temporary file, default =', temp_file_dic)
            print('')
            print('  Filter:')
            print('    -sc [INT]       max softclip length, filter if softclip length is longer than [INT], default =', max_soft_clip)
            print('    -el [INT]       max number of nucleotides between deletion signal and terminal, filter if number is bigger than [INT], default =', end_length)
            print('')
            print('  Other:')
            print('    -h / --help     help')
            exit(0)


    # make a dictionary to store temp file
    ran_str = ''.join(random.sample(string.ascii_letters + string.digits, 8))
    temp_file_dic = temp_file_dic + '_' + ran_str + '_' + realignment_forward.last_file_name(input_bam)

    os.mkdir(temp_file_dic)

    res_multi_proc(input_bam, output_bam, max_soft_clip, end_length, temp_file_dic, core_num)




