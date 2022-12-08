##########################################################################################################
# to remove multiple mapping in realignment results, using "SC" tag
# Author: Zhe Jiang
# created on 2021/12/16
##########################################################################################################


##############################################################
# import packages
##############################################################

import os
import random
import string
import sys
import shutil
import pysam
import realignment_forward
from multiprocessing import Pool, Manager


##############################################################
# define functions
##############################################################

# define a function to read results into a dictionary
def load_bam(input_bam_file_load):

    input_bam_load = pysam.AlignmentFile(input_bam_file_load, "rb")   # input bam
    bam_dic_load = {}                                            # define a dictionary to store result
    # {"read name" : [[score1, score2, ...], [number of best score]]}

    all_reads_load = input_bam_load.fetch(until_eof = True)

    # add score into bam_dic
    for read_load in all_reads_load:
        try:                                                # if read name exists, add directly
            bam_dic_load[read_load.query_name][0].append(read_load.get_tag("SC"))
        except:                                             # if read name does not exist, create a new key
            bam_dic_load[read_load.query_name] = [[read_load.get_tag("SC")],[]]

    # stat the number of best score and add into bam_dic
    for key_load in bam_dic_load:
        score_list_rm = bam_dic_load[key_load][0]
        score_list_rm.sort(reverse = True)
        num_be_sc_rm = 0                               # store the number of best score
        best_score_rm = score_list_rm[0]                  # store the best score

        for i in range(len(score_list_rm)):
            if score_list_rm[i] == best_score_rm:
                num_be_sc_rm += 1
            else:
                break

        bam_dic_load[key_load][1].append(num_be_sc_rm)

    input_bam_load.close()

    return bam_dic_load


# define a function to output result after removing multiple mapping, keep same best score
def rm_mp_keep(input_bam_file_before_rm, output_bam_file_rm):

    input_bam_file_rm = input_bam_file_before_rm[:-4] + 'score_change.bam'

    bam_score_cha(input_bam_file_before_rm, input_bam_file_rm)

    input_bam_rm = pysam.AlignmentFile(input_bam_file_rm, "rb")
    bam_header_rm = input_bam_rm.header

    output_read_list_rm = []

    bam_dic_rm = load_bam(input_bam_file_rm)

    all_reads_rm = input_bam_rm.fetch(until_eof = True)

    # only add non-multiple mapping results in to output bam file
    for read_rm in all_reads_rm:
        if read_rm.get_tag("NH") == 1:                                                         # no multiple mapping
            output_read_list_rm.append(read_rm)
        else:                                                                               # with multiple mapping
            best_score = bam_dic_rm[read_rm.query_name][0][0]                                     # store best score
            if read_rm.get_tag("SC") == best_score and bam_dic_rm[read_rm.query_name][1][0] == 1:    # the only best alignment than do next
                out_read_rm = pysam.AlignedSegment()
                out_read_rm.query_name = read_rm.query_name                           # name
                out_read_rm.flag = flag_change(read_rm.flag)                          # flag
                out_read_rm.reference_id = read_rm.reference_id                       # reference name
                out_read_rm.reference_start = read_rm.reference_start                 # reference position
                out_read_rm.mapping_quality = 60                                      # mapping quality
                out_read_rm.cigar = read_rm.cigar                                     # cigar
                out_read_rm.next_reference_id = read_rm.next_reference_id             # reference name of paired read (mate)
                out_read_rm.next_reference_start = read_rm.next_reference_start       # reference position of paired read (mate)
                out_read_rm.template_length = read_rm.template_length                 # length between read and mate
                out_read_rm.query_sequence = read_rm.query_sequence                   # squence
                out_read_rm.query_qualities = read_rm.query_qualities                 # quality
                out_read_rm.tags = (("NH", 1),                                        # set "NH" tag to store mapping times
                                 ("SC", read_rm.get_tag("SC")))                       # create a tag "SC" to store best result

                output_read_list_rm.append(out_read_rm)
            elif read_rm.get_tag("SC") == best_score and bam_dic_rm[read_rm.query_name][1][0] > 1:
                out_read_rm = pysam.AlignedSegment()
                out_read_rm.query_name = read_rm.query_name                           # name
                out_read_rm.flag = read_rm.flag                                       # flag
                out_read_rm.reference_id = read_rm.reference_id                       # reference name
                out_read_rm.reference_start = read_rm.reference_start                 # reference position
                out_read_rm.mapping_quality = 1                                       # mapping quality
                out_read_rm.cigar = read_rm.cigar                                     # cigar
                out_read_rm.next_reference_id = read_rm.next_reference_id             # reference name of paired read (mate)
                out_read_rm.next_reference_start = read_rm.next_reference_start       # reference position of paired read (mate)
                out_read_rm.template_length = read_rm.template_length                 # length between read and mate
                out_read_rm.query_sequence = read_rm.query_sequence                   # squence
                out_read_rm.query_qualities = read_rm.query_qualities                 # quality
                out_read_rm.tags = (("NH", bam_dic_rm[read_rm.query_name][1][0]),     # set "NH" tag to store mapping times
                                 ("SC", read_rm.get_tag("SC")))                       # create a tag "SC" to store best result

                output_read_list_rm.append(out_read_rm)

    realignment_ver11.write_reads(output_read_list_rm, output_bam_file_rm, bam_header_rm)

    input_bam_rm.close()



# define a function to output result after removing multiple mapping, remove if there are scores equal
def rm_mp_remove(input_bam_file_before_rm, output_bam_file_rm):

    input_bam_file_rm = input_bam_file_before_rm[:-4] + 'score_change.bam'

    bam_score_cha(input_bam_file_before_rm, input_bam_file_rm)

    input_bam_rm = pysam.AlignmentFile(input_bam_file_rm, "rb")
    bam_header_rm = input_bam_rm.header

    output_read_list_rm = []

    bam_dic_rm = load_bam(input_bam_file_rm)

    all_reads_rm = input_bam_rm.fetch(until_eof = True)

    # only add non-multiple mapping results in to output bam file
    for read_rm in all_reads_rm:
        if read_rm.get_tag("NH") == 1:                                                         # no multiple mapping
            output_read_list_rm.append(read_rm)
        else:                                                                               # with multiple mapping
            best_score_rm = bam_dic_rm[read_rm.query_name][0][0]                                     # store best score
            if read_rm.get_tag("SC") == best_score_rm and bam_dic_rm[read_rm.query_name][1][0] == 1:    # the only best alignment than do next
                out_read_rm = pysam.AlignedSegment()
                out_read_rm.query_name = read_rm.query_name                           # name
                out_read_rm.flag = flag_change(read_rm.flag)                          # flag
                out_read_rm.reference_id = read_rm.reference_id                       # reference name
                out_read_rm.reference_start = read_rm.reference_start                 # reference position
                out_read_rm.mapping_quality = 60                                   # mapping quality
                out_read_rm.cigar = read_rm.cigar                                     # cigar
                out_read_rm.next_reference_id = read_rm.next_reference_id             # reference name of paired read (mate)
                out_read_rm.next_reference_start = read_rm.next_reference_start       # reference position of paired read (mate)
                out_read_rm.template_length = read_rm.template_length                 # length between read and mate
                out_read_rm.query_sequence = read_rm.query_sequence                   # squence
                out_read_rm.query_qualities = read_rm.query_qualities                 # quality
                out_read_rm.tags = (("NH", 1),                                     # set "NH" tag to store mapping times
                                 ("SC", read_rm.get_tag("SC")))                    # create a tag "SC" to store best result

                output_read_list_rm.append(out_read_rm)

    realignment_forward.write_reads(output_read_list_rm, output_bam_file_rm, bam_header_rm)

    input_bam_rm.close()


# define a function to change flag
def flag_change(flag_input):

    if flag_input > 255:
        return (flag_input - 256)
    else:
        return (flag_input)


# define a function to sort bam and create index
def sort_bam_index_rm(input_bam_unsort, output_bam_sorted, cores_num_sort):

    core_str_sort = str(cores_num_sort)

    pysam.sort("-@", core_str_sort, "-m", "2G", "-O", "BAM", "-o", output_bam_sorted, input_bam_unsort)

    output_bam_bai = output_bam_sorted + '.bai'
    pysam.index(output_bam_sorted, output_bam_bai)


# define a function to change best score by cigar
def score_change(score_cha, cigar_cha):

    # if there is one deletion in cigar, score - 10
    del_num_cha = 0
    cig_num_cha = ''
    len_cig_cha = len(cigar_cha)

    for i in range(len_cig_cha):
        if cigar_cha[i].isdigit():
            cig_num_cha += cigar_cha[i]
        else:
            if cigar_cha[i] == 'D':
                del_num_cha += int(cig_num_cha)
                cig_num_cha = ''
            else:
                cig_num_cha= ''

    return (score_cha - 20 * del_num_cha)


# define a function to change best score in the bam
def bam_score_cha(input_bam_file_cha, output_bam_file_cha):

    input_bam_cha = pysam.AlignmentFile(input_bam_file_cha, 'rb')
    bam_header_cha = input_bam_cha.header

    output_read_ls_cha = []

    all_reads_cha = input_bam_cha.fetch(until_eof = True)

    for read_cha in all_reads_cha:

        out_read_cha = pysam.AlignedSegment()
        out_read_cha.query_name = read_cha.query_name  # name
        out_read_cha.flag = read_cha.flag  # flag
        out_read_cha.reference_id = read_cha.reference_id  # reference name
        out_read_cha.reference_start = read_cha.reference_start  # reference position
        out_read_cha.mapping_quality = read_cha.mapping_quality  # mapping quality
        out_read_cha.cigar = read_cha.cigar  # cigar
        out_read_cha.next_reference_id = read_cha.next_reference_id  # reference name of paired read (mate)
        out_read_cha.next_reference_start = read_cha.next_reference_start  # reference position of paired read (mate)
        out_read_cha.template_length = read_cha.template_length  # length between read and mate
        out_read_cha.query_sequence = read_cha.query_sequence  # squence
        out_read_cha.query_qualities = read_cha.query_qualities  # quality
        out_read_cha.tags = (("NH", 1),  # set "NH" tag to store mapping times
                            ("SC", score_change(read_cha.get_tag("SC"), read_cha.cigarstring)))  # create a tag "SC" to store best result

        output_read_ls_cha.append(out_read_cha)

    realignment_forward.write_reads(output_read_ls_cha, output_bam_file_cha, bam_header_cha)


# define a function to do all process
def remove_mp(rm_input_bam_file, rm_output_bam_file, rm_mode, cores_num_rm, tmp_file_dic_rm):

    ran_str = ''.join(random.sample(string.ascii_letters + string.digits, 8))
    rm_bam_output_tmp_dic = tmp_file_dic_rm + '_' + ran_str + "_" + realignment_forward.last_file_name(rm_input_bam_file) + 'rm_multi'

    os.mkdir(rm_bam_output_tmp_dic)

    # create tmp file
    mk_tmp_result = realignment_forward.mk_temp_file(rm_input_bam_file, rm_bam_output_tmp_dic, 1000000)
    tmp_file_ls_rm = mk_tmp_result[0]

    # create output tmp file
    tmp_out_file_ls_rm = []

    for i in range(len(tmp_file_ls_rm)):

        tmp_out_file_ls_rm.append(tmp_file_ls_rm[i][:-4] + '.rm_out.bam')

    # choose mode
    if rm_mode == "keep":
        rm_mp = rm_mp_keep
    elif rm_mode == "remove":
        rm_mp = rm_mp_remove
    else:
        sys.stderr.write("Fatal error: The remove mode must be keep or remove!")
        exit(1)

    # run multiple processing
    po = Pool(cores_num_rm)
    for ii in range(len(tmp_file_ls_rm)):
        po.apply_async(rm_mp, (tmp_file_ls_rm[ii], tmp_out_file_ls_rm[ii]))

    po.close()
    po.join()

    # output file
    bam_input_all_nsort_rm = pysam.AlignmentFile(rm_input_bam_file, "rb")
    bam_header_rm = bam_input_all_nsort_rm.header

    # to create a temp file to store unsorted bam
    ran_str = ''.join(random.sample(string.ascii_letters + string.digits, 8))
    output_bam_file_temp_all_rm = rm_bam_output_tmp_dic + '/' + 'output_temp_' + ran_str + '.bam'

    bam_output_all_rm = pysam.AlignmentFile(output_bam_file_temp_all_rm, "wb", header=bam_header_rm)

    bam_input_all_nsort_rm.close()

    for ii in range(len(tmp_out_file_ls_rm)):
        bam_output_all_rm = realignment_forward.bam_merge(tmp_out_file_ls_rm[ii], bam_output_all_rm)

    bam_output_all_rm.close()

    # sort by coordinate and create index
    sort_bam_index_rm(output_bam_file_temp_all_rm, rm_output_bam_file, cores_num_rm)

    # remove all tmp files
    shutil.rmtree(rm_bam_output_tmp_dic)




###################################################################
# Main program
###################################################################


if __name__ == "__main__":

    ################################################################
    # Input
    ################################################################


    bam_input_file = ""
    bam_output_file = ""

    cores_num = 1
    tmp_dic = 'tmp'

    remove_mode = "keep"


    for i in range(len(sys.argv)):
        input_para = sys.argv[i]
        if input_para == '-i':
            bam_input_file = sys.argv[i + 1]
        elif input_para == '-o':
            bam_output_file = sys.argv[i + 1]
        elif input_para == '-rm':
            remove_mode = sys.argv[i + 1]
        elif input_para == '-h' or input_para == '--help':
            print('-----------------------------------------------------Welcome to use PRAISE realignment multiple mapping removing-----------------------------------------------------')
            print('Usage:')
            print('    python remove_multi-mapping.py [Options]* {-i <bam_in>} [-o <bam_out>]')
            print('')
            print('    <bam_in>        Input bam file, need to be sorted and have index (.bai file)')
            print('    <bam_out>       Output sorted bam file, along with index (.bai), removing multiple mapping result')
            print('')
            print('[Options]')
            print('')
            print('  Input:')
            print('    -t [INT]        cores number, default =', cores_num)
            print('    -T [Directory]  directory of temporary file, default =', tmp_dic)
            print('')
            print('  Remove multiple mapping mode:')
            print('    -rm keep        keep all mapping results if their scores all equal to the best score, default mode')
            print("    -rm remove      remove all mappling results if there are more than one read's scores equal to the best score")
            print('')
            print('  Other:')
            print('    -h / --help     help')
            exit(0)


    ###############################################################
    # Processing
    ###############################################################

    remove_mp(bam_input_file, bam_output_file, remove_mode, cores_num, tmp_dic)







