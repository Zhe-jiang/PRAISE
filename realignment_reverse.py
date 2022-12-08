##########################################################################################################
# realignment for PRAISE no-spliced alignment result, reverse sequence
# Author: Zhe Jiang
# created on 2021/11/21
##########################################################################################################


##############################################################
# import packages
##############################################################

from biopython.Bio import pairwise2
import substitution_matrices
import pysam
import sys
import os
from multiprocessing import Pool, Manager
import shutil
import random
import string
import time
import remove_multi_mapping


##############################################################
# define functions
##############################################################

# define gap function in reference (insert in reads)
def gap_function_ref(open_site, open_length, site_base):
    if open_site == 0:
        return (-1 * open_length)
    else:
        if open_length == 0:
            return 0
        else:
            return (-9.3 - (open_length - 1) * 3.4)


# define gap function in reads
def gap_function_reads(open_site, open_length, site_base):
    if open_site == 0:
        return (-1 * open_length)
    else:
        if open_length == 0:
            return 0
        else:
            if site_base == 'A':
                return (-3.1 - (open_length - 1) * 2.9)
            else:
                return (-5.2 - (open_length - 1) * 4.3)


# define re-score function and output the cigar & start site
def ali_score(alis, start_site_ini):

    start_bool = True   # record whether it is start gap in seq_read
    best_score = -10000000      # store current best score
    best_cigar = ''     # store cigar of best score
    best_ali = alis[0]  # initiate
    best_ali_start_site = start_site_ini    # initiate


    for ali in alis:
        ali_ref = ali.seqA
        ali_read = ali.seqB

        # to reverse and make sure in re-score system, it is 5' to 3'
        ali_ref = ali_ref[::-1]
        ali_read = ali_read[::-1]

        current_cig = ''  # store current cigar 'M', 'I', 'D', 'S'
        current_len = 0  # store length in current cigar, int
        current_cigar_string = ''   # record cigar of this alignment
        current_score = 0   # store current score

        start_bool = True               # initiate start_bool
        start_site = start_site_ini     # initiate start_site

        for i in range(len(ali_ref)):

            # start gap in ali_read: just move start_site forward
            if ali_ref[i] != '-' and ali_read[i] == '-' and start_bool:
                start_site += 1

            # start gap in ali_gap: record as softclip
            elif ali_ref[i] == '-' and ali_read[i] != '-' and start_bool:
                if current_cig == 'S':      # judge if it is extension or start
                    current_len += 1
                else:
                    current_cig = 'S'
                    current_len = 1
                    # because it is start gap, so there is no need to record previous current_len and current_cigar into current_cigar_string

            # match & mismatch
            elif ali_read[i] != '-' and ali_ref[i] != '-':
                # match
                if ali_read[i] == ali_ref[i]:
                    current_score += 5  # score
                    # cigar
                    if current_cig == 'M':  # judge if it is extension or start
                        current_len += 1
                    elif current_cig == '':  # start of alignment
                        current_cig = 'M'
                        current_len = 1
                    else:  # not start of alignment
                        current_cigar_string += str(current_len)
                        current_cigar_string += current_cig
                        current_cig = 'M'
                        current_len = 1
                    start_bool = False  # it not start ever
                # mismatch
                else:
                    if start_bool or judge_end_mismatch(ali_ref, ali_read, i):  # mismatch at start or end: deal as softclip
                        if current_cig == 'S':      # judge if it is extension or start
                            current_len += 1
                            if start_bool:          # if softclip caused by start mismatch, move start site forward
                                start_site += 1
                        elif current_cig == '':     # start of alignment
                            current_cig = 'S'
                            current_len = 1
                            if start_bool:          # if softclip caused by start mismatch, move start site forward
                                start_site += 1
                        else:   # not start of alignment
                            current_cigar_string += str(current_len)
                            current_cigar_string += current_cig
                            current_cig = 'S'
                            current_len = 1
                    else:   # mismatch in the middle
                        start_bool = False  # it not start ever
                        current_score -= 5
                        # cigar
                        if current_cig == 'M':  # judge if it is extension or start
                            current_len += 1
                        else:
                            current_cigar_string += str(current_len)
                            current_cigar_string += current_cig
                            current_cig = 'M'
                            current_len = 1

            # end gap in ali_read: just stop
            elif ali_ref[i] != '-' and ali_read[i] == '-' and judge_end_gap(ali_read, i):
                break

            # end gap in ali_ref: deal as softclip
            elif ali_ref[i] == '-' and ali_read[i] != '-' and judge_end_gap(ali_ref, i):
                start_bool = False  # it not start ever
                if current_cig == 'S':  # judge if it is extension or start
                    current_len += 1
                elif current_cig == '':  # start of alignment
                    current_cig = 'S'
                    current_len = 1
                else:  # not start of alignment
                    current_cigar_string += str(current_len)
                    current_cigar_string += current_cig
                    current_cig = 'S'
                    current_len = 1

            # gap in ali_read: deletion
            elif ali_ref[i] != '-' and ali_read[i] == '-' and (not start_bool):
                start_bool = False  # it not start ever
                if current_cig == 'D':  # judge if it is extension or start
                    current_len += 1
                    current_score -= 3  # extension score -3
                else:  # start
                    current_cigar_string += str(current_len)
                    current_cigar_string += current_cig
                    current_cig = 'D'
                    current_len = 1
                    if ali_ref[i] == 'A':   # T del
                        current_score += 5  # start T del score 5
                    else:   # not T del
                        current_score -= 5  # start not T del score -5

            # gap in ali_ref: insert
            elif ali_ref[i] == '-' and ali_read[i] != '-' and (not start_bool):
                start_bool = False  # it not start ever
                if current_cig == 'I':  # judge if it is extension or start
                    current_len += 1
                    current_score -= 3  # extension score -3
                else:  # start
                    current_cigar_string += str(current_len)
                    current_cigar_string += current_cig
                    current_cig = 'I'
                    current_len = 1
                    current_score -= 5    # start score -5


        # add final cigar
        current_cigar_string += str(current_len)
        current_cigar_string += current_cig

        if current_score > best_score:
            best_score = current_score
            best_cigar = current_cigar_string
            best_ali = ali
            best_ali_start_site = start_site

    start_site_dif_sco = start_site_ini - best_ali_start_site

    return best_ali, best_score, best_cigar, best_ali_start_site, start_site_dif_sco


# define a function that judge whether it is end gap
def judge_end_gap(sequence, cur_site):

    jud = True

    for j in range(cur_site, len(sequence)):
        if sequence[j] != '-':
            jud = False
            break

    return jud


# define a function that judge whether it is end mismatch
def judge_end_mismatch(sequence_ref, sequence_read, cur_site):

    jud = True

    for j in range(cur_site, len(sequence_ref)):
        if sequence_read[j] == sequence_ref[j] or (sequence_ref[j] == 'G' and sequence_read[j] == 'A'):
            jud = False
            break

    return jud


# define a function change string cigar into tuples cigar
def cigar_change(cig):

    cig_list = []
    num = ''
    len_cig = len(cig)

    for i in range(len_cig):
        if cig[i].isdigit():
            num += cig[i]
        elif cig[i] == 'M':
            cig_num = int(num)      # store the number of this cigar
            num = ''
            cig_op = 0              # store of operation of this cigar: M means 0
            cig_list.append((cig_op, cig_num))
        elif cig[i] == 'I':
            cig_num = int(num)  # store the number of this cigar
            num = ''
            cig_op = 1  # store of operation of this cigar: I means 1
            cig_list.append((cig_op, cig_num))
        elif cig[i] == 'D':
            cig_num = int(num)  # store the number of this cigar
            num = ''
            cig_op = 2  # store of operation of this cigar: D means 1
            cig_list.append((cig_op, cig_num))
        elif cig[i] == 'N':
            cig_num = int(num)  # store the number of this cigar
            num = ''
            cig_op = 3  # store of operation of this cigar: N means 3
            cig_list.append((cig_op, cig_num))
        elif cig[i] == 'S':
            cig_num = int(num)  # store the number of this cigar
            num = ''
            cig_op = 4  # store of operation of this cigar: S means 4
            cig_list.append((cig_op, cig_num))
        elif cig[i] == 'H':
            cig_num = int(num)  # store the number of this cigar
            num = ''
            cig_op = 5  # store of operation of this cigar: H means 5
            cig_list.append((cig_op, cig_num))
        elif cig[i] == 'P':
            cig_num = int(num)  # store the number of this cigar
            num = ''
            cig_op = 6  # store of operation of this cigar: P means 6
            cig_list.append((cig_op, cig_num))
        elif cig[i] == '=':
            cig_num = int(num)  # store the number of this cigar
            num = ''
            cig_op = 7  # store of operation of this cigar: = means 7
            cig_list.append((cig_op, cig_num))
        elif cig[i] == 'X':
            cig_num = int(num)  # store the number of this cigar
            num = ''
            cig_op = 8  # store of operation of this cigar: I means 8
            cig_list.append((cig_op, cig_num))
        elif cig[i] == 'B':
            cig_num = int(num)  # store the number of this cigar
            num = ''
            cig_op = 9  # store of operation of this cigar: I means 9
            cig_list.append((cig_op, cig_num))

    return tuple(cig_list)


# define a function to load reference into a dictionary
def load_ref(ref_path):

    ref_dic = {}    # dictionary to store reference {ref_name : ref_seq}
    ref_name = ''
    ref_seq = ''

    for line in open(ref_path, 'r'):
        if line[0] == '>':  # name line
            if ref_name != '':
                ref_dic[ref_name] = ref_seq     # not start line, add last result into dic
            whole_line = line[1:-1]
            ref_name = whole_line.split(" ")[0]
            ref_seq = ''
        else:               # sequence line
            ref_seq += line[:-1]

    ref_dic[ref_name] = ref_seq                 # add last record

    return ref_dic


# define a function to load seq_ref
def load_seq_ref(ref_seq, start_posi, cig, len_read):

    len_ref_seq = len(ref_seq)
    len_cig = len(cig)
    cig_num = ''

    start_extend = 3        # even if there is no softclip, extend forward 3 bp at the start
    end_extend = 3          # extend forward 3 bp at the end

    for i in range(len_cig):                # add len of start softclip to start_extend
        if cig[i].isdigit():
            cig_num += cig[i]
        else:
            if cig[i] == 'S':
                start_extend += int(cig_num)
                break
            else:
                break


    start_posi_out = max(0, start_posi - start_extend)
    end_posi_out = min(len_ref_seq, start_posi + len_read + end_extend)

    return start_posi_out, ref_seq[start_posi_out : end_posi_out]


# define a function to do filter
def align_filter(cig, score, softclip_len, min_score):
    # cig is the cigar of alignment
    # score is the total score of whole alignment
    # softclip_len is the longest soft clip that will be reversed : int
    # min_score is the minimiun score per base that will be reversed : float

    if not judge_softclip(cig, softclip_len):
        return False                                        # return false if softclip is too long
    else:
        if score_calc(cig, score) < min_score:
            return False                                    # return false if score per base is low
        else:
            return True


# define a funtion to judge if the length of softclip is qualified
def judge_softclip(cig, len_softclip):
    # if softclip in cig > len_softclip, return False; else return True

    result = True
    len_cig = len(cig)
    cig_num = ''

    # softclip in the start
    for i in range(len_cig):
        if cig[i].isdigit():
            cig_num += cig[i]
        else:
            if cig[i] == 'S' and int(cig_num) > len_softclip:
                result = False
            break

    cig_num = ''
    # softclip in the end
    if cig[-1] == 'S':
        for i in range(len_cig - 1):
            if cig[len_cig - i - 2].isdigit():
                cig_num += cig[len_cig - i - 2]
            else:
                if int(cig_num[::-1]) > len_softclip:
                    result = False
                break

    return result


# define a function to calculate score per base
def score_calc(cig, score):

    base_num = 0
    cig_num = ''
    len_cig = len(cig)

    for i in range(len_cig):
        if cig[i].isdigit():
            cig_num += cig[i]
        else:
            if cig[i] == 'S':
                cig_num = ''
            elif cig[i] == 'M':
                base_num += int(cig_num)
                cig_num = ''
            elif cig[i] == 'D':
                base_num += int(cig_num)
                cig_num = ''
            elif cig[i] == 'I':
                base_num += int(cig_num)
                cig_num = ''
            else:
                cig_num = ''

    return (score/base_num)


# define a function to split file to make temp file in order to do multiprocessing
def mk_temp_file(input_file_name, temp_dic, reads_num_mk):

    temp_file_list = []     # store the name of temp_file

    last_input_file_name = last_file_name(input_file_name)

    if temp_dic[-1] == '/':
        temp_dic = temp_dic[:-1]

    file_index = 1
    temp_file_name = temp_dic + '/' + last_input_file_name[:-4] + '.temp.' + str(file_index) + '.bam'
    temp_file_list.append(temp_file_name)

    # store whole input bam
    save_mk = pysam.set_verbosity(0)        # help suppressing HTSlib's message
    bam_input_all = pysam.AlignmentFile(input_file_name, "rb")
    pysam.set_verbosity(save_mk)
    temp_bam_header = bam_input_all.header  # record the header

    bam_temp = pysam.AlignmentFile(temp_file_name, "wb", header = temp_bam_header)

    all_bam_reads = bam_input_all.fetch(until_eof = True)

    last_read_name = ''                 # to store last read name, make sure paired read will be in the same file

    count = 1

    # read by read, write into a new temp file, 50000 reads per file
    for read in all_bam_reads:
        if count > reads_num_mk and read.query_name != last_read_name:
            # close this temp file
            count = 1
            bam_temp.close()
            # create the next temp file
            file_index += 1
            temp_file_name = temp_dic + '/' + last_input_file_name[:-4] + '.temp.' + str(file_index) + '.bam'
            temp_file_list.append(temp_file_name)
            bam_temp = pysam.AlignmentFile(temp_file_name, "wb", header = temp_bam_header)
            bam_temp.write(read)
        else:
            count += 1
            bam_temp.write(read)
            last_read_name = read.query_name

    # close the last bam_temp
    bam_temp.close()


    return temp_file_list, file_index


# define a function that run realignment
def realign(input_temp_bam, ref_dic, max_softclip_temp, min_score_per_base_temp, max_softclip_input_temp, matrix_temp, q_num, q_bool, tot_tasks, lock, rea_mode_rea):
    # input bam file
    save_rea = pysam.set_verbosity(0)         # help suppressingHTTSlib's warning message
    bam_input_temp = pysam.AlignmentFile(input_temp_bam, "rb")
    pysam.set_verbosity(save_rea)

    # load all reads in allreads
    allreads = bam_input_temp.fetch(until_eof = True)

    # create a bam file to store reads after realignment
    bam_header_temp_out = bam_input_temp.header  # record the header

    output_bam_file_temp = input_temp_bam[:-4] + '.realign.bam'
    output_bam_filter_file_temp = input_temp_bam[:-4] + '.realign.filter.bam'


    reads_output_list = []
    reads_filter_list = []

    last_reads_name_rea = ""            # store last reads name
    output_dic_rea = {}                 # store realign output of one reads as dictionary
        # key: seq_ref_ini
        # output_dic_rea[key] = []
        # output_dic_rea[key][0] : best score
        # output_dic_rea[key][1] : cigar of best result
        # output_dic_rea[key][2] : difference of start site initial and best aligment start site

    output_fil_dic_rea = {}  # store filter output of one reads as dictionary
        # key: seq_ref_ini
        # output_dic_rea[key] = []
        # output_dic_rea[key][0] : best score
        # output_dic_rea[key][1] : cigar of best result
        # output_dic_rea[key][2] : difference of start site initial and best aligment start site

    # use a boolean to judge the mode
    if rea_mode_rea == 'fast':
        fast_bool = True
    else:
        fast_bool = False


    for read in allreads:

        # filter if softclip is longer than max_softclip_input or unmapped or mapping to reverse
        if (not judge_softclip(read.cigarstring, max_softclip_input_temp)) or read.is_unmapped:
            reads_filter_list.append(read)
            continue

        else:

            # judge whether it is fast mode or not
            if fast_bool:
                if cigar_judge(read.cigarstring):

                    # get ref and read
                    seq_read_ini = read.query_sequence
                    seq_read_ini = seq_read_ini.upper()

                    seq_ref_ini = read.get_reference_sequence()
                    seq_ref_ini = seq_ref_ini.upper()

                    # calculate the score
                    score_fas_rea = fast_score(seq_ref_ini, seq_read_ini)
                        # score_fas_rea[0]: total score
                        # score_fas_rea[1]: score per base

                    out_read = pysam.AlignedSegment()
                    out_read.query_name = read.query_name  # name
                    out_read.flag = read.flag  # flag
                    out_read.reference_id = read.reference_id  # reference name
                    out_read.reference_start = read.reference_start  # reference position
                    out_read.mapping_quality = 60  # mapping quality
                    out_read.cigar = read.cigar  # cigar
                    out_read.next_reference_id = read.next_reference_id  # reference name of paired read (mate)
                    out_read.next_reference_start = read.next_reference_start  # reference position of paired read (mate)
                    out_read.template_length = read.template_length  # length between read and mate
                    out_read.query_sequence = read.query_sequence  # squence
                    out_read.query_qualities = read.query_qualities  # quality
                    out_read.tags = (("NH", read.get_tag("NH")),  # set "NH" tag to store mapping times
                                     ("SC", score_fas_rea[0]))  # create a tag "SC" to store best result

                    if score_fas_rea[1] < min_score_per_base_temp:
                        reads_filter_list.append(out_read)
                    else:
                        reads_output_list.append(out_read)

                    continue


            read_name_rea = read.query_name

            # if there is a new read, initialize
            if last_reads_name_rea != read_name_rea:
                output_dic_rea = {}
                output_fil_dic_rea = {}

            last_reads_name_rea = read_name_rea

            # load reference and read sequence
            seq_read_ini = read.query_sequence  # load read sequence

            reference_full_seq = ref_dic[read.reference_name]  # load ref sequence
            result_load_seq_ref = load_seq_ref(reference_full_seq,
                                               read.reference_start,
                                               read.cigarstring,
                                               len(seq_read_ini)
                                               )

            seq_ref_ini = result_load_seq_ref[1]  # reference's sequence
            seq_ref_ini = seq_ref_ini.upper()


            output_bool_rea = output_dic_rea.__contains__(seq_ref_ini)              # if seq ref in output file list
            output_fil_bool_rea = output_fil_dic_rea.__contains__(seq_ref_ini)      # if seq ref in filtered file list


            # if the reference sequence is realigned
            if output_bool_rea:

                out_read = pysam.AlignedSegment()
                out_read.query_name = read_name_rea  # name
                out_read.flag = read.flag  # flag
                out_read.reference_id = read.reference_id  # reference name
                out_read.reference_start = result_load_seq_ref[0] - output_dic_rea[seq_ref_ini][2]  # reference position
                out_read.mapping_quality = 60  # mapping quality
                out_read.cigar = cigar_change(output_dic_rea[seq_ref_ini][1])  # cigar
                out_read.next_reference_id = read.next_reference_id  # reference name of paired read (mate)
                out_read.next_reference_start = read.next_reference_start  # reference position of paired read (mate)
                out_read.template_length = read.template_length  # length between read and mate
                out_read.query_sequence = read.query_sequence  # squence
                out_read.query_qualities = read.query_qualities  # quality
                out_read.tags = (("NH", read.get_tag("NH")),  # set "NH" tag to store mapping times
                                 ("SC", output_dic_rea[seq_ref_ini][0]))  # create a tag "SC" to store best result

                reads_output_list.append(out_read)

                continue
            # if the reference sequence is filtered
            elif output_fil_bool_rea:

                out_read = pysam.AlignedSegment()
                out_read.query_name = read_name_rea  # name
                out_read.flag = read.flag  # flag
                out_read.reference_id = read.reference_id  # reference name
                out_read.reference_start = result_load_seq_ref[0] - output_fil_dic_rea[seq_ref_ini][2]  # reference position
                out_read.mapping_quality = 60  # mapping quality
                out_read.cigar = cigar_change(output_fil_dic_rea[seq_ref_ini][1])  # cigar
                out_read.next_reference_id = read.next_reference_id  # reference name of paired read (mate)
                out_read.next_reference_start = read.next_reference_start  # reference position of paired read (mate)
                out_read.template_length = read.template_length  # length between read and mate
                out_read.query_sequence = read.query_sequence  # squence
                out_read.query_qualities = read.query_qualities  # quality
                out_read.tags = (("NH", read.get_tag("NH")),  # set "NH" tag to store mapping times
                                 ("SC", output_fil_dic_rea[seq_ref_ini][0]))  # create a tag "SC" to store best result

                reads_filter_list.append(out_read)

                continue


            # make sure that sequence from 3' to 5'
            seq_ref = seq_ref_ini[::-1]
            seq_read = seq_read_ini[::-1]

            ref_start_posi = result_load_seq_ref[0]  # start position in the reference

            alignments = pairwise2.align.globaldc(seq_ref,
                                                  seq_read,
                                                  matrix_temp,
                                                  gap_function_ref,
                                                  gap_function_reads
                                                  )

            if len(alignments) == 0:
                reads_filter_list.append(read)
                continue

            else:

                best_result = ali_score(alignments, ref_start_posi)
                # best_result[0]: alignment
                # best_result[1]: best score
                # best_result[2]: cigar of best result
                # best_result[3]: start position of best result
                # best_result[4]: difference of start site initial and best aligment start site

            # output the reads after realignment into bam_output
            out_read = pysam.AlignedSegment()
            out_read.query_name = read_name_rea                        # name
            out_read.flag = read.flag                                  # flag
            out_read.reference_id = read.reference_id                  # reference name
            out_read.reference_start = best_result[3]                  # reference position
            out_read.mapping_quality = 60                              # mapping quality
            out_read.cigar = cigar_change(best_result[2])              # cigar
            out_read.next_reference_id = read.next_reference_id        # reference name of paired read (mate)
            out_read.next_reference_start = read.next_reference_start  # reference position of paired read (mate)
            out_read.template_length = read.template_length            # length between read and mate
            out_read.query_sequence = read.query_sequence              # squence
            out_read.query_qualities = read.query_qualities            # quality
            out_read.tags = (("NH",read.get_tag("NH")),                # set "NH" tag to store mapping times
                             ("SC", best_result[1]))                   # create a tag "SC" to store best result

            if align_filter(best_result[2], best_result[1], max_softclip_temp,
                            min_score_per_base_temp):           # if filter = True, add to output
                reads_output_list.append(out_read)
                output_dic_rea[seq_ref_ini] = [best_result[1], best_result[2], best_result[4]]
                    # write best alignment into the dic
            else:                                               # if filter = False, add to output_filter
                reads_filter_list.append(out_read)
                output_fil_dic_rea[seq_ref_ini] = [best_result[1], best_result[2], best_result[4]]
                    # write best alignment into the dic

    write_reads(reads_filter_list, output_bam_filter_file_temp, bam_header_temp_out)
    write_reads(reads_output_list, output_bam_file_temp, bam_header_temp_out)

    # count finished processing number and output LOG, need a lock
    lock.acquire()  # set a lock

    finished_num = q_num.get()
    current_progress_bool = q_bool.get()

    finished_num += 1
    current_progress_ratio = finished_num / tot_tasks
    updated_progress_bool = log_pa(current_progress_ratio, current_progress_bool)

    q_num.put(finished_num)
    q_bool.put(updated_progress_bool)

    lock.release()  # release the lock


# define a function to judge whether the cigar is XXX M
def cigar_judge(cigar_jud):

    if cigar_jud[-1] == 'M':
        if cigar_jud[:-1].isdigit():
            return True
        else:
            return False
    else:
        return False


# define a score calculate function that work for fast mode
def fast_score(ref_seq_fas, ref_read_fas):

    score_fas = 0
    len_fas = len(ref_seq_fas)
    for i in range(len_fas):
        if ref_seq_fas[i] == ref_read_fas[i]:
            score_fas += 5
        else:
            score_fas -= 5

    return score_fas, score_fas/len_fas


# define a function to merge output temp bam
def bam_merge(temp_bam, all_bam):

    save_mer = pysam.set_verbosity(0)       # set suppressing HTSlib's warning message
    temp_bam_input = pysam.AlignmentFile(temp_bam, "rb")
    pysam.set_verbosity(save_mer)

    temp_allreads = temp_bam_input.fetch(until_eof = True)

    for read_temp in temp_allreads:
        all_bam.write(read_temp)

    temp_bam_input.close()

    return all_bam


# define a function to get last file name in input file path
def last_file_name(input_file):

    output_file = ''
    for i in range(len(input_file)):
        if input_file[len(input_file) - i -1] == '/':
            break
        else:
            output_file += input_file[len(input_file) - i -1]

    return output_file[::-1]


# define a function to sort bam and create index
def sort_bam_index(input_bam_unsort, output_bam_sorted, core_num_sort):

    core_str_sort = str(core_num_sort)
    pysam.sort("-@", core_str_sort, "-m", "2G", "-O", "BAM", "-n", "-o", output_bam_sorted, input_bam_unsort)


# define a function to write reads into a bam file
def write_reads(reads_ls_wri, bam_file_write_reads_wri, header_of_bam_wri):

    bam_output = pysam.AlignmentFile(bam_file_write_reads_wri, "wb", header = header_of_bam_wri)

    for read_write in reads_ls_wri:
        bam_output.write(read_write)

    bam_output.close()


# define a funtion to output log file during pairwise alignment
def log_pa(current_p, progress_b):

    if current_p >= 0.05 and progress_b["5%"] == False:
        sys.stderr.write("        %s\tPairwise alignment: 5%% Done! \n" % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))
        progress_b["5%"] = True
    elif current_p >= 0.1 and progress_b["10%"] == False:
        sys.stderr.write("        %s\tPairwise alignment: 10%% Done! \n" % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))
        progress_b["10%"] = True
    elif current_p >= 0.15 and progress_b["15%"] == False:
        sys.stderr.write("        %s\tPairwise alignment: 15%% Done! \n" % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))
        progress_b["15%"] = True
    elif current_p >= 0.2 and progress_b["20%"] == False:
        sys.stderr.write("        %s\tPairwise alignment: 20%% Done! \n" % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))
        progress_b["20%"] = True
    elif current_p >= 0.25 and progress_b["25%"] == False:
        sys.stderr.write("        %s\tPairwise alignment: 25%% Done! \n" % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))
        progress_b["25%"] = True
    elif current_p >= 0.3 and progress_b["30%"] == False:
        sys.stderr.write("        %s\tPairwise alignment: 30%% Done! \n" % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))
        progress_b["30%"] = True
    elif current_p >= 0.35 and progress_b["35%"] == False:
        sys.stderr.write("        %s\tPairwise alignment: 35%% Done! \n" % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))
        progress_b["35%"] = True
    elif current_p >= 0.4 and progress_b["40%"] == False:
        sys.stderr.write("        %s\tPairwise alignment: 40%% Done! \n" % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))
        progress_b["40%"] = True
    elif current_p >= 0.45 and progress_b["45%"] == False:
        sys.stderr.write("        %s\tPairwise alignment: 45%% Done! \n" % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))
        progress_b["45%"] = True
    elif current_p >= 0.5 and progress_b["50%"] == False:
        sys.stderr.write("        %s\tPairwise alignment: 50%% Done! \n" % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))
        progress_b["50%"] = True
    elif current_p >= 0.55 and progress_b["55%"] == False:
        sys.stderr.write("        %s\tPairwise alignment: 55%% Done! \n" % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))
        progress_b["55%"] = True
    elif current_p >= 0.6 and progress_b["60%"] == False:
        sys.stderr.write("        %s\tPairwise alignment: 60%% Done! \n" % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))
        progress_b["60%"] = True
    elif current_p >= 0.65 and progress_b["65%"] == False:
        sys.stderr.write("        %s\tPairwise alignment: 65%% Done! \n" % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))
        progress_b["65%"] = True
    elif current_p >= 0.7 and progress_b["70%"] == False:
        sys.stderr.write("        %s\tPairwise alignment: 70%% Done! \n" % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))
        progress_b["70%"] = True
    elif current_p >= 0.75 and progress_b["75%"] == False:
        sys.stderr.write("        %s\tPairwise alignment: 75%% Done! \n" % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))
        progress_b["75%"] = True
    elif current_p >= 0.8 and progress_b["80%"] == False:
        sys.stderr.write("        %s\tPairwise alignment: 80%% Done! \n" % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))
        progress_b["80%"] = True
    elif current_p >= 0.85 and progress_b["85%"] == False:
        sys.stderr.write("        %s\tPairwise alignment: 85%% Done! \n" % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))
        progress_b["85%"] = True
    elif current_p >= 0.9 and progress_b["90%"] == False:
        sys.stderr.write("        %s\tPairwise alignment: 90%% Done! \n" % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))
        progress_b["90%"] = True
    elif current_p >= 0.95 and progress_b["95%"] == False:
        sys.stderr.write("        %s\tPairwise alignment: 95%% Done! \n" % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))
        progress_b["95%"] = True
    elif current_p >= 1 and progress_b["100%"] == False:
        sys.stderr.write("        %s\tPairwise alignment: 100%% Done! \n" % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))
        progress_b["100%"] = True

    return progress_b


# Convert seconds into days-hours-minutes-seconds format
def convert_seconds(tot_seconds):
    out_seconds = tot_seconds % 60
    tot_minutes = int(tot_seconds // 60)
    out_minutes = int(tot_minutes % 60)
    tot_hours = int(tot_minutes // 60)
    out_hours = int(tot_hours % 24)
    tot_days = int(tot_hours // 24)

    return ("%d days  %d hours  %d minutes  %.1f seconds." %(tot_days, out_hours, out_minutes, out_seconds))


###################################################################
# Main program
###################################################################

if __name__ == "__main__":

    ################################################################
    # Input
    ################################################################


    # make sure the squences are input from 3' to 5'
    start_time = time.time()    # for the use of stat total time

    # file name and path
    input_bam_file = ""                                 # path of input bam
    output_bam_file = ""                                # path of output bam
    output_bam_filter_file = ""                         # path of bam contains reads to be filtered
    input_reference_file = ""                           # path of reference file
    temp_file_dic = 'tmp'                               # path of temp file directory
    # parameter
    max_softclip = 3                                    # defaulted max_softclip : 3
    min_score_per_base = 4.8                            # defaulted min score per base : 4.8
    max_softclip_input = 8                              # defaulted max_softclip input : 8
    core_num = 1                                        # defaulted core used : 1
    matrix_file = "test_8"                              # defaulted matrix file
    rm_multi_mode = "keep"                              # the mode to remove multiple mapping, default: "keep"
    rm_multi_switch = "off"                             # a signal weather need remove multiple mapping or not, default: off
        # keep: keep all mapping results if their scores all equal to the best score
        # remove: remove all mappling results if there are more than one read's scores equal to the best score
    rea_mode = "normal"
        # normal: just do realignment
        # fast: skip all XXX M reads to accelerate

    command_string = "python"  # store command

    for i in range(len(sys.argv)):
        input_para = sys.argv[i]
        command_string = command_string + " " + sys.argv[i]
        if input_para == '-i':
            input_bam_file = sys.argv[i + 1]
        elif input_para == '-o':
            output_bam_file = sys.argv[i + 1]
        elif input_para == '-f':
            output_bam_filter_file = sys.argv[i + 1]
        elif input_para == '-x':
            input_reference_file = sys.argv[i + 1]
        elif input_para == '-sc':
            max_softclip = int(sys.argv[i + 1])
        elif input_para == '-sc-i':
            max_softclip_input = int(sys.argv[i + 1])
        elif input_para == '-ms':
            min_score_per_base = float(sys.argv[i + 1])
        elif input_para == '-t':
            core_num = int(sys.argv[i + 1])
        elif input_para == '-m':
            matrix_file = sys.argv[i + 1]
        elif input_para == '-T':
            temp_file_dic = sys.argv[i + 1]
        elif input_para == '-rm-mode':
            rm_multi_mode = sys.argv[i + 1]
        elif input_para == '-rm':
            rm_multi_switch = sys.argv[i + 1]
        elif input_para == '--fast':
            rea_mode = "fast"
        elif input_para == '--normal':
            rea_mode = "normal"
        elif input_para == '-h' or input_para == '--help':
            print('-----------------------------------------------------Welcome to use PRAISE realignment-----------------------------------------------------')
            print('Usage:')
            print('    python realignment_reverse.py [Options]* -x <reference> {-i <bam_in>} [-o <bam_out> -f <bam_fil_out>]')
            print('')
            print('    <reference>     Genome reference (.fa file)')
            print('    <bam_in>        Input bam file, need to be sorted by name')
            print('    <bam_out>       Output sorted bam file with its index (.bai file)')
            print('    <bam_fil_out>   Output sorted bam file containing filtered reads, with its index (.bai file)')
            print('')
            print('[Options]')
            print('')
            print('  Input:')
            print('    -t [INT]        cores number, default =', core_num)
            print('    -m [MATRIX]     pairwise matrix, default =', matrix_file)
            print('    -T [Directory]  directory of temporary file, default =', temp_file_dic)
            print('')
            print('  Realignment mode:')
            print('    --normal        normal realignment mode, relatively slow')
            print('    --fast          fast realignment mode, skip reads with cigar XXX M, relatively fast')
            print('')
            print('  Filter:')
            print('    -sc [INT]       max softclip length, filter if softclip length is longer than [INT], default =', max_softclip)
            print('    -sc-i [INT]     max input softclip length, filter before realignment if softclip length is longer than [INT], default =', max_softclip_input)
            print('    -ms [FLOAT]     min score per base of reads, filter if score per base is lower than [FLOAT], default =', min_score_per_base)
            print('')
            print('  Remove multiple mapping:')
            print('    -rm on/off      Whether need remove multiple mapping after realignment or not, default =', rm_multi_switch)
            print('    -rm-mode keep   keep all mapping results if their scores all equal to the best score, default mode')
            print("    -rm-mode remove remove all mappling results if there are more than one read's scores equal to the best score")
            print('')
            print('  Other:')
            print('    -h / --help     help')
            exit(0)

    # LOG
    sys.stderr.write("%s\tRealignment Starts \n \n" % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

    # make a dictionary to store temp file
    ran_str = ''.join(random.sample(string.ascii_letters + string.digits, 8))
    temp_file_dic = temp_file_dic + '_' + ran_str + '_' + last_file_name(input_bam_file)

    os.mkdir(temp_file_dic)

    # LOG
    sys.stderr.write(
        "Input information: \n"
        "    Command: %s \n"
        "    Reference: %s \n"
        "    Input file: %s \n"
        "    Output file: %s \n"
        "    Output filtered file: %s \n"
        "    Directory of temporary file: %s \n"
        "    Number of cores used: %d \n"
        "    Matrix file used: %s \n \n"
        %(command_string, input_reference_file, input_bam_file, output_bam_file, output_bam_filter_file,
          temp_file_dic, core_num, matrix_file)
    )
    sys.stderr.write(
        "Filter information: \n"
        "    Max input softclip: %d \n"
        "    Max output softclip: %d \n"
        "    Min score per base of reads: %.2f/%.2f \n \n"
        %(max_softclip_input, max_softclip, min_score_per_base, 5.0)
    )
    sys.stderr.write(
        "Realignment mode: %s\n \n" %(rea_mode)
    )

    if rm_multi_switch == "on":
        sys.stderr.write(
            "Multiple mapping removal information: \n"
            "    State: On \n"
            "    Mode: %s \n \n"
            % (rm_multi_mode)
        )
    elif rm_multi_switch == "off":
        sys.stderr.write(
            "Multiple mapping removal information: \n"
            "    State: Off \n \n"
        )


    ###############################################################
    # Processing
    ###############################################################


    # load matrix
    substitution_matrices.load()
    matrix = substitution_matrices.load(matrix_file)


    # load reference
    reference_dic = load_ref(input_reference_file)

    # LOG
    sys.stderr.write("Processing information: \n")

    ##########
    #  making temp files part
    ##########

    # LOG
    sys.stderr.write("    %s\tMaking temporary files: Starts \n" % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

    mk_temp_result = mk_temp_file(input_bam_file, temp_file_dic, 20000)
    temp_file = mk_temp_result[0]

    temp_out_file = []
    temp_out_filter_file = []


    # create name list of temp output bam
    for ii in range(len(temp_file)):

        temp_out_file.append(temp_file[ii][:-4] + '.realign.bam')
        temp_out_filter_file.append(temp_file[ii][:-4] + '.realign.filter.bam')

    # LOG
    sys.stderr.write("    %s\tMaking temporary files: Done! \n"
                     "        %d temporary files made \n"
                     % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), mk_temp_result[1]))


    ##########
    # Pairwise alignment part
    ##########

    # Create a dictionary to store the progress
    progress_bool = {"10%": False,
                     "20%": False,
                     "30%": False,
                     "40%": False,
                     "50%": False,
                     "60%": False,
                     "70%": False,
                     "80%": False,
                     "90%": False,
                     "100%": False,
                     "5%": False,
                     "15%": False,
                     "25%": False,
                     "35%": False,
                     "45%": False,
                     "55%": False,
                     "65%": False,
                     "75%": False,
                     "85%": False,
                     "95%": False}

    # use Manager().Queue() to achieve communication between processing
    q_num = Manager().Queue()   # store number of finished processing
    q_bool = Manager().Queue()  # store progress_bool

    q_num.put(0)
    q_bool.put(progress_bool)

    # create a lock
    lock = Manager().Lock()

    # LOG
    sys.stderr.write("    %s\tPairwise alignment: Starts \n" % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

    # multiprocessing part
    po = Pool(core_num)
    for ii in range(len(temp_file)):
        po.apply_async(realign, (temp_file[ii], reference_dic, max_softclip, min_score_per_base, max_softclip_input, matrix, q_num, q_bool, len(temp_file), lock, rea_mode))

    po.close()
    po.join()


    # LOG
    sys.stderr.write("    %s\tPairwise alignment: Done! \n" % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

    ##########
    # Merge part
    ##########

    # LOG
    sys.stderr.write("    %s\tBam merge: Starts \n" % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

    # output file
    bam_input_all_outer = pysam.AlignmentFile(input_bam_file, "rb")
    bam_header = bam_input_all_outer.header

    # to create a temp file to store unsorted bam
    ran_str = ''.join(random.sample(string.ascii_letters + string.digits, 8))
    output_bam_file_temp_all = temp_file_dic + '/' + 'output_temp_' + ran_str + '.bam'
    output_bam_filter_file_temp_all = temp_file_dic + '/' + 'output_filter_temp_' + ran_str + '.bam'

    bam_output_all = pysam.AlignmentFile(output_bam_file_temp_all, "wb", header = bam_header)
    bam_output_filter_all = pysam.AlignmentFile(output_bam_filter_file_temp_all, "wb", header = bam_header)

    bam_input_all_outer.close()


    for ii in range(len(temp_out_file)):
        bam_output_all = bam_merge(temp_out_file[ii], bam_output_all)

    for ii in range(len(temp_out_filter_file)):
        bam_output_filter_all = bam_merge(temp_out_filter_file[ii], bam_output_filter_all)

    bam_output_all.close()
    bam_output_filter_all.close()

    if rm_multi_switch == 'on':
        output_bam_file_temp_all_sort = output_bam_file_temp_all[:-4] + '.sort.bam'

        # output bam sort and make index
        sort_bam_index(output_bam_file_temp_all, output_bam_file_temp_all_sort, core_num)
        sort_bam_index(output_bam_filter_file_temp_all, output_bam_filter_file, core_num)

        # LOG
        sys.stderr.write("    %s\tBam merge: Done! \n" % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

        ##########
        # Remove multiple mapping part
        ##########

        # LOG
        sys.stderr.write("    %s\tRemove multiple mapping: Starts \n" % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

        remove_multi_mapping.remove_mp(output_bam_file_temp_all_sort, output_bam_file, rm_multi_mode)

        # LOG
        sys.stderr.write("    %s\tRemove multiple mapping: Done! \n" % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

    elif rm_multi_switch == 'off':
        # output bam sort and make index
        sort_bam_index(output_bam_file_temp_all, output_bam_file, core_num)
        sort_bam_index(output_bam_filter_file_temp_all, output_bam_filter_file, core_num)

        # LOG
        sys.stderr.write("    %s\tBam merge: Done! \n" % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

    ##########
    # Remove temp files part
    ##########

    # LOG
    sys.stderr.write("    %s\tRemove temporary files: Starts \n" % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

    # remove temp directory and un sort bam
    shutil.rmtree(temp_file_dic)

    # LOG
    sys.stderr.write("    %s\tRemove temporary files: Done! \n" % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

    end_time = time.time()
    total_time = end_time - start_time

    sys.stderr.write("\n%s\tRealignment: Done! \n"
                     "    Total time used: %s \n"
                     % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), convert_seconds(total_time)))




