##########################################################################################################
# to split bam into several parts
# Author: Zhe Jiang
# created on 2021/12/17

##########################################################################################################


##############################################################
# import packages
##############################################################

import pysam
import sys


##############################################################
# define functions
##############################################################

# define a function to split bam file
def bam_split(input_file_name, temp_dic, file_num):


    last_input_file_name = last_file_name(input_file_name)

    if temp_dic[-1] == '/':
        temp_dic = temp_dic[:-1]

    # determine the reads number in every bam file
    tot_reads = int(pysam.view("-c", input_file_name)[:-1])
    reads_num = (tot_reads // file_num) + 1

    # store whole input bam
    bam_input_all = pysam.AlignmentFile(input_file_name, "rb")
    temp_bam_header = bam_input_all.header  # record the header

    # create file first
    for file_index in range(1, file_num + 1):
        temp_file_name = temp_dic + '/' + "hisat2_R2.sort" + '.split.' + str(file_index) + '.bam'
        bam_temp = pysam.AlignmentFile(temp_file_name, "wb", header = temp_bam_header)
        bam_temp.close()

    file_index = 1
    temp_file_name = temp_dic + '/' + "hisat2_R2.sort" + '.split.' + str(file_index) + '.bam'
    bam_temp = pysam.AlignmentFile(temp_file_name, "wb", header = temp_bam_header)

    all_bam_reads = bam_input_all.fetch(until_eof = True)

    last_read_name = ''                 # to store last read name, make sure paired read will be in the same file

    count = 1

    # read by read, write into a new temp file, 50000 reads per file
    for read in all_bam_reads:
        if count > reads_num and read.query_name != last_read_name:
            # close this temp file
            count = 1
            bam_temp.close()
            # create the next temp file
            file_index += 1
            temp_file_name = temp_dic + '/' + "hisat2_R2.sort" + '.split.' + str(file_index) + '.bam'
            bam_temp = pysam.AlignmentFile(temp_file_name, "wb", header = temp_bam_header)
            bam_temp.write(read)
        else:
            count += 1
            bam_temp.write(read)
            last_read_name = read.query_name

    # close the last bam_temp
    bam_temp.close()


def last_file_name(input_file):

    output_file = ''
    for i in range(len(input_file)):
        if input_file[len(input_file) - i -1] == '/':
            break
        else:
            output_file += input_file[len(input_file) - i -1]

    return output_file[::-1]




##############################################################
# Input
##############################################################

bam_input_file = ""
output_dir = ""
bam_num = 10                        # default number of bam file

for i in range(len(sys.argv)):
    input_para = sys.argv[i]
    if input_para == '-i':
        bam_input_file = sys.argv[i + 1]
    elif input_para == '-o':
        output_dir = sys.argv[i + 1]
    elif input_para == '-n':
        bam_num = int(sys.argv[i + 1])
    elif input_para == '-h' or input_para == '--help':
        print('-----------------------------------------------------Welcome to use PRAISE alignment BAM split-----------------------------------------------------')
        print('Usage:')
        print('    python bam_split.py [Options]* {-i <bam_in>} [-o <out_bam_dir>]')
        print('')
        print('    <bam_in>        Input bam file, need to be sorted by name')
        print('    <bam_out_dir>   Directory of output bam file')
        print('')
        print('[Options]')
        print('')
        print('  Output:')
        print('    -n [INT]        Number of output bam file, default = ', bam_num)
        print('')
        print('  Other:')
        print('    -h / --help     help')
        exit(0)





##############################################################
# Processing
##############################################################

bam_split(bam_input_file, output_dir, bam_num)
