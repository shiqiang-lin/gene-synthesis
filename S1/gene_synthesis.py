#
#usage:
#python3.10 gene_synthesis.py gene_sequence.fasta
#


import sys
import os,shutil
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
from Bio import SeqIO
from operator import itemgetter



#Here are the functions used for running the py file.

"""
The input is a DNA fragment string and a set Tm, and the output is the conbination of
overlap sequence, Tm of overlap sequence, and the |Tm - Tm_set| with the |Tm - Tm_set|
the lowest among all the conbinations.
"""
def get_overlap_by_Tm(DNA_fragment_str,Tm_set):
    slice_list = []
    slices_list = []
    min_overlap_len = 10
    max_overlap_len = 25
    
    DNA_fragment_len = len(DNA_fragment_str)
    for i in range(DNA_fragment_len-min_overlap_len,DNA_fragment_len-max_overlap_len-1,-1):
        overlap_slice = DNA_fragment_str[i:]
        overlap_slice_seq = Seq(overlap_slice)
        mt_of_overlap_slice_seq = mt.Tm_NN(overlap_slice_seq)
        mt_distance = abs(mt_of_overlap_slice_seq - Tm_set)

        slice_list.append(overlap_slice)
        slice_list.append(mt_of_overlap_slice_seq)
        slice_list.append(mt_distance)

        slices_list.append(slice_list)
        slice_list = []

    slices_list_sorted = sorted(slices_list,key=itemgetter(2))
    result_slices_list = slices_list_sorted[0]

    return result_slices_list[0],result_slices_list[1],result_slices_list[2]



"""
The input is overlap_str, remaining_gene_str, set_fragment_len, and the
output is fragment_Li_str, new_remaining_gene_str.
The running condition is len(overlap_str)+len(remaining_gene_str) >= set_fragment_len.
We see the condition from outside the function.
"""
def get_next_fragment(overlap_str,remaining_gene_str,set_fragment_len):
    fragment_Li_str = overlap_str + remaining_gene_str[0:set_fragment_len-len(overlap_str)]
    new_remaining_gene_str = remaining_gene_str[set_fragment_len-len(overlap_str):]

    return fragment_Li_str,new_remaining_gene_str



"""
The input is overlap_str and remaining_gene_str. The function gets the tail
fragment(s) and stores in the fragments_list.
The running condition is 0<len(fragments_list[-1][1])+len(remaining_gene_str)
<set_fragment_len.
"""
def get_tail_fragments(overlap_str,remaining_gene_str):
    global fragments_list

    min_remaining_gene_str_len = 20
    fetch_bases_from_Ln_1_len = 20
    if len(remaining_gene_str) <= min_remaining_gene_str_len:
        fetch_bases_from_Ln_1 = fragments_list[-1][0][-fetch_bases_from_Ln_1_len:]
        new_fragment_Ln_1 = fragments_list[-1][0][0:len(fragments_list[-1][0])-fetch_bases_from_Ln_1_len]

        new_Ln_1_slice = get_overlap_by_Tm(new_fragment_Ln_1,Tm_set)

        Li_fragment_seq_in_func_tmp_list = []
        Li_fragment_seq_in_func_tmp_list.append(new_fragment_Ln_1)

        
        new_Ln_1 = Li_fragment_seq_in_func_tmp_list + list(new_Ln_1_slice)
        #Here is the new Ln-1.
        fragments_list[-1] = new_Ln_1
        
        Ln_fragment_str = new_Ln_1[1] + fetch_bases_from_Ln_1 + remaining_gene_str

        Lj_fragment_seq_in_func_tmp_list = []
        Lj_fragment_seq_in_func_tmp_list.append(Ln_fragment_str)
        
        fragments_list.append(Lj_fragment_seq_in_func_tmp_list)
        #Now, we have obtained all the needed sense-strand fragments.

    else:
        Ln_fragment_str = fragments_list[-1][1] + remaining_gene_str

        Lk_fragment_seq_in_func_tmp_list = []
        Lk_fragment_seq_in_func_tmp_list.append(Ln_fragment_str)
        
        fragments_list.append(Lk_fragment_seq_in_func_tmp_list)
        #Now, we have obtained all the needed sense-strand fragments.



"""
The input is gene_str and Tm_set, and the output is the sense primer and
antisense primer (including Tm and |Tm - Tm_set|) which have Tms that are
nearest to the Tm_set.
""" 
def get_gene_primers(gene_str,Tm_set):
    min_primer_bases = 15
    max_primer_bases = 25

    gene_str_sense_seq = Seq(gene_str)
    gene_str_antisense_seq = gene_str_sense_seq.reverse_complement()
    gene_antisense_str = str(gene_str_antisense_seq)

    sense_primers_list = []
    antisense_primers_list = []
    sense_primer_tmp_list = []
    antisense_primer_tmp_list = []
    

    for i in range(min_primer_bases-1,max_primer_bases-1):
        sense_primer_str = gene_str[0:i+1]
        sense_primer_str_seq = Seq(sense_primer_str)
        sense_primer_Tm = mt.Tm_NN(sense_primer_str_seq)
        sense_primer_tmp_list.append(sense_primer_str)
        sense_primer_tmp_list.append(sense_primer_Tm)
        sense_primer_tmp_list.append(abs(sense_primer_Tm-Tm_set))
        sense_primers_list.append(sense_primer_tmp_list)
        sense_primer_tmp_list = []

        antisense_primer_str = gene_antisense_str[0:i+1]
        antisense_primer_str_seq = Seq(antisense_primer_str)
        antisense_primer_Tm = mt.Tm_NN(antisense_primer_str_seq)
        antisense_primer_tmp_list.append(antisense_primer_str)
        antisense_primer_tmp_list.append(antisense_primer_Tm)
        antisense_primer_tmp_list.append(abs(antisense_primer_Tm-Tm_set))
        antisense_primers_list.append(antisense_primer_tmp_list)
        antisense_primer_tmp_list = []

    sense_primers_list_sorted = sorted(sense_primers_list,key=itemgetter(2))
    antisense_primers_list_sorted = sorted(antisense_primers_list,key=itemgetter(2))

    return sense_primers_list_sorted[0],antisense_primers_list_sorted[0]



#Here is the main part of py file.
#get the DNA sequence
gene_file = sys.argv[1]
gene_str = str(SeqIO.read(gene_file, "fasta").seq)

if len(gene_str) < 61:
    print("The DNA is less than 61 bases.")
    print("Please synthesize the sense and antisense strands, anneal them,")
    print("and ligate the gene to the vector.")
    sys.exit()


print("Here is the gene sequence:")
print(gene_str)
print("The length of the gene sequence is %s." % len(gene_str))
print("\n", end = '')

"""
This parameter can be set according to the DNA synthesis company.
It can be at least 50 bps, and can reach 70 bps or longer.
You may need to consult the company to decide the parameter.
"""
set_fragment_len = 59   

#list containing fragments L1, L2, L3,..., Ln-1, and their relevant
#objects, which are overlap, Tm, delta Tm, for each Tm_set.
fragments_list = []

"""
This list contains the fragments for all Tm_sets. Each element includes
Tm_set, framents_list, gene_primers, mean Tm, variance for the Tms of all
overlaps and antisense_gene_primer.
"""
fragments_for_all_Tm_sets_list = []

fragments_for_Tm_set_tmp_list = []

#Get the fragments for each Tm_set.
"""
The range of Tm_set is from 54 to 68, including 68. However, this
range can be wider, for example from 52 to 71, or other values of
lower bound and upper bound.
"""
for Tm_set in range(54,69):
    #get the L1 fragment
    L1_fragment_list = gene_str[0:set_fragment_len]
    L1_fragment_str =''.join(L1_fragment_list)
    L1_slice = get_overlap_by_Tm(L1_fragment_str,Tm_set)
    Li_fragment_seq_tmp_list = []
    Li_fragment_seq_tmp_list.append(L1_fragment_str)

    L1 = Li_fragment_seq_tmp_list + list(L1_slice)
    fragments_list.append(L1)

    #get the rest fragments
    remaining_gene_str = gene_str[set_fragment_len:]
    overlap_str_plus_remaining_gene_str_len = len(fragments_list[-1][1]) + len(remaining_gene_str)

    while(overlap_str_plus_remaining_gene_str_len >= set_fragment_len):
        overlap_str = fragments_list[-1][1]
        fragment_L_tmp_str,remaining_gene_str = get_next_fragment(overlap_str,remaining_gene_str,set_fragment_len)
        L_tmp_slice = get_overlap_by_Tm(fragment_L_tmp_str,Tm_set)
        
        Li_fragment_seq_in_func_tmp_list = []
        Li_fragment_seq_in_func_tmp_list.append(fragment_L_tmp_str)
        
        L_tmp = Li_fragment_seq_in_func_tmp_list + list(L_tmp_slice)
        fragments_list.append(L_tmp)

        overlap_str_plus_remaining_gene_str_len = len(fragments_list[-1][1]) + len(remaining_gene_str)

    #get the tail fragments
    if len(remaining_gene_str) == 0:
        fragments_list[-1] = fragments_list[-1][0]
       
    if len(remaining_gene_str) > 0:
        overlap_str = fragments_list[-1][1]
        get_tail_fragments(overlap_str,remaining_gene_str)

    """
    Next, we need to get the optimal Tm for the needed sense fragments. The mean Tm
    of the needed sense fragments is used as the optimal Tm.
    """
    Tm_recommend = 0

    for i in range(len(fragments_list)-1):
        Tm_recommend = Tm_recommend + fragments_list[i][2]
        
    Tm_recommend = Tm_recommend / (len(fragments_list)-1)


    #get the gene primers according to the Tm recommend
    gene_primers = []
    sense_gene_primer,antisense_gene_primer = get_gene_primers(gene_str,Tm_recommend)
    gene_primers.append(sense_gene_primer)
    gene_primers.append(antisense_gene_primer)
    ###Now, we have obtained all the needed fragments and the gene primers.
        

    fragments_for_Tm_set_tmp_list = []
    fragments_for_Tm_set_tmp_list.append(Tm_set)
    fragments_for_Tm_set_tmp_list.append(fragments_list)
    fragments_for_Tm_set_tmp_list.append(gene_primers)

    variance_of_overlaps_and_antisense_gene_primer = 0
    for i in range(len(fragments_list)-1):
        variance_of_overlaps_and_antisense_gene_primer = variance_of_overlaps_and_antisense_gene_primer + \
                                                         pow(fragments_list[i][2] - Tm_recommend,2)

    variance_of_overlaps_and_antisense_gene_primer = variance_of_overlaps_and_antisense_gene_primer + \
                                                     pow(gene_primers[1][1],2)

    variance_of_overlaps_and_antisense_gene_primer = variance_of_overlaps_and_antisense_gene_primer / \
                                                     len(fragments_list)

    fragments_for_Tm_set_tmp_list.append(Tm_recommend)
    fragments_for_Tm_set_tmp_list.append(variance_of_overlaps_and_antisense_gene_primer)

    fragments_for_all_Tm_sets_list.append(fragments_for_Tm_set_tmp_list)
    fragments_for_Tm_set_tmp_list = []
    #empty the container
    fragments_list = []

#get the fragments for the Tm_set that have the least variance
fragments_for_all_Tm_sets_list_sorted = sorted(fragments_for_all_Tm_sets_list,key=itemgetter(4))

print("The optimal Tm_set is %s Celsius degree." % fragments_for_all_Tm_sets_list_sorted[0][0])
print("The recommend Tm for gene synthesis by PCR is %s Celsius degree." % \
      round(fragments_for_all_Tm_sets_list_sorted[0][3],2))



#plot the Tms for overlaps
import matplotlib.pyplot as plt

x = [i for i in range(len(fragments_for_all_Tm_sets_list_sorted[0][1])-1)]

Tm_for_each_overlap_list = []
for i in range(len(fragments_for_all_Tm_sets_list_sorted[0][1])-1):
    Tm_for_each_overlap_list.append(fragments_for_all_Tm_sets_list_sorted[0][1][i][2])

plt.scatter(x,Tm_for_each_overlap_list,c='royalblue',s=20)

plt.axhline(fragments_for_all_Tm_sets_list_sorted[0][3])

plt.xlabel("overlap")
plt.ylabel("Tm(Celsius degree)")
plt.show()



#print the needed fragments to screen
result_fragments_list = fragments_for_all_Tm_sets_list_sorted[0][1]


#add the start position and end position for each overlap
end_position_of_overlap = len(result_fragments_list[0][0]) - 1
start_position_of_overlap = end_position_of_overlap + 1 - len(result_fragments_list[0][1])

#Please notice that the start position and end position of each overlap are calculated with the first
#base as 1. We need only to add 1 to the positions of first overlap.
result_fragments_list[0].append(start_position_of_overlap + 1)
result_fragments_list[0].append(end_position_of_overlap + 1)

for i in range(1,len(result_fragments_list)-1):
    end_position_of_overlap = result_fragments_list[i-1][5] + len(result_fragments_list[i][0]) \
                              - len(result_fragments_list[i-1][1])
    start_position_of_overlap = end_position_of_overlap + 1 - len(result_fragments_list[i][1])

    result_fragments_list[i].append(start_position_of_overlap)
    result_fragments_list[i].append(end_position_of_overlap)


print("\n", end = '')
print("Here are the needed fragments:")
for i in range(len(result_fragments_list)):
    print("L" + str(i+1) + ":" + "  " + result_fragments_list[i][0])

#print the gene primers
result_gene_primers_list = fragments_for_all_Tm_sets_list_sorted[0][2]

print("\n", end = '')
print("Here are the gene primers:")

print("sense gene primer:   " + result_fragments_list[0][0])
print("antisense gene primer:   " + result_fragments_list[1][0])

#store the result to txt files
#get current path
current_directory = os.getcwd()
print("\n", end = '')
print("Current directory is %s." % current_directory)
print("\n", end='')


#make a directory for storing result txt files, if the folder has alreadly existed, first
#delete, then make the directory. 
current_directory = os.getcwd()
foldername = 'gene_synthesis'
dirs = os.listdir(current_directory)

if(foldername not in dirs):
    os.mkdir(foldername)
else:
    shutil.rmtree(foldername)
    os.mkdir(foldername)

#go to the gene_synthesis folder
os.chdir(foldername)

#open the all_sense_fragments.txt file
all_sense_fragments_file = open("all_sense_fragments.txt",'w')

#overlap".ljust(25,' ') is the same as max_primer_bases
print('\r\n',end='',file=all_sense_fragments_file)
print("fragment".ljust(set_fragment_len+10,' '),\
      "overlap".ljust(25+10,' '),\
      "Tm".ljust(10,' '),\
      "Tm-Tm_mean".ljust(20,' '),\
      "fragment_len".ljust(20,' '),\
      "overlap_len".ljust(20,' '),\
      "start_pos-end_pos".ljust(20,' '),\
      sep="",\
      file=all_sense_fragments_file
      )
print('\r\n',end='',file=all_sense_fragments_file)


Tm_mean = fragments_for_all_Tm_sets_list_sorted[0][3]


for i in range(len(result_fragments_list)-1):
    C1_print_line_str = "L" + str(i+1) + ":" + "  " + str(result_fragments_list[i][0])
    C2_print_line_str = str(result_fragments_list[i][1])
    C3_print_line_str = str(round(result_fragments_list[i][2],2))
    C4_print_line_str = str(round(result_fragments_list[i][2] - Tm_mean,2))
    C5_print_line_str = str(len(result_fragments_list[i][0]))
    C6_print_line_str = str(len(result_fragments_list[i][1]))
    #overlap_start_end_pos
    C7_print_line_str = str(result_fragments_list[i][4]) + '-' + str(result_fragments_list[i][5])    


    
    print(C1_print_line_str.ljust(set_fragment_len+10,' '),\
          C2_print_line_str.ljust(25+10,' '),\
          C3_print_line_str.ljust(10,' '),\
          C4_print_line_str.ljust(20,' '),\
          C5_print_line_str.ljust(20,' '),\
          C6_print_line_str.ljust(20,' '),\
          C7_print_line_str.ljust(20,' '),\
          sep="",\
          file=all_sense_fragments_file
          )

    
C1_print_line_str = "L" + str(len(result_fragments_list)) + ":" + "  " + \
                    str(result_fragments_list[len(result_fragments_list)-1][0])
C2_print_line_str = ''
C3_print_line_str = ''
C4_print_line_str = ''
C5_print_line_str = str(len(result_fragments_list[len(result_fragments_list)-1][0]))
C6_print_line_str = ''
C7_print_line_str = ''

print(C1_print_line_str.ljust(set_fragment_len+10,' '),\
      C2_print_line_str.ljust(25+10,' '),\
      C3_print_line_str.ljust(10,' '),\
      C4_print_line_str.ljust(20,' '),\
      C5_print_line_str.ljust(20,' '),\
      C6_print_line_str.ljust(20,' '),\
      C7_print_line_str.ljust(20,' '),\
      sep="",\
      file=all_sense_fragments_file
      )
all_sense_fragments_file.close()


#open the sense_antisense_fragments.txt file
sense_antisense_fragments_file = open("sense_antisense_fragments.txt",'w')

#overlap".ljust(25,' ') is the same as max_primer_bases
print('\r\n',end='',file=sense_antisense_fragments_file)
print("fragment".ljust(set_fragment_len+10,' '),\
      "overlap".ljust(25+10,' '),\
      "Tm".ljust(10,' '),\
      "Tm-Tm_mean".ljust(20,' '),\
      "fragment_len".ljust(20,' '),\
      "overlap_len".ljust(20,' '),\
      "start_pos-end_pos".ljust(20,' '),\
      sep="",\
      file=sense_antisense_fragments_file
      )
print('\r\n',end='',file=sense_antisense_fragments_file)


Tm_mean = fragments_for_all_Tm_sets_list_sorted[0][3]


for i in range(len(result_fragments_list)-1):
    #The following two lines reverse complement the even fragments in the sense_antisense_fragments.txt.
    #We count the fragments in the sense_antisense_fragments.txt from 1.
    if i % 2 == 1:
        C1_print_line_str = "L" + str(i+1) + ":" + "  " + \
                            str(Seq(result_fragments_list[i][0]).reverse_complement())
    else:
        C1_print_line_str = "L" + str(i+1) + ":" + "  " + str(result_fragments_list[i][0])
                   
    C2_print_line_str = str(result_fragments_list[i][1])
    C3_print_line_str = str(round(result_fragments_list[i][2],2))
    C4_print_line_str = str(round(result_fragments_list[i][2] - Tm_mean,2))
    C5_print_line_str = str(len(result_fragments_list[i][0]))
    C6_print_line_str = str(len(result_fragments_list[i][1]))
    #overlap_start_end_pos
    C7_print_line_str = str(result_fragments_list[i][4]) + '-' + str(result_fragments_list[i][5])    


    
    print(C1_print_line_str.ljust(set_fragment_len+10,' '),\
          C2_print_line_str.ljust(25+10,' '),\
          C3_print_line_str.ljust(10,' '),\
          C4_print_line_str.ljust(20,' '),\
          C5_print_line_str.ljust(20,' '),\
          C6_print_line_str.ljust(20,' '),\
          C7_print_line_str.ljust(20,' '),\
          sep="",\
          file=sense_antisense_fragments_file
          )

    
C1_print_line_str = "L" + str(len(result_fragments_list)) + ":" + "  " + \
                    str(result_fragments_list[len(result_fragments_list)-1][0])
C2_print_line_str = ''
C3_print_line_str = ''
C4_print_line_str = ''
C5_print_line_str = str(len(result_fragments_list[len(result_fragments_list)-1][0]))
C6_print_line_str = ''
C7_print_line_str = ''

print(C1_print_line_str.ljust(set_fragment_len+10,' '),\
      C2_print_line_str.ljust(25+10,' '),\
      C3_print_line_str.ljust(10,' '),\
      C4_print_line_str.ljust(20,' '),\
      C5_print_line_str.ljust(20,' '),\
      C6_print_line_str.ljust(20,' '),\
      C7_print_line_str.ljust(20,' '),\
      sep="",\
      file=sense_antisense_fragments_file
      )
sense_antisense_fragments_file.close()

print("Please note that in the sense_antisense_fragments.txt, we list the overlap sequences \
in the second columm with all sense sequences.")
print('\r\n',end='')

#open the gene_primers.txt file
gene_primers_file = open("gene_primers.txt",'w')

print('\r\n',end='',file=gene_primers_file)
print("gene_primer".ljust(70,' '),\
      "Tm".ljust(25,' '),\
      "Tm-Tm_mean".ljust(20,' '),\
      sep="",\
      file=gene_primers_file
      )
print('\r\n',end='',file=gene_primers_file)

for i in range(len(result_gene_primers_list)):
    if i % 2 == 0:
        C1_print_line_str = "sense_gene_primer:" + "   " + str(result_gene_primers_list[i][0])
    else:
        C1_print_line_str = "antisense_gene_primer:" + "   " + str(result_gene_primers_list[i][0])

    C2_print_line_str = str(round(result_gene_primers_list[i][1],2))
    C3_print_line_str = str(round(result_gene_primers_list[i][2],2))
      
    print(C1_print_line_str.ljust(70,' '),\
          C2_print_line_str.ljust(25,' '),\
          C3_print_line_str.ljust(20,' '),\
          sep="",\
          file=gene_primers_file
          )

gene_primers_file.close()
print("The result txt files have been stored in the folder 'gene_synthesis'.")
print("\n", end='')

#change to current_directory
os.chdir(current_directory)



