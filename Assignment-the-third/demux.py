#!/usr/bin/env python

#barcodes = "/projects/bgmp/shared/2017_sequencing/indexes.txt"

#./demux.py -R1 ../TEST-input_FASTQ/test_R1.fq.gz -R2 ../TEST-input_FASTQ/test_R2.fq.gz -R3 ../TEST-input_FASTQ/test_R3.fq.gz -R4 ../TEST-input_FASTQ/test_R4.fq.gz -bar /projects/bgmp/shared/2017_sequencing/indexes.txt

import gzip
import bioinfo
import itertools
import argparse
import math
import matplotlib.pyplot as plt

def get_args():
    parser = argparse.ArgumentParser(description="A program to introduce yourself")
    parser.add_argument("-R1", help = "Read 1: Biological Read", required = True, type = str)
    parser.add_argument("-R2", help = "Read 2: Index 1", required = True, type = str)
    parser.add_argument("-R3", help = "Read 3: Index 2", required = True, type = str)
    parser.add_argument("-R4", help = "Read 4: Biological Read", required = True, type = str)
    parser.add_argument("-bar", help = "Known Barcodes/Indexes", required = True, type = str)
    return parser.parse_args()





args = get_args()
R1 = args.R1
R2 = args.R2
R3 = args.R3
R4 = args.R4
barcodes = args.bar


#ADD LATER:

# R1 = "R1_test.fq.gz" #"/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz" 
# R2 = "R2_test.fq.gz" #"/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz"
# R3 = "R3_test.fq.gz" #"/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz"
# R4 = "R4_test.fq.gz" #"/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz"


#Establish Dictionary
bases = {"A":"T", 
         "C":"G", 
         "T":"A", 
         "G":"C",
         "N":"N"}


def rev_comp(seq):
        rev_comp = ""
        seq = seq.strip()
        for nuc in seq: #need add new line to match input
            rev_comp += bases[nuc]
        rev_comp = rev_comp[::-1] 
        return(rev_comp)



known_indexes = set() #establish an empty set 
with open(barcodes,"r") as fh:
    next(fh)
    for line in fh: 
        new_line = line.strip().split("\t") #strips the new line, splits by tab
        known_indexes.add(new_line[4])
#print(line)


#Open all the files, use the R1 file handles 
R1_fh = gzip.open(R1, "rt") #change the R1-R4 variables to be adapted to argparse
R2_fh = gzip.open(R2, "rt")
R3_fh = gzip.open(R3, "rt")
R4_fh = gzip.open(R4, "rt")



def open_files():
    all_files = {}
    for barcodes in known_indexes:
          file_1 = f'outfiles/{barcodes}-{barcodes}_R1.fq' #KEY IS PATH TO WHERE STUFF GOES 
          file_2 = f'outfiles/{barcodes}-{barcodes}_R2.fq'
          all_files[barcodes+"R1"] = open(file_1,"w") #VALUE IS IMPORTANT, KEY DOESN'T MATTER 
          all_files[barcodes+"R4"] = open(file_2,"w") #VALUE ALLOWS US TO WRITE TO THIS PATH
    all_files["outfiles/unk_R1.fq"] = open("outfiles/unk_R1.fq","w") #HARDCODED FOR UNK AND HOPPED 
    all_files["outfiles/unk_R2.fq"] = open("outfiles/unk_R2.fq","w")
    all_files["outfiles/hopped_R1.fq"] = open("outfiles/hopped_R1.fq","w")
    all_files["outfiles/hopped_R2.fq"] = open("outfiles/hopped_R2.fq","w")
    return(all_files)

all_files = open_files()
#Isolate the Sequence Line in all the Files 
#ISOLATE ALL RECORDS 
record_4_lines_R1= [] #empty string where we append every record from R1 in to 
record_4_lines_R4 = []  #empty string where we append every record from R4 in to 
counter=0 #to keep track of the lines to write into a file
i=0 #to count the line number
tot_matched = 0
tot_hopped = 0
tot_unknown = 0 
#ESTABLISH PERMUTATIONS AND DICT TO KEEP TRACK OF INDEX PAIRS
index_dict = {}
for permut in itertools.product(list(known_indexes), repeat=2):
    index_dict[permut] = 0
index_dict[("unknown", "unknown")] = 0  

    
for line_R1,line_R2,line_R3,line_R4 in zip(R1_fh, R2_fh, R3_fh, R4_fh):
    i+=1 #ensure that modulus works later on 
    line_R1 = line_R1.strip('\n') #literally stripping the new line 
    line_R2 = line_R2.strip('\n')
    line_R3 = line_R3.strip('\n')
    line_R4 = line_R4.strip('\n')
    #print(i)
    if i%4==2: #Sequence Line
        index1 = line_R2
        index2 = line_R3
    if i%4==0: #pulling phred score 
        phred1= line_R2
        phred2 = line_R3
    if counter<4: #Because the record is 4 lines, if less than 4 we know that it is still the record
        record_4_lines_R1.append(line_R1) #we append the line to make a full record 
        record_4_lines_R4.append(line_R4)
        counter+=1  
        if counter==4: 
            counter=0 #reset everything
            rev_index2 = rev_comp(index2)
            add_to_header = str(index1) +  "-" + rev_index2
            new_head_R1 = record_4_lines_R1[0]+' '+ add_to_header
            new_head_R4 = record_4_lines_R4[0]+' '+ add_to_header
            record_4_lines_R1[0] = new_head_R1
            record_4_lines_R4[0] = new_head_R1
#ADDING HEADERS:
            # record_copy_R1 = record_4_lines_R1 
            # record_copy_R4 = record_4_lines_R4
            #THIS ELIMINATES UNKNOWNS
            if bioinfo.qual_score(phred1) < 26 or bioinfo.qual_score(phred2) < 26 or index1 not in known_indexes or rev_index2 not in known_indexes: #looks at qscore of index files ONLY:
                all_files["outfiles/unk_R1.fq"].write(f'{record_4_lines_R1[0]}\n{record_4_lines_R1[1]}\n{record_4_lines_R1[2]}\n{record_4_lines_R1[3]}\n') 
                all_files["outfiles/unk_R2.fq"].write(f'{record_4_lines_R4[0]}\n{record_4_lines_R4[1]}\n{record_4_lines_R4[2]}\n{record_4_lines_R4[3]}\n')
                index_dict[("unknown","unknown")] += 1
                tot_unknown += 1
            elif index1 != rev_index2:
                all_files["outfiles/hopped_R1.fq"].write(f'{record_4_lines_R1[0]}\n{record_4_lines_R1[1]}\n{record_4_lines_R1[2]}\n{record_4_lines_R1[3]}\n')
                all_files["outfiles/hopped_R2.fq"].write(f'{record_4_lines_R4[0]}\n{record_4_lines_R4[1]}\n{record_4_lines_R4[2]}\n{record_4_lines_R4[3]}\n')
                index_dict[(index1,rev_index2)] += 1
                tot_hopped += 1
            #EVERYTHING THAT WAS MATCHED
            else: 
                all_files[f'{index1}R1'].write(f'{record_4_lines_R1[0]}\n{record_4_lines_R1[1]}\n{record_4_lines_R1[2]}\n{record_4_lines_R1[3]}\n') #VALUE IS IMPORTANT, KEY DOESN'T MATTER 
                all_files[f'{index1}R4'].write(f'{record_4_lines_R4[0]}\n{record_4_lines_R4[1]}\n{record_4_lines_R4[2]}\n{record_4_lines_R4[3]}\n')
                index_dict[(index1,rev_index2)] += 1
                tot_matched += 1
            record_4_lines_R1=[]
            record_4_lines_R4 = []
total_records = i/4
#print(total_records)

title = f'index pairs\t\tNumber of Reads\tPercent Mapped'
print(title)
for key in index_dict:
    percent = 100*index_dict[key]/total_records
    out = f'{key[0]}\t{key[1]}\t{index_dict[key]}\t{percent}' 
    print(out)     

#WRITE STATISTICS TO FILES
all_percentages = [] #store all the percentages as a list outside the loop
matched_indexes = [] #list of matching indexes for graph
total_counts = []
with open("output.txt", "w") as out:
    new_tot_unk = 100*tot_unknown/total_records
    new_tot_matched = 100*tot_matched/total_records
    new_tot_hopped = 100*tot_hopped/total_records
    out.write(f'Summary of Statistics\n')
    out.write(f'Total percent matched reads: {new_tot_matched}\n')
    out.write(f'Total percent hopped reads: {new_tot_hopped}\n')
    out.write(f'Total percent unknown reads: {new_tot_unk}\n')
    for keys in index_dict:
        percent = 100*index_dict[keys]/total_records
        if keys == ("unknown","unknown"): 
            out.write(f'The percent of mapped reads for {keys[0]}-{keys[1]} is: {percent}% ({index_dict[keys]})\n')
        elif keys[0] == keys[1]:
            out.write(f'The percent of mapped reads for {keys[0]}-{keys[1]} is: {percent}% ({index_dict[keys]})\n')
            matched_indexes.append(keys[0]) #write matched keys to a new list for a plot
            all_percentages.append(percent)
            total_counts.append(index_dict[keys])
        elif keys[0] != keys[1]:
            out.write(f'The percent of mapped reads for {keys[0]}-{keys[1]} is: {percent}% ({index_dict[keys]})\n')


with open("matched_stats.tsv", "w") as matched:
    matched.write(f'Index\tBarcodes\tTotal Counts\tPercent Mapped\n')
    for i, values in enumerate(matched_indexes):
        matched.write(f'{i}\t{matched_indexes[i]}\t{total_counts[i]}\t{all_percentages[i]}%\n')


#OUTPUT A HORIZONTAL GRAPH

x_value = all_percentages
y_value = matched_indexes
plt.barh(y_value,x_value, color = 'orchid')
plt.xlabel('Percent Mapped')
plt.ylabel('Indexes')
plt.title("Percent Mapped Reads for Matched Indexes")
plt.tight_layout() #fixes axis cut off 
plt.savefig("output_fig.png")






# for value in index_dict:
#     percent = 100*index_dict[keys]/total_records
#     x = f'{value[0]}-{value[1]}
#     y =







    
        
# index_dict[keys]
# for keys,values in index_dict: 
#     mapped_reads = values/total_records
# print(mapped_reads)

# values = index_dict.values()
# print(values)





#table = open("index_dict","w")
# index_keys = index_dict.keys()
# index_values = index_dict.values()
# print(f'{index_keys}\t{index_values}')

# for key in index_dict:
    
#table = table.write(f'{index_keys}\t{index_values}')

    

    #             #print("record_copy")
    # if i%4==2: #Sequence Line
    #     index1 = line_R2
    #     index2 = line_R3
    # if i%4==0: #pulling phred score 
    #     phred1= line_R2
    #     phred2 = line_R3
         
#STATISTICS
#print(index_dict)
         
       

#Reverse Complement R3 Index (sequence line) in the Header and Reverse the Quality Too 

#Check if the indexes match or not: barcode or unknown NOT hopped 

#ESTABLISH COMBINATIONS 



    


# file_handles = all_files.keys()
# print(file_handles)

# index_dict_R1 = {}
# index_dict_R2 = {}
# counts_1 = 0 
# counts_2 = 0 
# for key in known_indexes:
#     fh_1 = f'{key}-{key}_R1.fq'
#     fh_2 = f'{key}-{key}_R2.fq'
#     index_dict_R1[fh_1] = counts_1
#     index_dict_R2[fh_2] = counts_2
# print(index_dict_R1)
# print(index_dict_R2)


#     index_dict[indices] =
    
#permut = list(itertools.product(['A', 'C', 'T','G'], repeat=8)) #gives us the number of possible ACTG index combinations, set to 8 bc indexes are 8 char
#print(len(permut))
# permut = ', '.join([' '.join(sub) for sub in permut]) #should make this into a string


all_files=open_files()
def close_files():
     for file in all_files.values():
          file.close()
          

 



