#!/usr/bin/env python

barcodes_path = "/projects/bgmp/shared/2017_sequencing/indexes.txt"

import gzip
import bioinfo
import itertools
import argparse

def get_args():
    parser = argparse.ArgumentParser(description="A program to introduce yourself")
    parser.add_argument("-R1", help = "Read 1: Biological Read", required = True)
    parser.add_argument("-R2", help = "Read 2: Index 1", required = True)
    parser.add_argument("-R3", help = "Read 3: Index 2", required = True)
    parser.add_argument("-R4", help = "Read 4: Biological Read", required = True)
    parser.add_argument("-index", help = "Known Indexes", required = True)
    return parser.parse_args()



args = get_args()
R1 = args.r1
R2 = args.r2
R3 = args.r3
R4 = args.r4
index = args.index


#ADD LATER:

R1 = "R1_test.fq.gz" #"/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz" 
R2 = "R2_test.fq.gz" #"/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz"
R3 = "R3_test.fq.gz" #"/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz"
R4 = "R4_test.fq.gz" #"/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz"






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


 #Make Set of Indexes: 

# #barcodes = {"GTAGCGTA","CGATCGAT","GATCAAGG"
# "AACAGCGA", "TAGCCATG", "CGGTAATC"
# "CTCTGGAT",    	"TACCGGAT",    	"CTAGCTCA"
# "CACTTCAC",    	"GCTACTCT"   	"ACGATCAG"
# "TATGGCAC"   	"TGTTCCGT"   	"GTCCTAAG"
# "TCGACAAG", "TCTTCGAC"    	"ATCATGCG"
# "ATCGTGGT"    "TCGAGAGT"    	"TCGGATTC"
# "GATCTTGC",    "AGAGTCCA"   	"AGGATAGC"}


known_indexes = set() #establish an empty set 
with open(barcodes_path,"r") as fh:
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
          file_1 = f'outfiles/{barcodes}-{barcodes}_R1.fq' #KEY IS PATH TO WHERE SHIT GOES 
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
            elif index1 != rev_index2:
                all_files["outfiles/hopped_R1.fq"].write(f'{record_4_lines_R1[0]}\n{record_4_lines_R1[1]}\n{record_4_lines_R1[2]}\n{record_4_lines_R1[3]}\n')
                all_files["outfiles/hopped_R2.fq"].write(f'{record_4_lines_R4[0]}\n{record_4_lines_R4[1]}\n{record_4_lines_R4[2]}\n{record_4_lines_R4[3]}\n')
            #EVERYTHING THAT WAS MATCHED
            else: 
                all_files[f'{index1}R1'].write(f'{record_4_lines_R1[0]}\n{record_4_lines_R1[1]}\n{record_4_lines_R1[2]}\n{record_4_lines_R1[3]}\n') #VALUE IS IMPORTANT, KEY DOESN'T MATTER 
                all_files[f'{index1}R4'].write(f'{record_4_lines_R4[0]}\n{record_4_lines_R4[1]}\n{record_4_lines_R4[2]}\n{record_4_lines_R4[3]}\n')
                
            record_4_lines_R1=[]
            record_4_lines_R4 = []
    
    #             #print("record_copy")
    # if i%4==2: #Sequence Line
    #     index1 = line_R2
    #     index2 = line_R3
    # if i%4==0: #pulling phred score 
    #     phred1= line_R2
    #     phred2 = line_R3
         

         
       
    

    



#Reverse Complement R3 Index (sequence line) in the Header and Reverse the Quality Too 

#Check if the indexes match or not: barcode or unknown NOT hopped 

#if the indexes don't match: hopped or unknown 

#Write Function that Adds the Reverse Complement to the Header 


permut = list(itertools.product(['A', 'C', 'T','G'], repeat=8)) #gives us the number of possible ACTG index combinations, set to 8 bc indexes are 8 char
print(len(permut))
# permut = ', '.join([' '.join(sub) for sub in permut]) #should make this into a string


all_files=open_files()
def close_files():
     for file in all_files.values():
          file.close()
          

 



