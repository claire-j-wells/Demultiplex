#!/usr/bin/env python

barcodes_path = "/projects/bgmp/shared/2017_sequencing/indexes.txt"


#ADD LATER:

R1 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz"
R2 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz"
R3 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz"
R4 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz"



#import numpy as np
import gzip
import bioinfo

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
        new_line = line.strip.split("\t") #strips the new line, splits by tab
        known_indexes.add(new_line[4])


#Open all the files, use the R1 file handles 
R1_fh = gzip.open(R1, "r") #change the R1-R4 variables to be adapted to argparse
R2_fh = gzip.open(R2, "r")
R3_fh = gzip.open(R3, "r")
R4_fh = gzip.open(R4, "r")


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
i = 0
record_4_lines_R1= ""
record_4_lines_R4 = ""
counter=0 #to keep track of the lines to write into a file
i=0
for line_R1,line_R2,line_R3,line_R4 in zip(R1_fh, R2_fh, R3_fh, R4_fh):
    i+=1
    line_R1 = line.strip('\n')
    line_R2 = line.strip('\n')
    line_R3 = line.strip('\n')
    line_R4 = line.strip('\n')
    if counter<4: #Because the record is 4 lines, if less than 4 we know that it is still the record
        record_4_lines_R1 += line_R1 + "/n" #we append the line
        record_4_lines_R4 += line_R4 + "/n"

        counter+=1
        if counter==4: 
            counter=0 #reset everything
            rev_index2 = rev_comp(index2)
            add_to_header = str(index1) +  "-" + rev_index2
#ADDING HEADERS:




            # record_copy_R1 = record_4_lines_R1 
            # record_copy_R4 = record_4_lines_R4
        
    
            #THIS ELIMINATES UNKNOWNS
            if bioinfo.qual_score(phred1) < 26 or bioinfo.qual_score(phred2) < 26 or index1 not in known_indexes or rev_index2 not in known_indexes #looks at qscore of index files ONLY 
                all_files["outfiles/unk_R1.fq"].write(record_4_lines_R1) 
                all_files["outfiles/unk_R2.fq"].write(record_4_lines_R4)
            elif index1 != rev_index2:
                all_files["outfiles/hopped_R1.fq"].write(record_4_lines_R1)
                all_files["outfiles/hopped_R2.fq"].write(record_4_lines_R4)
            else: 
                all_files[f'{index1}+"R1"'].write(record_4_lines_R1) #VALUE IS IMPORTANT, KEY DOESN'T MATTER 
                all_files[f'{index1}+"R4"'].write(record_4_lines_R4)
                
            record_4_lines_R1=""
            record_4_lines_R4 = ""
    
                #print("record_copy")
        if i%4==2: #Sequence Line
            index1 = line_R2
            index2 = line_R3
    
    if i%4==4: #pulling phred score 
         phred1= line_R2
         phred2 = line_R3
         

         
       
    

    



#Reverse Complement R3 Index (sequence line) in the Header and Reverse the Quality Too 

#Check if the indexes match or not: barcode or unknown NOT hopped 

#if the indexes don't match: hopped or unknown 

#Write Function that Adds the Reverse Complement to the Header 






all_files=open_files()
def close_files():
     for file in all_files.values():
          file.close()
          

 



