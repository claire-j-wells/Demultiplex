#!/usr/bin/env python

#write out test commands here: 

#./demux.py -R1 ../TEST-input_FASTQ/test_R1.fq.gz -R2 ../TEST-input_FASTQ/test_R2.fq.gz -R3 ../TEST-input_FASTQ/test_R3.fq.gz -R4 ../TEST-input_FASTQ/test_R4.fq.gz -bar /projects/bgmp/shared/2017_sequencing/indexes.txt

#IMPORT ALL NECESSARY TOOLS! 
import gzip
import bioinfo
import itertools
import argparse
import math
import matplotlib.pyplot as plt

#DEFINE ARGS
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

#Establish a dictionary for the rev_comp function to work. 
bases = {"A":"T", 
         "C":"G", 
         "T":"A", 
         "G":"C",
         "N":"N"}

def rev_comp(seq:str) -> str:
        '''This function takes in a string and reverse complements using
        a dictionary. It returns a new reverse complemented string'''
        rev_comp = ""
        seq = seq.strip()
        for nuc in seq: #need add new line to match input
            rev_comp += bases[nuc]
        rev_comp = rev_comp[::-1] 
        return(rev_comp)


#This takes in the barcodes file and cleans it and puts it into a set. 
#We want a set because it is unordered and much faster than a list
known_indexes = set() #establish an empty set 
with open(barcodes,"r") as fh:
    next(fh) #skips the weird designators in indexes (B1, B9, etc stuff)
    for line in fh: 
        new_line = line.strip().split("\t") #strips the new line, splits by tab
        known_indexes.add(new_line[4])


#Open all the files. Accessing files like this so it is NOT in a loop
R1_fh = gzip.open(R1, "rt") 
R2_fh = gzip.open(R2, "rt")
R3_fh = gzip.open(R3, "rt")
R4_fh = gzip.open(R4, "rt")


#Function to open files in the format of a dictionary 
def open_files():
    '''This function establishes a dictionary with file handles as keys and opening 
    and writing as values'''
    all_files = {} #establish all files dictionary 
    for barcodes in known_indexes:
          file_1 = f'outfiles/{barcodes}-{barcodes}_R1.fq' #These two lines establish variables for what we want output files 
          file_2 = f'outfiles/{barcodes}-{barcodes}_R2.fq' #to look like as we write to them
          all_files[barcodes+"R1"] = open(file_1,"w") #Inserting barcodes+"R1" as key in dictionary with opening and writing to file as value
          all_files[barcodes+"R4"] = open(file_2,"w") #Key doesn't matter as much, we just need to know what to call downstream when we want to write 
    all_files["outfiles/unk_R1.fq"] = open("outfiles/unk_R1.fq","w") #Repeated the same thing here but this is hardcoded since we are only outputting 
    all_files["outfiles/unk_R2.fq"] = open("outfiles/unk_R2.fq","w") #an R1/R2 file each for unknown and hopped reads. 
    all_files["outfiles/hopped_R1.fq"] = open("outfiles/hopped_R1.fq","w")
    all_files["outfiles/hopped_R2.fq"] = open("outfiles/hopped_R2.fq","w")
    return(all_files)

all_files = open_files() #call on this to open all the files!

#Isolate the Sequence Line in all the Files 
#ISOLATE ALL RECORDS 
record_4_lines_R1= [] #empty list where we append every record from R1 in to 
record_4_lines_R4 = []  #empty list where we append every record from R4 in to 
counter=0 #to keep track of the lines to write into a file
i=0 #to count the line number
tot_matched = 0 #establish this counter for overall matched statistic
tot_hopped = 0 #establish this counter for overall hopped statistic
tot_unknown = 0 #establish this counter for overall unknown statistic

#ESTABLISH PERMUTATIONS AND DICT TO KEEP TRACK OF INDEX PAIRS
index_dict = {} #dictionary for all the permutations for statistics downstream
for permut in itertools.product(list(known_indexes), repeat=2): #.product produces a list of known indexes that are permuted against each other and themselves
    index_dict[permut] = 0 #set those unique permutations as keys with a value of 0
index_dict[("unknown", "unknown")] = 0  #"unknown" is not a category in the known_indexes list so add this in

    
for line_R1,line_R2,line_R3,line_R4 in zip(R1_fh, R2_fh, R3_fh, R4_fh): #use zip to connect line_R1 to R1.fh etc., be able to read first lines of all files simulatenously
    i+=1 #ensure that modulus works later on. i counts the number of LINES 
    line_R1 = line_R1.strip('\n') #literally stripping the new line character
    line_R2 = line_R2.strip('\n')
    line_R3 = line_R3.strip('\n')
    line_R4 = line_R4.strip('\n')
    #print(i) #left over from checking process. 
    if i%4==2: #Sequence Line
        index1 = line_R2 #assign variables to index lines for R2 and R3 respectively 
        index2 = line_R3
    if i%4==0: #pulling phred score 
        phred1= line_R2 #assign variables to phred score lines for R2 and R3 respectively 
        phred2 = line_R3 #by doing this, it allows us to use these variables later on
    if counter<4: #Because the record is 4 lines, if less than 4 we know that it is still the record and KEEP GOING
        record_4_lines_R1.append(line_R1) #we append the line to make a full record 
        record_4_lines_R4.append(line_R4)
        counter+=1  #counter counts records. Add 1 when less than 4, not a complete record yet. 
        if counter==4: 
            counter=0 #reset everything so we know to STOP appending and restart from scratch. 
            rev_index2 = rev_comp(index2) #reverse complement the second index so that way we can test matching against index 1
            add_to_header = str(index1) +  "-" + rev_index2 #set up the new portion of the header to be added
            new_head_R1 = record_4_lines_R1[0]+' '+ add_to_header #builder new header by calling on the [0] index of list containing individual record
            new_head_R4 = record_4_lines_R4[0]+' '+ add_to_header 
            record_4_lines_R1[0] = new_head_R1 #resetting the older header to just be the previously defined new header 
            record_4_lines_R4[0] = new_head_R1
#ADDING HEADERS:
            #THIS ELIMINATES UNKNOWNS
            if bioinfo.qual_score(phred1) < 26 or bioinfo.qual_score(phred2) < 26 or index1 not in known_indexes or rev_index2 not in known_indexes: #looks at qscore of index files ONLY:
                all_files["outfiles/unk_R1.fq"].write(f'{record_4_lines_R1[0]}\n{record_4_lines_R1[1]}\n{record_4_lines_R1[2]}\n{record_4_lines_R1[3]}\n') 
                all_files["outfiles/unk_R2.fq"].write(f'{record_4_lines_R4[0]}\n{record_4_lines_R4[1]}\n{record_4_lines_R4[2]}\n{record_4_lines_R4[3]}\n')
                index_dict[("unknown","unknown")] += 1 #this counter adds the record to the dictionary made for statistics. It's a record counter
                tot_unknown += 1 #total matched for unknowns
            #THIS LOOKS FOR HOPPED READS BY CHECKING INDEX1 AGAINST THE REVERSE COMPLEMENT OF INDEX 2
            elif index1 != rev_index2:
                all_files["outfiles/hopped_R1.fq"].write(f'{record_4_lines_R1[0]}\n{record_4_lines_R1[1]}\n{record_4_lines_R1[2]}\n{record_4_lines_R1[3]}\n')
                all_files["outfiles/hopped_R2.fq"].write(f'{record_4_lines_R4[0]}\n{record_4_lines_R4[1]}\n{record_4_lines_R4[2]}\n{record_4_lines_R4[3]}\n')
                index_dict[(index1,rev_index2)] += 1 #this counter adds the record to the dictionary made for statistics. It's a record counter
                tot_hopped += 1 #total matched for hopped
            #EVERYTHING THAT WAS MATCHED AKA EVERYTHING ELSE 
            else: 
                all_files[f'{index1}R1'].write(f'{record_4_lines_R1[0]}\n{record_4_lines_R1[1]}\n{record_4_lines_R1[2]}\n{record_4_lines_R1[3]}\n') #VALUE IS IMPORTANT, KEY DOESN'T MATTER 
                all_files[f'{index1}R4'].write(f'{record_4_lines_R4[0]}\n{record_4_lines_R4[1]}\n{record_4_lines_R4[2]}\n{record_4_lines_R4[3]}\n')
                index_dict[(index1,rev_index2)] += 1 #this counter adds the record to the dictionary made for statistics. It's a record counter
                tot_matched += 1 #total matched, regardless of index
            record_4_lines_R1=[] #empty the record at the end of the loop
            record_4_lines_R4 = [] #empty the record at the end of the loop
total_records = i/4 #create a variable for total records outside the loop using i as a line counter for total number of records 
#print(total_records)

#MAKE A PRINTED OUTPUT 
title = f'index pairs\t\tNumber of Reads\tPercent Mapped'
print(title)
for key in index_dict:
    percent = 100*index_dict[key]/total_records
    out = f'{key[0]}\t{key[1]}\t{index_dict[key]}\t{percent}' 
    print(out)     

#WRITE STATISTICS TO FILES
all_percentages = [] #store list of matched percentages as a list outside the loop
matched_indexes = [] #store a list of matching indexes for graph
total_counts = [] #store a list of total raw individual counts outside loop
with open("outputs/output.txt", "w") as out:
    new_tot_unk = 100*tot_unknown/total_records #math for totals 
    new_tot_matched = 100*tot_matched/total_records #math for totals
    new_tot_hopped = 100*tot_hopped/total_records #math for totals
    out.write(f'Summary of Statistics\n') #Write stuff to output file
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
            all_percentages.append(percent) #append the percentages to a list for the plot since math is accessed within a dictionary, don't want to make a new for loop
            total_counts.append(index_dict[keys]) #putting the value of the matched keys into total counts 
        elif keys[0] != keys[1]:
            out.write(f'The percent of mapped reads for {keys[0]}-{keys[1]} is: {percent}% ({index_dict[keys]})\n')

#WRITE TO TSV FILE 
with open("outputs/matched_stats.tsv", "w") as matched:
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
plt.savefig("outputs/output_fig.png")



#CLOSE ALL FILES 
all_files=open_files()
def close_files():
     '''This function closes all files'''
     for file in all_files.values():
          file.close()
          

 



