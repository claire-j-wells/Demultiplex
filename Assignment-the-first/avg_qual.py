#!/usr/bin/env python

#USER Notes:

#abspath for files: /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz

#./avg_qual.py -f /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz -len 101 -o Average_Quality_Scores_Read_1

#GOAL: For the file being input, take the average quality score at each given position and generate a distribution
import argparse
import gzip
import bioinfo 
def get_args():
    parser = argparse.ArgumentParser(description="A program to introduce yourself")
    parser.add_argument("-len", help = "len_of_list", type = int)
    parser.add_argument("-f", help = "filename", type = str, required = False)
    parser.add_argument("-o", help = "output_file", required = False)
    return parser.parse_args()
args = get_args()
filename = args.f
len_list = args.len
output = args.o 

#Initialize an empty list of elements: 
#Depending on the number of individual characters in each file, you can initialize your list 
#to be whatever size you want using len_list arg

def init_list(lst: list, value: float=0.0) -> list:
    '''This function takes an empty list and will populate it with
    the value passed in "value". If no value is passed, initializes list
    with 101 values of 0.0.'''
    for inconsequential_value in range(len_list):
        lst.append(0.0)     
    return lst

my_list = []
my_list = init_list(my_list)
#print(my_list)

# def convert_phred(letter: str) -> int:
#     '''Converts a single character into a phred score'''
#     return ord(letter)-33


#print(my_list)
def populate_list(file: str) -> tuple[list, int]:
    with gzip.open(filename, "rt") as fh: 
        i = 0
        for record in fh:
            i+=1
            record = record.strip('\n')
            if i%4 == 0:
                for index, score in enumerate(record):
                    #if index==101:
                        #print(score)    
                    variable = bioinfo.convert_phred(score) #this converts each individual score and assigns it to the variable "variable"
                    my_list[index]+=variable  #this takes the variable recovered and adds it by index from my list
        #print(my_list)
        return(my_list,i)

index_list=[]
mean_value_list = []

my_list, num_lines = populate_list(filename)

for index, my_sum in enumerate(my_list):
    mean_value = my_sum/(num_lines/4)
    my_list[index] = mean_value
    print(f"{index}\t{mean_value}")
    index_list.append(index)
    mean_value_list.append(mean_value)

import matplotlib.pyplot as plt

plt_title = f'{output}'

plt.bar(index_list, mean_value_list,width=0.5, color = 'skyblue')
plt.xlabel("Base Pair (bp)")
plt.ylabel("Mean Value")
plt.title(plt_title)
plt.savefig("test_3.png")



