#!/usr/bin/env python

# Author: Claire Wells <optional@email.address>

# Check out some Python module resources:
#   - https://docs.python.org/3/tutorial/modules.html
#   - https://python101.pythonlibrary.org/chapter36_creating_modules_and_packages.html
#   - and many more: https://www.google.com/search?q=how+to+write+a+python+module

'''This module is a collection of useful bioinformatics functions
written during the Bioinformatics and Genomics Program coursework.
You should update this docstring to reflect what you would like it to say'''

__version__ = "0.3"         # Read way more about versioning here:
                            # https://en.wikipedia.org/wiki/Software_versioning

DNA_bases = {'A', 'T', 'C','G','N'}
RNA_bases =  {'A', 'U', 'G', 'C','N'}

def convert_phred(letter: str) -> int:
    '''Converts a single character into a phred score'''
    return ord(letter)-33

def qual_score(phred_score: str) -> float:
    """Calculate avg qual score of whole phred string.#score was a letter, 
    convert_phred turns it into an integer, 
    the output was a number but the number wasn't stored anywhere"""
    sum = 0
    for score in phred_score:
        score = convert_phred(score) 
        sum = score+sum
        average = sum/len(phred_score)
    return average  
        
def validate_base_seq(seq: str, RNAflag: bool=False) -> bool:
    '''This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive.'''
    seq = seq.upper()     #make the sequence all uppercase
    seq = set(seq)        #make the sequence a set so that we can compare a set to a set 

    if RNAflag:                       #If it passes the RNA flag (so if it's RNA), validate base_seq to RNA 
        return set.issubset(seq,RNA_bases) #.issubset reports whether another set contains this set 
                                        #is seq contained within RNA is what set.issubset(seq,RNA) is saying 
    else: 
        return set.issubset(seq,DNA_bases)    


    '''This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive.'''
    seq = seq.upper()
    return len(seq) == seq.count("A") + seq.count("U" if RNAflag else "T") + seq.count("G") + seq.count("C")
    
def gc_content(DNA:str) -> float:
    assert validate_base_seq(DNA), "String contains invalid characters!"
    '''Returns GC content of a DNA or RNA sequence as a decimal between 0 and 1.'''
    DNA = DNA.upper()  #change to all upper 
    return (DNA.count('G')  + DNA.count('C'))/len(DNA)
    

def calc_median(lst: list) -> float:
    '''Given a sorted list, returns the median value of the list'''
    if len(lst) %2 == 0:
        a = lst[int(len(lst)/2)]
        b = lst[int(((len(lst)/2))-1)]
        return (a+b)/2
    else:
        return lst[int(((len(lst)+1)/2)-1)]


def oneline_fasta(filename) -> str:
    '''This function makes multiple fasta sequence lines into one really long line'''
    with open(filename, "r") as fh, open("first_output.fa", "w") as wf:
        first_time =  True
    for line in fh:
      if first_time:
            wf.write(line)
            first_time = False
      elif ">" in line:
            wf.write(f'\n{line}')
      else:
           wf.write(line.strip("\n"))
    return "first_output.fa"



if __name__ == "__main__":
    # write tests for functions above, Leslie has already populated some tests for convert_phred
    # These tests are run when you execute this file directly (instead of importing it)
    assert convert_phred("I") == 40, "wrong phred score for 'I'"
    assert convert_phred("C") == 34, "wrong phred score for 'C'"
    assert convert_phred("2") == 17, "wrong phred score for '2'"
    assert convert_phred("@") == 31, "wrong phred score for '@'"
    assert convert_phred("$") == 3, "wrong phred score for '$'"
    print("Your convert_phred function is working! Nice job")

#Qual Score Asserts

    assert qual_score("EEE") == 36
    assert qual_score("#I") == 21
    assert qual_score("EJ") == 38.5
    assert qual_score("IC") == 37, "wrong average phred score"
    assert qual_score("8:") == 24
    print("You calcluated the correct average phred score")

#Validate Base Seq Asserts: 
    assert validate_base_seq("AATAGAT", False) == True, "Validate base seq does not work on DNA"
    assert validate_base_seq("AAUAGAU", True) == True, "Validate base seq does not work on RNA"
    assert validate_base_seq("TATUC",False) == False
    assert validate_base_seq("UCUGCU", False) == False
    print("Passed DNA and RNA tests")

#Calc_Median
    assert calc_median([1,2,100]) == 2, "calc_median function does not work for odd length list"
    assert calc_median([1,2]) == 1.5, "calc_median function does not work for even length list"
    assert calc_median([2,4,5]) == 4, "calc_median function does not work for odd length list"
    assert calc_median([12,54,70,120]) == 62, "calc_median function does not work for even length list"
    print("Median successfully calculated")
    
#GC Content 
    assert gc_content("GCGCGC") == 1
    assert gc_content("AATTATA") == 0
    assert gc_content("GCATCGAT") == 0.5
    print("correctly calculated GC content")

