Note: Refer to written out flowchart of pseudo-code for visual reference
Given: 
---
- Four files and a list of all the known indexes:  

Files and Labels:
| File name | label | 
|---|---|
| 1294_S1_L008_R1_001.fastq.gz | read 1 | 
| 1294_S1_L008_R2_001.fastq.gz | index 1 | 
| 1294_S1_L008_R3_001.fastq.gz | index 2 | 
| 1294_S1_L008_R4_001.fastq.gz | read 2 | 

List of Known Indexes: 

```
B1	GTAGCGTA    A5	CGATCGAT    C1	GATCAAGG
B9	AACAGCGA    C9	TAGCCATG    C3	CGGTAATC
B3	CTCTGGAT    C4	TACCGGAT    A11	CTAGCTCA
C7	CACTTCAC    B2	GCTACTCT    A1	ACGATCAG
B7	TATGGCAC    A3	TGTTCCGT    B4	GTCCTAAG
A12	TCGACAAG    C10	TCTTCGAC    A2	ATCATGCG
C2	ATCGTGGT    A10	TCGAGAGT    B8	TCGGATTC
A7	GATCTTGC    B10	AGAGTCCA    A8	AGGATAGC
```
1. We have R1, R2, R3 and R4. These are the primary files that we are working with. R1 and R4 contain reads 1 and reads 2 and R2 and R3 contain the indexes. 

```
with gzip.open(file, "rt") as fh 
    for line in fh 
    *use modulus to isolate specific sequence lines
```


2. The first thing we need to do is reverse complement the R3 indexes. If the indexes of R2 and R3 are matching, they **will** be the reverse complement of one another so in order to match them and prove they are matching we must write a function that will make these reverse complemented reads "un-reverse complemented". 

```
def reverse_complement_bases()
    input: str of rv complemented bases 
    output: str of bases 
```
   
2a. Between these two steps we also need to consider how the code would work to actually make headers and then attach them to the headers of every file

```
header = print(f'{index} - {rev_comp_index}')
```


3. Add Headers. Once the R3 indexes have been un-reverse complemented, we need to add headers to R1 and R4. The headers might look something like this: `AAA-AAA` The second AAA is the converted form of the original reverse complemented TTT that is originally displayed in the R3 file. 

```
#Goal: Add headers to everything that way it doesn't become an issue downstream!

def index_seq_to_header(): 
    input: f-string header? 
    output: new files with new headers 
```


4. We will now check our indexes against the given available list of indexes using the headers we added to R1 and R4 in the step above. 

```
if header == given_index:
    keep going to check qual_score 
elif header != given_index:
    goes to unknown
elif....(not sure how I would form this statement atm)
    *goes of this statement would be to see if index hopping is occurring

```

5. If an index contains an N or does not match the given list of indexes, the record will be considered unknown. These records will be stored in a file named: `Unk_R1.fq` and `Unk_R2.fq`

6. If R2 index and R3 index (converted) match one another and are in the given list of indexes, it will move on to be checked for its quality score. More on this continued in Step 8. 

7. If the record is in the given list of indexes but DOES NOT match its complement, then this is considered index hopping and will go into a designated folder called `hopped_R1.fq` and `hopped_R2.fq`

8. The records that have matching indexes must be checked for their quality score. If they pass the quality threshold, they will be placed in designated folders called `index_R1.fq` and `index_R2.fq`. These files will be labeled appropriately with the matching index. 

Notes: No QS threshold has been established yet. 

Potential functions to be used/written: 

def index_seq_to_header() <br>
def write_record_to_file() <br>
def calc_qual_score() <br>
    call convert_phred() <br>



