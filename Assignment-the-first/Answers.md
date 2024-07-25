# Assignment the First

## Part 1
1. Be sure to upload your Python script. Provide a link to it here:

| File name | label | Read length | Phred encoding |
|---|---|---|---|
| 1294_S1_L008_R1_001.fastq.gz |  |  |  |
| 1294_S1_L008_R2_001.fastq.gz |  |  |  |
| 1294_S1_L008_R3_001.fastq.gz |  |  |  |
| 1294_S1_L008_R4_001.fastq.gz |  |  |  |

2. Per-base NT distribution
    1. Use markdown to insert your 4 histograms here.
    2. **YOUR ANSWER HERE**
    3. **YOUR ANSWER HERE**
    
## Part 2
1. Define the problem
There are a few problems. We are given four fastq files (R1,R2,R3,R4). R1 and R4 are Read 1 and Read 2 and R2 and R3 are Index 1 and Index 2 respectively. These indexes are used to match the reads in each sample. We are also given a list of 24 indexes that are known and used to generate the RNA sequence data. 

The end goal of this is to demultiplex the data. This means we need to get the reads from R1 and R4 and get those reads into their own designated files labeled by matching index (ex: AAAA-AAAA_R2.fq). R1 will go into an index_R1.fq file, this will be a collection of records from R1 only and R2 will go into another seperate file of index_R2.fq. There will be 24 different files times 2 because there is both R1 and R4 files to parse through for a total of 48 files. 

In addition, we will also have a unk_R1.fq and a unk_R2.fq folder. In this folder, we will have records that do NOT match against the set of 24 given indexes. This includes anything that has an "N" in it. IF the indexes are identical, then they will go into the index_R2.fq and index_R1.fq files mentioned above. IF the indexes are not identical BUT did match with the set of 24 given indexes in some way, then they will go into it's own designated folder of index hopping called hopped_R1.fq and hopped_R2.fq. 

2. Describe output <br>
For the output, we will have 52 total files. <br>
For the files that have matching indexes to the given list of indexes, they will have their own files designated as `index_R1.fq` and `index_R2.fq`. In place of index, the matched up indexes will be the title of these FQ files. 

We will have another set of files for indexes that are present in our list of 24 but do *NOT* match oone another. This will be designated as index hopping and will be in a designated `hopped_R1.fq` and `hopped_R2.fq` folder. Finally, we have files for unknown 

3. Upload your [4 input FASTQ files](../TEST-input_FASTQ) and your [>=6 expected output FASTQ files](../TEST-output_FASTQ).



4. Pseudocode
5. High level functions. For each function, be sure to include:
    1. Description/doc string
    2. Function headers (name and parameters)
    3. Test examples for individual functions
    4. Return statement
