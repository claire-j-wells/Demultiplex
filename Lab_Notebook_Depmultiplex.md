
Thursday, July 25th, 2024: I think I'm gonna be sick
--- 

- Starting part 1 of the Demultiplex assignment. I put made a copy of the demultiplex assignment and put it in my own github repo. From there, I git cloned into `Bi622/PS/` folder
- DO NOT COPY OR MOVE ANY OF THE FILES! THEY ARE HUGE! 

- Now performing some preliminary data cleaning 

- `ls -lah' in order to see the general size of the files of interest. Results pasted below: 

`-rw-r-xr--+  1 coonrod  is.racs.pirg.bgmp  20G Jul 30  2018 1294_S1_L008_R1_001.fastq.gz <br>
-rw-r-xr--+  1 coonrod  is.racs.pirg.bgmp 2.6G Jul 30  2018 1294_S1_L008_R2_001.fastq.gz <br>
-rw-r-xr--+  1 coonrod  is.racs.pirg.bgmp 2.8G Jul 30  2018 1294_S1_L008_R3_001.fastq.gz <br>
-rw-r-xr--+  1 coonrod  is.racs.pirg.bgmp  21G Jul 30  2018 1294_S1_L008_R4_001.fastq.gz <br>`

- `zcat 1294_S1_L008_R1_001.fastq.gz | wc -l ` in order to check the number of lines: 

- output:`1452986940`

This number matches Leslie's number! I am choosing not to run the other files at the moment for the sake of time  

`zcat 1294_S1_L008_R4_001.fastq.gz | sed -n '2~4p' | awk '{print length($0)}' | head -4`

output: 101

`zcat 1294_S1_L008_R1_001.fastq.gz | sed -n '2~4p' | awk '{print length($0)}' | head -4`

output: 101

`zcat 1294_S1_L008_R2_001.fastq.gz | sed -n '2~4p' | awk '{print length($0)}' | head -4`

output: 8 

`zcat 1294_S1_L008_R2_001.fastq.gz | sed -n '2~4p' | awk '{print length($0)}' | head -4`

output: 8 


| File name | label | Read length | Phred encoding |
|---|---|---|---|
| 1294_S1_L008_R1_001.fastq.gz | read 1 | 101 | phred+33 |
| 1294_S1_L008_R2_001.fastq.gz | index 1 | 8 | phred+33 |
| 1294_S1_L008_R3_001.fastq.gz | index 2 | 8 | phred+33 |
| 1294_S1_L008_R4_001.fastq.gz | read 2 | 101 | phred+33 |


phred+33 because we can see symbols in the phred scores. If there's letters ONLY then it's phred+64 and if there are symbols then we know it's phred+33. 


July 26th, 2024
----
sbatch avg_qual_R1.sh <br>
Submitted batch job 7641698 <br>
Submitted new batch job 7641704 <br>
sbatch avg_qual_R2.sh <br>
Submitted batch job 7639401 <br>

sbatch avg_qual_R3.sh <br>
Submitted batch job 7641697 <br>
^this batch job above failed - this failed following changing the "savefig" command to be the same name as the title to prevent overwriting of png outputs. Now rerunning JUST R1 and will later rerun R3 after I modified savefig to output a hardcoded image name.  
Submitted batch job 7641710
Submitted batch job 7651636
sbatch avg_qual_R4.sh <br>
Submitted batch job 7639403 <br>

All of the submitted runs returned with a exit status of 0. Unfortunately I forgot to add the command that counts time on slurm scripts, however, from memory and based on email notifs I received, Read 1 and Read 4 took the longest being the sequence reads and were about 2-3 hours and Read 2 and 3 were significantly shorter at approximately 20 minutes or so. 

July 27th, 2024
---

So I was having some issues with my bioinfo module. For some reason, pylance is really just not like that I'm importing the module. I tried to run my code `avg_qual.py` using the module and it didn't work but when I just paste in the convert_phred() function and run the code it does work just fine. Likely because it just can't access what it needs to. 

Resolved: Literally didn't `chmod` the function. Yikes. Won't forget next time! Ran code again on test and it successfully ran! Hooray! 

August 1st, 2024 - GRIND! Start Part 3!
---
So today we started the actual `demux.py` script. First thing I did was set up an area to import any necessary packages. In the final script I did use `argparse` but for now I'm just hardcoding everything for the sake of ease when running things. 

I made some test files from the original files using : 

```
zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz | head -100000 > R1_test.fq
```
For each file, I modified the command from R1 to be R2, R3, R4 etc. 

I had been working on this previously but I was trying to design my reverse complement function. This function was super close to being complete but got Leslie's help with formating indents so that the loop wouldn't flip flop the base as it was being appended to the line. This function returns the new string of DNA reverse complemented. 

Next we needed to make the list of indexes a set because sets are much easier to work with than lists. 

Open all files using `gzip.open`. We need to us `gzip` because the actual files are gzipped and we used it with the `"rt"` function because we need the file to be readable `r` and in text mode `t`. 

We are NOT doing this in a for loop because it takes up too much to open the file and read through the entire thing every time we run this script. Therefore we need to think of an efficient way to go about this that doesn't involve a for-loop. The solution to this is using a dictionary with some sort of file handle as a key and the value being the opening and writing to that file if given that key. We decided to write this as a function and we were able to call this function at the beginning of the script. It's important to remember to CLOSE the files though at the end. 

The general workflow for this entire script was to establish all my variables and then once those variable were established, I could call on them later on in my if-statements. I set up a for loop that runs through all the lines from all my files  ***concurrently*** by using zip. I wanted to work with an individual record EVERY TIME and not the entire file. This helps save time and memory. 

See actual script for in-detail annotations for every line. 

Saturday August 3rd, 2024
---

Once I got the bulk of the sorting completed, we needed to do statistics! Yay....

The main goal was to get the percent of mapped reads for matching indexes. First off, start doing this by making a dictionary of all the possible permutations of the known indexes. I used `permutations.products` because not only did it list possible permutations across indexes but it also listed possible matched permutations. I set these combinations at keys in a dictionary and set the value to 0. I also added an "(unknown, unknown)" key set to 0 for the unknown category. Listed as a tuple because the output of itertools is a tuple object. 

See notes on actual script for what I did to actually count all the records. 

Monday August 5th, 2024
---

Had to do some finishing touches on the outputs. I made a horizontal bar plot of matched indexes using `plt.barh`. 
Note: `plt.tight_layout()` was a miracle in preventing cutoff of indexes on the horizontal x-axis. Not exactly sure why it works but noted for future horizontal bar graphs. 

Also for outputs I decided to make a .tsv summarizing actual raw counts of records for each of the matched indexes as well as their percentages. 

`output.txt` file also has this information as well as a summary of total counts. 

Tuesday August 6th, 2024
---

So I did a final check just to make sure that all the finishing touches were in in terms of actually running the code as a whole so I wanted to run my unit tests from Part 1 to make sure everything was working okay and as expected. 

I ran into my first problem right off the bat when I realized that I made up my own indexes and as a result I would either need to hardcode a new set of indexes and adapt it to be similar to the given txt file ORRRRRR I could just modify my unit tests to have new barcodes that did actually match the ones that Leslie gave us. I decided to change the barcodes. Turns out this was a helpful processb because there were some errors in my unit tests that were a result of being confused back in part one so by redoing it things made a lot more sense. 

I ran my unit tests and for some reason all of my files actually went into the unknown folder. After getting some help from Varsheni, we realized that it came from a bunch of inviisble spaces that were left over from copy pasting and moving things around. Would not have figured that one out. It's good to keep in mind for the future. 

I ran them again thinking that they would be working soundly and well I was wrong. I got confused with the R3 index because I should have put in reverse complement of an item in the given list so that way once it was reverse complemented, the script would recognize it as a thing in the list but NOT matching index 1 which would result in that record being hopped. 

Once I fixed this bug, the unit tests worked out as planned. For some reason it was very difficult for me to wrap my head around the unit test aspect. It was like thinking backwards and from many different perspectives since we have so many files and pieces moving around. 

I did a final run through and wrote my slurm script. Also added in argparse so I could use the real files. 

Slurm ID: 7996977

```
Command being timed: "./demux.py -R1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz -R2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz -R3 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz -R4 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz -bar /projects/bgmp/shared/2017_sequencing/indexes.txt"
User time (seconds): 4471.00
System time (seconds): 60.35
Percent of CPU this job got: 88%
Elapsed (wall clock) time (h:mm:ss or m:ss): 1:24:55
```

Three outputs of `demux.py` script:

output_fig.png (This is an output graph visualizing the percent of mapped reads across matching indexes)
matched_stats.tsv (This is a tsv file that has the matched index, the raw counts for records and the actual percent)
output.txt (This text file has a summary of stats)




