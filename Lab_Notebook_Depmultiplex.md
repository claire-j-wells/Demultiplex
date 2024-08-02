
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





