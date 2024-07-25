
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

This number matches Leslie's number! I am choosing not to run the other files at the moment for the sake of time but all files should match. 

`zcat 1294_S1_L008_R4_001.fastq.gz | sed -n '2~4p' | awk '{print length($0)}' | head -4`

output: 101

`zcat 1294_S1_L008_R2_001.fastq.gz | sed -n '2~4p' | awk '{print length($0)}' | head -4`

output: 8 

- Because the index reads and the sequence reads are the same length, we did not run every single file, only one each


| File name | label | Read length | Phred encoding |
|---|---|---|---|
| 1294_S1_L008_R1_001.fastq.gz | read 1 | 101 | phred+33 |
| 1294_S1_L008_R2_001.fastq.gz | index 1 | 8 | phred+33 |
| 1294_S1_L008_R3_001.fastq.gz | index 2 | 8 | phred+33 |
| 1294_S1_L008_R4_001.fastq.gz | read 2 | 101 | phred+33 |


phred+33 because we can see symbols in the phred scores. If there's letters ONLY then it's phred+64 and if there are symbols then we know it's phred+33. 





