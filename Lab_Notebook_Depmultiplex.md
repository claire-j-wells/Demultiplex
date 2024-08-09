
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


Conda Environment Information:
# packages in environment at /projects/bgmp/cwell/miniforge3/envs/bgmp_py312:
#
# Name                    Version                   Build  Channel
_libgcc_mutex             0.1                 conda_forge    conda-forge
_openmp_mutex             4.5                       2_gnu    conda-forge
alsa-lib                  1.2.12               h4ab18f5_0    conda-forge
attr                      2.5.1                h166bdaf_1    conda-forge
brotli                    1.1.0                hd590300_1    conda-forge
brotli-bin                1.1.0                hd590300_1    conda-forge
bzip2                     1.0.8                hd590300_5    conda-forge
ca-certificates           2024.7.4             hbcca054_0    conda-forge
cairo                     1.18.0               hbb29018_2    conda-forge
certifi                   2024.7.4           pyhd8ed1ab_0    conda-forge
contourpy                 1.2.1           py312h8572e83_0    conda-forge
cycler                    0.12.1             pyhd8ed1ab_0    conda-forge
dbus                      1.13.6               h5008d03_3    conda-forge
expat                     2.6.2                h59595ed_0    conda-forge
font-ttf-dejavu-sans-mono 2.37                 hab24e00_0    conda-forge
font-ttf-inconsolata      3.000                h77eed37_0    conda-forge
font-ttf-source-code-pro  2.038                h77eed37_0    conda-forge
font-ttf-ubuntu           0.83                 h77eed37_2    conda-forge
fontconfig                2.14.2               h14ed4e7_0    conda-forge
fonts-conda-ecosystem     1                             0    conda-forge
fonts-conda-forge         1                             0    conda-forge
fonttools                 4.53.1          py312h41a817b_0    conda-forge
freetype                  2.12.1               h267a509_2    conda-forge
gettext                   0.22.5               h59595ed_2    conda-forge
gettext-tools             0.22.5               h59595ed_2    conda-forge
glib                      2.80.3               h8a4344b_1    conda-forge
glib-tools                2.80.3               h73ef956_1    conda-forge
graphite2                 1.3.13            h59595ed_1003    conda-forge
gst-plugins-base          1.24.5               hbaaba92_0    conda-forge
gstreamer                 1.24.5               haf2f30d_0    conda-forge
harfbuzz                  8.5.0                hfac3d4d_0    conda-forge
icu                       73.2                 h59595ed_0    conda-forge
keyutils                  1.6.1                h166bdaf_0    conda-forge
kiwisolver                1.4.5           py312h8572e83_1    conda-forge
krb5                      1.21.3               h659f571_0    conda-forge
lame                      3.100             h166bdaf_1003    conda-forge
lcms2                     2.16                 hb7c19ff_0    conda-forge
ld_impl_linux-64          2.40                 hf3520f5_7    conda-forge
lerc                      4.0.0                h27087fc_0    conda-forge
libasprintf               0.22.5               h661eb56_2    conda-forge
libasprintf-devel         0.22.5               h661eb56_2    conda-forge
libblas                   3.9.0           22_linux64_openblas    conda-forge
libbrotlicommon           1.1.0                hd590300_1    conda-forge
libbrotlidec              1.1.0                hd590300_1    conda-forge
libbrotlienc              1.1.0                hd590300_1    conda-forge
libcap                    2.69                 h0f662aa_0    conda-forge
libcblas                  3.9.0           22_linux64_openblas    conda-forge
libclang-cpp15            15.0.7          default_h127d8a8_5    conda-forge
libclang13                18.1.8          default_h6ae225f_0    conda-forge
libcups                   2.3.3                h4637d8d_4    conda-forge
libdeflate                1.20                 hd590300_0    conda-forge
libedit                   3.1.20191231         he28a2e2_2    conda-forge
libevent                  2.1.12               hf998b51_1    conda-forge
libexpat                  2.6.2                h59595ed_0    conda-forge
libffi                    3.4.2                h7f98852_5    conda-forge
libflac                   1.4.3                h59595ed_0    conda-forge
libgcc-ng                 14.1.0               h77fa898_0    conda-forge
libgcrypt                 1.11.0               h4ab18f5_0    conda-forge
libgettextpo              0.22.5               h59595ed_2    conda-forge
libgettextpo-devel        0.22.5               h59595ed_2    conda-forge
libgfortran-ng            14.1.0               h69a702a_0    conda-forge
libgfortran5              14.1.0               hc5f4f2c_0    conda-forge
libglib                   2.80.3               h8a4344b_1    conda-forge
libgomp                   14.1.0               h77fa898_0    conda-forge
libgpg-error              1.50                 h4f305b6_0    conda-forge
libiconv                  1.17                 hd590300_2    conda-forge
libjpeg-turbo             3.0.0                hd590300_1    conda-forge
liblapack                 3.9.0           22_linux64_openblas    conda-forge
libllvm15                 15.0.7               hb3ce162_4    conda-forge
libllvm18                 18.1.8               hc9dba70_0    conda-forge
libnsl                    2.0.1                hd590300_0    conda-forge
libogg                    1.3.5                h4ab18f5_0    conda-forge
libopenblas               0.3.27          pthreads_hac2b453_1    conda-forge
libopus                   1.3.1                h7f98852_1    conda-forge
libpng                    1.6.43               h2797004_0    conda-forge
libpq                     16.3                 ha72fbe1_0    conda-forge
libsndfile                1.2.2                hc60ed4a_1    conda-forge
libsqlite                 3.46.0               hde9e2c9_0    conda-forge
libstdcxx-ng              14.1.0               hc0a3c3a_0    conda-forge
libsystemd0               255                  h3516f8a_1    conda-forge
libtiff                   4.6.0                h1dd3fc0_3    conda-forge
libuuid                   2.38.1               h0b41bf4_0    conda-forge
libvorbis                 1.3.7                h9c3ff4c_0    conda-forge
libwebp-base              1.4.0                hd590300_0    conda-forge
libxcb                    1.16                 hd590300_0    conda-forge
libxcrypt                 4.4.36               hd590300_1    conda-forge
libxkbcommon              1.7.0                h2c5496b_1    conda-forge
libxml2                   2.12.7               h4c95cb1_3    conda-forge
libzlib                   1.3.1                h4ab18f5_1    conda-forge
lz4-c                     1.9.4                hcb278e6_0    conda-forge
matplotlib                3.9.1           py312h7900ff3_0    conda-forge
matplotlib-base           3.9.1           py312h9201f00_0    conda-forge
mpg123                    1.32.6               h59595ed_0    conda-forge
munkres                   1.1.4              pyh9f0ad1d_0    conda-forge
mysql-common              8.3.0                hf1915f5_4    conda-forge
mysql-libs                8.3.0                hca2cd23_4    conda-forge
ncurses                   6.5                  h59595ed_0    conda-forge
nspr                      4.35                 h27087fc_0    conda-forge
nss                       3.102                h593d115_0    conda-forge
numpy                     2.0.0           py312h22e1c76_0    conda-forge
openjpeg                  2.5.2                h488ebb8_0    conda-forge
openssl                   3.3.1                h4ab18f5_1    conda-forge
packaging                 24.1               pyhd8ed1ab_0    conda-forge
pcre2                     10.44                h0f59acf_0    conda-forge
pillow                    10.4.0          py312h287a98d_0    conda-forge
pip                       24.0               pyhd8ed1ab_0    conda-forge
pixman                    0.43.2               h59595ed_0    conda-forge
ply                       3.11               pyhd8ed1ab_2    conda-forge
pthread-stubs             0.4               h36c2ea0_1001    conda-forge
pulseaudio-client         17.0                 hb77b528_0    conda-forge
pyparsing                 3.1.2              pyhd8ed1ab_0    conda-forge
pyqt                      5.15.9          py312h949fe66_5    conda-forge
pyqt5-sip                 12.12.2         py312h30efb56_5    conda-forge
python                    3.12.4          h194c7f8_0_cpython    conda-forge
python-dateutil           2.9.0              pyhd8ed1ab_0    conda-forge
python_abi                3.12                    4_cp312    conda-forge
qhull                     2020.2               h434a139_5    conda-forge
qt-main                   5.15.8              ha2b5568_22    conda-forge
readline                  8.2                  h8228510_1    conda-forge
setuptools                70.1.1             pyhd8ed1ab_0    conda-forge
sip                       6.7.12          py312h30efb56_0    conda-forge
six                       1.16.0             pyh6c4a22f_0    conda-forge
tk                        8.6.13          noxft_h4845f30_101    conda-forge
toml                      0.10.2             pyhd8ed1ab_0    conda-forge
tomli                     2.0.1              pyhd8ed1ab_0    conda-forge
tornado                   6.4.1           py312h9a8786e_0    conda-forge
tzdata                    2024a                h0c530f3_0    conda-forge
wheel                     0.43.0             pyhd8ed1ab_1    conda-forge
xcb-util                  0.4.1                hb711507_2    conda-forge
xcb-util-image            0.4.0                hb711507_2    conda-forge
xcb-util-keysyms          0.4.1                hb711507_0    conda-forge
xcb-util-renderutil       0.3.10               hb711507_0    conda-forge
xcb-util-wm               0.4.2                hb711507_0    conda-forge
xkeyboard-config          2.42                 h4ab18f5_0    conda-forge
xorg-kbproto              1.0.7             h7f98852_1002    conda-forge
xorg-libice               1.1.1                hd590300_0    conda-forge
xorg-libsm                1.2.4                h7391055_0    conda-forge
xorg-libx11               1.8.9                hb711507_1    conda-forge
xorg-libxau               1.0.11               hd590300_0    conda-forge
xorg-libxdmcp             1.1.3                h7f98852_0    conda-forge
xorg-libxext              1.3.4                h0b41bf4_2    conda-forge
xorg-libxrender           0.9.11               hd590300_0    conda-forge
xorg-renderproto          0.11.1            h7f98852_1002    conda-forge
xorg-xextproto            7.3.0             h0b41bf4_1003    conda-forge
xorg-xf86vidmodeproto     2.3.1             h7f98852_1002    conda-forge
xorg-xproto               7.0.31            h7f98852_1007    conda-forge
xz                        5.2.6                h166bdaf_0    conda-forge
zlib                      1.3.1                h4ab18f5_1    conda-forge
zstd                      1.5.6                ha6fb4c9_0    conda-forge


