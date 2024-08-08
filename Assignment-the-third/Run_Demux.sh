#!/bin/bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH -c 8
#SBATCH --mem=100G
#SBATCH --time=0-3
#SBATCH --job-name=Demux
#SBATCH --mail-user=cwell@uoregon.edu     #optional: if you'd like email
#SBATCH --mail-type=ALL                   #optional: must set email first, what type of email you wan
#SBATCH --output=out/slurm-special-name-%j.out
#SBATCH --error=out_err/slurm-special-name-%j.err



conda activate bgmp_py312


/usr/bin/time -v ./demux.py -R1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz \
 -R2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz \
 -R3 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz \
 -R4 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz \
 -bar /projects/bgmp/shared/2017_sequencing/indexes.txt

exit 