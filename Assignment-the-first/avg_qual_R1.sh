#!/bin/bash

#SBATCH --job-name=R1
#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --cpus-per-task=4
#SBATCH -c 8
#SBATCH --mem=100G
#SBATCH --nodes=1
#SBATCH --output=slurm-R1-%j.out
#SBATCH --error=slurm-R1-%j.err
#SBATCH --mail-user=cwell@uoregon.edu
#SBATCH --mail-type=ALL


conda activate bgmp_py312

./avg_qual.py -f /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz -len 101 -o Average_Quality_Scores_Read_1

exit 