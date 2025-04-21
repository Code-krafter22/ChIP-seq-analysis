#!/bin/bash

#SBATCH --chdir=/home/k/kashvichirag/midterm
#SBATCH --job-name=bigwig
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=10000
#SBATCH -t 02:00:00
#SBATCH -o run.out
#SBATCH -e run.err
#SBATCH --mail-user=kashvichirag@usf.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

source ~/.bash_profile
module purge
conda activate deeptools_env

bamCoverage --bam SRR1810900.bam --normalizeUsing RPKM --extendReads 200 -o SRR1810900.bw
bamCoverage --bam SRR1810907.bam --normalizeUsing RPKM --extendReads 200 -o SRR1810907.bw
bamCoverage --bam SRR1810912.bam --normalizeUsing RPKM --extendReads 200 -o SRR1810912.bw
bamCoverage --bam SRR1810918.bam --normalizeUsing RPKM --extendReads 200 -o SRR1810918.bw
