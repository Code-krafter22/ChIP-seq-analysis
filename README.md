##ChIP-seq Analysis Pipeline

This repository contains SLURM job scripts for a basic ChIP-seq analysis pipeline, including:

Quality control using FastQC

Read alignment using Bowtie2 and Samtools

Peak calling using MACS3

📁 Contents
File	Description
quality_check.sh	Downloads raw FASTQ files and runs FastQC
chipseq_mapping.sh	Aligns reads with Bowtie2 and processes with Samtools
macs3_peak_calling.sh	Calls peaks using MACS3 for YAP1, TAZ, TEAD4
🔧 Requirements
SLURM-managed HPC environment

Software modules or conda environments for:

fastqc

bowtie2

samtools

macs3

🧪 Step-by-Step Workflow
1. 🔬 Quality Control (FastQC) (run : sbatch quality_check.sh)
Script: quality_check.sh
Description: Downloads .fastq.gz files and performs quality check using FastQC.

Output: .html and .zip FastQC reports for each sample in the midterm/ directory.

2. 🧬 Read Alignment (Bowtie2 + Samtools) (run : sbatch chipseq_mapping.sh)
Script: chipseq_mapping.sh
Description: Aligns reads to the reference genome and converts .sam files to sorted and indexed .bam files.

FASTQ files: SRR1810900.fastq.gz, SRR1810907.fastq.gz, SRR1810912.fastq.gz, SRR1810918.fastq.gz
Reference genome: new_refgenome/GRCh38_noalt_as (pre-built Bowtie2 index)

Output:
Aligned .sam and .bam files
Indexed .bam.bai files

3. ⛰️ Peak Calling (MACS3) (run : sbatch macs3_peak_calling.sh)
Script: macs3_peak_calling.sh
Description: Calls ChIP-seq peaks for YAP1, TAZ, and TEAD4 using MACS3 with a shared input control.


Input BAMs:

ChIP: SRR1810900.bam (YAP1), SRR1810907.bam (TAZ), SRR1810918.bam (TEAD4)

Input: SRR1810912.bam

Output:

YAP1_macs3_peak/, TAZ_macs3_peak/, TEAD4_macs3_peak/

Each directory contains peak files (.narrowPeak, .xls, etc.)

📂 Directory Structure
python
Copy
Edit
midterm/
├── *.fastq.gz            # Raw reads
├── *.sam                 # Alignments
├── *.bam, *.bai          # Sorted, indexed BAM files
├── *_macs3_peak/         # Peak calling output for each TF
├── *.html, *.zip         # FastQC reports
✍️ Author

Kashvi Shah
Graduate Student in Bioinformatics
Email: kashvichirag@usf.edu

🪪 License
This project is licensed under the MIT License.
