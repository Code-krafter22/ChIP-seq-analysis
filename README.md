# 🔬 ChIP-seq Pipeline: From Raw Reads to Peaks

Welcome to my ChIP-seq analysis pipeline! This project processes ChIP-seq data for **YAP1**, **TAZ**, and **TEAD4**, moving from raw sequencing files to peak discovery using standard bioinformatics tools on a SLURM HPC cluster.

> 🧠 Mapping molecules to meaning  
> ☕ Fueled by science, shell scripts, and curiosity

---

## 🧪 Pipeline Overview

- 🔍 **Quality Check**: Assess read quality with `FastQC`
- 🧷 **Alignment**: Map reads to the human genome using `Bowtie2` & sort/index with `Samtools`
- 🏔️ **Peak Calling**: Identify enriched regions using `MACS3`

---

## 📁 Included Scripts

| Script                    | Function              |
|---------------------------|-----------------------|
| `quality_check.sh`        | Downloads FASTQ files and runs FastQC |
| `chipseq_mapping.sh`      | Aligns reads and prepares sorted/indexed BAMs |
| `macs3_peak_calling.sh`   | Calls peaks for YAP1, TAZ, and TEAD4 using MACS3 |

---

## 🚀 Quickstart

```bash
# Submit jobs to SLURM
sbatch quality_check.sh
sbatch chipseq_mapping.sh
sbatch macs3_peak_calling.sh
