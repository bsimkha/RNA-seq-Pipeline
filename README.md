# RNA-seq pipeline
This pipeline was built in **Nextflow DSL2**, and performs end-to-end processing of pair-ended **RNA-seq FASTQ files**. Processing include trimming, rRNA filtering, QC, alignment, and gene quantification

## Overview
This pipeline automates the following steps:

**Read trimming** — `fastp`    
**rRNA removal** — `bbmap`   
**Quality control** — `FastQC`   
**Alignment** — `STAR`   
**Quantification** — `featureCounts`   
**Count merging** — combined expression matrix generation   

## Requirements   
Install required tools and ensure they are available in your system `PATH`, or update their paths in `config.yaml`.   
### Core Requirements
* Nextflow
* Java (Version 17 -25)

### Bioinformatics Tools
* fastp
* bbmap (bbduk.sh)
* FastQC
* STAR
* Subread (featureCounts)

### Optional
* R (for preflight check)

## Installation
```bash
git clone https://github.com/bsimkha/RNA-seq-Pipeline.git
cd RNA-seq-Pipeline
```

## Configuration
Edit `config.yaml` before running the pipeline to include proper `input` and `output` directory paths.  

## Input Requirements
* FASTQ files
* File pattern should be edited in `config.yaml`

## Running the Pipeline

### Preflight Check (Optional)
Validate all FASTQ have proper mates
```
Rscript Scripts/Preflight.R
```
### Nextflow Pipeline
This command will run the actual pipeline
```
nextflow run Scripts/main.nf
```

