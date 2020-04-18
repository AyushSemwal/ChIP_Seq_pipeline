# ChIP_Seq_pipeline

This pipeline creating script is written in python. For each sample it creates a directory where different pbs scripts are created and stored. Each pbs script is the part of the whole pipeline. 

The script takes in following parameters (in order) as input:

1) bowtie2 built reference index directory
2) fastq files directory
3) output directory
4) samples metadata file
5) controls metadata file
6) SNPs file

You just need to submit bowtie2_align.pbs script to the cluster for each sample. This is done in such a manner so that the pipleine can be run for eacxh sample in parallel.

Depending upon whether allele specific analysis has to be done or not, SNPsplit is run. The output folder contains Bowtie2 Bam file, de-duplicated bam files, SNP splitted bam files (if allele-specific) and bigwig files.
