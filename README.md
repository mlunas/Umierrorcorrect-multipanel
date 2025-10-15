# Umierrorcorrect-multipanel
Alternative umierrorcorrect workflow to filter reads belonging to overlapping panels sequenced under the same index.

Reference
---------

UMIErrorCorrect has been published in Clinical Chemistry.

[Link to the Umierrorcorrect paper](https://doi.org/10.1093/clinchem/hvac136)

Osterlund T., Filges S., Johansson G., Stahlberg A. UMIErrorCorrect and UMIAnalyzer: Software for Consensus Read Generation, Error Correction, and Visualization Using Unique Molecular Identifiers, Clinical Chemistry, 2022;, hvac136

Installation
------------
WIP
Dependencies
------------

Umi-errorcorrect runs using Python 3 and requires the following programs/libraries to be installed (if you run through docker all dependencies are already handled):

Python-libraries (should be installed automatically):

    pysam (v 0.8.4 or greater)

External programs:

    bwa (bwa mem command is used)
    Either of gzip or pigz (parallel gzip)

Install the external programs and add them to the path.

Usage
-----

This modified pipeline is built upon Umierrorcorrect v0.29 and it shares all functionalities with the original pipeline. A complete guide on umierrorcorrect usage can be found under the [Stahlberg group github Umierrorcorrect repository](https://github.com/stahlberggroup/umierrorcorrect/tree/master). This customized pipeline includes the split_bam.py script. The pipeline can be run in the alternative mode by including a second .bed file with the -b2 flag:

```
run_umierrorcorrect.py -r1 read1.fastq.gz -r2 read2.fastq.gz -ul umi_length -sl spacer_length -r reference_fasta_file.fasta -o output_directory -b panel1.bed -b2 panel2.bed
```

