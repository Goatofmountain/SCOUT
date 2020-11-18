## SCOUT ##
A genotyper designed for single-cell whole genome sequencing data.

## Introduction ##
SCOUT is a local-smoothing mixture generative model based genotyper that makes
accurate single-cell SNV detection by using the other SNV information from the
local genome. Different from the other single-cell genotyper, our method requires
the single-cell WGS data of interest when genotyping.

## Dependencies ##
- Python version >= 3.7
- numpy version >= 1.19.1
- pandas version >= 0.24.2
- pysam version >= 0.15.3
- sklearn version >= 0.20.3

## Installation ##
Clone the SCOUT repository:
```
git clone git@github.com:Goatofmountain/SCOUT.git
```

This repository contains the following files which illustrate the methods described in our paper.

The following files pertain to this:

- source/Calculate/Candidate.py (functions and python class used in the software)
- source/Calculate/SomaticSNV.py
- bin/SCOUT.py (the main program for small genomic intervals)
- bin/SCOUT_WholeGenome.py (the main program for the whole genome genotyping)

## Usage ##
The program requires bam file of the target single-cell WGS. The bam file should be sorted by coordinates.


Assuming indexed reference genome file to be ref.fa, the aim single-cell WGS bam file to be sc_bam.fa,
the output name to be output_Name, the result would be placed on OUTDIR, and the entire calculation
process would use NumCPU threads, to run SCOUT for the whole genome:

```
python SCOUT/bin/SCOUT_WholeGenome.py -N <output_Name> -r <ref.fa> -i <sc_bam.fa> -o <OUTDIR> -P <NumCPU> -m WholeGenome
```
run SCOUT for specific genome region (aimchrom:Start-End):
```
python SCOUT/bin/SCOUT_WholeGenome.py -N <output_Name> -r <ref.fa> -i <sc_bam.fa> -o <OUTDIR> -P <NumCPU> -c <aimchrom> -S <Start> -E <End>
```
The arguments of Monovar are as follows:
```
-N: The name of output dir and files.
-r: Reference genome file.
-i: The input scWGS data.
-o: The output dir. (Default value: "./")
-P: Number of threads to use in multiprocessing (Default value: 1)
-c: Aim chromosome for genotyping. (Default value: chr1)
-S: Start point for genotyping. (Default value: 100000000)
-E: End point for genotyping. (Default value: 103000000)
-M: Weight decay parameter associated with the half life of wight. (Default value: np.log(3))
-W: Window size for the range of adjacent SNVs considered for the calculation. (Default value: 30000)
```

We recommand to preprocess the data with the MarkDuplicate and ApplyBQSR pipeline of GATK before using our
software.

Detecting the somatic mutation for single-cell data with SCOUT/source/Calculate/SomaticSNV.py.
Assuming the germline vcf file to be germ.vcf, the SCOUT output vcf to be sc.vcf, and the output
somatic snv to be somatic.vcf, calculate the somatic SNV:
```
python SCOUT/source/Calculate/SomaticSNV.py -i <sc.vcf> -g <germ.vcf> -o <somatic.vcf>
```
