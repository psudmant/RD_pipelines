## Read-Depth CNV calling pipelines

[![Snakemake](https://img.shields.io/badge/snakemake-≥3.5.4-brightgreen.svg?style=flat-square)](http://snakemake.bitbucket.org)

pipelines implemented for various read-depth CNV calling tasks, including:

* RD_setup: performing QC, GC-normalization, parameter fitting, browser track generation and DenseTrackSet windowed copy number estimation
* GC_analysis: asssessing GC_associated coverage biases
* dCGH: calling CNVs using dCGH

Quick setup:
```
git clone https://github.com/eichlerlab/RD_pipelines --recursive
cd RD_pipelines
. env.conf
make
```
