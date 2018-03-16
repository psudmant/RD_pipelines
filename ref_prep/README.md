# Pipeline for preparing a new reference

This pipeline creates the files necessary for read depth analysis on a new reference.
These steps include:
  * Reference masking
  * Mask track creation
  * GC window analysis
  * SUNK identification
  * DTS window creation

## Quick start
Setup `config.yaml` with a reference name and paths to a masked
and unmasked reference fasta files and their fai index files.
By default, dgv, gaps, and segdups bed files are empty, but if the data exist
for your reference, they can be included. Run the pipeline with some sane parameters:

```code bash
snakemake --drmaa " -V -cwd -w n -e ./log -o ./log {params.sge_opts} -S /bin/bash" -j 6 -w 60 -kT
```

## Control locations
Control locations must be specified for the read depth to copy number regression in `RD_setup`.
A common way to do this is to map Illumina data for your reference against the closest existing reference,
then identify large (>50 kbp) regions of fixed copy. Another approach is to use control regions for a
closely related species (e.g. within great apes). This pipeline includes optional steps for lifting over
these regions using a UCSC chain file. To run it, specify `input_control_locs` and `chain_file` in `config.yaml`
and run rule `get_control_locations`:

```code bash
snakemake --drmaa " -V -cwd -w n -e ./log -o ./log {params.sge_opts} -S /bin/bash" -j 6 -w 60 -kT get_control_locations
```

## RD_setup configuration
Correctly matching output files with the necessary files for RD_setup
can be difficult. To simplify this process, once you have control locations
for your reference, you can generate text that can be appended to `config.yaml`
in `RD_setup` using `snakemake get_rd_setup_text`.

## Notes for future maintainers
Many of the scripts in this analysis can't be updated to use the python modules
in our wssd_sunk repository due to a dependence on a particular version 
of pygr. These scripts use Peter Sudmant's old environment.
