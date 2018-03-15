# Pipeline for preparing a new reference

This pipeline creates the files necessary for read depth analysis on a new reference.
These steps include:
  * Reference masking
  * Mask track creation
  * GC window analysis
  * SUNK identification
  * DTS window creation

## Notes for future maintainers
Many of the scripts in this analysis can't be updated to use the python modules
in our wssd_sunk repository due to a dependence on a particular version 
of pygr. These scripts use Peter Sudmant's old environment.
