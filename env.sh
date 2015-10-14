# Setup Python environment.
export PYTHONPATH=/net/eichler/vol7/home/psudmant/EEE_Lab/1000G/1000genomesScripts/windowed_analysis/DTS_window_analysis:$PYTHONPATH
export PYTHONPATH=/net/eichler/vol7/home/psudmant/EEE_Lab/1000G/1000genomesScripts:$PYTHONPATH
export PYTHONPATH=/net/gs/vol2/home/psudmant/EEE_Lab/projects/common_code:$PYTHONPATH
export PYTHONPATH=/net/eichler/vol7/home/psudmant/EEE_Lab/projects/common_code/ssf_DTS_caller:$PYTHONPATH
export PYTHONPATH=/net/eichler/vol7/home/psudmant/EEE_Lab/1000G/1000genomesScripts/get_gc_correction:$PYTHONPATH

module purge

module load modules modules-init modules-gs
module load modules-eichler
module load python/2.7.2
module load hdf5/1.8.8
module load pytables/2.3.1_hdf5-1.8.8
module load numpy/1.6.1
module load scipy/0.10.0
module load lzo/2.06
module load zlib/1.2.6
module load R/2.15.0
