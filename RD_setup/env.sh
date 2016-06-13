# Setup Python environment.
#export PYTHONPATH=$PYTHONPATH:/net/eichler/vol7/home/psudmant/EEE_Lab/1000G/1000genomesScripts/windowed_analysis/DTS_window_analysis
#export PYTHONPATH=$PYTHONPATH:/net/eichler/vol7/home/psudmant/EEE_Lab/1000G/1000genomesScripts
#export PYTHONPATH=$PYTHONPATH:/net/gs/vol2/home/psudmant/EEE_Lab/projects/common_code
#export PYTHONPATH=$PYTHONPATH:/net/eichler/vol7/home/psudmant/EEE_Lab/projects/common_code/ssf_DTS_caller
#export PYTHONPATH=$PYTHONPATH:/net/eichler/vol7/home/psudmant/EEE_Lab/1000G/1000genomesScripts/get_gc_correction

module purge

module load modules modules-init modules-gs
module load modules-eichler
module load anaconda/2.3.0
module load lzo/2.06
module load zlib/1.2.6
module load R/2.15.0
