#!/bin/bash
module load acemd/3212.pro
# work around bad libmipwrapper
export LD_LIBRARY_PATH=/home/mi/fab/local/lib
export OMP_NUM_THREADS=2
export OMP_DYNAMIC=TRUE
nice -20 acemd input
