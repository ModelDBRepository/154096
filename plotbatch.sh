#!/bin/sh

# Use this script to generate a certain plot (or plots) for a whole batch of
# saved data.
# USAGE: ./plotbatch.sh <parentdirectory> <plottypes>
# WHERE: parentdirectory contains subdirectories, each containing spks and vars
# WHERE: plottypes in format "raster 0-1000" or "scale noinhib: activity noinhib"
# (including quotes)

dir=$1
plottypes=$2

for subdir in $( ls $dir ); do
  plotpath=$dir/$subdir
  
  export plotpath
  export plottypes

  echo "msub -v $plotpath -v $plottypes -o $plotpath -e $plotpath plotbatchcluster.sh"
  msub -v plotpath -v plottypes -o $plotpath -e $plotpath plotbatchcluster.sh
done
