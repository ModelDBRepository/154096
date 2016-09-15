#!/bin/sh
#$ -M m.s.rowan@cs.bham.ac.uk
#MOAB -q bbtest
# Use the current working directory when looking for plot.py
#$ -cwd
#$ -m e
#$ -N alzplot

cd  "$PBS_O_WORKDIR"
echo "python plot.py $plotpath $plottypes" # activity noinhib: scale noinhib: raster: power"
#python plot.py $plotpath activity noinhib: scale noinhib: power: raster
python plot.py $plotpath $plottypes

