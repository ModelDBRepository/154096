#!/bin/sh
#MOAB -l walltime=48:0:0,nodes=1:ppn=1,pmem=2gb
#$ -M m.s.rowan@cs.bham.ac.uk
# Use the current working directory when looking for runsim script
#$ -cwd
#$ -m e

###########################
# BlueBEAR cluster script #
###########################
cd  "$PBS_O_WORKDIR"
echo "runsim.sh $savepath $var=$val $args"
./runsim.sh $savepath "$var=$val" $args
python plotavg.py `dirname $savepath` # Make average-plots for all 
