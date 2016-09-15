#!/bin/sh
#
# Usage:
# runsim.sh [savepath] [args]

# EXAMPLES:

# Basic run with default parameters from alz.hoc (saved in data/<time>):
#    ./runsim.sh

# Basic run with default alz.hoc params, but saved to a defined path:
#    ./runsim.sh data/testname

# Run with one alternative param:
#    ./runsim.sh data/testname infotriallength=8000

# Run with multiple alternative params (must be wrapped in "{ and }" ):
#    ./runsim.sh data/testname "{infotriallength=8000 segmentlength=1600e3}"

echo "Starting runsim.sh"
echo "savepath = $1"
echo "variable = $2"
echo "args = $@"

MODL_INCLUDE="./mod"

if [ $# -lt 1 ]; then
  savepath="data/`date +%Y-%m-%d_%H-%M-%S`" # Default to saving in 'data/<current_time>';
else
  savepath=$1 # Take first argument as the save path
  variable=$2 # Take second argument as the variable / list of variables to be pre-set
  shift # Remove savepath argument from the list
  shift # Remove variable argument from the list
fi

echo "Saving to $savepath"
# If savepath doesn't exist, create it
if [ ! -d $savepath ]; then mkdir -p $savepath; fi


# Run simulation.
# To pass extra args to nrniv, enter each after a '-c' flag, in the form:
#   ./runsim data/test -c "\"filepath=\"hello\"\" [-c ...]
echo nrniv -dll mod/`arch`/.libs/libnrnmech.so $@ -c \"$variable\" -c \"strdef filepath\" -c \"filepath=\\\"$savepath\\\"\" sim.hoc
nrniv -dll mod/`arch`/.libs/libnrnmech.so $@ -c "$variable" -c "strdef filepath" -c "filepath=\"$savepath\"" sim.hoc


# Make graphs
python plot.py $savepath activity noinhib: scale noinhib: deletionscale noinhib #: raster: power  # Don't plot raster, power or info by default as they take a long time
#python plot.py $savepath all noinhib
