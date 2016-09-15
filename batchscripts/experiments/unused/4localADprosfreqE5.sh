#!/bin/sh
# Deletion, and default electrostimulation, but only for E5

# Set save path
datadir="rowanms-data/neurostim/4localADprosfreqE5/"

# Set job name
jobname="4localfreqE5"

# Set name of parameter to vary
var="prosfreq"

# Set list of values to try
#vals="1 2 3 4 5 6 7 8 9 10"
vals="0.5 1 2 3 4 5 10 20 30"

# Set random seeds
seeds=""

# Set non-varying parameters with '-c' prefix
# e.g. args="-c {activitybeta=10e-7} -c {activitytau=100e3}"
args="-c {useprosthesis=1} -c {deleting=1} -c {nproscellpops=2}"

./batchcommon $datadir $var "$vals" "$args" "$jobname" "$seeds"
