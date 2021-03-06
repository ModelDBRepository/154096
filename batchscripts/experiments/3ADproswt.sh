#!/bin/sh
# Deletion, and default electrostimulation

# Set save path
datadir="rowanms-data/neurostim/3ADproswt/"

# Set job name
jobname="3adproswt"

# Set name of parameter to vary
var="proswt"

# Set list of values to try
#vals="1 2 3 4 5 6 7 8 9 10"
vals="0.5 1 2 3 4 5 6"

# Set random seeds
seeds=""

# Set non-varying parameters with '-c' prefix
# e.g. args="-c {activitybeta=10e-7} -c {activitytau=100e3}"
args="-c {useprosthesis=1} -c {deleting=1}"

./batchcommon $datadir $var "$vals" "$args" "$jobname" "$seeds"
