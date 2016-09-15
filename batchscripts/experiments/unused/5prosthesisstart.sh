#!/bin/sh
# Deletion, and 'best' electrostimulation, with varying prosthesisstarts

# Set save path
datadir="rowanms-data/neurostim/5prosthesisstart"

# Set job name
jobname="5prosstart"

# Set name of parameter to vary
var="prosthesisstart"

# Set list of values to try
#vals="1 2 3 4 5 6 7 8 9 10"
vals="16000e3 24000e3 40000e3 64000e3 90000e3 120000e3 160000e3"

# Set random seeds
seeds=""

# Set non-varying parameters with '-c' prefix
# e.g. args="-c {activitybeta=10e-7} -c {activitytau=100e3}"
args="-c {useprosthesis=1} -c {deleting=1}"

./batchcommon $datadir $var "$vals" "$args" "$jobname" "$seeds"
