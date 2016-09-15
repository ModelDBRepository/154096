#!/bin/sh
# No stimulation, just deletion

# Set save path
datadir="rowanms-data/neurostim/1AD/"

# Set job name
jobname="1AD"

# Set name of parameter to vary
var="abc"

# Set list of values to try
#vals="1 2 3 4 5 6 7 8 9 10"
vals="0"

# Set random seeds
seeds=""

# Set non-varying parameters with '-c' prefix
# e.g. args="-c {activitybeta=10e-7} -c {activitytau=100e3}"
args="-c {useprosthesis=0}"

./batchcommon $datadir $var "$vals" "$args" "$jobname" "$seeds"
