#!/bin/sh
# Regular deletion without scaling, to show necessity of scaling

# Set save path
datadir="rowanms-data/neurostim/0deletionNoscaling/"

# Set job name
jobname="0deletionNoscaling"

# Set name of parameter to vary
var="abc"

# Set list of values to try
#vals="1 2 3 4 5 6 7 8 9 10"
vals="0"

# Set random seeds
seeds=""

# Set non-varying parameters with '-c' prefix
# e.g. args="-c {activitybeta=10e-7} -c {activitytau=100e3}"
args="-c {useprosthesis=0} -c {scaling=0} -c {dynamicdelete=0}"

./batchcommon $datadir $var "$vals" "$args" "$jobname" "$seeds"
