# plotavg.py
# Mark Rowan, School of Computer Science, University of Birmingham, UK
# December 2012

# For plotting 'average of averages' of activity and scale when a number
# of runs with identical parameters, but different random seeds, have been
# done.

# Given a filepath, this simple script checks all subdirectories for any
# 'activity.npz' and 'scale.npz' files (which contain activity and scale
# factor data as plotted in the PDFs).

# It calculates the mean and std of the various runs' data and presents
# this on a similar-looking graph in the directory in which the script
# was called, as scale.pdf and activity.pdf


import sys
import os
import matplotlib
matplotlib.use('agg') # Prevent pyplot from trying to open a graphical front-end on headless clients
from matplotlib import pyplot
import numpy as np

filepath = sys.argv[1] # argv[0] is the name of the script, argv[1] is the filepath
print "\nLoading data from %s" % filepath


# Create list of sub-directories
#os.chdir(filepath) # Change to given directory
dirlist = [o for o in os.listdir(filepath)] # if os.path.isdir(o)]
print dirlist

# Initialise arrays for y-axes and one x-axis
scale = np.array([])
activity = np.array([])
x = np.array([])

# For each sub-directory
for dir in dirlist:
    # obtain activity.npz and scale.npz if present
    scalefile = "%s/%s/scale.npz" % (filepath, dir)
    activityfile = "%s/%s/activity.npz" % (filepath, dir)
    if not os.path.isfile(scalefile) or not os.path.isfile(activityfile):
        print "Missing scale.npz or activity.npz in %s" % dir
        # print error if not present, but continue
    else:
        # Load scale file
        print scalefile
        scaledata = np.load(scalefile)
        npscaledata = scaledata['y']
        print "Scale size: %d" % np.size(npscaledata)
        # Check if this data is shorter than it should be
        # Can't check scale.shape[1] if scale is currently only 1-D
        if np.size(scale.shape) > 1 and np.size(npscaledata) < scale.shape[1]:
            sizediff = scale.shape[1] - np.size(npscaledata)
            # Pad the array to bring it up to correct size
            npscaledata = np.hstack((npscaledata, np.zeros(sizediff)))
            print "Padded by %d" % sizediff

        # Append y-axis data to 'scale' (or assign scale=y if y is empty)
        if np.size(scale) < 1:
            scale = npscaledata
        else:
            scale = np.vstack((scale, npscaledata))
       
        # Grab this file's x-axis if it's bigger than the current one
        if np.size(x) < 1 or (np.size(x.shape) > 1 and x.shape[1] < np.size(scaledata['x'])):
            x = scaledata['x']


        # Load activity file
        print activityfile
        activitydata = np.load(activityfile)
        npactivitydata = activitydata['y']
        print "Activity size: %d" % np.size(npactivitydata)
        # Check if this data is shorter than it should be
        # Can't check activity.shape[1] if activity is currently only 1-D
        if np.size(activity.shape) > 1 and np.size(npactivitydata) < activity.shape[1]:
            sizediff = activity.shape[1] - np.size(npactivitydata)
            # Pad the array to bring it up to correct size
            npactivitydata = np.hstack((npactivitydata, np.zeros(sizediff)))
            print "Padded by %d" % sizediff

        # Append y-axis data to 'activity' (or assign activity=y if y is empty)
        if np.size(activity) < 1:
            activity = npactivitydata
        else:
            activity = np.vstack((activity, npactivitydata))

# Remove NaNs from data
print "Removing NaNs"
print "%d from scale" % np.sum(np.isnan(scale))
scale = np.nan_to_num(scale)
print "%d from activity" % np.sum(np.isnan(activity))
activity = np.nan_to_num(activity)

# Plot scale
pyplot.errorbar(x, np.mean(scale,0), np.std(scale,0), ecolor='grey', linestyle='-', marker='.', markersize=1.0)
# Draw labels
pyplot.xlabel("Time (days)")
pyplot.ylabel("Scale factor")
# Save at 300 DPI as 'filepath/scale.pdf'
pyplot.savefig("%s/scale.pdf" % filepath, dpi=300, format="pdf")
np.savez("%s/scale" % filepath, x=x, y=np.mean(scale,0), err=np.std(scale,0))

pyplot.clf() # Clear plot

# Plot activity
pyplot.errorbar(x, np.mean(activity,0), np.std(activity,0), ecolor='grey', linestyle='-', marker='.', markersize=1.0)
# Draw labels
pyplot.xlabel("Time (days)")
pyplot.ylabel("Activity (Hz)")
# Save at 300 DPI as 'filepath/activity.pdf'
pyplot.savefig("%s/activity.pdf" % filepath, dpi=300, format="pdf")
np.savez("%s/activity" % filepath, x=x, y=np.mean(activity,0), err=np.std(activity,0))

