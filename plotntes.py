# plotntes.py
# Mark Rowan, School of Computer Science, University of Birmingham, UK
# Aug 2013

# For a given directory containing experiments, containing multiple runs,
# and a supplied list of data segments, obtain and plot the nTE per population
# for each of the data segments (averaged over all runs)

# E.g. for experiment 'proswt3', containing runs '1', '2'... ,
# and for the requested segments '3', '12', and '50',
# plot the nTE of each population for each segment across all runs

# Usage: python plotntes.py <path_to_experiment_variable> <data_segments>

# ----------------------------------------------------------------------
# loadspks(filepath)
# takes name of file to be read and calls grvec read procedure to load file
def loadspks(filepath):
    filename = filepath + "/spks"

    if os.path.exists(filename):
        h.grv_.read_vfile(filename)  # open file using grvec
        print "There are %d segments in this file" % numsegs()
    else:
        print "ERROR: No spks file found at path %s!" % filename
        sys.exit()


# numsegs()
# returns number of data segments in the file
def numsegs():
    return int(h.grv_.llist.count())


# read(index)
# takes index of file segment to be read (1 to numsegs)
# returns vec and tvec spiking data in a 2-D array
# data[0] is the spike time vector
# data[1] is the vector of cell IDs which spike at times in data[0]
def read(index):
    print "Reading grvec data"
    data =  np.zeros( (1, 1) ) # Preassign zero-length data variable

    # Check we're not requesting an index beyond the extent of the file
    if index > numsegs():
        print "Index must be <= %d" % numsegs()

    else:
        h.grv_.rv_readvec(index, h.tvec, h.vec) # Read segment 'index' into vec/tvec

        # Convert h.tvec and h.vec into Python / numpy arrays.
        data = np.zeros( (2, h.vec.size()) ) # Preallocate for speed

        for i in range(int(h.vec.size())):
            data[0,i] = h.tvec.x[i]
            data[1,i] = h.vec.x[i]

        #print data.shape
        return data

# loadvars(filepath)
# Reads variables saved in the vars file into memory
def loadvars(filepath):
    filename = filepath + "/vars"

    if os.path.exists(filename):
        # Read the file into 'auxvars'
        import imp
        f = open(filename)
        auxvars = imp.load_source('auxvars', '', f) # Read parameters (e.g. numcells) from aux file
        f.close()
    else:
        print "ERROR: No vars file found at path %s! Ignoring this run.\n" % filename
        #sys.exit()
        auxvars = 0
    return auxvars

# ----------------------------------------------------------------------

# Imports
import sys
import readline
import numpy as np
import string
import os
import linecache
import matplotlib
matplotlib.use('agg') # Prevent pyplot from trying to open a graphical front-end on headless clients
from matplotlib import pyplot
import copy


# Check filename was supplied
if len(sys.argv) < 3:
    print "Usage:\npython plotntes.py <data dir path> <data segment>"
    sys.exit()


# Handle NEURON imports
print "Loading NEURON/grvec routines... \n\n"
from neuron import h

# Load grvec / intfzip simulation
h('xopen("setup.hoc")')
h('xopen("nrnoc.hoc")')
h('load_file("init.hoc")')
h('load_file("grvec.hoc")')
h('load_file("labels.hoc")')
h('load_file("infot.hoc")')

# Load file
global auxvars # Vars from the vars file are held here, to be accessible by all methods
# Set index numbers for the aux file data lines (format: CELL,TYPE,ACTIVITY(firing rate),POSTSYNACTIVITY,TARGET,SCALE,DEAD)
global CELLID
CELLID = 0
global TYPE
TYPE = 1
global FIRINGRATE
FIRINGRATE = 3
global ACTIVITY
ACTIVITY = 2
global TARGET
TARGET = 4
global SCALE
SCALE = 5
global DEAD 
DEAD = 6

# MAIN LOOP
filepath = sys.argv[1] # Describes the top-level experiment directory
segnum = int(sys.argv[2]) # Data segment to be plotted

h.usetable_infot = 0 # Turn off log lookup tables
binsize = 10.0 # ms for bin size
poplabels = np.array(['E6', 'I6', 'I6L', 'E5B', 'E5R', 'I5', 'I5L', 'E4', 'I4', 'I4L', 'E2', 'I2', 'I2L'])
base_popsizes = np.array([59, 25, 13, 17, 65, 25, 13, 30, 20, 14, 150, 25, 13]) # First element was 60, but cell 0 always seems to be missing

# Get list of sub-directories for each experimental parameter value
dirlist = [o for o in os.listdir(filepath) if os.path.isdir(os.path.join(filepath,o))]
print dirlist
allnormtes = np.zeros( (len(dirlist), len(base_popsizes)) ) # Pre-allocate x-axis array

for rundir in dirlist:
    rundirpath = "%s/%s" % (filepath, rundir)

    print "\nrundir %s" % rundir
    print "rundirpath %s" % rundirpath
    runnum = dirlist.index(rundir) # Find indexof(rundir)
    print "Loading data for run %s of %s from dir '%s'" % (runnum, len(dirlist), rundir)
    auxvars = loadvars(rundirpath) # Get vars (for wt/freq)
    netscale = round(auxvars.numcells / 470)
    popsizes = base_popsizes * netscale # Scale up if we have > 470 cells
    # Plot normalised transfer entropy (nTE)
    # Uses a time-binned MUA vector over two given segments of data, for comparison
    print "\nLoading data from %s" % rundirpath
    loadspks(rundirpath)
    print "Reading segment %d" % segnum
    data = read(segnum) # Read NEURON data for segment i

    # Create empty numpy array for number of spikes-per-bin
    firstspk = data[0,0]
    lastspk = data[0,len(data[1])-1]
    timecovered = lastspk-firstspk
    numbins = int(timecovered / binsize)

    spikesperbin = np.zeros( [auxvars.numcells, numbins] ) # x,y matrix (== cellid, spikes per bin)

    # Make MUA time series vector by counting all population spikes during each 'binsize' ms
    print "Binning file segment %d into %d bins (window size = %f ms)" % (segnum, numbins, binsize)
    for j in range(len(data[1])):
        spiketime = data[0,j]
        cellid = data[1,j]
        bin = int(float((spiketime % (numbins * binsize)) / binsize)) # Find bin number into which this spike should go
        spikesperbin[cellid,bin] += 1

    # Combine cells from each population together
    allpopMUAs = []
    currentpopcell = 0
    for popnum in range(len(popsizes)):
        finalpopcell = currentpopcell + popsizes[popnum] - 1
        popMUA = np.sum(spikesperbin[currentpopcell:finalpopcell, :], axis=0)
        currentpopcell += popsizes[popnum]
        allpopMUAs.append(popMUA) # append all popMUAs into one structure

    #print allpopMUAs
    #print "\n"
    normte = np.zeros(len(popsizes))
    # For each population, obtain its MUA against all other populations
    for popnum in range(len(allpopMUAs)):
        allotherpops = copy.deepcopy(allpopMUAs)
        thispop = allotherpops.pop(popnum) # Gives us this population, AND removes it from 'remaining' list
        thispopvec = h.Vector(thispop) # Convert thispop to a hoc Vector

        # Now, for each 'other' population, compare it to 'thispop'
        for otherpop in allotherpops:
            # Convert otherpop to a hoc Vector
            otherpopvec = h.Vector(otherpop)
            normtevec = h.normte(thispopvec, otherpopvec, 30) # Try 30 shuffles
            #normtevec.printf("%8.4f\n")
            #print "using %f, popnum %d" % (normtevec.x[2], popnum)
            normte[popnum] += normtevec.x[2]
            #print normte
            #print "\n"

    allnormtes[runnum] = normte
    print "allnormtes"
    print allnormtes

xaxis = range(len(base_popsizes))
meanperpop = np.mean(allnormtes, axis=0)
stdperpop = np.std(allnormtes, axis=0)
print meanperpop
print stdperpop
pyplot.errorbar(xaxis, meanperpop, yerr=stdperpop, linestyle='-', marker='x', markersize='4')
# Draw labels
xlocations=range(len(popsizes)) # Allow plot to be spaced equally on x-axis, independent of the value     
pyplot.xticks(xlocations,poplabels) # Display frequency/weight values over the x tick locations 
pyplot.xlabel("Population")
pyplot.ylabel("nTE")

# Save at 300 DPI as 'filepath/deletionscale.pdf'
matplotlib.rcParams.update({'font.size': 16})
pyplot.savefig("%s/nte%s.pdf" % (filepath, segnum), dpi=300, format="pdf")
# Save processed plot data to a file for later analysis
np.savez("%s/nte%s" % (filepath, segnum), x=xaxis, y=meanperpop, yerr=stdperpop, xlabels=poplabels, xlocations=xlocations)
print "Done"
