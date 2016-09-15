# plotdeletions.py
# Mark Rowan, School of Computer Science, University of Birmingham, UK
# May 2013

# For a given directory containing experiments, each containing multiple runs,
# obtain the proportion of deleted cells (I or E) and plot an error graph

# E.g. for experiment 'ADproswt', containing directories 'proswt0.5', 'proswt2',
# each of which contains runs 1-20: find the proportion of dead E/I cells from 'aux' file
# and plot against the neurostimulation weight value, obtained from 'vars' file.

# Usage: python plotdeletions.py <path_to_experiment>

# ----------------------------------------------------------------------

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


# readlastauxdata(filepath)
# Reads the last "t = ..." segment from a given aux file and returns it as a list of lists
def readlastauxdata(filepath):
    filename = filepath + "/aux"
    datalist = [] # Create empty list

    if os.path.exists(filename):
        # Search to end of file to find line number (nasty hack, but easy to code)
        with open(filename) as f:
            for numlines, l in enumerate(f):
                pass
        numlines +=2 # Account for EOF (1 line) and "t = ..." (1 line)
        #print "Number of lines in aux file = %d" % numlines

        curlineno = numlines - auxvars.numcells # Number of lines to backup
        while curlineno < numlines:
            curline = linecache.getline(filename, curlineno)
            splitline = curline.split(",") # Read this line from the file
            datalist.append(splitline) # datalist contains cells 0->numcells
            curlineno += 1
    else:
        print "ERROR: No aux file found at path %s!\n" % filename
        sys.exit()

    #print datalist
    return datalist
    # There is no sanity-checking here. Ideally we should at least ensure that each
    # line we read from the file is actually a data-line with the correct format,
    # rather than the EOF or "t = ..." markers (in case of a truncated file, or similar)


# getIcellIDs(auxdata)
# Returns a list of all cell IDs for I-cells (so they can be explicitly plotted / ignored)
def getIcellIDs(auxdata):
    cell_list = []
    for i in range(auxvars.numcells):
        celltype = int(auxdata[i][TYPE])
        if (h.strm(h.CTYP.o(celltype).s,"^I")):
            # Cell is I type
            cell_list.append(int(auxdata[i][CELLID]))
    return cell_list


# getEcellIDs(auxdata)
# Returns a list of all cell IDs for E-cells (so they can be explicitly plotted / ignored)
def getEcellIDs(auxdata):
    cell_list = []
    for i in range(auxvars.numcells):
        celltype = int(auxdata[i][TYPE])
        if (h.strm(h.CTYP.o(celltype).s,"^I")) == False:
            # Cell is E type
            cell_list.append(int(auxdata[i][CELLID]))
    return cell_list


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


# Check filename was supplied
if len(sys.argv) < 2:
    print "Usage:\npython plotdeletions.py <data dir path>"
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
# Determine whether we are plotting frequency or weight
plotfreq = 0
plotwt = 0
if "freq" in filepath:
    plotfreq = 1
else:
    if "wt" in filepath:
        plotwt = 1
    # Otherwise, we are plotting stimulation start time

# Get list of sub-directories for each experimental parameter value
dirlist = [o for o in os.listdir(filepath) if os.path.isdir(os.path.join(filepath,o))]
print dirlist

if "baseline" in dirlist:
    dirlist.remove("baseline")
else:
    print "'baseline' dir not found. Please add a symlink to a baseline comparison dir'"
    sys.exit

# Pre-allocate global graph x/y axes
yaxisI = np.zeros([2,len(dirlist)]) # Dimension 1 is mean, dimension 2 is std
yaxisE = np.zeros([2,len(dirlist)]) # Dimension 1 is mean, dimension 2 is std
xaxis = np.zeros([1,len(dirlist)])

# For each experiment's sub-directory (i.e. wt/freq value on the x-axis)
for exptdir in dirlist:
    exptdirpath = "%s/%s" % (filepath, exptdir)

    rundirlist = [o for o in os.listdir(exptdirpath) if os.path.isdir(os.path.join(exptdirpath,o))]
    print rundirlist
    exptnum = dirlist.index(exptdir) # Find indexof(exptdir)
    # Pre-allocate y axis arrays for *this run*
    thisyaxisI = np.array([])
    thisyaxisE = np.array([])

    # For each experimental run (seeds 1-20), to get error bars
    for rundir in rundirlist:
        print "\nrundir %s" % rundir
        rundirpath = "%s/%s" % (exptdirpath, rundir)
        print "rundirpath %s" % rundirpath
        runnum = rundirlist.index(rundir) # Find indexof(rundir)
        print "Loading data for run %s of %s from dir '%s'" % (runnum, len(rundirlist), rundir)
        auxvars = loadvars(rundirpath) # Get vars (for wt/freq)
        if auxvars:
            # Only if the 'vars' file was successfully loaded...
    
            # Obtain the wt/freq/time value from 'vars' -> xaxis
            if plotfreq:
                xaxis[0,exptnum]=auxvars.prosfreq
            else:
                if plotwt:
                    xaxis[0,exptnum]=auxvars.proswt
                else:
                    xaxis[0,exptnum]=float(auxvars.prosthesisstart)/1000/60/60

            # Read auxdata
            auxdata = readlastauxdata(rundirpath)

            # Get list of cells
            listofIcells = getIcellIDs(auxdata)
            listofEcells = getEcellIDs(auxdata)

            # For each I cell
            deadIproportion = 0
            for Icell in listofIcells:
                # Obtain the 'dead' flag value from 'aux' at end of run
                if int(auxdata[Icell][DEAD]):
                    # If cell is dead, add to proportion of dead cells for this run
                    deadIproportion += 1 # float(1.0/len(listofIcells))
                
            # For each E cell
            deadEproportion = 0
            for Ecell in listofEcells:
                # Obtain the 'dead' flag value from 'aux' at end of run
                if int(auxdata[Ecell][DEAD]):
                    # If cell is dead, add to proportion of dead cells for this run
                    deadEproportion += 1 # float(1.0/len(listofEcells))

            # Compare to baseline for this particular run
            # Find baseline level of deletion for this particular run
            BLrundirpath = "%s/baseline/%s" % (filepath, rundir)
            print "Loading baseline data from %s" % BLrundirpath
            
            # Read auxdata
            BLauxdata = readlastauxdata(BLrundirpath)

            # For each I cell
            BLdeadIproportion = 0
            for Icell in listofIcells:
                # Obtain the 'dead' flag value from 'aux' at end of run
                if int(BLauxdata[Icell][DEAD]):
                    # If cell is dead, add to proportion of dead cells for this run
                    BLdeadIproportion += 1 # float(1.0/len(listofIcells))
                
            # For each E cell
            BLdeadEproportion = 0
            for Ecell in listofEcells:
                # Obtain the 'dead' flag value from 'aux' at end of run
                if int(BLauxdata[Ecell][DEAD]):
                    # If cell is dead, add to proportion of dead cells for this run
                    BLdeadEproportion += 1 # float(1.0/len(listofEcells))

            
            thisyaxisI = np.hstack([thisyaxisI,deadIproportion - BLdeadIproportion])
            thisyaxisE = np.hstack([thisyaxisE,deadEproportion - BLdeadEproportion])
            print thisyaxisE
            print thisyaxisI

    # Stack this run's y axis values onto the global graph axes
    yaxisI[0,exptnum] = np.mean(thisyaxisI)
    yaxisI[1,exptnum] = np.std(thisyaxisI)
    yaxisE[0,exptnum] = np.mean(thisyaxisE)
    yaxisE[1,exptnum] = np.std(thisyaxisE)

print xaxis
print yaxisI
print yaxisE

# Sort data
xsorted = xaxis[0][np.argsort(xaxis[0])]
yisorted = yaxisI[0][np.argsort(xaxis[0])]
yierrsorted = yaxisI[1][np.argsort(xaxis[0])]
yesorted = yaxisE[0][np.argsort(xaxis[0])]
yeerrsorted = yaxisE[1][np.argsort(xaxis[0])]

# Save processed plot data to a file for later analysis
#np.savez("%s/avgdeadnum.npz" % filepath, x=xsorted, yi=yisorted, ye=yesorted, yierr=yierrsorted, yeerr=yeerrsorted)

xlocations=np.arange(0,len(xsorted)) # Allow plot to be spaced equally on x-axis, independent of the value

pyplot.axhline(y=0, linewidth=2, color='0.5', linestyle='dashed') # Draw dashed line at y=0
pyplot.errorbar(xlocations, yisorted, yerr=yierrsorted, ecolor='b', elinewidth=2, linestyle='None', marker='x', markersize=10, markeredgewidth=3, markeredgecolor='b')
pyplot.errorbar(xlocations, yesorted, yerr=yeerrsorted, ecolor='r', elinewidth=2, linestyle='None', marker='x', markersize=10, markeredgewidth=3, markeredgecolor='r')

pyplot.xticks(xlocations,xsorted) # Display frequency/weight values over the x tick locations 
x1,x2,y1,y2 = pyplot.axis() 
pyplot.xlim(xmin = x1-1, xmax = x2+1) # Give a bit of space either side of the data

#pyplot.ylabel('Change in proportion of dead cells')
pyplot.ylabel('Change in # of dead cells')
if plotfreq:
    pyplot.xlabel('Stimulation frequency (Hz)')
else:
    if plotwt:
        pyplot.xlabel('Stimulation weight')
    else:
        pyplot.xlabel('Stimulation start time (hours)')

matplotlib.rcParams.update({'font.size': 16})
pyplot.savefig("%s/avgdeadnum.pdf" % filepath, dpi=300, format="pdf")

