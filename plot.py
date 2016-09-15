# plot.py
# Mark Rowan, School of Computer Science, University of Birmingham, UK
# April 2012

# Use NEURON / neurosim's grvec hoc routines to read saved spike data
# and then to allow manipulation and plotting in a Python+pyplot environment

# Run interactively, where <filepath> is a directory containing spks, .spks, aux and vars files:
#       python plot.py <filepath>
# or in batch:
#       python plot.py <filepath> [options] [plots]

# Options and plot types are listed at the bottom of this file, and can be
# viewed by running interactively (see above). A good default is 'all [noinhib]'.
# Each option and plot type should be separated by a colon, e.g:
#       python plot.py <filepath> detail 4: activity 1-15: info noinhib



############## DEFINE ROUTINES ###############

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
        print "ERROR: No vars file found at path %s!\n" % filename
        sys.exit()
    return auxvars


# getauxdata(filepath, cellID)
# Returns a list of all aux data items for a given cell ID, read via random access from the given aux file
# Note: does not sanity-check that requested cell ID < numcells!
# Note: all cells must be present and listed in order at each save step, with one 't = ...' line separator
# e.g. cell 0 must always be on line 1, numcells + 1 + 1, numcells * 2 + 1 + 1
def getauxdata(filepath, cell):
    global detail # Declare that we want to use the global 'detail' value

    filename = filepath + "/aux"
    datalist = [] # Create empty list

    if os.path.exists(filename):
        lineno = 2 + cell # Start at line 2 (after "t = ..." line)
        curline = "-1"

        while curline != "":
            # linecache.getline returns "" on read errors, e.g. eof

            curline = linecache.getline(filename, lineno)
            lineno += (auxvars.numcells + 1) * detail  # +1 accounts for the "t = ..." lines between segments
            splitline = curline.split(",")

            # Check that this is correct cell before adding to list; exit otherwise!
            if len(splitline) > 1:

                if int(splitline[CELLID]) == cell:
                    datalist.append(splitline)
                else:
                    print "ERROR: Expected cell %d, found cell %s at line %d!" % (cell, splitline[CELLID], lineno)
                    sys.exit()
    else:
        print "ERROR: No aux file found at path %s!\n" % filename
        sys.exit()

    return datalist;


# getIcellIDs(filepath)
# Returns a list of all cell IDs for I-cells (so they can be explicitly plotted / ignored)
def getIcellIDs(filepath):
    cell_list = []
    for i in range(auxvars.numcells):
        celldata = getauxdata(filepath, i) # Read data
        celltype = int(celldata[0][TYPE])
        if (h.strm(h.CTYP.o(celltype).s,"^I")):
            # Cell is I type
            cell_list.append(int(celldata[0][CELLID]))
    return cell_list


# getEcellIDs(filepath)
# Returns a list of all cell IDs for E-cells (so they can be explicitly plotted / ignored)
def getEcellIDs(filepath):
    cell_list = []
    for i in range(auxvars.numcells):
        celldata = getauxdata(filepath, i) # Read data
        celltype = int(celldata[0][TYPE])
        if (h.strm(h.CTYP.o(celltype).s,"^I")) == False:
            # Cell is E type
            cell_list.append(int(celldata[0][CELLID]))
    return cell_list


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


# Depending on the length of the time series of scale and activity graphs,
# return an x-axis delimited in seconds, hours, or days.
# Returns two values: a list with the x axis points, and a string for the units
def makexaxis(start, numpoints, interval, detail):
    # Default to s
    units = 0.001
    unitstring = "seconds"

    # If > 0.8e7 ms (8000s), use 'hours'
    if start + numpoints * interval * detail > 0.8e7:
	units = 0.001 / 60 / 60
	unitstring = "hours"

    # If > 0.8e8 ms (80 000s, 22 hrs), use 'days'
    if start + numpoints * interval * detail > 0.8e8:
	units = 0.001 / 60 / 60 / 24
	unitstring = "days"

    return np.arange(float(start*units), float((start + numpoints * interval * detail)*units), float((interval * detail)*units)), unitstring


############# DEFINE DIFFERENT PLOT TYPES #################

# raster
# basic spike raster (cell index vs time)
def raster(filepath, rasterfrom, rasterto):
    global titleson # Declare that we want to use the global 'titleson' variable

    rasterfrom = float(rasterfrom)
    rasterto = float(rasterto)
    if rasterto == rasterfrom: 
        # By default (when rasterto = 0), draw a separate plot for each data segment.
        # If rasterfrom != rasterto, draw only a single plot for the data between rasterfrom and rasterto

        # Create save directory if not already present
        if (not os.path.isdir("%s/rasters" % filepath)):
            os.mkdir("%s/rasters" % filepath)

        for i in range(numsegs()):
            print i
            data = read(i) # Read data for segment i

            pyplot.plot(data[0], data[1], '.', markersize=1, color='black')
            pyplot.xlabel("Time (ms)")
            pyplot.ylabel("Cell ID")

            if titleson:
                pyplot.suptitle(filepath)

            pyplot.savefig("%s/rasters/%d.png" % (filepath, i), dpi=300, format="png") # Save at 300 DPI as 'filepath/raster/segment.png'
            pyplot.hold(False) # Don't add next segment's data to current figure
        pyplot.clf() # Clear figure for next plot

    else:
        segmentstart = int(rasterfrom / auxvars.buffertime)
        segmentend = int(rasterto / auxvars.buffertime)
        print "%d - %d" % (segmentstart, segmentend)
        plotx = []
        ploty = []

        for i in range (segmentstart, segmentend + 1):
            print "Reading segment %d" % i
            data = read(i) # Read data for segment i
            for j in range(len(data[0])):
                if data[0,j] >= rasterfrom and data[0,j] <= rasterto:
                    # If data time is within requested bounds, add to plotdata
                    plotx.append(data[0,j])
                    ploty.append(data[1,j])

        print "Plotting..."
        pyplot.plot(plotx, ploty, '.', markersize=1, color='black')
        pyplot.xlabel("Time (ms)")
        pyplot.ylabel("Cell ID")

        if titleson:
            pyplot.suptitle(filepath)

        pyplot.savefig("%s/raster.png" % filepath, dpi=300, format="png") # Save at 300 DPI
        pyplot.clf() # Clear figure for next plot



# info
# Information contribution of each neuron in the network, per data segment
# Using mlabwrap to call Crumiller et al.'s MATLAB Fourier information calculation code
def info(filepath, cells, noinhib):
    # Spike data is saved in alternating 'unique' and 'repeat' stimulation blocks
    # Crumiller's code requires a N-by-M matrix (N = number of trials, M = number of cells),
    # in which each element is a vector containing the spike times of cell M in trial N.
    # This procedure takes the first n saved data segments (up to the auxvars.deletionstart time)
    # and packs all the unique/repeat trials' spike times into separate CSV files
    print "Processing spike data into CSV for information contribution calculations..."

    stimlength = auxvars.infotriallength # Length of 'repeat' or 'unique' stimulation block (ms)
    buffertime = auxvars.buffertime # Length of whole data segment (ms)
    ignoretime = auxvars.GetInfoDur - (auxvars.numinfotrials * stimlength * 2) # How much of the start of the data is not supposed to be part of the info trials?

    if cells == "":
        # Make a list of cells whose data we want to write out
        # (Either we requested all cells, or all-but-I cells)
        if noinhib:
            cells = getEcellIDs(filepath)
        else:
            cells = range(auxvars.numcells)

    if numsegs() * buffertime < auxvars.GetInfoDur:
        print "Not enough spike segments to obtain pre-deletion information (needs at least %d)" % int(auxvars.GetInfoDur / buffertime)
    else:   
        uniquesfile = open("%s/info-uniques.csv" % filepath, 'w')
        repeatsfile = open("%s/info-repeats.csv" % filepath, 'w')
        
        for i in range(int( ignoretime/buffertime ), int( auxvars.GetInfoDur / buffertime)):
            print "Scanning segment %d (up to %d) to add to CSV file" % (i, int(auxvars.GetInfoDur/buffertime)-1)
        
            # For each data segment in the file
            data = read(i) # Read NEURON data for segment i
            for index in range(len(data[0])):
                spiketime = data[0, index] # What time was this cell spike?
                cellID = data[1, index] # Which cell just fired?
                if cellID in cells:
                    trialnum = int((spiketime-ignoretime) / (stimlength * 2)) # Obtain trial number.
                    # Force int, as in Python 2 and 3 '/' operator means one or other of integer or float div...
                    # Find out whether we are in 1st stimlength ms (repeats), or last stimlength ms (uniques)
                    uniquestrial = ((spiketime % float(stimlength * 2)) / float(stimlength * 2)) > 0.5
                    # Tests (stimlength = 4000):
                    # spiketime = 0 => trialnum = 0, uniquestrial = false
                    # spiketime = 3999 => trialnum = 0, uniquestrial = false
                    # spiketime = 4000 => trialnum = 0, uniquestrial = true
                    # spiketime = 7999 => trialnum = 0, uniquestrial = true
                    # spiketime = 8000 => trialnum = 1, uniquestrial = false
                    # spiketime = 11999 => trialnum = 1, uniquestrial = false
                    # spiketime = 12000 => trialnum = 1, uniquestrial = true
                    # spiketime = 15999 => trialnum = 1, uniquestrial = true
                    # spiketime = 16000 => trialnum = 2, uniquestrial = false

                    # Remove offset so that each unique/repeat block has spikes beginning at 0ms
                    # Also div by 1000 to convert ms to s
                    spiketime = (spiketime % stimlength) / 1000
               
                    # If spike was recorded in second half of stimlength, then write to file
                    # (we ignore first-half spikes from each trial, so the network can safely
                    # transition from the previous trial and settle into a stable state).
                    if spiketime >= ((stimlength/1000) / float(2)):
                        # Dump line to file
                        if uniquestrial:
                            uniquesfile.write("%d,%d,%f\n" % (cellID,trialnum,spiketime))
                        else:
                            repeatsfile.write("%d,%d,%f\n" % (cellID,trialnum,spiketime)) 
    
        uniquesfile.close()
        repeatsfile.close()

        # Now find out how long each cell lived for, and write to a third CSV file
        print "Finding cell death times..."
        celldeathtimes = open("%s/info-celldeaths.csv" % filepath, 'w') # Store time-of-death for each cell

        for cellID in cells:
            # For each cell in the requested list
            celldata = getauxdata(filepath, cellID) # Get a list of all entries for this particular cell
            # Find out at which point this cell died (if at all), and its scalefactor at this time
            for entrynum in range(len(celldata)):
                if int(celldata[entrynum][DEAD]) or entrynum >= len(celldata)-1:
                    # If cell has just died, or we've reached the end of the data...
                    break # Stop processing

            cellscale = float(celldata[entrynum][SCALE]) # Get scalefactor
            deathtime = auxvars.t_start + (entrynum * auxvars.recording_interval) # Time of death
            celldeathtimes.write("%d,%d,%f\n" % (cellID,deathtime,cellscale))
            
        celldeathtimes.close()
        print "Done"


# scale
# Scalefactors over time, read from the aux file
# 'cells' is a list of cells which should be plotted (or leave blank to plot all)
def scale(filepath, cells, noinhib):
    print "Plotting scale factors..."
    drawlegend = False
    if len(cells) < 8 and len(cells) > 1:
        drawlegend = True # Only draw a legend when examining 2-7 specific cells

    plotavg = False
    if cells == "":
        plotavg = True
        # Make a list of cells whose data we want to write out
        # (Either we requested all cells, or all-but-I cells)
        if noinhib:
            cells = getEcellIDs(filepath)
        else:
            cells = range(auxvars.numcells)

    celldata = getauxdata(filepath, 0) # Get data for an arbitrary cell, to obtain number of entries for pre-allocating array
    numentries = len(celldata)
    yaxis = np.zeros([len(cells), numentries]) # Collect cell data for plotting average values
    #xaxis = range(auxvars.t_start, (auxvars.t_start + numentries * auxvars.recording_interval * detail), auxvars.recording_interval * detail) # Create x-axis
    xaxis, xaxisunits = makexaxis(auxvars.t_start, numentries, auxvars.recording_interval, detail)

    for cellnum in range(len(cells)):
        # For each cell in the requested list
        celldata = getauxdata(filepath, int(cells[cellnum])) # Get a list of all entries for this particular cell
        # If the cell is not dead, add the data point to the y-axis
        for entrynum in range(len(celldata)):
            if int(celldata[entrynum][DEAD]):
                yaxis[cellnum, entrynum] = -1
            else:
                yaxis[cellnum, entrynum] = float(celldata[entrynum][SCALE])

    if plotavg:
        yavg = np.zeros(numentries)
        stdy = np.zeros(numentries)
        for i in range(numentries):
            yavg[i] = np.mean(yaxis[ yaxis[:,i]!=-1, i ] ) # Get avg scale for each non-dead time point
            stdy[i] = np.std(yaxis[ yaxis[:,i]!=-1, i ] ) # Get std dev for each non-dead time point
        pyplot.errorbar(xaxis, yavg, yerr=stdy, ecolor='grey', linestyle='-', marker='.', markersize=1.0)
    else:
        pyplot.plot(xaxis, yaxis) # Plot each cell individually

    # Set axes limits to appropriate values
    pyplot.ylim(ymin=0)
    pyplot.xlim(xmax=xaxis[len(xaxis)-1], xmin=0) # Cut graph off at end of data, rather than leaving a margin until the next 'significant point'

    # Draw labels
    pyplot.xlabel("Time (%s)" % xaxisunits)
    pyplot.ylabel("Scale factor")
    if titleson:
        pyplot.suptitle(filepath)

    # Draw legend if only a small number of cells
    if drawlegend:
        pyplot.legend(cells, loc='upper right') # Draw legend

    # Save at 300 DPI as 'filepath/scale.pdf'
    pyplot.savefig("%s/scale.pdf" % filepath, dpi=300, format="pdf")
    # Save processed plot data to a file for later analysis
    np.savez("%s/scale" % filepath, x=xaxis, y=yavg, err=stdy)
    pyplot.clf() # Clear figure for next plot
    print "Done"


# activity
# Activity values over time, read from the aux file
# 'cells' is a list of cells which should be plotted (or leave blank to plot all)
def activity(filepath, cells, noinhib):
    global titleson # Declare that we want to use the global 'titleson' variable

    print "Plotting activity values..."
    drawlegend = False
    if len(cells) < 8 and len(cells) > 1:
        drawlegend = True # Only draw a legend when examining 2-7 specific cells

    plotavg = False
    if cells == "":
        plotavg = True
        # Make a list of cells whose data we want to write out
        # (Either we requested all cells, or all-but-I cells)
        if noinhib:
            cells = getEcellIDs(filepath)
        else:
            cells = range(auxvars.numcells)

    celldata = getauxdata(filepath, 0) # Get data for an arbitrary cell, to obtain number of entries for pre-allocating array
    numentries = len(celldata)
    yaxis = np.zeros([len(cells), numentries]) # Collect cell data for plotting average values
    #xaxis = range(auxvars.t_start, (auxvars.t_start + numentries * auxvars.recording_interval * detail), auxvars.recording_interval * detail) # Create x-axis
    xaxis, xaxisunits = makexaxis(auxvars.t_start, numentries, auxvars.recording_interval, detail)

    for cellnum in range(len(cells)):
        # For each cell in the requested list
        celldata = getauxdata(filepath, int(cells[cellnum])) # Get a list of all entries for this particular cell
        # If the cell is not dead, add the data point to the y-axis
        for entrynum in range(len(celldata)):
            if int(celldata[entrynum][DEAD]):
                yaxis[cellnum, entrynum] = -1
            else:
                yaxis[cellnum, entrynum] = float(celldata[entrynum][ACTIVITY]) * 1000

    if plotavg:
        yavg = np.zeros(numentries)
        stdy = np.zeros(numentries)
        for i in range(numentries):
            yavg[i] = np.mean(yaxis[ yaxis[:,i]!=-1, i ] ) # Get avg activity for each non-dead time point
            stdy[i] = np.std(yaxis[ yaxis[:,i]!=-1, i ] ) # Get std dev for each non-dead time point
        pyplot.errorbar(xaxis, yavg, yerr=stdy, ecolor='grey', linestyle='-', marker='.', markersize=1.0)
    else:
        pyplot.plot(xaxis, yaxis) # Plot each cell individually

    # Set axes limits to appropriate values
    pyplot.ylim(ymin=0)
    pyplot.xlim(xmax=xaxis[len(xaxis)-1], xmin=0) # Cut graph off at end of data, rather than leaving a margin until the next 'significant point'

    # Draw labels
    pyplot.xlabel("Time (%s)" % xaxisunits)
    pyplot.ylabel("Activity (Hz)")
    if titleson:
        pyplot.suptitle(filepath)

    # Draw legend if only a small number of cells
    if drawlegend:
        pyplot.legend(cells, loc='upper right') # Draw legend

    # Save at 300 DPI as 'filepath/activity.pdf'
    pyplot.savefig("%s/activity.pdf" % filepath, dpi=300, format="pdf")
    # Save processed plot data to a file for later analysis
    np.savez("%s/activity" % filepath, x=xaxis, y=yavg, err=stdy)
    pyplot.clf() # Clear figure for next plot
    print "Done"


# deletionscale
# Plots relative scaling factors of cells (compared to mean) at time-of-deletion, vs total time
def deletionscale(filepath, noinhib):
    global titleson # Declare that we want to use the global 'titleson' variable
    print "Plotting deletion scale values..."
    if noinhib:
        cells = getEcellIDs(filepath)
    else:
        cells = range(auxvars.numcells)

    celldata = getauxdata(filepath, 0) # Get data for an arbitrary cell, to obtain number of entries for pre-allocating array
    numentries = len(celldata)

    yaxis = [] # Individual cell scale values at time x
    ymeanscale = [] # Mean scale values over all cells at time x

    # Get xaxisunits
    # But throw away the actual xaxis vector, as we don't want fixed-size intervals
    xaxis, xaxisunits = makexaxis(auxvars.t_start, numentries, auxvars.recording_interval, detail)
    xaxis = []


    # Get mean scale of all the population at each entry in the x-axis
    print "Finding mean scale values..."
    # Go through all data once for every entrynum and calculate mean population scale at that point
    allscales = np.zeros([len(cells), numentries]) # Collect cell data for plotting average values

    for cellnum in range(len(cells)):
        # For each cell in the requested list
        celldata = getauxdata(filepath, int(cells[cellnum])) # Get cell data
        # If the cell is not dead, add the data point to the y-axis
        for entrynum in range(len(celldata)):
            if int(celldata[entrynum][DEAD]):
                allscales[cellnum, entrynum] = -1
            else:
                allscales[cellnum, entrynum] = float(celldata[entrynum][SCALE])
    meanscales = np.zeros(numentries)
    for i in range(numentries):
        meanscales[i] = np.mean(allscales[ allscales[:,i]!=-1, i ] ) # Get avg scale for each time point


    # Get individual scale factor for each cell at time of deletion
    print "Finding individual scale values..."
    for cellnum in range(len(cells)):
        # For each cell in the requested list
        celldata = getauxdata(filepath, int(cells[cellnum])) # Get a list of all entries for this particular cell
        # Find out at which point this cell died (if at all), and its scalefactor at this time
        for entrynum in range(len(celldata)):
            if int(celldata[entrynum][DEAD]) or entrynum >= len(celldata)-1:
                # If cell has just died, or we've reached the end of the data...
                cellscale = float(celldata[entrynum][SCALE]) # Get scalefactor
                yaxis.append(cellscale) # Add scalefactor to y-axis
                ymeanscale.append(meanscales[entrynum])  # Obtain mean scalefactor and add to y
                xaxis.append(auxvars.t_start + (entrynum * auxvars.recording_interval)) # Add time to x
                break # Stop processing


    # Data is in arbitrary order (x-axis time depends on cell ID, rather than x-axis incrementing in time)
    # This makes it impossible to draw a line from 'start' to 'finish'. So sort the x,y,mean(y) data first.
    sortedmeanscales = sorted(zip(xaxis, ymeanscale))
    sortedscales = sorted(zip(xaxis, yaxis))
    # Now unzip
    xaxis = [point[0] for point in sortedscales]
    yaxis = [point[1] for point in sortedscales]
    ymeanscale = [point[1] for point in sortedmeanscales]
    
    # Draw plots
    pyplot.plot(xaxis, ymeanscale, '', linestyle='-', color='r') # Plot line for mean scale
    pyplot.plot(xaxis, yaxis, '.', linestyle='', color='b') # Plot each cell's individual scale


    # Set axes limits to appropriate values
    pyplot.ylim(ymin=0)
    pyplot.xlim(xmax=xaxis[len(xaxis)-1], xmin=0) # Cut graph off at end of data, rather than leaving a margin until the next 'significant point'

    # Draw labels
    pyplot.xlabel("Time (%s)" % xaxisunits)
    pyplot.ylabel("Scale factor")
    if titleson:
        pyplot.suptitle(filepath)

    # Save at 300 DPI as 'filepath/deletionscale.pdf'
    pyplot.savefig("%s/deletionscale.pdf" % filepath, dpi=300, format="pdf")
    # Save processed plot data to a file for later analysis
    np.savez("%s/deletionscale" % filepath, x=xaxis, y=yaxis)
    pyplot.clf() # Clear figure for next plot
    print "Done"

# getmeanpopscale
# Helper function for deletionscale
# Returns mean scalefactor of the currently-alive population at a given time
def getmeanpopscale(cells,entrynum):
    totalscale = 0
    livingcells = 0
    for cellnum in range(len(cells)):
        celldata = getauxdata(filepath, int(cells[cellnum])) # Obtain each cell's data
        if not int(celldata[entrynum][DEAD]):
            # If we have found a cell which is still alive
            livingcells = livingcells + 1 # Increase counter of living cells
            totalscale = totalscale + float(celldata[entrynum][SCALE]) # Add to scalefactor sum
    return totalscale / livingcells


# power
# Plots power spectra of the network separately for E and I cells, obtained from the neural spike trains
def power(filepath):
    binsize = 5.0 # ms for bin size (5ms = 200Hz sampling rate)
    Icells = getIcellIDs(filepath) # Make list of I cells

    # Load each segment of saved spike data
    for i in range(numsegs()):
        # For each data segment in the file
        data = read(i) # Read NEURON data for segment i

        # Create empty numpy array for number of spikes-per-bin
        firstspk = data[0,0]
        lastspk = data[0,len(data[1])-1]
        timecovered = lastspk-firstspk
        numbins = int(timecovered / binsize)

        spikesperbinI = np.zeros( numbins )
        spikesperbinE = np.zeros( numbins )

        # Make MUA time series vector by counting all population spikes during each 'binsize' ms
        print "Binning file segment %d into %d bins (window size = %f ms)" % (i, numbins, binsize)

        for j in range(len(data[1])):
            spiketime = data[0,j]
            bin = int(float((spiketime % (numbins * binsize)) / binsize)) # Find bin number into which this spike should go

            cellID = data[1,j]
            if cellID in Icells:
                # Check if this is an I or E cell
                spikesperbinI[bin] = spikesperbinI[bin] + 1
            else:
                spikesperbinE[bin] = spikesperbinE[bin] + 1 

        # Find and subtract mean from the MUA vector
        meanmuaI = float(np.sum(spikesperbinI) / len(spikesperbinI))
        spikesperbinI = spikesperbinI - meanmuaI
        meanmuaE = float(np.sum(spikesperbinE) / len(spikesperbinE))
        spikesperbinE = spikesperbinE - meanmuaE

        # Run MTSpec's spectral analysis on the MUA time series
        # MTSpec documentation: http://krischer.github.com/mtspec/mtspec.multitaper.mtspec.html#mtspec.multitaper.mtspec
        delta = 1.0 / float(1000.0/binsize) # 1 / sample rate (Hz)
        time_bandwidth = 4
        nfft = None
        number_of_tapers = None
        quadratic = False
        adaptive = True
        verbose = False
        optional_output = False
        statistics = False
        rshape = False
        fcrit = False
        
        print "Running mtspec for I cells"
        #spec, freq = multitaper.mtspec(spikesperbin, delta, time_bandwidth, nfft, quadratic, adaptive, verbose, optional_output, statistics, rshape, fcrit)
        specI, freqI = mtspec(spikesperbinI, delta, time_bandwidth)

        print "Running mtspec for E cells"
        specE, freqE = mtspec(spikesperbinE, delta, time_bandwidth)

        # Save to file for further examination
        #numpy.savetxt("%s/power-%d-muaE.txt" % (filepath, i), spikesperbinE)
        #numpy.savetxt("%s/power-%d-specE.txt" % (filepath, i), specE)
        #numpy.savetxt("%s/power-%d-freqE.txt" % (filepath, i), freqE)

        print "Generating plots\n"
        # Create save directory if not already present
        if (not os.path.isdir("%s/power" % filepath)):
            os.mkdir("%s/power" % filepath)

        # I power
        # Draw labels
        pyplot.xlabel("Frequency (Hz)")
        pyplot.ylabel("Normalised power")
        if titleson:
            pyplot.suptitle(filepath)
        pyplot.plot(freqI,specI)
        # Set axes limits to appropriate values
        pyplot.ylim(ymax=0.08, ymin=0)
        pyplot.xlim(xmax=100, xmin=0)
        pyplot.savefig("%s/power/I-%d.pdf" % (filepath, i), dpi=300, format="pdf")
        pyplot.clf()
        
        # E power
        # Draw labels
        pyplot.xlabel("Frequency (Hz)")
        pyplot.ylabel("Normalised power")
        if titleson:
            pyplot.suptitle(filepath)
        pyplot.plot(freqE,specE)
        # Set axes limits to appropriate values
        pyplot.ylim(ymax=0.08, ymin=0)
        pyplot.xlim(xmax=100, xmin=0)
        pyplot.savefig("%s/power/E-%d.pdf" % (filepath, i), dpi=300, format="pdf")
        pyplot.clf()
        
        # Save processed plot data to a file for later analysis
        np.savez("%s/power/I-%d.npz" % (filepath, i), x=freqI, y=specI)
        np.savez("%s/power/E-%d.npz" % (filepath, i), x=freqE, y=specE)

    print "Done"


def nte(filepath,segs):
    # Plot normalised transfer entropy (nTE)
    # Uses a time-binned MUA vector over two given segments of data, for comparison
    # Input requires IDs of two data segments
    # Output as a pair of line plots comparing nTEs per population, for the two segments
    h.usetable_infot = 0 # Turn off log lookup tables
    binsize = 10.0 # ms for bin size
    netscale = round(auxvars.numcells / 470)
    poplabels = ['E6' 'I6' 'I6L' 'E5B' 'E5R' 'I5' 'I5L' 'E4' 'I4' 'I4L' 'E2' 'I2' 'I2L']
    popsizes = np.array([59, 25, 13, 17, 65, 25, 13, 30, 20, 14, 150, 25, 13]) # First element was 60, but cell 0 always seems to be missing
    popsizes *= netscale # Scale-up if we have > 470 cells

    # Load each segment of spike data -- each segment will be plotted as a separate line
    for istring in segs:
        i = int(istring) # Convert char to int
        print "Reading segment %d" % i
        # For each data segment in the list
        data = read(i) # Read NEURON data for segment i

        # Create empty numpy array for number of spikes-per-bin
        firstspk = data[0,0]
        lastspk = data[0,len(data[1])-1]
        timecovered = lastspk-firstspk
        numbins = int(timecovered / binsize)

        spikesperbin = np.zeros( [auxvars.numcells, numbins] ) # x,y matrix (== cellid, spikes per bin)

        # Make MUA time series vector by counting all population spikes during each 'binsize' ms
        print "Binning file segment %d into %d bins (window size = %f ms)" % (i, numbins, binsize)
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

        print allpopMUAs
        print "\n"
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
                print normte
                print "\n"

        pyplot.plot(range(len(popsizes)), normte, linestyle='-', marker='x')

    # Draw labels
    xlocations=range(len(popsizes)) # Allow plot to be spaced equally on x-axis, independent of the value     
    pyplot.xticks(xlocations,poplabels) # Display frequency/weight values over the x tick locations 
    pyplot.xlabel("Population")
    pyplot.ylabel("nTE")
    if titleson:
        pyplot.suptitle(filepath)

    # Save at 300 DPI as 'filepath/deletionscale.pdf'
    pyplot.savefig("%s/nte.pdf" % filepath, dpi=300, format="pdf")
    # Save processed plot data to a file for later analysis
    #np.savez("%s/deletionscale" % filepath, x=xaxis, y=yaxis)
    pyplot.clf() # Clear figure for next plot
    print "Done"

################# COMMAND LINE PROCESSING #################

def processargs(args):
    # Checks list of arguments, and returns a list of cells from either a comma-separated list or
    # a hyphenated range (e.g. 0,1,3,7,9 or 0-10)
    if len(args.split()) > 3:
        print "\n\nERROR: Too many arguments! Permitted: 'noinhib', 'cell1,cell2,..,celln', or '0-100'"

    if len(args.split()) == 1 or (len(args.split()) == 2 and 'noinhib' in args):
        # No arguments (or just 'noinhib') supplied
        return "" # Return an empty list
    else:
        # A cell / cell range / time range argument was supplied...
        args = args.split()
        args.pop(0) # Remove command name from the front of the list of args
        if "-" in args[0]:
            # We're dealing with a request for a range of cells
            args = args[0].split("-") # Split dash-separated range
            cellslist = []
            for i in range(int(args[0]), int(args[1])+1):
                cellslist.append(i)
            return cellslist
        else:
            args = args[0].split(",") # Split comma-separated list
            return args


def makegraph(graphtype):
    if "noinhib" in graphtype:
        noinhib = True
    else:
        noinhib = False

    command = graphtype.split()[0]

    if command == "quit":
        sys.exit()

    elif command == "raster":
        if len(graphtype.split()) == 2:
            rasterrange = graphtype.split()[1]
            rasterrange = rasterrange.split('-')
            rasterfrom = rasterrange[0]
            rasterto = rasterrange[1]
        else:
            rasterfrom = 0
            rasterto = 0
        raster(filepath, rasterfrom, rasterto)

    elif command == "info":
        info(filepath, processargs(graphtype), noinhib)

    elif command == "scale":
        scale(filepath, processargs(graphtype), noinhib)

    elif command == "activity":
        activity(filepath, processargs(graphtype), noinhib)

    elif command == "power":
        power(filepath)

    elif command == "deletionscale":
        deletionscale(filepath, noinhib)

    elif command == "all":
        scale(filepath, "", noinhib)
        activity(filepath, "", noinhib)
        deletionscale(filepath, noinhib)
        raster(filepath, 0, 0)
        info(filepath, "", noinhib)
        power(filepath)

    elif command == "titles":
        global titleson
        titleson = not titleson
        print "\n\n\nTitles on: %s" % titleson

    elif command == "detail":
        global detail
        detail = int(graphtype.split()[1])
        print "\nSetting detail level to %d" % detail

    elif command == "nte":
        nte(filepath, processargs(graphtype))

    else:
        print "\nUnknown graph type or option '%s'" % command



############### MAIN METHOD ################

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
from mtspec import * # Multitaper spectral analysis (MTSpec library from http://pypi.python.org/pypi/mtspec)
import copy


# Check filename was supplied
if len(sys.argv) < 2:
    print "Usage:\npython plot.py <data dir path> [option arg:] [plot_type args:] [plot_type args:] [...:]"
    print "For a list of plot types and arguments, run interactively (i.e. python plot.py <data dir>)"
    sys.exit()

if len(sys.argv) >= 3:
    plot_as_arg = True
    args = sys.argv[2:] # Get all plot type args
    graphtypelist = ""
    for item in args:
        graphtypelist += item + " " # Concatenate the args string back together (including a space between each item)
else:
    plot_as_arg = False

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
global titleson
titleson = False
global detail
detail = 1 # Take every nth point for plotting (set == 1 to plot every point, or > 1 to reduce points)

filepath = sys.argv[1]
print "\nLoading data from %s" % filepath
loadspks(filepath)
auxvars = loadvars(filepath)
print "Done."


if plot_as_arg:
    for graphtype in graphtypelist.split(":"):
        # Make all the requested graph types, then quit
        makegraph(graphtype)
    sys.exit()
else:
    while 1:
        # Get user's graph preference interactively
        print "\n\n\n"
        print filepath
        print "titles                                   toggle graph titles (use 'off' for printing; default off)"
        print "detail                                   set level of detail for plots (every 'n'th point; default=1)"
        print "\n"
        print "all [noinhib]                            draw all plots [excluding inhibitory cells]"
        print "raster [starttime-endtime]               basic spike raster [between given ms times; default=ALL]"
        print "info [cells | noinhib]                   information contribution [comma-separated list of individual cells or range of cells e.g. 0-20 | without inhibitory cells]"
        print "scale [cell1s | noinhib]                 scale factor values [comma-separated list of individual cells or range of cells e.g. 0-20 | without inhibitory cells]"
        print "activity [cells | noinhib]               activity values [comma-separated list of individual cells or range of cells e.g. 0-20 | without inhibitory cells]"
        print "deletionscale                            scalefactors of cells at time-of-deletion, plotted against total time"
        print "power                                    oscillatory power spectra"
        print "quit                                     exit program"

        graphtype = raw_input("\ne.g. \"raster 100000-150000\", \"scale 5,23,190\", \"info 0-49\", or \"activity no-inhib\"\n>")
        makegraph(graphtype)
