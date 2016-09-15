"""
PYHOC

Python module for using other Python modules via the hoc interpreter.
You will need to begin with this:

	objref p
	p = new PythonObject()
	nrnpython("import pyhoc")

Many of the usage examples expect an input consisting of one or two
long vectors representing time series, here called A and B. The following
code generates such vectors. Simply copy and paste this code into the
NEURON interpreter, then copy and paste the example code for a given
module, and the script should run.

	objref A, B, r1, r2 // Initialize variables
	r1 = r2 = new Random() // Create random variables
	r1.normal(0,1) r2.normal(0,1) // Set random variables to normal distribution
	A = new Vector(20000) // Initialize the first time series
	A.indgen(0.01) A.sin(5,0) A.addrand(r1) // Populate it
	B = new Vector(20000) // Initialize the second time series
	B.indgen(0.01) B.sin(5,15) B.addrand(r2) // Populate it

Alternatively, if you're working in an intfcol-based environment,
a more realistic pair of time series (i.e. actual LFP time series)
can be obtained as follows:

	objref A, B
    A=nqLFP[0].v[0]
    B=nqLFP[0].v[1]

Version: 2011apr28
"""

# BSMART
def bsmart(nqx1,nqx2,ntrls=1,npts=-1,p=12,fs=200,freq=100): # Set defaults for everything except the data x
    """
    This is the wrapper for the BSMART code.
    
    Usage is similar to bsmart.py:
    	grangernqs=pyhoc.bsmart(x1,x2,[ntrls,npts,p,fs,freq]);
    where
    	x1 = vector representing first time series
        x2 = vector representing second time series
    	ntrls = number of trials in the time series (best set to 1)
    	npts = length of input (if set to -1, is calculated automatically)
    	p = polynomial order for fitting (lower = smoother fit)
    	fs = sampling rate for the time series (in Hz)
    	freq = maximum frequency to be returned (usually fs/2)
    
    where grangernqs has the following fields:
        F -- vector of frequencies for each of the following
        pp1 -- power spectrum for first time series
        pp2 -- power spectrum for second time series
        cohe -- coherence between the two time series
        Fx2y -- causality from first to second time series
        Fy2x -- causality from second to first time series
        Fxy -- nondirectional causality
        
    Example usage from NEURON is as follows: 
        objref output
        output=p.pyhoc.bsmart(A,B)
        
        output.gr("Fx2y","F") // Strong causality
        output.gr("Fy2x","F") // No causality

    Version: 2011apr21
    """
## Import packages
    from numpy import array, zeros, size, shape # Shorten useful functions
    from bsmart import timefreq, pwcausalr
    from neuron import h
    
## Initialize data vectors
    tmp1=array(nqx1) # Convert NQS table to Numpy arrays
    tmp2=array(nqx2)
    if npts==-1: npts=size(tmp1,0) # Reset npts if needed
    x=array(zeros((2,npts))) # Store both time series in one matrix
    x[0,]=tmp1
    x[1,]=tmp2
    
## Do the analysis
    F,pp,cohe,Fx2y,Fy2x,Fxy=pwcausalr(x,int(ntrls),int(npts),int(p),fs,int(freq)); # Do the analysis
    
## Initialize hoc objects
    h('objref grangernqs') # Initialize NQS object
    h('objref F, pp1, pp2, cohe, Fx2y, Fy2x, Fxy, tmp') # Initialize vectors
    h('F   =new Vector()')
    h('pp1 =new Vector()')
    h('pp2 =new Vector()')
    h('cohe=new Vector()')
    h('Fx2y=new Vector()')
    h('Fy2x=new Vector()')
    h('Fxy=new Vector()')
    
## Convert from Python to hoc
    h.tmp=F        ; h('F=F.from_python(tmp)')
    h.tmp=pp[0,:]  ; h('pp1=pp1.from_python(tmp)')
    h.tmp=pp[1,:]  ; h('pp2=pp2.from_python(tmp)')
    h.tmp=cohe[0,:]; h('cohe=cohe.from_python(tmp)')
    h.tmp=Fx2y[0,:]; h('Fx2y=Fx2y.from_python(tmp)')
    h.tmp=Fy2x[0,:]; h('Fy2x=Fy2x.from_python(tmp)')
    h.tmp=Fxy[0,:] ; h('Fxy=Fxy.from_python(tmp)')
    
## Convert from hoc to Python
    h('grangernqs=new NQS("F","pp1","pp2","cohe","Fx2y","Fy2x","Fxy")')
    h('grangernqs.setcol("F",F)') # Save the data to the NQS table
    h('grangernqs.setcol("pp1",pp1)')
    h('grangernqs.setcol("pp2",pp2)') 
    h('grangernqs.setcol("cohe",cohe)')
    h('grangernqs.setcol("Fx2y",Fx2y)')
    h('grangernqs.setcol("Fy2x",Fy2x)')
    h('grangernqs.setcol("Fxy",Fxy)')
    grangernqs=h.grangernqs
    
    return grangernqs
    
    
    
    
    
    
# DOWNSAMPLE
def downsample(olddata,oldrate=10000,newrate=200): # Too different from the original code to even call.
    """
    This function downsamples a given vector or matrix.
    
    Usage:
        newdata=pyhoc.downsample(olddata,origrate,newrate)
    where:
        newdata = downsampled data
        olddata = data at original sampling rate
        origrate = original sampling rate (default 10 kHz)
        newrate = desired sampling rate (default 200 Hz)
    
    If olddata has multiple columns, these are assumed to be different time
    series. Thus, an original matrix of N rows by M columns will be downsampled
    to a matrix of N' rows and M columns, where N' = N*origrate/newrate.
    
    Example usage from NEURON is as follows:
        objref output1, output2
        output1=p.pyhoc.downsample(A) // 
        output2=p.pyhoc.downsample(A,10000,1000)
        
        A.size() // = 20000 -- original vector size
        output1.size() // = 400 -- downsampled by a factor of 50
        output2.size() // = 2000 -- downsampled by a factor of 10
        
    Version: 2011apr28
    """
## Load packages
    from scipy import array, shape, size, reshape, zeros
    from neuron import h
    
## Convert data
    ratio=oldrate/float(newrate) # Calculate ratio of sampling rates
    olddata=array(olddata) # Make sure it's an array
    if olddata.ndim==1: olddata=reshape(olddata,(size(olddata,0),1)) # Turn vector into an array
    rows,cols=shape(olddata) # Find out how many rows and columns there are
    newrows=int(rows/ratio); # Calculate how many rows the new file will have
    
## Perform downsampling
    newdata=zeros((newrows,cols)); # Initialize new array
    for i in range(cols): # Loop over time series
        for j in range(newrows): # Loop over new time points
            tstart=int(j*ratio) # This is the starting time of the average
            tfinish=int((j+1)*ratio) # This is the finishing time of the average
            newdata[j,i]=olddata[tstart:tfinish,i].mean(); # Calculate mean across the original time points
        
## Convert from PythonObjet to hoc array
    h('objref tmpinput, tmpoutput')
    h('tmpoutput = new Vector()')
    h.tmpinput=newdata
    h('tmpoutput=tmpoutput.from_python(tmpinput)')
    output=h.tmpoutput
    return output



# SPKLFP
def spklfp():
    """
    This function takes the data structures generated by an
    intfcol-based simulation and uses them to plot every
    quantity of general interest: a spike raster, per-cell
    firing rate histogram, population firing rates, raw LFP
    time series, and LFP spectra.
    
    Usage is as follows:
        p.pyhoc.spklfp()
    
    It requires the following to run:
    	- Matplotlib 1.0.1 or later
    	- NQS table storing LFPs called "nqLFP"
    	- NQS table storing spikes called "snq"
    
    Version: 2011apr28
    """
    print 'Converting data...'
    flattendata() # Convert data
    import spklfp # Run code
    return 0



# SPECTROGRAM
def spectrogram(ts,fs=200,window=2,maxfreq=50,tsmooth=2,fsmooth=2):
    """
    This function takes a given time series and turns it into a spectrogram
    (i.e. a 3-D plot where one axis is time, one is frequency, and one is
    amplitude).
    
    Usage:
        pyhoc.spectrogram(ts,[fs,window,maxfreq,tsmooth,fsmooth])
    where:
        ts = time series to be spectrogrammed
        fs = sampling rate (in Hz)
        window = length of window for computing spectra (in s)
        maxfreq = maximum frequency to plot (in Hz)
        tsmooth = amount of smoothing to do along time axis
        fsmooth = amount of smoothing to do along frequency axis
    
    Example usage from NEURON is as follows:
        p.pyhoc.spectrogram(A)
        
    Version: 2011apr28
    """
    from spectrogram import plotspectrogram
    from pylab import array
    ts=array(ts) 
    plotspectrogram(ts,fs,window,maxfreq,int(tsmooth),int(fsmooth))
    return 0








# VIEWLFPS
def viewlfps(ncols=1,trimdata=1,fs=200,tmax=0,fmax=50,fftsmooth=50,mtpar=4,order=12):
    """
    This function provides an interative way of visualizing LFPs for a particular
    simulation -- it allows you to visualize LFP time series or spectra, the latter
    calculated in one of three ways (plain FFT, multitaper spectrum, or auto-
    regressive fitting via BSMART). 
    
    The non-hoc version allows for the comparison of multiple columns and multiple 
    simulations; however, due to the limitations of the hoc interpreter, only one
    simulation can be viewed at a time in this version. The simulation must have 
    been run in such a way as to generate "nqLFP" (i.e. it must be an intfcol-based
    simulation) , which is then read in by this script.
    
    Usage:
        pyhoc.viewlfps([ncols,trimdata,fs,tmax,fmax,fftsmooth,mtpar,order])
    where:
        ncols = number of columns; M = ncols * number of layers (default 1)
        trimdata = number of seconds' worth of data to trim off each end of the LFP (default 1)
        fs = data sampling rate (default 200 Hz)
        tmax = maximum time to display; set to 0 for all (default 0)
        fmax = maximum frequency to display; set to 0 for all (default 50)
        fftsmooth = the amount of smoothing to do on the FFT (default 50)
        mtpar = the window size for the multitaper method (default 4)
        order = polynomial order for BSMART (default 12)
    
    Example usage from NEURON is as follows:
        p.pyhoc.viewlfps()
    
    Requires:
    	- NQS table storing LFPs called "nqLFP"
    
    Version: 2011apr28
    """
    from pylab import loadtxt, shape
    from viewlfps import plotlfps
    
## Define options
    fs=200.0 # Sampling rate in Hz
    tmax=0 # Maximum time in s
    fmax=50 # Maximum frequency in Hz
    toplot=0 # Which to plot -- 0=time series, 1=plain FFT, 2=multitaper, 3=BSMART
    fftsmooth=50 # How much to smooth the raw FFT -- 50-100 is good
    mtpar=3.5 # The parameter for the multitaper method -- 2-4 is good
    order=10 # The polynomial order for BSMART -- 10-30 is good
    
## Convert data
    print 'Converting data...'
    flattendata() # Convert data
    
## Import data
    print 'Importing data...'
    filenames=['/tmp/pyhoc-lfp.txt'] # Originally /home/cliffk/bill/ckintf/data/1102/juemo/lfp; alternative '/home/cliffk/bill/ckintf/data/1104/07-chrislfp1.txt'
    killdata=trimdata*fs # How much data to cut off each end
    alldata=[]
    alldata.append(loadtxt(filenames[0]))
    npts=shape(alldata[0])[0]
    if npts<=killdata*3: # If killdata is too big for the length of the data, make it smaller
        print 'Warning: trimming data would have result in nothing left!'
        killdata=int(npts/3.)
    alldata[0]=alldata[0][killdata:-killdata-1,:] # Remove bad data
	
## Plot LFPs
    plotlfps(alldata,ncols,fs,tmax,fmax,int(fftsmooth),mtpar,int(order))
    print '...done.'
    return 0






# FLATTENDATA
def flattendata():
    """
    This function is a hoc script to convert NQS tables generated by an
    intfcol-based simulation to a form readable by Python. Not to be used 
    directly by the user. Based on $ckintf/batch.hoc.
    
    Version: 2011apr28
    """
    from neuron import h
    from subprocess import call
    h('oldhz=nqLFP.cob.v.size/tstop*1000 // Original sampling rate; *1000 because tstop is in ms')
    h('newhz=200 // The new frequency to sample at, in Hz')
    h('ratio=oldhz/newhz // Calculate the ratio betwen the old and new sampling rates')
    h('npts=tstop/1000*newhz // Number of points in the resampled time seris')
    h('nlayers=nqLFP.m // Number of layers (usually 5 -- 2/3, 4, 5, 6, all)')
    h('objref tempvec // Temporary place to store NQS column as a vector')
    h('objref tempstr // Name of the NQS column being selected')
    h('objref storelfp // Create matrix to store results in')
    h('storelfp = new Matrix(npts, nlayers*numcols) // Combine layers/columns into one dimension')
    h('count=-1 // Set column of storelfp to zero')
    h('for i=0,numcols-1 { for j=0,nlayers-1 { count+=1 tempstr=nqLFP[i].s[j] tempvec=nqLFP[i].getcol(tempstr.s) for k=0,npts-1 {storelfp.x[k][count]=tempvec.mean(k*ratio,(k+1)*ratio-1)}}}')
    h('objref fobj')
    h('fobj = new File("/tmp/pyhoc-lfp.txt")')
    h('fobj.wopen()')
    h('storelfp.fprint(fobj,"%10.1f") // Its usually in the thousands so one d.p. should do')
    h('fobj.close()')
    h('skipsnq=0 // flag to create NQS with spike times, one per column')
    h('initAllMyNQs() // setup of NQS objects with spike/other information')
    h('objref storespikes, tmpt, tmpid, tmptype, tmpcol // Initialize vectors and matrices -- the tmp vectors are for storing parts of the NQS arrays')
    h('totalnumberofspikes=0 // Calculate the total number of spikes generated across all columns')
    h('for i=0,numcols-1 totalnumberofspikes+=snq[i].cob.v.size')
    h('storespikes = new Matrix(totalnumberofspikes, 4) // Four columns: spike time, cell ID, cell type, and spike time')
    h('count=-1 // Initialize row count')
    h('for i=0,numcols-1 { tmpt=snq[i].getcol("t") tmpid=snq[i].getcol("id") tmptype=snq[i].getcol("type") tmpcol=snq[i].getcol("col") for j=0,snq[i].cob.v.size-1 { count+=1 storespikes.x[count][0]=tmpt.x[j] storespikes.x[count][1]=tmpid.x[j] storespikes.x[count][2]=tmptype.x[j] storespikes.x[count][3]=tmpcol.x[j]}}')
    h('objref fobj2')
    h('fobj2 = new File("/tmp/pyhoc-spk.txt")')
    h('fobj2.wopen()')
    h('storespikes.fprint(fobj2,"%6.0f") // All quantities are integers, so this should be fine')
    h('fobj2.close()')
    call(['sed','-i','1d','/tmp/pyhoc-spk.txt'])
    call(['sed','-i','1d','/tmp/pyhoc-lfp.txt'])
