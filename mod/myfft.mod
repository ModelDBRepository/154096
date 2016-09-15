: $Id: myfft.mod,v 1.3 2009/02/16 15:56:24 samn Exp $ 



NEURON {
 SUFFIX myfft
 GLOBAL INSTALLED
 GLOBAL verbose
}

PARAMETER {
  INSTALLED=0
  verbose=0
}

VERBATIM

/*
	Windowed Sinc FIR Generator
	Bob Maling (BobM.DSP@gmail.com)
	Contributed to musicdsp.org Source Code Archive
	Last Updated: April 12, 2005

	Usage:
		Lowpass:	wsfirLP(h, N, WINDOW, fc)
		Highpass:	wsfirHP(h, N, WINDOW, fc)
		Bandpass:	wsfirBP(h, N, WINDOW, fc1, fc2)
		Bandstop:	wsfirBS(h, N, WINDOW, fc1, fc2)

	where:
		h (double[N]):	filter coefficients will be written to this array
		N (int):		number of taps
		WINDOW (int):	Window (W_BLACKMAN, W_HANNING, or W_HAMMING)
		fc (double):	cutoff (0 < fc < 0.5, fc = f/fs)
						--> for fs=48kHz and cutoff f=12kHz, fc = 12k/48k => 0.25
		
		fc1 (double):	low cutoff (0 < fc < 0.5, fc = f/fs)
		fc2 (double):	high cutoff (0 < fc < 0.5, fc = f/fs)


	Windows included here are Blackman, Hanning, and Hamming. Other windows	can be
	added by doing the following:
		1. "Window type constants" section: Define a global constant for the new window.
		2. wsfirLP() function: Add the new window as a case in the switch() statement.
		3. Create the function for the window
		   
		   For windows with a design parameter, such as Kaiser, some modification
		   will be needed for each function in order to pass the parameter.
*/

#include <math.h>
#include <stdlib.h>

// Function prototypes
void wsfirLP(double h[], int N, int WINDOW, double fc, int Vsz);
void wsfirHP(double h[], int N, int WINDOW, double fc, int Vsz);
void wsfirBS(double h[], int N, int WINDOW, double fc1, double fc2, int Vsz);
void wsfirBP(double h[], int N, int WINDOW, double fc1, double fc2, int Vsz);
void genSinc(double sinc[], int N, double fc, int Vsz);
void wBlackman(double w[], int N, int Vsz);
void wHanning(double w[], int N, int Vsz);
void wHamming(double w[], int N, int Vsz);

// Window type contstants
#define W_BLACKMAN 1
#define W_HANNING  2
#define W_HAMMING  3

// Generate lowpass filter
// 
// This is done by generating a sinc function and then windowing it
void wsfirLP(double h[],		// h[] will be written with the filter coefficients
			 int N,		// size of the filter (number of taps)
			 int WINDOW,	// window function (W_BLACKMAN, W_HANNING, etc.)
			 double fc,	// cutoff frequency
                         int Vsz)
{
	int i;
	double *w = (double*) calloc(N,sizeof(double));		// window function
        double *sinc = (double*) calloc(N,sizeof(double));	// sinc function
    
	// 1. Generate Sinc function
	genSinc(sinc, N, fc, N);
    
	// 2. Generate Window function
	switch (WINDOW) {
		case W_BLACKMAN:	// W_BLACKMAN
			wBlackman(w, N, N);
			break;
		case W_HANNING:		// W_HANNING
			wHanning(w, N, N);
			break;
		case W_HAMMING:		// W_HAMMING
			wHamming(w, N, N);
			break;
		default:
			break;
	}

	// 3. Make lowpass filter
	for (i = 0; i < N; i++) {
		h[i] = sinc[i] * w[i];
	}

        for ( ; i < Vsz; i++) {
                h[i] = 0.;
        }

	// Delete dynamic storage
        free(w);
        free(sinc);

	return;
}

void lpwin(void* vv) {
  double* x, fc;
  int fsz,win,vsz;
  vsz = vector_instance_px(vv,&x); 
  fc = *getarg(1);
  fsz = ifarg(2) ? (int)*getarg(2) : vsz;
  win= ifarg(3) ? (int)*getarg(3) : W_BLACKMAN;
  wsfirLP(x,fsz,win,fc,vsz);
}

// Generate highpass filter
//
// This is done by generating a lowpass filter and then spectrally inverting it
void wsfirHP(double h[],		// h[] will be written with the filter coefficients
			 int N,		// size of the filter
			 int WINDOW,	// window function (W_BLACKMAN, W_HANNING, etc.)
			 double fc,	// cutoff frequency
                         int Vsz)
{
	int i;

	// 1. Generate lowpass filter
	wsfirLP(h, N, WINDOW, fc, Vsz);

	// 2. Spectrally invert (negate all samples and add 1 to center sample) lowpass filter
	// = delta[n-((N-1)/2)] - h[n], in the time domain
	for (i = 0; i < N; i++) {
		h[i] *= -1.0;	// = 0 - h[i]
	}
	h[(N-1)/2] += 1.0;	// = 1 - h[(N-1)/2]

	return;
}

void hpwin(void* vv) {
  double* x, fc;
  int fsz,win,vsz;
  vsz = vector_instance_px(vv,&x); 
  fc = *getarg(1);
  fsz = ifarg(2) ? (int)*getarg(2) : vsz;
  win = ifarg(3) ? (int)*getarg(3) : W_BLACKMAN;
  wsfirHP(x,fsz,win,fc,vsz);
}


// Generate bandstop filter
//
// This is done by generating a lowpass and highpass filter and adding them
void wsfirBS(double h[],		// h[] will be written with the filter taps
			 int N,		// size of the filter
			 int WINDOW,	// window function (W_BLACKMAN, W_HANNING, etc.)
			 double fc1,	// low cutoff frequency
			 double fc2,	// high cutoff frequency
                         int Vsz)
{
	int i;
	double *h1 = (double*) calloc(N,sizeof(double));
	double *h2 = (double*) calloc(N,sizeof(double));

	// 1. Generate lowpass filter at first (low) cutoff frequency
	wsfirLP(h1, N, WINDOW, fc1, N);

	// 2. Generate highpass filter at second (high) cutoff frequency
	wsfirHP(h2, N, WINDOW, fc2, N);

	// 3. Add the 2 filters together
	for (i = 0; i < N; i++) {
		h[i] = h1[i] + h2[i];
	}

        for (; i < Vsz; i++) {
               h[i] = 0.;
        }

	// Delete dynamic memory
        free(h1);
        free(h2);

	return;
}

void bswin(void* vv) {
  double* x, fc1, fc2;
  int fsz,win,vsz;
  vsz = vector_instance_px(vv,&x); 
  fc1 = *getarg(1);
  fc2 = *getarg(2);
  fsz = ifarg(3) ? (int)*getarg(3) : vsz;
  win = ifarg(4) ? (int)*getarg(4) : W_BLACKMAN;
  wsfirBS(x,fsz,win,fc1,fc2,vsz);
}


// Generate bandpass filter
//
// This is done by generating a bandstop filter and spectrally inverting it
void wsfirBP(double h[],		// h[] will be written with the filter taps
			 int N,		// size of the filter
			 int WINDOW,	// window function (W_BLACKMAN, W_HANNING, etc.)
			 double fc1,	// low cutoff frequency
			 double fc2,	// high cutoff frequency
                         int Vsz)
{
	int i;

	// 1. Generate bandstop filter
	wsfirBS(h, N, WINDOW, fc1, fc2, Vsz);

	// 2. Spectrally invert bandstop (negate all, and add 1 to center sample)
	for (i = 0; i < N; i++) {
		h[i] *= -1.0;	// = 0 - h[i]
	}
	h[(N-1)/2] += 1.0;	// = 1 - h[(N-1)/2]

	return;
}

void bpwin(void* vv) {
  double* x, fc1, fc2;
  int fsz,win,vsz;
  vsz = vector_instance_px(vv,&x); 
  fc1 = *getarg(1);
  fc2 = *getarg(2);
  fsz = ifarg(3) ? (int)*getarg(3) : vsz;
  win = ifarg(4) ? (int)*getarg(4) : W_BLACKMAN;
  wsfirBP(x,fsz,win,fc1,fc2,vsz);
}


// Generate sinc function - used for making lowpass filter
void genSinc(double sinc[],		// sinc[] will be written with the sinc function
			 int N,		// size (number of taps)
			 double fc,	// cutoff frequency
                         int Vsz)
{
	int i;
	double M = N-1;
	double n;

	// Constants
	double PI = 3.14159265358979323846;

	// Generate sinc delayed by (N-1)/2
	for (i = 0; i < N; i++) {
		if (i == M/2.0) {
			sinc[i] = 2.0 * fc;
		}
		else {
			n = (double)i - M/2.0;
			sinc[i] = sin(2.0*PI*fc*n) / (PI*n);
		}
	}        

        for (; i < Vsz; i++) {
               sinc[i] = 0.;
        }

	return;
}

void sincwin(void* vv) {
  double* x, fc;
  int fsz,vsz;
  vsz = vector_instance_px(vv,&x); 
  fc = *getarg(1);
  fsz = ifarg(2) ? (int) *getarg(2) : vsz;
  genSinc(x,fsz,fc,vsz);
}


// Generate window function (Blackman)
void wBlackman(double w[],		// w[] will be written with the Blackman window
			   int N,	// size of the window
                           int Vsz)
{
	int i;
	double M = N-1;
	double PI = 3.14159265358979323846;

	for (i = 0; i < N; i++) {
		w[i] = 0.42 - (0.5 * cos(2.0*PI*(double)i/M)) + (0.08*cos(4.0*PI*(double)i/M));
	}

        for (; i < Vsz; i++) {
                w[i] = 0.;
        }

	return;
}

void blackmanwin(void* vv) {
  double* x;
  int fsz,vsz;
  vsz = vector_instance_px(vv,&x); 
  fsz = ifarg(1) ? (int) *getarg(1) : vsz;
  wBlackman(x,fsz,vsz);
}


// Generate window function (Hanning)
void wHanning(double w[],		// w[] will be written with the Hanning window
			  int N,		// size of the window
                          int Vsz)      //size of output buffer (extra space for zero padding)
{
	int i;
	double M = N-1;
	double PI = 3.14159265358979323846;

	for (i = 0; i < N; i++) {
		w[i] = 0.5 * (1.0 - cos(2.0*PI*(double)i/M));
	}

        for (; i < Vsz; i++) {
                w[i] = 0.;
        }

	return;
}

void hanningwin(void* vv) {
  double* x;
  int fsz,vsz;
  vsz = vector_instance_px(vv,&x); 
  fsz = ifarg(1) ? (int) *getarg(1) : vsz;
  wHanning(x,fsz,vsz);
}


// Generate window function (Hamming)
void wHamming(double w[],		// w[] will be written with the Hamming window
			  int N,		// size of the window
                          int VSz)       //size of output buffer (extra space for zero padding)
{
	int i;
	double M = N-1;
	double PI = 3.14159265358979323846;

	for (i = 0; i < N; i++) {
		w[i] = 0.54 - (0.46*cos(2.0*PI*(double)i/M));
	}

        for (; i < VSz; i++) {
                w[i] = 0.;
        }

	return;
}

void hammingwin(void* vv) {
  double* x;
  int fsz,vsz;
  vsz = vector_instance_px(vv,&x); 
  fsz = ifarg(1) ? (int) *getarg(1) : vsz;
  wHamming(x,fsz,vsz);
}

double* wrap(double* x,int n,int flen){
  double* y = (double*) calloc(n,sizeof(double));
  int i,j=0;
  for(i=flen/2+1;i<flen;i++)    y[j++]=x[i];
  j=n-flen/2-1;
  for(i=0;i<=flen/2;i++)   y[j++]=x[i];
  return y;
}

void wraparound(void* vv) {
  double* x,*y;
  int vsz,fsz,i;
  vsz = vector_instance_px(vv,&x);
  fsz = (int) *getarg(1);
  if(fsz > vsz) {
    printf("wraparound ERRA: invalid filter len %d > vector len %d!\n",fsz,vsz);
    return;
  }
  y = wrap(x,vsz,fsz);
  for(i=0;i<vsz;i++) x[i]=y[i];
  free(y);
}


/*************************************************************************
 *                                                                       *
 *               ROUTINES IN THIS FILE:                                  *
 *                                                                       *
 *                      fft_float(): calling routine for complex fft     *
 *                                   of a real sequence                  *
 *                                                                       *
 *                      fft_pow(): calling routine for complex fft       *
 *                                   of a real sequence that concludes   *
 *                                   by computing power spectrum         *
 *                                                                       *
 *                      FAST(): actual fft routine                       *
 *                                                                       *
 * 			FR2TR(): radix 2 transform                       *
 *                                                                       *
 * 			FR4TR(): radix 4 transform                       *
 *                                                                       *
 * 			FORD1(): re-ordering routine                     *
 *                                                                       *
 * 			FORD2(): other re-ordering routine               *
 *                                                                       *
 *  			fastlog2(): just what it sounds like             *
 *                                                                       *
 ************************************************************************/

/*
** Discrete Fourier analysis routine
** from IEEE Programs for Digital Signal Processing
** G. D. Bergland and M. T. Dolan, original authors
** Translated from the FORTRAN with some changes by Paul Kube
**
** Modified to return the power spectrum by Chuck Wooters
**
** Modified again by Tony Robinson (ajr@eng.cam.ac.uk) Dec 92

** Slight naming mods by N. Morgan, July 1993
	(fft_chuck -> fft_pow)
	( calling args long ll -> long winlength)
			(long m -> long log2length)
*/
#include <math.h>
#include <stdio.h>

typedef double real;

extern void four1(double mdata[], unsigned long nn, int isign);
extern void twofft(double data1[], double data2[], double fft1[], double fft2[],
	    unsigned long n);
extern void realft(double mdata[], unsigned long n, int isign);
extern void convlv(double mdata[], unsigned long n, double respns[], unsigned long m,
	    int isign, double ans[]);
extern void correl(double data1[], double data2[], unsigned long n, double ans[]);
extern void spctrm(double mdata[], double p[], int m, int k);

/*
  convlv()  -- convolution of a read data set with a response function

  the respns function must be stored in wrap-around order

  N.R.C  p. 430
*/

void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	hoc_execerror("Numerical Recipes run-time error...\n", error_text);
}

/*
  twofft()  -- discrete FFT of two real functions simultaneously
  N.R.C  p. 414
*/

void twofft(double data1[], double data2[], double fft1[], double fft2[],
	unsigned long n)
{
	void nrfour1(double mdata[], unsigned long nn, int isign);
	unsigned long nn3,nn2,jj,j;
	double rep,rem,aip,aim;

	nn3=1+(nn2=2+n+n);
	for (j=1,jj=2;j<=n;j++,jj+=2) {
		fft1[jj-1]=data1[j];
		fft1[jj]=data2[j];
	}
	nrfour1(fft1,n,1);
	fft2[1]=fft1[2];
	fft1[2]=fft2[2]=0.0;
	for (j=3;j<=n+1;j+=2) {
		rep=0.5*(fft1[j]+fft1[nn2-j]);
		rem=0.5*(fft1[j]-fft1[nn2-j]);
		aip=0.5*(fft1[j+1]+fft1[nn3-j]);
		aim=0.5*(fft1[j+1]-fft1[nn3-j]);
		fft1[j]=rep;
		fft1[j+1]=aim;
		fft1[nn2-j]=rep;
		fft1[nn3-j] = -aim;
		fft2[j]=aip;
		fft2[j+1] = -rem;
		fft2[nn2-j]=aip;
		fft2[nn3-j]=rem;
	}
}

/*
  realft()  -- discrete FFT of a real function with 2n data pts
  N.R.C  p. 417
*/

void realft(double mdata[], unsigned long n, int isign)
{
	void nrfour1(double mdata[], unsigned long nn, int isign);
	unsigned long i,i1,i2,i3,i4,np3;
	double c1=0.5,c2,h1r,h1i,h2r,h2i;
	double wr,wi,wpr,wpi,wtemp,theta;

	theta=3.141592653589793/(double) (n>>1);
	if (isign == 1) {
		c2 = -0.5;
		nrfour1(mdata,n>>1,1);
	} else {
		c2=0.5;
		theta = -theta;
	}
	wtemp=sin(0.5*theta);
	wpr = -2.0*wtemp*wtemp;
	wpi=sin(theta);
	wr=1.0+wpr;
	wi=wpi;
	np3=n+3;
	for (i=2;i<=(n>>2);i++) {
		i4=1+(i3=np3-(i2=1+(i1=i+i-1)));
		h1r=c1*(mdata[i1]+mdata[i3]);
		h1i=c1*(mdata[i2]-mdata[i4]);
		h2r = -c2*(mdata[i2]+mdata[i4]);
		h2i=c2*(mdata[i1]-mdata[i3]);
		mdata[i1]=h1r+wr*h2r-wi*h2i;
		mdata[i2]=h1i+wr*h2i+wi*h2r;
		mdata[i3]=h1r-wr*h2r+wi*h2i;
		mdata[i4] = -h1i+wr*h2i+wi*h2r;
		wr=(wtemp=wr)*wpr-wi*wpi+wr;
		wi=wi*wpr+wtemp*wpi+wi;
	}
	if (isign == 1) {
		mdata[1] = (h1r=mdata[1])+mdata[2];
		mdata[2] = h1r-mdata[2];
	} else {
		mdata[1]=c1*((h1r=mdata[1])+mdata[2]);
		mdata[2]=c1*(h1r-mdata[2]);
		nrfour1(mdata,n>>1,-1);
	}
}


static double nrsqrarg;
#define SQR(a) ((nrsqrarg=(a)) == 0.0 ? 0.0 : nrsqrarg*nrsqrarg)
void convlv(double mdata[], unsigned long n, double respns[], unsigned long m,
	int isign, double ans[])
{
	void realft(double mdata[], unsigned long n, int isign);
	void twofft(double data1[], double data2[], double fft1[], double fft2[],
		unsigned long n);
	unsigned long i,no2;
	double dum,mag2,*fft;

	fft=nrvector(1,n<<1);
	for (i=1;i<=(m-1)/2;i++)
		respns[n+1-i]=respns[m+1-i];
	for (i=(m+3)/2;i<=n-(m-1)/2;i++)
		respns[i]=0.0;
	twofft(mdata,respns,fft,ans,n);
	no2=n>>1;
	for (i=2;i<=n+2;i+=2) {
		if (isign == 1) {
			ans[i-1]=(fft[i-1]*(dum=ans[i-1])-fft[i]*ans[i])/no2;
			ans[i]=(fft[i]*dum+fft[i-1]*ans[i])/no2;
		} else if (isign == -1) {
			if ((mag2=SQR(ans[i-1])+SQR(ans[i])) == 0.0)
				nrerror("Deconvolving at response zero in convlv");
			ans[i-1]=(fft[i-1]*(dum=ans[i-1])+fft[i]*ans[i])/mag2/no2;
			ans[i]=(fft[i]*dum-fft[i-1]*ans[i])/mag2/no2;
		} else nrerror("No meaning for isign in convlv");
	}
	ans[2]=ans[n+1];
	realft(ans,n,-1);
	nrfree_vector(fft,1,n<<1);
}



double rconvlv(double* mdata,int n,double* respns,int m,double* ans) {
//  Vect* v3 = (Vect*)v;

//  Vect* v1, *v2;

  // data set
//  v1 = vector_arg(1);


  // filter
//  v2 = vector_arg(2);

  // convolve unless isign is -1, then deconvolve!
  int isign = 1;
//  if (ifarg(3)) isign = (int)(*getarg(3)); else isign = 1;

  // make both data sets equal integer power of 2
//  int v1n = v1->capacity();
//  int v2n = v2->capacity();
//  int m = (nIN > mIN) ? nIN : mIN;
//  int n = 1;
//  while(n < m) n*=2;

//  double *data = (double *)calloc(n,(unsigned)sizeof(double));
//  int i;
//  for (i=0;i<v1n;++i) data[i] = v1->elem(i);

  // assume respns is given in "wrap-around" order,
  // with countup t=0..t=n/2 followed by countdown t=n..t=n/2
  // v2n should be an odd <= n

//  double *respns = (double *)calloc(n,(unsigned)sizeof(double));
//  for (i=0;i<v2n;i++) respns[i] = v2->elem(i);

//  double *ans = (double *)calloc(2*n,(unsigned)sizeof(double));

//  convlv(&data[0]-1,n,&respns[0]-1,v2n,isign,&ans[0]-1);
  convlv(&mdata[0]-1,n,&respns[0]-1,m,isign,&ans[0]-1);

//  if (v3->capacity() != n) v3->resize(n);
//  for (i=0;i<n;++i) v3->elem(i)=ans[i];

//  free((char *)data);
//  free((char *)respns);
//  free((char *)ans);

//  return v3->temp_objvar();

  return 1;
}


#define PI  3.1415926535897932
#define PI8 0.392699081698724 /* PI / 8.0 */
#define RT2 1.4142135623731  /* sqrt(2.0) */
#define IRT2 0.707106781186548  /* 1.0/sqrt(2.0) */

#define signum(i) (i < 0 ? -1 : i == 0 ? 0 : 1)

int  FAST(real*, int);
void FR2TR(int, real*, real*);
void FR4TR(int, int, real*, real*, real*, real*);
void FORD1(int, real*);
void FORD2(int, real*);
int  fastlog2(int);

double two_to_the(double N){
  return (int)(pow(2.0,(N))+0.5);
}

int fft_real(real *orig, real *fftd, int npoint) {
  int i;

  for(i = 0; i< npoint; i++) fftd[i] = orig[i];

  if(FAST(fftd, npoint) == 0 ){
    fprintf(stderr, "Error calculating fft.\n");
    return 0;
  }
  return 1;
}

/*
double bandfilt(double* p,int sz,int lohz,int hihz){

  double dret=0.0;

  int inv = 1;

  // make data set integer power of 2
  int n = 1;
  while(n < sz) n*=2;

  double *mdata = (double *)calloc(n,(unsigned)sizeof(double));
  int i;
  for (i=0;i<sz;++i) mdata[i] = p[i];

  realft(&mdata[0]-1,n,inv);

  if (v3->capacity() != n) v3->resize(n);
  for (i=0;i<n;++i) v3->elem(i)=mdata[i];

  free((char *)mdata);

  return v3->temp_objvar();



  dret=1.0;
  return dret;
}
*/
double fftband(void* v){
  double* x; //input vector
  int n = vector_instance_px(v,&x); 
  double dret = 0.;
  double* pout = 0;
  int outlen = vector_arg_px(1,&pout); //size of filtered output
  int lowhz = (int)*getarg(2); //lowest Hz we want to keep
  int hihz = (int)*getarg(3); //highest Hz we want to keep
  int log2length = ceil(log((double)(n))/log(2.0));
  int fftlen = two_to_the((double)log2length);
  double* fftd = (double*)calloc(fftlen,sizeof(double)); //stores fft results
  if(outlen  < fftlen){
    printf("fftband ERRA: outlen %d < %d",outlen,fftlen);
    goto FBCLEANUP;
  }
  if(!fft_real(x,fftd,n)){
    printf("fftband ERRB: couldn't perform fft\n");
    goto FBCLEANUP;
  }
  int hz = 0 , di = 0;
  int npoints  = (int) (pow(2.0,(real) log2length) + 0.5);
  int npoints2 = npoints / 2;
//  fftd[0]=0.; fftd[npoints2]=0.;
fftd[0]=fftd[1]=0.;
  int i,j,k;
  for(i=1;k<npoints2;i++){
    j=2*i;
    k=2*i+1;
    if(i<lowhz){
      di = lowhz - i;
      if(di==1){
        fftd[j] /= 2.; fftd[k] /= 2.;
      } else if(di==2){
        fftd[j] /= 4.; fftd[k] /= 4.;
      //} else if(di==3){
       // fftd[j] /= 8.; fftd[k] /= 8.;
      } else {
        fftd[j] = 0.; fftd[k] = 0.;
      }
    } else if(i>hihz){
      di = i - hihz;
      if(di==1){
        fftd[j] /= 2.; fftd[k] /= 2.;
      } else if(di==2){
        fftd[j] /= 4.; fftd[k] /= 4.;
      //} else if(di==3){
       // fftd[j] /= 8.; fftd[k] /= 8.;
      } else {
        fftd[j] = 0.; fftd[k] = 0.;
      }
    }
  }
//  if(!fft_real(fftd,pout,fftlen)){
  if(!fft_real(fftd,pout,n)){
    printf("fftbadn ERRC: couldn't perform fft\n");
    goto FBCLEANUP;
  }
  dret=1.;
FBCLEANUP:
  if(fftd) free(fftd);
  return dret;
}

int fft_pow(real *orig, real *power, long winlength, long log2length) {
  int i, j, k;
  real *temp = NULL;
  int npoints, npoints2;

  npoints  = (int) (pow(2.0,(real) log2length) + 0.5);
  npoints2 = npoints / 2;
  temp = (real*) malloc(npoints * sizeof(real));
  if(temp == (real*) NULL) {
    fprintf(stderr, "Error mallocing memory in fft_pow()\n");
    return 0;
  }
  for(i=0;i<winlength;i++) temp[i] = (real) orig[i];
  for(i = winlength; i < npoints; i++) temp[i] = 0.0; //0-padding at end

  if(FAST(temp, npoints) == 0 ){
    fprintf(stderr,"Error calculating fft.\n");
    free(temp);
    return 0;
  }
  // convert the complex data to power 
  power[0] = temp[0]*temp[0];
  power[npoints2] = temp[1]*temp[1];
  //  Only the first half of the power[] array is filled with data.
  //  The second half would just be a mirror image of the first half.
  for(i=1;i<npoints2;i++){
    j=2*i;
    k=2*i+1;
    power[i] = temp[j]*temp[j]+temp[k]*temp[k];
  }
  free(temp);
  return 1;
}

int fastlog2(int n);

double tlog2(double d){
  static double l2 = 0.693147181; 
  return log(d)/l2;
}

double dfftpow(double* x,int n,double* ppow,int powlen,int* fftlen){
  int log2length = ceil(log((double)(n))/log(2.0));
  *fftlen = two_to_the((double)log2length);
  if(powlen < *fftlen/2 + 1){
    printf("dfftpow ERRA: powlen=%d < fftlength/2+1 = %d\n",powlen,*fftlen/2+1);
    return 0.0;
  }
  double dret =  (double) fft_pow(x,ppow,n,log2length);
  return dret;
}

double fftpow(void* v){
  double* x;
  int n = vector_instance_px(v,&x);
  double* ppow;
  int powlen = vector_arg_px(1,&ppow);
  int fftlen=0;
  if(dfftpow(x,n,ppow,powlen,&fftlen)){
    vector_resize(vector_arg(1),fftlen/2+1);
    return 1.0;
  }
  return 0.0;
}

/*
** FAST(b,n)
** This routine replaces the real vector b
** of length n with its finite discrete fourier transform.
** DC term is returned in b[0];
** n/2th harmonic real part in b[1].
** jth harmonic is returned as complex number stored as
** b[2*j] + i b[2*j + 1]
** (i.e., remaining coefficients are as a DPCOMPLEX vector).
**
*/
int FAST(real *b, int n) {
  real fn;
  int i, in, nn, n2pow, n4pow, nthpo;

  n2pow = fastlog2(n);
  if(n2pow <= 0) return 0;
  nthpo = n;
  fn = nthpo;
  n4pow = n2pow / 2;

  /* radix 2 iteration required; do it now */
  if(n2pow % 2) {
    nn = 2;
    in = n / nn;
    FR2TR(in, b, b + in);
  }
  else nn = 1;

  /* perform radix 4 iterations */
  for(i = 1; i <= n4pow; i++) {
    nn *= 4;
    in = n / nn;
    FR4TR(in, nn, b, b + in, b + 2 * in, b + 3 * in);
  }

  /* perform inplace reordering */
  FORD1(n2pow, b);
  FORD2(n2pow, b);

  /* take conjugates */
  for(i = 3; i < n; i += 2) b[i] = -b[i];

  return 1;
}

/* radix 2 subroutine */
void FR2TR(int in, real *b0, real *b1) {
  int k;
  real mt;
  for(k = 0; k < in; k++) {
    mt = b0[k] + b1[k];
    b1[k] = b0[k] - b1[k];
    b0[k] = mt;
  }
}

/* radix 4 subroutine */
void FR4TR(int in, int nn, real *b0, real *b1, real *b2, real* b3) {
  real arg, piovn, th2;
  real *b4 = b0, *b5 = b1, *b6 = b2, *b7 = b3;
  real t0, t1, t2, t3, t4, t5, t6, t7;
  real r1, r5, pr, pi;
  real c1, c2, c3, s1, s2, s3;

  int j, k, jj, kk, jthet, jlast, ji, jl, jr, int4;
  int L[16], L1, L2, L3, L4, L5, L6, L7, L8, L9, L10, L11, L12, L13, L14, L15;
  int j0, j1, j2, j3, j4, j5, j6, j7, j8, j9, j10, j11, j12, j13, j14;
  int k0, kl;

  L[1] = nn / 4;
  for(k = 2; k < 16; k++) {  /* set up L's */
    switch (signum(L[k-1] - 2)) {
    case -1:
      L[k-1]=2;
    case 0:
      L[k]=2;
      break;
    case 1:
      L[k]=L[k-1]/2;
    }
  }

  L15=L[1]; L14=L[2]; L13=L[3]; L12=L[4]; L11=L[5]; L10=L[6]; L9=L[7];
  L8=L[8];  L7=L[9];  L6=L[10]; L5=L[11]; L4=L[12]; L3=L[13]; L2=L[14];
  L1=L[15];

  piovn = PI / nn;
  ji=3;
  jl=2;
  jr=2;

  for(j1=2;j1<=L1;j1+=2)
  for(j2=j1;j2<=L2;j2+=L1)
  for(j3=j2;j3<=L3;j3+=L2)
  for(j4=j3;j4<=L4;j4+=L3)
  for(j5=j4;j5<=L5;j5+=L4)
  for(j6=j5;j6<=L6;j6+=L5)
  for(j7=j6;j7<=L7;j7+=L6)
  for(j8=j7;j8<=L8;j8+=L7)
  for(j9=j8;j9<=L9;j9+=L8)
  for(j10=j9;j10<=L10;j10+=L9)
  for(j11=j10;j11<=L11;j11+=L10)
  for(j12=j11;j12<=L12;j12+=L11)
  for(j13=j12;j13<=L13;j13+=L12)
  for(j14=j13;j14<=L14;j14+=L13)
  for(jthet=j14;jthet<=L15;jthet+=L14)
    {
      th2 = jthet - 2;
      if(th2<=0.0)
	{
	  for(k=0;k<in;k++)
	    {
	      t0 = b0[k] + b2[k];
	      t1 = b1[k] + b3[k];
	      b2[k] = b0[k] - b2[k];
	      b3[k] = b1[k] - b3[k];
	      b0[k] = t0 + t1;
	      b1[k] = t0 - t1;
	    }
	  if(nn-4>0)
	    {
	      k0 = in*4 + 1;
	      kl = k0 + in - 1;
	      for (k=k0;k<=kl;k++)
		{
		  kk = k-1;
		  pr = IRT2 * (b1[kk]-b3[kk]);
		  pi = IRT2 * (b1[kk]+b3[kk]);
		  b3[kk] = b2[kk] + pi;
		  b1[kk] = pi - b2[kk];
		  b2[kk] = b0[kk] - pr;
		  b0[kk] = b0[kk] + pr;
		}
	    }
	}
      else
	{
	  arg = th2*piovn;
	  c1 = cos(arg);
	  s1 = sin(arg);
	  c2 = c1*c1 - s1*s1;
	  s2 = c1*s1 + c1*s1;
	  c3 = c1*c2 - s1*s2;
	  s3 = c2*s1 + s2*c1;

	  int4 = in*4;
	  j0=jr*int4 + 1;
	  k0=ji*int4 + 1;
	  jlast = j0+in-1;
	  for(j=j0;j<=jlast;j++)
	    {
	      k = k0 + j - j0;
	      kk = k-1; jj = j-1;
	      r1 = b1[jj]*c1 - b5[kk]*s1;
	      r5 = b1[jj]*s1 + b5[kk]*c1;
	      t2 = b2[jj]*c2 - b6[kk]*s2;
	      t6 = b2[jj]*s2 + b6[kk]*c2;
	      t3 = b3[jj]*c3 - b7[kk]*s3;
	      t7 = b3[jj]*s3 + b7[kk]*c3;
	      t0 = b0[jj] + t2;
	      t4 = b4[kk] + t6;
	      t2 = b0[jj] - t2;
	      t6 = b4[kk] - t6;
	      t1 = r1 + t3;
	      t5 = r5 + t7;
	      t3 = r1 - t3;
	      t7 = r5 - t7;
	      b0[jj] = t0 + t1;
	      b7[kk] = t4 + t5;
	      b6[kk] = t0 - t1;
	      b1[jj] = t5 - t4;
	      b2[jj] = t2 - t7;
	      b5[kk] = t6 + t3;
	      b4[kk] = t2 + t7;
	      b3[jj] = t3 - t6;
	    }
	  jr += 2;
	  ji -= 2;
	  if(ji-jl <= 0) {
	    ji = 2*jr - 1;
	    jl = jr;
	  }
	}
    }
}

/* an inplace reordering subroutine */
void FORD1(int m, real *b) {
  int j, k = 4, kl = 2, n = 0x1 << m;
  real mt;

  for(j = 4; j <= n; j += 2) {
    if(k - j>0) {
      mt = b[j-1];
      b[j - 1] = b[k - 1];
      b[k - 1] = mt;
    }
    k -= 2;
    if(k - kl <= 0) {
      k = 2*j;
      kl = j;
    }
  }
}

/*  the other inplace reordering subroutine */
void FORD2(int m, real *b) {
  real mt;

  int n = 0x1<<m, k, ij, ji, ij1, ji1;

  int l[16], l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11, l12, l13, l14, l15;
  int j1, j2, j3, j4, j5, j6, j7, j8, j9, j10, j11, j12, j13, j14;


  l[1] = n;
  for(k=2;k<=m;k++) l[k]=l[k-1]/2;
  for(k=m;k<=14;k++) l[k+1]=2;

  l15=l[1];l14=l[2];l13=l[3];l12=l[4];l11=l[5];l10=l[6];l9=l[7];
  l8=l[8];l7=l[9];l6=l[10];l5=l[11];l4=l[12];l3=l[13];l2=l[14];l1=l[15];

  ij = 2;

  for(j1=2;j1<=l1;j1+=2)
  for(j2=j1;j2<=l2;j2+=l1)
  for(j3=j2;j3<=l3;j3+=l2)
  for(j4=j3;j4<=l4;j4+=l3)
  for(j5=j4;j5<=l5;j5+=l4)
  for(j6=j5;j6<=l6;j6+=l5)
  for(j7=j6;j7<=l7;j7+=l6)
  for(j8=j7;j8<=l8;j8+=l7)
  for(j9=j8;j9<=l9;j9+=l8)
  for(j10=j9;j10<=l10;j10+=l9)
  for(j11=j10;j11<=l11;j11+=l10)
  for(j12=j11;j12<=l12;j12+=l11)
  for(j13=j12;j13<=l13;j13+=l12)
  for(j14=j13;j14<=l14;j14+=l13)
  for(ji=j14;ji<=l15;ji+=l14) {
    ij1 = ij-1; ji1 = ji - 1;
    if(ij-ji<0) {
      mt = b[ij1-1];
      b[ij1-1]=b[ji1-1];
      b[ji1-1] = mt;

      mt = b[ij1];
      b[ij1]=b[ji1];
      b[ji1] = mt;
    }
    ij += 2;
  }
}

int fastlog2(int n) {
  int num_bits, power = 0;

  if((n < 2) || (n % 2 != 0)) return(0);
  num_bits = sizeof(int) * 8;   /* How big are ints on this machine? */

  while(power <= num_bits) {
    n >>= 1;
    power += 1;
    if(n & 0x01) {
      if(n > 1)	return(0);
      else return(power);
    }
  }
  return(0);
}
ENDVERBATIM

FUNCTION hfastlog2(){
  VERBATIM
  return fastlog2((int)*getarg(1));
  ENDVERBATIM
}

FUNCTION log2(){
  VERBATIM
  static double l2 = 0.693147181; 
  return log(*getarg(1))/l2;
  ENDVERBATIM
}

PROCEDURE install () {
  if (INSTALLED==1) {
    printf("already installed myfft.mod")
  } else {
    INSTALLED=1
    VERBATIM
    install_vector_method("fftpow", fftpow);
    install_vector_method("fftband",fftband);
    install_vector_method("lpwin",lpwin);
    install_vector_method("hpwin",hpwin);
    install_vector_method("bswin",bswin);
    install_vector_method("bpwin",bpwin);
    install_vector_method("sincwin",sincwin);
    install_vector_method("blackmanwin",blackmanwin);
    install_vector_method("hanningwin",hanningwin);
    install_vector_method("hammingwin",hammingwin);
    install_vector_method("wraparound",wraparound);
    ENDVERBATIM
  }
}
