: $Id: tsa.mod,v 1.18 2011/10/25 21:45:39 billl Exp $ 


NEURON {
 SUFFIX tsa
 GLOBAL INSTALLED
 GLOBAL verbose
}

PARAMETER {
  INSTALLED=0
  verbose=0
}

VERBATIM
#include "misc.h"
static double *x1x, *y1y, *z1z;
//void spctrm(double data[], double p[], int m, int k);
extern double dfftpow(double* x,int n,double* ppow,int powlen,int* fftlen);

#define WINDOW(j,a,b) (1.0-fabs((((j)-1)-(a))*(b)))       /* Bartlett */

double *nrvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+1)*sizeof(double)));
	if (!v) { nrerror("allocation failure in vector()"); return 0x0;}
	return v-nl+1;
}

void nrfree_vector(double *v, long nl, long nh)
/* free a double vector allocated with vector() */
{
  free((char*) (v+nl-1));
}

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

static double nrsqrarg;
#define SQR(a) ((nrsqrarg=(a)) == 0.0 ? 0.0 : nrsqrarg*nrsqrarg)


void nrfour1(double mdata[], unsigned long nn, int isign)
{
	unsigned long n,mmax,m,j,istep,i;
	double wtemp,wr,wpr,wpi,wi,theta;
	double tempr,tempi;

	n=nn << 1;
	j=1;
	for (i=1;i<n;i+=2) {
		if (j > i) {
			SWAP(mdata[j],mdata[i]);
			SWAP(mdata[j+1],mdata[i+1]);
		}
		m=n >> 1;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}
	mmax=2;
	while (n > mmax) {
		istep=mmax << 1;
		theta=isign*(6.28318530717959/mmax);
		wtemp=sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi=sin(theta);
		wr=1.0;
		wi=0.0;
		for (m=1;m<mmax;m+=2) {
			for (i=m;i<=n;i+=istep) {
				j=i+mmax;
				tempr=wr*mdata[j]-wi*mdata[j+1];
				tempi=wr*mdata[j+1]+wi*mdata[j];
				mdata[j]=mdata[i]-tempr;
				mdata[j+1]=mdata[i+1]-tempi;
				mdata[i] += tempr;
				mdata[i+1] += tempi;
			}
			wr=(wtemp=wr)*wpr-wi*wpi+wr;
			wi=wi*wpr+wtemp*wpi+wi;
		}
		mmax=istep;
	}
}

void nrspctrm(double mdata[], double p[], int m, int k)
{
  // assume overlap = 1
	void nrfour1(double mdata[], unsigned long nn, int isign);
	int mm,m44,m43,m4,kk,joffn,joff,j2,j,c=0;
	double w,facp,facm,*w1,*w2,sumw=0.0,den=0.0;

	mm=m+m;
	m43=(m4=mm+mm)+3;
	m44=m43+1;
	w1=nrvector(1,m4);
	w2=nrvector(1,m);
	facm=m;
	facp=1.0/m;
	for (j=1;j<=mm;j++) sumw += SQR(WINDOW(j,facm,facp));
	for (j=1;j<=m;j++) p[j]=0.0;
	for (j=1;j<=m;j++) w2[j] = mdata[c++];
	for (kk=1;kk<=k;kk++) {
		for (joff = -1;joff<=0;joff++) {
			for (j=1;j<=m;j++) w1[joff+j+j]=w2[j];
			for (j=1;j<=m;j++) w2[j] = mdata[c++];
			joffn=joff+mm;
			for (j=1;j<=m;j++) w1[joffn+j+j]=w2[j];
		}
		for (j=1;j<=mm;j++) {
			j2=j+j;
			w=WINDOW(j,facm,facp);
			w1[j2] *= w;
			w1[j2-1] *= w;
		}
		nrfour1(w1,mm,1);
		p[1] += (SQR(w1[1])+SQR(w1[2]));
		for (j=2;j<=m;j++) {
			j2=j+j;
			p[j] += (SQR(w1[j2])+SQR(w1[j2-1])
				+SQR(w1[m44-j2])+SQR(w1[m43-j2]));
		}
		den += sumw;
	}
	den *= m4;
	for (j=1;j<=m;j++) p[j] /= den;
	nrfree_vector(w2,1,m);
	nrfree_vector(w1,1,m4);
}

double mymax(double x, double y){
  return x > y ? x : y;
}

/* ********************************************************************* */

void mysvd(int m, int n, double** u, double w[], double** v, int* ierr)
/*
 *   This subroutine is a translation of the Algol procedure svd,
 *   Num. Math. 14, 403-420(1970) by Golub and Reinsch.
 *   Handbook for Auto. Comp., Vol II-Linear Algebra, 134-151(1971).
 *
 *   This subroutine determines the singular value decomposition
 *        t
 *   A=usv  of a real m by n rectangular matrix, where m is greater
 *   then or equal to n.  Householder bidiagonalization and a variant
 *   of the QR algorithm are used.
 *  
 *
 *   On input.
 *
 *      m is the number of rows of A (and u).
 *
 *      n is the number of columns of A (and u) and the order of v.
 *
 *      u contains the rectangular input matrix A to be decomposed.
 *
 *   On output.
 *
 *      w contains the n (non-negative) singular values of a (the
 *        diagonal elements of s).  they are unordered.  if an
 *        error exit is made, the singular values should be correct
 *        for indices ierr+1,ierr+2,...,n.
 *
 *      u contains the matrix u (orthogonal column vectors) of the
 *        decomposition.
 *        if an error exit is made, the columns of u corresponding
 *        to indices of correct singular values should be correct.
 *
 *      v contains the matrix v (orthogonal) of the decomposition.
 *        if an error exit is made, the columns of v corresponding
 *        to indices of correct singular values should be correct.
 *
 *      ierr is set to
 *        zero       for normal return,
 *        k          if the k-th singular value has not been
 *                   determined after 30 iterations,
 *        -1         if memory allocation fails.
 *
 *   Questions and comments should be directed to B. S. Garbow,
 *   Applied Mathematics division, Argonne National Laboratory
 *
 *   Modified to eliminate machep
 *
 *   Translated to C by Michiel de Hoon, Human Genome Center,
 *   University of Tokyo, for inclusion in the C Clustering Library.
 *   This routine is less general than the original svd routine, as
 *   it focuses on the singular value decomposition as needed for
 *   clustering. In particular,
 *     - We require m >= n
 *     - We calculate both u and v in all cases
 *     - We pass the input array A via u; this array is subsequently
 *       overwritten.
 *     - We allocate for the array rv1, used as a working space,
 *       internally in this routine, instead of passing it as an
 *       argument. If the allocation fails, svd sets *ierr to -1
 *       and returns.
 *   2003.06.05
 */
{ int i, j, k, i1, k1, l1, its;
  double c,f,h,s,x,y,z;
  int l = 0;
  double g = 0.0;
  double scale = 0.0;
  double anorm = 0.0;
  double* rv1 = (double*)malloc(n*sizeof(double));
  if (!rv1)
  { *ierr = -1;
    return;
  }
  *ierr = 0;
  /* Householder reduction to bidiagonal form */
  for (i = 0; i < n; i++)
  { l = i + 1;
    rv1[i] = scale * g;
    g = 0.0;
    s = 0.0;
    scale = 0.0;
    for (k = i; k < m; k++) scale += fabs(u[k][i]);
    if (scale != 0.0)
    { for (k = i; k < m; k++)
      { u[k][i] /= scale;
        s += u[k][i]*u[k][i];
      }
      f = u[i][i];
      g = (f >= 0) ? -sqrt(s) : sqrt(s);
      h = f * g - s;
      u[i][i] = f - g;
      if (i < n-1)
      { for (j = l; j < n; j++)
        { s = 0.0;
          for (k = i; k < m; k++) s += u[k][i] * u[k][j];
          f = s / h;
          for (k = i; k < m; k++) u[k][j] += f * u[k][i];
        }
      }
      for (k = i; k < m; k++) u[k][i] *= scale;
    }
    w[i] = scale * g;
    g = 0.0;
    s = 0.0;
    scale = 0.0;
    if (i<n-1)
    { for (k = l; k < n; k++) scale += fabs(u[i][k]);
      if (scale != 0.0)
      { for (k = l; k < n; k++)
        { u[i][k] /= scale;
          s += u[i][k] * u[i][k];
        }
        f = u[i][l];
        g = (f >= 0) ? -sqrt(s) : sqrt(s);
        h = f * g - s;
        u[i][l] = f - g;
        for (k = l; k < n; k++) rv1[k] = u[i][k] / h;
        for (j = l; j < m; j++)
        { s = 0.0;
          for (k = l; k < n; k++) s += u[j][k] * u[i][k];
          for (k = l; k < n; k++) u[j][k] += s * rv1[k];
        }
        for (k = l; k < n; k++)  u[i][k] *= scale;
      }
    }
    anorm = mymax(anorm,fabs(w[i])+fabs(rv1[i]));
  }
  /* accumulation of right-hand transformations */
  for (i = n-1; i>=0; i--)
  { if (i < n-1)
    { if (g != 0.0)
      { for (j = l; j < n; j++) v[j][i] = (u[i][j] / u[i][l]) / g;
        /* double division avoids possible underflow */
        for (j = l; j < n; j++)
        { s = 0.0;
          for (k = l; k < n; k++) s += u[i][k] * v[k][j];
          for (k = l; k < n; k++) v[k][j] += s * v[k][i];
        }
      }
    }
    for (j = l; j < n; j++)
    { v[i][j] = 0.0;
      v[j][i] = 0.0;
    }
    v[i][i] = 1.0;
    g = rv1[i];
    l = i;
  }
  /* accumulation of left-hand transformations */
  for (i = n-1; i >= 0; i--)
  { l = i + 1;
    g = w[i];
    if (i!=n-1)
      for (j = l; j < n; j++) u[i][j] = 0.0;
    if (g!=0.0)
    { if (i!=n-1)
      { for (j = l; j < n; j++)
        { s = 0.0;
          for (k = l; k < m; k++) s += u[k][i] * u[k][j];
          /* double division avoids possible underflow */
          f = (s / u[i][i]) / g;
          for (k = i; k < m; k++) u[k][j] += f * u[k][i];
        }
      }
      for (j = i; j < m; j++) u[j][i] /= g;
    }
    else
      for (j = i; j < m; j++) u[j][i] = 0.0;
    u[i][i] += 1.0;
  }
  /* diagonalization of the bidiagonal form */
  for (k = n-1; k >= 0; k--)
  { k1 = k-1;
    its = 0;
    while(1)
    /* test for splitting */
    { for (l = k; l >= 0; l--)
      { l1 = l-1;
        if (fabs(rv1[l]) + anorm == anorm) break;
        /* rv1[0] is always zero, so there is no exit
         * through the bottom of the loop */
        if (fabs(w[l1]) + anorm == anorm)
        /* cancellation of rv1[l] if l greater than 0 */
        { c = 0.0;
          s = 1.0;
          for (i = l; i <= k; i++)
          { f = s * rv1[i];
            rv1[i] *= c;
            if (fabs(f) + anorm == anorm) break;
            g = w[i];
            h = sqrt(f*f+g*g);
            w[i] = h;
            c = g / h;
            s = -f / h;
            for (j = 0; j < m; j++)
            { y = u[j][l1];
              z = u[j][i];
              u[j][l1] = y * c + z * s;
              u[j][i] = -y * s + z * c;
            }
          }
          break;
        }
      }
      /* test for convergence */
      z = w[k];
      if (l==k) /* convergence */
      { if (z < 0.0)
        /*  w[k] is made non-negative */
        { w[k] = -z;
          for (j = 0; j < n; j++) v[j][k] = -v[j][k];
        }
        break;
      }
      else if (its==30)
      { *ierr = k;
        break;
      }
      else
      /* shift from bottom 2 by 2 minor */
      { its++;
        x = w[l];
        y = w[k1];
        g = rv1[k1];
        h = rv1[k];
        f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
        g = sqrt(f*f+1.0);
        f = ((x - z) * (x + z) + h * (y / (f + (f >= 0 ? g : -g)) - h)) / x;
        /* next qr transformation */
        c = 1.0;
        s = 1.0;
        for (i1 = l; i1 <= k1; i1++)
        { i = i1 + 1;
          g = rv1[i];
          y = w[i];
          h = s * g;
          g = c * g;
          z = sqrt(f*f+h*h);
          rv1[i1] = z;
          c = f / z;
          s = h / z;
          f = x * c + g * s;
          g = -x * s + g * c;
          h = y * s;
          y = y * c;
          for (j = 0; j < n; j++)
          { x = v[j][i1];
            z = v[j][i];
            v[j][i1] = x * c + z * s;
            v[j][i] = -x * s + z * c;
          }
          z = sqrt(f*f+h*h);
          w[i1] = z;
          /* rotation can be arbitrary if z is zero */
          if (z!=0.0)
          { c = f / z;
            s = h / z;
          }
          f = c * g + s * y;
          x = -s * g + c * y;
          for (j = 0; j < m; j++)
          { y = u[j][i1];
            z = u[j][i];
            u[j][i1] = y * c + z * s;
            u[j][i] = -y * s + z * c;
          }
        }
        rv1[l] = 0.0;
        rv1[k] = f;
        w[k] = x;
      }
    }
  }
  free(rv1);
  return;
}

int mycompare( const void* a, const void* b ) {
  double* arg1 = (double*) a;
  double* arg2 = (double*) b;
  if( *arg1 < *arg2 ) return -1;
  else if( *arg1 == *arg2 ) return 0;
  else return 1;
}  

double* myspectrum(double* v1,int dc,int* anslen){
  // n data pts will be divided into k partitions of size m
  // the spectrum will have m values from 0 to m/2 cycles/dt.
  //  n = (2*k+1)*m

  // data set
  //  Vect* v1 = vector_arg(1);

  //  int dc = v1->capacity();
  int mr;
  //  if (ifarg(2)) mr = (int)(*getarg(2)); else mr = dc/8;
  mr = dc / 8;

  // make sure the length of partitions is integer power of 2
  int m = 1;
  while(m<mr) m*=2;
  *anslen = m;

  int k = (int) (ceil(((double)(dc)/m-1.)/2.));
  int n = (2*k+1)*m;

  int i;

  double *mdata = (double *)calloc(n,(unsigned)sizeof(double));
  for (i=0;i<dc;++i) mdata[i] = v1[i]; //v1->elem(i);

  //if (ans->capacity() < m) ans->resize(m);
  double* ans = (double *)calloc(m,(unsigned)sizeof(double));
  //spctrm(&mdata[0], &ans->elem(0)-1, m, k);
  nrspctrm(&mdata[0], &ans[0], m, k);

  free((char *)mdata);

  return ans;
}

double myspct(void* v){
  double* x;
  int n = vector_instance_px(v,&x) , anslen = 0 , i = 0;
  double* ans = myspectrum(x,n,&anslen);
  for(i=0;i<anslen;i++) printf("%g ",ans[i]);
  printf("\n");
  double* p;
  int outlen = vector_arg_px(1,&p);
  for(i=0;i<outlen && i<anslen;i++) p[i]=ans[i];
  free(ans);
  return 1.0;
}

int iinrange(double* x,int sz,double dmin,double dmax){
  int i = 0, cnt = 0;
  for(i=0;i<sz;i++) if(x[i]>=dmin && x[i]<=dmax) cnt++;
  return cnt;
}

double inrange(void* v){
  double* x , dmin, dmax;
  int sz = vector_instance_px(v,&x) , i , cnt = 0;
  return (double) iinrange(x,sz,*getarg(1),*getarg(2));
}

extern double dfftpow(double* x,int n,double* ppow,int powlen,int* fftlen);

ENDVERBATIM

: get averagecorrelations between channels, but first bandpass filter them
: tscoravgband(output vector,vec-list of eeg time series,window size,increment,vec:60/120Hz power,vec:saturation,
: vec:delta power,vec:theta,vec:beta,vec:gamma,vec:ripple,sampling rate,vec:low-cutoffs,vec:high-cutoffs, [optional
: output list of vectors for filtered time-series]
:  --
: low and high cutoffs are different bands that will be kept, so can do theta+gamma, i.e. lo=4,31 , hi=12,100
: or can just do theta : lo=4 , hi=12
FUNCTION tscoravgband () {
  VERBATIM

  int i;

  double dRet = 0.0;
  ListVec* pTS = 0;

  double* pcvin = 0;
  double** pfilter = 0x0;

  double** pTSN = 0; // normalized time series

  double* vTCor = 0; //correlation matrix abs sum avg for a time point
  int corsz = vector_arg_px(1,&vTCor);

  double* pTMP = 0x0;

  ListVec* pTSOut = 0;

  pTS = AllocListVec(*hoc_objgetarg(2)); //list of eeg time series
  
  if(!pTS){
    printf("tscoravg ERRA: problem initializing 1st 2 args!\n");
    goto CLEANUP;
  }
  int iChans = pTS->isz; // # of channels
  if(iChans < 2) {
    printf("tscoravg ERRB: must have at least 2 EEG channels!\n");
    goto CLEANUP;
  }
  int iWinSz = (int) *getarg(3); //window size for analysis
  int iSamples = pTS->plen[0]; //make sure all time series have same length
  int iINC = ifarg(4)?(int)*getarg(4):1; 
  double *phz = 0, *psat = 0, *pfft = 0; int sz = 0 , iNoiseTest = 0;
  if(ifarg(5) && ifarg(6)){ //do saturation and 60/120 Hz Frequency test
    if((sz=vector_arg_px(5,&phz))!=corsz){
      printf("tscoravg ERRG: invalid size for phz should be %d, is %d!\n",corsz,sz);
      goto CLEANUP;
    } else if((sz=vector_arg_px(6,&psat))!=corsz){
      printf("tscoravg ERRH: invalid size for psat should be %d, is %d!\n",corsz,sz);
      goto CLEANUP;
    } else {	      
      iNoiseTest = 1;
      pfft = (double*)calloc(iWinSz+1,sizeof(double));
    }
  }
  int iOutSz = (iSamples - iWinSz + 1) / iINC;//output List size

  double *pdelta=0,*ptheta=0,*pbeta=0,*pgamma=0,*pripple=0;
  if(ifarg(7) && vector_arg_px(7,&pdelta)<iOutSz){ printf("bad pdelta size\n"); goto CLEANUP; } //delta 0-3
  if(ifarg(8) && vector_arg_px(8,&ptheta)<iOutSz){ printf("bad ptheta size\n"); goto CLEANUP; } //theta 4-12
  if(ifarg(9)&& vector_arg_px(9,&pbeta)<iOutSz){ printf("bad pbeta size\n"); goto CLEANUP; } //beta 13-29
  if(ifarg(10)&& vector_arg_px(10,&pgamma)<iOutSz){ printf("bad pgamma size\n"); goto CLEANUP; } //gamma 30-100
  if(ifarg(11)&& vector_arg_px(11,&pripple)<iOutSz){ printf("bad pripple size\n"); goto CLEANUP; } //ripple 101-300

  double sampr = *getarg(12);

  double *plohz=0x0,*phihz=0x0;
  int iFilters = 1;
  if((iFilters=vector_arg_px(13,&plohz))!=(i=vector_arg_px(14,&phihz))) {
    printf("invalid arg 13,14 lohz/hihz filter sizes: %d %d\n",iFilters,i); goto CLEANUP;
  } else if(iFilters<1) { printf("no filters specified!\n"); goto CLEANUP; }
  int N = 1, M = 1025;
  while(N < iWinSz) N*=2;
  if(N>M) M=N+1;
  pfilter = getdouble2D(iFilters,N);
  //normalize low/high cutoff frequencies by sampling rate (should be between 0 and 0.5)
  for(i=0;i<iFilters;i++){ plohz[i]/=sampr; phihz[i]/=sampr; }
  for(i=0;i<iFilters;i++){ wsfirBP(pfilter[i], M, 1, plohz[i], phihz[i], N); wrap(pfilter[i],N,M); }
  pcvin = (double*) calloc(N,sizeof(double));

  if(verbose) printf("iWinSz=%d iINC=%d iSamples=%d\n",iWinSz,iINC,iSamples);
  for(i=1;i<iChans;i++) {
    if(pTS->plen[i] != iSamples){
      printf("tscoravg ERRC: time series of unequal size %d %d!\n",iSamples,pTS->plen[i]);
      goto CLEANUP;
    }
  }

  if(corsz < iOutSz){
    printf("tscoravg ERRD: need output vector of at least %d as arg 1, have only %d!\n",iOutSz,corsz);
    goto CLEANUP;
  }
  if(verbose) printf("iOutSz=%d\n",iOutSz);

  pTSN = getdouble2D(iChans,2*N); if(verbose) printf("got pTSN\n");
  if(!pTSN){
    printf("tscor ERRD: out of memory!\n");
    goto CLEANUP;
  }

  if(ifarg(15) && (pTSOut = AllocListVec(*hoc_objgetarg(15)))) { //optional - copy filtered output
    if(pTSOut->isz!=iChans) { printf("incorrect output vec-list size: %d %d\n",pTSOut->isz,iChans); goto CLEANUP; }
    for(i=0;i<iChans;i++) if(pTSOut->plen[i] < iSamples) { printf("incorrect output vec-size %d %d\n",pTSOut->plen[i],iSamples);
      goto CLEANUP; }}

  pTMP = (double*) calloc(2*N,sizeof(double));

  double* psums = (double*)calloc(iChans,sizeof(double)); //stores sum of values in channel's time window
  double* psums2 = (double*)calloc(iChans,sizeof(double));//stores sum of squared values in channel's time window

  double dVal = 0.0;
  int iStartIDX = 0 , ierr = 0 , iOutIDX = 0 , iFilt = 0;
  for(iOutIDX=0;iOutIDX<corsz;iOutIDX++) vTCor[iOutIDX] = 0.;
  if( phz ) for(iOutIDX=0;iOutIDX<corsz;iOutIDX++) phz[iOutIDX]=0.;
  if( psat ) for(iOutIDX=0;iOutIDX<corsz;iOutIDX++) psat[iOutIDX]=0.;
  iOutIDX=0;
  double CP = iChans*(iChans-1.0)/2.0;
  for(iStartIDX=0;iOutIDX<iOutSz && iStartIDX+iWinSz<iSamples;iStartIDX+=iINC,iOutIDX++){
    if(iOutIDX%1000==0) printf("%d\n",iOutIDX);
    int IDX = iStartIDX , iChan = 0 , iEndIDX = iStartIDX + iWinSz;
    int iCurSz = iEndIDX - iStartIDX;
    if(phz) phz[iOutIDX]=0.; if(psat)psat[iOutIDX]=0.;
    if(pdelta) pdelta[iOutIDX]=0.; if(ptheta) ptheta[iOutIDX]=0.;
    if(pbeta) pbeta[iOutIDX]=0.; if(pgamma) pgamma[iOutIDX]=0.; if(pripple) pripple[iOutIDX]=0.;

    double dsatavg = 0., dhzavg = 0.;
    for(iChan=0;iChan<iChans;iChan++){       
      i=0;      
      int iSat = 0; 
      double* pc = &pTS->pv[iChan][iStartIDX]; //saturation test
      for(IDX=iStartIDX;IDX<iEndIDX;IDX++,pc++)if((*pc>=900. && *pc<=1010.) || (*pc>=-1010. && *pc<=-900.)) iSat++;
      dsatavg += (double)iSat/ (double)iWinSz;
      int fftlen = 0; memset(pfft,0,sizeof(double)*(iWinSz+1));
      if(dfftpow(&pTS->pv[iChan][iStartIDX],iWinSz,pfft,iWinSz+1,&fftlen)){       //60,120 Hz Frequency test
	int f = 0; pfft[0]=0.0; double fsum = 0.0; double* ppf = pfft;
	for(f=0;f<fftlen && f<=sampr/2;f++) fsum += *ppf++; //normalize power btwn 0-sampr/2Hz to 1.0
	for(f=55;f<=65 && f<fftlen;f++) dhzavg += pfft[f] / fsum; //60Hz contribution
	for(f=115;f<=125 && f<fftlen;f++) dhzavg += pfft[f] / fsum; //120Hz contribution
	if(pdelta) {
	  fsum=0.0;
	  for(f=0;f<fftlen;f++) fsum+=pfft[f];
	  for(f=0;f<=3;f++) pdelta[iOutIDX] += pfft[f]/fsum;
	  if(ptheta) for(f=4;f<=12;f++) ptheta[iOutIDX] += pfft[f]/fsum;
	  if(pbeta) for(f=13;f<=29;f++) pbeta[iOutIDX] += pfft[f]/fsum;
	  if(pgamma) for(f=30;f<=100;f++) pgamma[iOutIDX] += pfft[f]/fsum;
	  if(pripple) for(f=101;f<=300 && f<fftlen;f++) pripple[iOutIDX] += pfft[f]/fsum;
	}
      }	
    }
    phz[iOutIDX] = dhzavg / (double)iChans;
    psat[iOutIDX] = dsatavg / (double)iChans;
    if(pdelta) pdelta[iOutIDX] /= (double)iChans;
    if(ptheta) ptheta[iOutIDX] /= (double)iChans;
    if(pbeta) pbeta[iOutIDX] /= (double)iChans;
    if(pgamma) pgamma[iOutIDX] /= (double)iChans;
    if(pripple) pripple[iOutIDX] /= (double)iChans;

    //bandpass filter and normalize
    for(iChan=0;iChan<iChans;iChan++){
      memset(pTMP,0,sizeof(double)*2*N);
      double* pC = &pTS->pv[iChan][iStartIDX], *pIN = pcvin;
      if(0) for(IDX=iStartIDX;IDX<iEndIDX;IDX++) *pIN++ = *pC++; //copy input
      for(iFilt=0;iFilt<iFilters;iFilt++) {        
        if(1) {
          memset(pTSN[iChan],0,sizeof(double)*N);//zero it out first...
          pIN = pcvin; pC = &pTS->pv[iChan][iStartIDX];
          for(IDX=iStartIDX;IDX<iEndIDX;IDX++) *pIN++ = *pC++; //copy input
        }
        convlv(pcvin-1,N,pfilter[iFilt]-1,M,1,pTSN[iChan]-1);//do filtering
        for(IDX=0;IDX<N;IDX++) pTMP[IDX] += pTSN[iChan][IDX]; 
      }
      for(IDX=0;IDX<N;IDX++) pTSN[iChan][IDX] = pTMP[IDX]; //copy filtered output to pTSN
      if(pTSOut)for(IDX=0;IDX<N && IDX+iStartIDX<pTSOut->plen[iChan];IDX++) pTSOut->pv[iChan][iStartIDX+IDX]=pTSN[iChan][IDX];

      psums[iChan]=psums2[iChan]=0.;
      pC = pTSN[iChan]; //pointer to start of channel
      for(i=0;i<N;i++,pC++){ 
	psums[iChan] += *pC;//for average
	psums2[iChan] += (*pC * *pC); //for std-dev
      }
      double davg = psums[iChan] / N;
      double dstdev = psums2[iChan]/N - davg*davg; 
      if(dstdev > 0. ) dstdev = sqrt(dstdev); else dstdev=1.; //make sure no problems with nan
      if(dstdev <= 0.) dstdev = 1.0;
      pC = pTSN[iChan]; //pointer to start of channel
      for(i=0;i<N;i++,pC++) *pC = (*pC - davg)/dstdev; //normalization
    }
    //get average of pairwise correlations across channels
    int iC1,iC2;
    double dpsum = 0.0;
    for(iC1=0;iC1<iChans;iC1++){
      for(iC2=0;iC2<iC1;iC2++){
        double r = 0.0, *p1 = pTSN[iC1], *p2 = pTSN[iC2];
	for(i=0;i<N;i++) r += (*p1++ * *p2++);
        r /= (double) N;
        dpsum += r;
      }
    }
    dpsum /= CP;
    if(dpsum > 1.0 || dpsum < -1.0) printf("out of bounds = %g!\n",dpsum);
    vTCor[iOutIDX] = dpsum;
  }
  dRet = 1.0;
  printf("\n");

CLEANUP:
  if(pTS) FreeListVec(&pTS);                       if(pTSN) freedouble2D(&pTSN,iChans);
  free(psums);                                     free(psums2); 
  if(pfft) free(pfft);
  if(pfilter) freedouble2D(&pfilter,iFilters);     if(pcvin) free(pcvin);
  if(pTMP) free(pTMP);                             if(pTSOut) FreeListVec(&pTSOut);
  return dRet;
  ENDVERBATIM
}


: normalize values in list of vectors
FUNCTION tsnorm () {
  VERBATIM
  double dRet = 0.0;
  ListVec* pEIG = AllocListVec(*hoc_objgetarg(1));
  if(!pEIG || pEIG->isz < 1 || pEIG->plen[0]<1){
    printf("tsnorm ERRA: problem initializing 1st arg!\n");
    goto TSNCLEANUP;
  }
  double* dmean = 0, *dstdev = 0;
  int cols = vector_arg_px(2,&dmean);
  if(cols!=vector_arg_px(3,&dstdev) || cols!=pEIG->plen[0]){
    printf("tsnorm ERRB: problem initializing 2nd,3rd arg!\n");
    goto TSNCLEANUP;
  }
  int i,j;
  double val;
  for(i=0;i<pEIG->isz;i++)for(j=0;j<cols;j++)pEIG->pv[i][j]=(pEIG->pv[i][j]-dmean[j])/dstdev[j];
  dRet = 1.0;
TSNCLEANUP:
  if(pEIG) FreeListVec(&pEIG);
  return dRet;
  ENDVERBATIM
}

VERBATIM
int ddiff(double* pin,double* pout,int sz){
  int i;
  for(i=0;i<sz-1;i++) pout[i]=pin[i+1]-pin[i];
  return 1;
}
ENDVERBATIM

: get averages from tscorfull into Vectors
: tsgetmeans( 1 cor vec list, 2 phzlist, 3 psatlist, 4 pderivlist, 5 avg cor vec A, 6 avg cor vec B, 7 avg cor vec C, 8 phzt, 9 psatt, 10 pderivt , 11 numchannels , 12 phza factor, 13 vfid, 14 fid, 15 vprctbadch)
: arg 5 will have average correlation without getting rid of any channels on any epochs
: arg 6 will have average correlation with getting rid of saturated and flat channels
: arg 7 will have avg. cor. with getting rid of saturated,flat channels and channels with 60/120Hz power > phzt
: pfctr not used currently
: 13 vector of file ids so can selectively clean up particular files
: 14 file id to clean up -- iff == -1, do all
: 15 % of 'bad' channels for each time period -- channels excluded for "avgC" / # of channels
FUNCTION tsgetmeans () {
VERBATIM
  int iC1,iC2,iChans,idx,jdx;
  double dRet = 0.0;
  ListVec *plvcor = 0, *plvhz = 0, *plvsat = 0, *plvderiv = 0;
  double* pavgA=0,*pavgB=0,*pavgC=0,*pprctbadch=0,*pfid=0,*pbadc=0;
  int vsz = 0, fid=-1; 
  ////////////////////////////////////////////////////////////////////////////////////////////////
  //
  //                                     set up input args.
  //
  if(!(plvcor = AllocListVec(*hoc_objgetarg(1)))){
    printf("tsgetmeans ERRA: couldn't get cor vec list!\n");
    goto MCLEANUP;
  }
  if(!(plvhz = AllocListVec(*hoc_objgetarg(2)))){
    printf("tsgetmeans ERRB: couldn't get phz vec list!\n");
    goto MCLEANUP;
  }
  if(plvhz->isz < plvcor->isz){
    printf("tsgetmeans ERRC: phz vec list size < cor vec list size : %d %d\n",plvhz->isz,plvcor->isz);
    goto MCLEANUP;
  }
  if(!(plvsat = AllocListVec(*hoc_objgetarg(3)))){
    printf("tsgetmeans ERRD: couldn't get psat vec list!\n");
    goto MCLEANUP;
  }
  if(plvsat->isz < plvcor->isz){
    printf("tsgetmeans ERRE: psat vec list size < cor vec list size : %d %d\n",plvsat->isz,plvcor->isz);
    goto MCLEANUP;
  }  
  if(!(plvderiv = AllocListVec(*hoc_objgetarg(4)))){
    printf("tsgetmeans ERRF: couldn't get pderiv vec list!\n");
    goto MCLEANUP;
  }
  if(plvderiv->isz < plvcor->isz){
    printf("tsgetmeans ERRG: pderiv vec list size < cor vec list size : %d %d\n",plvderiv->isz,plvcor->isz);
    goto MCLEANUP;
  }  
  if(vector_arg_px(5,&pavgA) < plvcor->isz || 
     vector_arg_px(6,&pavgB) < plvcor->isz || 
     vector_arg_px(7,&pavgC) < plvcor->isz ||
     vector_arg_px(13,&pfid) < plvcor->isz ||
     vector_arg_px(15,&pprctbadch)< plvcor->isz){
    printf("tsgetmeans ERRH: output vecs must have size >= %d!\n",plvcor->isz);
    goto MCLEANUP;
  }
  double phzt = *getarg(8), psatt = *getarg(9) , pderivt = *getarg(10);
  iChans = (int)*getarg(11);
  double pfctr = *getarg(12);
  fid = (int)*getarg(14);
  pbadc = (double*) malloc(iChans*sizeof(double));
  //
  //
  ////////////////////////////////////////////////////////////////////////////////////////////////

  /////////
  // do the calculations  
  double dSumA,dSumB,dSumC;
  int cntA,cntB,cntC,cntBadCh;
  for(idx=0;idx<plvcor->isz;idx++){
    if(fid>=0 && pfid[idx]!=fid) continue;
    cntA=cntB=cntC=jdx=0;
    dSumA=dSumB=dSumC=0.;
    for(iC1=0;iC1<iChans;iC1++){ //mark bad channels
      if(plvsat->pv[idx][iC1] >= psatt ||
         plvderiv->pv[idx][iC1] >= pderivt)
        pbadc[iC1]=1;
      else
        pbadc[iC1]=0;
    }
    for(iC1=0;iC1<iChans;iC1++){
      if(pbadc[iC1])continue;
      for(iC2=0;iC2<iC1;iC2++,jdx++){
        dSumA += plvcor->pv[idx][jdx];
        cntA++;
        if(pbadc[iC2]) continue;
        dSumB += plvcor->pv[idx][jdx];
        cntB++;
      }
    }
    pavgA[idx] = dSumA / (double) cntA;
    if(cntB) pavgB[idx] = dSumB / (double) cntB; else pavgB[idx] = -666.;
    jdx=0;
    for(iC1=0;iC1<iChans;iC1++) if(plvhz->pv[idx][iC1] >= phzt) pbadc[iC1]=1; //mark bad channels
    for(iC1=0;iC1<iChans;iC1++){
      if(pbadc[iC1]) continue;
      for(iC2=0;iC2<iC1;iC2++,jdx++){
        if(pbadc[iC2]) continue;
        dSumC += plvcor->pv[idx][jdx];
        cntC++;
      }
    }
    if(cntC) pavgC[idx] = dSumC / (double) cntC; else pavgC[idx] = -666.;
    cntBadCh=0;
    for(iC1=0;iC1<iChans;iC1++) if(pbadc[iC1]) cntBadCh++;
    pprctbadch[idx]= (double)cntBadCh/(double)iChans;
  }

  dRet = 1.0;
MCLEANUP:
 if(pbadc) free(pbadc);
 return dRet;
ENDVERBATIM
}

: get full time series correlation matrix for each time window into list of vectors
: tscorfull( 1 cor vec list, 2 time series list, 3 winsz, 4 incsz, 5 phzlist, 6 psatlist, 7 pderivlist)
FUNCTION tscorfull () {
  VERBATIM

  int i;

  double dRet = 0.0;
  ListVec* pTS = 0, *plvcor = 0, *plvhz = 0, *plvsat = 0, *plvderiv = 0;

  double** pTSN = 0; // normalized time series
  double* pdiff = 0; // for straight line test
  double* psums = 0, *psums2 = 0, *pfft = 0;

  int corsz=0; //# of vectors in output correlation vector list
  plvcor = AllocListVec(*hoc_objgetarg(1)); //each vec is correlation mat in vector form for a time point
  if(!plvcor || (corsz=plvcor->isz)<1){
    printf("tscorfull ERRA: problem initializing arg1 plvcor!\n");
    goto FCLEANUP;
  }

  pTS = AllocListVec(*hoc_objgetarg(2)); //list of eeg time series  
  if(!pTS){
    printf("tscorfull ERRA: problem initializing arg2 pTS!\n");
    goto FCLEANUP;
  }

  int iChans = pTS->isz; // # of channels
  if(iChans < 2) {
    printf("tscorfull ERRB: must have at least 2 timeseries!\n");
    goto FCLEANUP;
  }

  int iMatSz = iChans*(iChans-1)/2;

  for(i=0;i<corsz;i++){   //make sure all output vectors have right size
    if(plvcor->plen[i]<iMatSz){
      printf("tscorfull ERRC: each cor vec must have sz >= %d, [%d] has only %d\n",iMatSz,i,plvcor->plen[i]);
      goto FCLEANUP;
    }
  }

  int iWinSz = (int) *getarg(3); //window size for analysis
  int iINC = ifarg(4)?(int)*getarg(4):1;  //increment size for analysis

  plvhz = AllocListVec(*hoc_objgetarg(5));
  plvsat = AllocListVec(*hoc_objgetarg(6));
  plvderiv = AllocListVec(*hoc_objgetarg(7));

  pfft = (double*)calloc(iWinSz+1,sizeof(double));
  pdiff = (double*)calloc(iWinSz,sizeof(double));

  int iSamples = pTS->plen[0]; //make sure all time series have same length
  int iOutSz = (iSamples - iWinSz + 1) / iINC;//output List size

  if(verbose) printf("iWinSz=%d iINC=%d iSamples=%d\n",iWinSz,iINC,iSamples);
  for(i=1;i<iChans;i++) {
    if(pTS->plen[i] != iSamples){
      printf("tscorfull ERRD: time series of unequal size %d %d!\n",iSamples,pTS->plen[i]);
      goto FCLEANUP;
    }
  }

  if(corsz < iOutSz){
    printf("tscorfull ERRE: need output list vector of at least sz %d as arg 1, have only %d!\n",iOutSz,corsz);
    goto FCLEANUP;
  }   if(verbose) printf("iOutSz=%d\n",iOutSz);
  if(plvhz->isz<iOutSz || plvsat->isz<iOutSz || plvderiv->isz<iOutSz){
    printf("tscorfull ERRF: need output plvhz,plvsat,plvderiv size of at least %d!\n",iOutSz);
    goto FCLEANUP;
  }
  for(i=0;i<iOutSz;i++){
    if(plvhz->plen[i]<iChans || plvsat->plen[i]<iChans || plvderiv->plen[i]<iChans){
      printf("tscorfull ERRG: each vec in plvhz,plvsat,plvderiv must be sz >= %d!\n",iChans);
      goto FCLEANUP;
    }
  }

  pTSN = getdouble2D(iChans,iWinSz); if(verbose) printf("got pTSN\n"); //normalized time series
  if(!pTSN){
    printf("tscor ERRD: out of memory!\n");
    goto FCLEANUP;
  }

  psums = (double*)calloc(iChans,sizeof(double));   //for sum(X)
  psums2 = (double*)calloc(iChans,sizeof(double));  //for sum(X^2)
  double* pc = 0;
  double dVal = 0.0;
  int iStartIDX = 0 , ierr = 0 , iOutIDX = 0;
  FillListVec(plvcor,0.0); FillListVec(plvhz,0.); FillListVec(plvsat,0.); FillListVec(plvderiv,0.);//zero results
  for(iStartIDX=0;iOutIDX<iOutSz && iStartIDX+iWinSz<iSamples;iStartIDX+=iINC,iOutIDX++){
    if(iOutIDX%1000==0) printf("%d\n",iOutIDX);
    int IDX = iStartIDX , iChan = 0 , iEndIDX = iStartIDX + iWinSz;
    int iCurSz = iEndIDX - iStartIDX;
    for(iChan=0;iChan<iChans;iChan++){ //get sum(X) , sum(X^2) for ALL channels
      psums[iChan]=psums2[iChan]=0.;
      pc = &pTS->pv[iChan][iStartIDX];
      for(IDX=iStartIDX;IDX<iEndIDX;IDX++,pc++){
        psums[iChan] += *pc; 
        psums2[iChan] += (*pc * *pc);
      }
    }
    for(iChan=0;iChan<iChans;iChan++){       //do noise checks and normalize for ALL channels
      i=0;
      double davg = psums[iChan] / iCurSz;
      double dstdev = psums2[iChan]/iCurSz - davg*davg; 
      if(dstdev > 0. ) dstdev = sqrt(dstdev); else dstdev=1.; //make sure no problems with nan
      if(dstdev <= 0.) dstdev = 1.0;
      //saturation test
      int iSat = 0;
      pc = &pTS->pv[iChan][iStartIDX];
      for(IDX=iStartIDX;IDX<iEndIDX;IDX++,i++,pc++){	
        if((*pc>=900. && *pc<=1010.) || (*pc>=-1010. && *pc<=-900.)) iSat++;
        pTSN[iChan][i] = (*pc - davg) / dstdev;	//normalization
      }
      plvsat->pv[iOutIDX][iChan] = (double)iSat/ (double)iWinSz;
      //60,120 Hz Frequency test
      int fftlen = 0; memset(pfft,0,sizeof(double)*(iWinSz+1));
      if(dfftpow(&pTS->pv[iChan][iStartIDX],iWinSz,pfft,iWinSz+1,&fftlen)){
        int f = 0; pfft[0]=0.0; double fsum = 0.0;
        for(f=0;f<fftlen && f<=1000;f++) fsum += pfft[f]; //normalize power btwn 0-1000Hz to 1.0
        pc = &plvhz->pv[iOutIDX][iChan];
        for(f=55;f<=65 && f<fftlen;f++)   *pc += pfft[f]; //60Hz contribution
        for(f=115;f<=125 && f<fftlen;f++) *pc += pfft[f]; //120Hz contribution
        *pc /= fsum;
      }
      //discrete derivative test to measure how flat of a line signal looks like
      ddiff(&pTS->pv[iChan][iStartIDX],pdiff,iWinSz);
      plvderiv->pv[iOutIDX][iChan] = (double)iinrange(pdiff, iWinSz-1, -1.0, 1.0) / (double)(iWinSz-1.0);
    }

    //form correlation vector
    int iC1,iC2; pc = &plvcor->pv[iOutIDX][0];
    for(iC1=0;iC1<iChans;iC1++){
      for(iC2=iC1+1;iC2<iChans;iC2++){
	for(i=0;i<iWinSz;i++) 
          *pc += pTSN[iC1][i]*pTSN[iC2][i];
        *pc /= (double) iWinSz;
        pc++;
      }
    }
  }
  dRet = 1.0;
  printf("\n");

FCLEANUP:
  if(pTS) FreeListVec(&pTS);
  if(pTSN) freedouble2D(&pTSN,iChans);
  if(psums) free(psums);
  if(psums2) free(psums2);
  if(pfft) free(pfft);
  if(pdiff) free(pdiff);
  return dRet;
  ENDVERBATIM
}

: doesnt use eigenvalues 
FUNCTION tscoravg () {
  VERBATIM

  int i;

  double dRet = 0.0;
  ListVec* pTS = 0;

  double** pTSN = 0; // normalized time series
  double** pCorrel = 0; //correlation matrix

  double* vTCor = 0; //correlation matrix abs sum avg for a time point
  int corsz = vector_arg_px(1,&vTCor);

  pTS = AllocListVec(*hoc_objgetarg(2)); //list of eeg time series
  
  if(!pTS){
    printf("tscoravg ERRA: problem initializing 1st 2 args!\n");
    goto CLEANUP;
  }

  int iChans = pTS->isz; // # of channels
  if(iChans < 2) {
    printf("tscoravg ERRB: must have at least 2 EEG channels!\n");
    goto CLEANUP;
  }
  int iWinSz = (int) *getarg(3); //window size for analysis


  int iSamples = pTS->plen[0]; //make sure all time series have same length

  int iINC = ifarg(4)?(int)*getarg(4):1; 

  double *phz = 0, *psat = 0, *pfft = 0; int sz = 0 , iNoiseTest = 0;
  if(ifarg(5) && ifarg(6)){ //do saturation and 60/120 Hz Frequency test
    if((sz=vector_arg_px(5,&phz))!=corsz){
      printf("tscoravg ERRG: invalid size for phz should be %d, is %d!\n",corsz,sz);
      goto CLEANUP;
    } else if((sz=vector_arg_px(6,&psat))!=corsz){
      printf("tscoravg ERRH: invalid size for psat should be %d, is %d!\n",corsz,sz);
      goto CLEANUP;
    } else {	      
      iNoiseTest = 1;
      pfft = (double*)calloc(iWinSz+1,sizeof(double));
    }
  }

  int idoabs = ifarg(7)?(int)*getarg(7):0; //do abs or not

  int iOutSz = (iSamples - iWinSz + 1) / iINC;//output List size

  double *pdelta=0,*ptheta=0,*pbeta=0,*pgamma=0,*pripple=0;
  if(ifarg(8) && vector_arg_px(8,&pdelta)<iOutSz){ printf("bad pdelta size\n"); goto CLEANUP; } //delta 0-3
  if(ifarg(9) && vector_arg_px(9,&ptheta)<iOutSz){ printf("bad ptheta size\n"); goto CLEANUP; } //theta 4-12
  if(ifarg(10)&& vector_arg_px(10,&pbeta)<iOutSz){ printf("bad pbeta size\n"); goto CLEANUP; } //beta 13-29
  if(ifarg(11)&& vector_arg_px(11,&pgamma)<iOutSz){ printf("bad pgamma size\n"); goto CLEANUP; } //gamma 30-100
  if(ifarg(12)&& vector_arg_px(12,&pripple)<iOutSz){ printf("bad pripple size\n"); goto CLEANUP; } //ripple 101-300

  if(verbose) printf("iWinSz=%d iINC=%d iSamples=%d\n",iWinSz,iINC,iSamples);
  for(i=1;i<iChans;i++) {
    if(pTS->plen[i] != iSamples){
      printf("tscoravg ERRC: time series of unequal size %d %d!\n",iSamples,pTS->plen[i]);
      goto CLEANUP;
    }
  }

  if(corsz < iOutSz){
    printf("tscoravg ERRD: need output vector of at least %d as arg 1, have only %d!\n",iOutSz,corsz);
    goto CLEANUP;
  }
  if(verbose) printf("iOutSz=%d\n",iOutSz);

  pCorrel = getdouble2D(iChans,iChans); if(verbose) printf("got pCorrel\n");
  pTSN = getdouble2D(iChans,iWinSz); if(verbose) printf("got pTSN\n");
  if(!pCorrel || !pTSN){
    printf("tscor ERRD: out of memory!\n");
    goto CLEANUP;
  }

  double* psums = (double*)calloc(iChans,sizeof(double));
  double* psums2 = (double*)calloc(iChans,sizeof(double));

  double dVal = 0.0;
  int iStartIDX = 0 , ierr = 0 , iOutIDX = 0;
  for(iOutIDX=0;iOutIDX<corsz;iOutIDX++) vTCor[iOutIDX] = 0.;
  if( phz ) for(iOutIDX=0;iOutIDX<corsz;iOutIDX++) phz[iOutIDX]=0.;
  if( psat ) for(iOutIDX=0;iOutIDX<corsz;iOutIDX++) psat[iOutIDX]=0.;
  iOutIDX=0;
  for(iStartIDX=0;iOutIDX<iOutSz && iStartIDX+iWinSz<iSamples;iStartIDX+=iINC,iOutIDX++){
    if(iOutIDX%1000==0) printf("%d\n",iOutIDX);
    int IDX = iStartIDX , iChan = 0 , iEndIDX = iStartIDX + iWinSz;
    int iCurSz = iEndIDX - iStartIDX;
    if(phz) phz[iOutIDX]=0.; if(psat)psat[iOutIDX]=0.;
    if(pdelta) pdelta[iOutIDX]=0.; if(ptheta) ptheta[iOutIDX]=0.;
    if(pbeta) pbeta[iOutIDX]=0.; if(pgamma) pgamma[iOutIDX]=0.; if(pripple) pripple[iOutIDX]=0.;
    if(iStartIDX==0 || iINC > 1){
      for(iChan=0;iChan<iChans;iChan++){ //normalize data
	psums[iChan]=psums2[iChan]=0.;
	for(IDX=iStartIDX;IDX<iEndIDX;IDX++){
	  dVal = pTS->pv[iChan][IDX];
	  psums[iChan] += dVal;
	  psums2[iChan] += dVal*dVal;
	}
      }
    } else {
      for(iChan=0;iChan<iChans;iChan++){
	dVal = pTS->pv[iChan][iStartIDX-1];
	psums[iChan] -= dVal;
	psums2[iChan] -= dVal*dVal;
	dVal = pTS->pv[iChan][iEndIDX-1];
	psums[iChan] += dVal;
	psums2[iChan] += dVal*dVal;
      }
    }
    double dsatavg = 0., dhzavg = 0.;
    for(iChan=0;iChan<iChans;iChan++){       
      i=0;
      double davg = psums[iChan] / iCurSz;
      double dstdev = psums2[iChan]/iCurSz - davg*davg; 
      if(dstdev > 0. ) dstdev = sqrt(dstdev); else dstdev=1.; //make sure no problems with nan
      if(dstdev <= 0.) dstdev = 1.0;
      if(iNoiseTest){	
	//saturation test
	int iSat = 0;
	for(IDX=iStartIDX;IDX<iEndIDX;IDX++,i++){	
	  dVal = pTS->pv[iChan][IDX];
	  if((dVal>=900. && dVal<=1010.) || (dVal>=-1010. && dVal<=-900.)) iSat++;
	  pTSN[iChan][i] = (dVal - davg) / dstdev;	
	}
	dsatavg += (double)iSat/ (double)iWinSz;
	//60,120 Hz Frequency test
	int fftlen = 0; memset(pfft,0,sizeof(double)*(iWinSz+1));
	if(dfftpow(&pTS->pv[iChan][iStartIDX],iWinSz,pfft,iWinSz+1,&fftlen)){
	  int f = 0; pfft[0]=0.0; double fsum = 0.0;
	  for(f=0;f<fftlen && f<=1000;f++) fsum += pfft[f]; //normalize power btwn 0-1000Hz to 1.0
	  for(f=55;f<=65 && f<fftlen;f++) dhzavg += pfft[f] / fsum; //60Hz contribution
	  for(f=115;f<=125 && f<fftlen;f++) dhzavg += pfft[f] / fsum; //120Hz contribution
	  if(pdelta) {
	    for(f=0;f<=3;f++) pdelta[iOutIDX] += pfft[f]/fsum;
	    if(ptheta) for(f=4;f<=12;f++) ptheta[iOutIDX] += pfft[f]/fsum;
	    if(pbeta) for(f=13;f<=29;f++) pbeta[iOutIDX] += pfft[f]/fsum;
	    if(pgamma) for(f=30;f<=100;f++) pgamma[iOutIDX] += pfft[f]/fsum;
	    if(pripple) for(f=101;f<=300 && f<fftlen;f++) pripple[iOutIDX] += pfft[f]/fsum;
	  }
	}	
      } else {
	for(IDX=iStartIDX;IDX<iEndIDX;IDX++,i++){	
	  dVal = pTS->pv[iChan][IDX];
	  pTSN[iChan][i] = (dVal - davg) / dstdev;	
	}
      }
    }
    if(iNoiseTest){
      phz[iOutIDX] = dhzavg / (double)iChans;
      psat[iOutIDX] = dsatavg / (double)iChans;
      if(pdelta) pdelta[iOutIDX] /= (double)iChans;
      if(ptheta) ptheta[iOutIDX] /= (double)iChans;
      if(pbeta) pbeta[iOutIDX] /= (double)iChans;
      if(pgamma) pgamma[iOutIDX] /= (double)iChans;
      if(pripple) pripple[iOutIDX] /= (double)iChans;
    }

    //form correlation matrix
    int iC1,iC2;
    for(iC1=0;iC1<iChans;iC1++){
      for(iC2=0;iC2<=iC1;iC2++){
	if(iC1==iC2){
	  pCorrel[iC1][iC2]=1.0; //diagonals == 1.0
	  continue;
	} else pCorrel[iC1][iC2]=0.; //init to 0.
	for(i=0;i<iWinSz;i++){
	  pCorrel[iC1][iC2] += pTSN[iC1][i]*pTSN[iC2][i];
	}
	pCorrel[iC1][iC2] /= (double)iWinSz;
	pCorrel[iC2][iC1] = pCorrel[iC1][iC2];
      }
    }
    double dpsum = 0.0; int icc = 0; //get avg of abs val of correlations
    if(idoabs){
      for(iC1=0;iC1<iChans;iC1++){
	for(iC2=0;iC2<iC1;iC2++){
	  dpsum += fabs(pCorrel[iC1][iC2]);
	  icc++;
	}
      }
    } else {
      for(iC1=0;iC1<iChans;iC1++){
	for(iC2=0;iC2<iC1;iC2++){
	  dpsum += pCorrel[iC1][iC2];
	  icc++;
	}
      }
    }
    dpsum /= (double) icc;
    if(dpsum > 1. || dpsum < -1.){
      printf("out of bounds = %g!\n",dpsum);
    }
    vTCor[iOutIDX] = dpsum;
  }
  dRet = 1.0;
  printf("\n");

CLEANUP:
  if(pTS) FreeListVec(&pTS);
  if(pCorrel) freedouble2D(&pCorrel,iChans);
  if(pTSN) freedouble2D(&pTSN,iChans);
  free(psums);
  free(psums2);
  if(pfft) free(pfft);
  return dRet;
  ENDVERBATIM
}


FUNCTION tscor () {
  VERBATIM

  int i;

  double dRet = 0.0;

  ListVec* pEIG = 0, *pTS = 0;

  double** pTSN = 0; // normalized time series
  double** pCorrel = 0; //correlation matrix
  double** pV = 0;
  double* pW = 0; // eigs

  pEIG = AllocListVec(*hoc_objgetarg(1)); //list of eigenvalue vectors, one for each time
  pTS = AllocListVec(*hoc_objgetarg(2)); //list of eeg time series
  
  if(!pEIG || !pTS){
    printf("tscor ERRA: problem initializing 1st 2 args!\n");
    goto CLEANUP;
  }

  int iChans = pTS->isz; // # of channels
  if(iChans < 2) {
    printf("tscor ERRB: must have at least 2 EEG channels!\n");
    goto CLEANUP;
  }
  int iWinSz = (int) *getarg(3); //window size for analysis
  int iINC = (int) *getarg(4); //increment for analysis

  int iSamples = pTS->plen[0]; //make sure all time series have same length

  if(verbose) printf("iWinSz=%d iINC=%d iSamples=%d\n",iWinSz,iINC,iSamples);
  for(i=1;i<iChans;i++) {
    if(pTS->plen[i] != iSamples){
      printf("tscor ERRC: time series of unequal size %d %d!\n",iSamples,pTS->plen[i]);
      goto CLEANUP;
    }
  }
  int iOutSz = (iSamples - iWinSz + 1) / iINC;//output List size
  if(pEIG->isz < iOutSz){
    printf("tscor ERRD: need List of size at least %d as arg 1, have only %d!\n",iOutSz,pEIG->isz);
    goto CLEANUP;
  }
  if(verbose) printf("iOutSz=%d\n",iOutSz);

  pCorrel = getdouble2D(iChans,iChans); if(verbose) printf("got pCorrel\n");
  pTSN = getdouble2D(iChans,iWinSz); if(verbose) printf("got pTSN\n");
  pV = getdouble2D(iChans,iChans);   if(verbose) printf("got pV\n");
  pW = (double*) malloc(sizeof(double)*iChans); if(verbose) printf("got pW\n");
  if(!pCorrel || !pTSN || !pV || !pW){
    printf("tscor ERRD: out of memory!\n");
    goto CLEANUP;
  }

  double dVal = 0.0;
  int iStartIDX = 0 , ierr = 0 , iTIDX = 0;
  for(iStartIDX=0;iStartIDX+iWinSz<iSamples;iStartIDX+=iINC,iTIDX++){
    if(iStartIDX%100==0) printf("%d\n",iStartIDX);
    int IDX = iStartIDX , iChan = 0 , iEndIDX = iStartIDX + iWinSz;
    int iCurSz = iEndIDX - iStartIDX;
    for(iChan=0;iChan<iChans;iChan++){ //normalize data
      double dSum = 0.0, dSum2 = 0.0;
      for(IDX=iStartIDX;IDX<iEndIDX;IDX++){
        dVal = pTS->pv[iChan][IDX];
	dSum += dVal;
	dSum2 += dVal*dVal;
      }
      dSum /= iCurSz; //mean
      dSum2 /= iCurSz; 
      dSum2 -= dSum*dSum;
      dSum2 = sqrt(dSum2); //standard deviation
      i=0;
      for(IDX=iStartIDX;IDX<iEndIDX;IDX++,i++){
	pTSN[iChan][i] = (pTS->pv[iChan][IDX] - dSum) / dSum2;
      }
    }
    //form correlation matrix
    int iC1,iC2;
    for(iC1=0;iC1<iChans;iC1++){
      for(iC2=0;iC2<=iC1;iC2++){
	if(iC1==iC2){
	  pCorrel[iC1][iC2]=1.0; //diagonals == 1.0
	  continue;
	} else pCorrel[iC1][iC2]=0.;
	for(i=0;i<iWinSz;i++){
	  pCorrel[iC1][iC2] += pTSN[iC1][i]*pTSN[iC2][i];
	}
	pCorrel[iC1][iC2] /= (double)iWinSz;
	pCorrel[iC2][iC1] = pCorrel[iC1][iC2];
      }
    }
    //perform SVD to get eigenvalues into pW
    mysvd(iChans,iChans,pCorrel,pW,pV,&ierr); 
    qsort(pW,iChans,sizeof(double),mycompare); //sort eigenvalues in ascending order
//    for(i=0;i<iChans;i++) pEIG->pv[iTIDX][i] = pW[i]; //store sorted eigenvalues
    memcpy(pEIG->pv[iTIDX],pW,sizeof(double)*iChans); //store sorted eigenvalues
  }
  dRet = 1.0;
  printf("\n");

CLEANUP:
  if(pEIG) FreeListVec(&pEIG);
  if(pTS) FreeListVec(&pTS);
  if(pCorrel) freedouble2D(&pCorrel,iChans);
  if(pTSN) freedouble2D(&pTSN,iChans);
  if(pW) free(pW);
  if(pV) freedouble2D(&pV,iChans);
  return dRet;
  ENDVERBATIM
}

VERBATIM

double getmean(double* p,int n) {
  if(n<1) return 0.0;
  int i = 0;
  double sum = 0.0, *pp=p;
  for(i=0;i<n;i++,pp++) sum += *pp;
  return sum / (double) n;
}

double getstd(double* p,int n,double X) {
  if(n<1) return 0.;
  int i;
  double X2=0., *pp=p;
  for(i=0;i<n;i++,pp++) X2 += (*pp * *pp);
  X2 = X2/(double)n - X*X;
  if(X2>0.) return sqrt(X2);
  return 0.;
}

double dnormv(double* p,int n) {
  if(n<1) return 0.0;
  double X = getmean(p,n);
  double s = getstd(p,n,X);
  double* pp = p;
  int i;
  if(s>0.)
    for(i=0;i<n;i++,pp++) *pp = (*pp - X)/s; 
  else
    for(i=0;i<n;i++,pp++) *pp -= X;
  return 1.0;
}

static double normv(void* v){ 
  double* x;
  int sz = vector_instance_px(v,&x);
  if(sz<1){
    printf("normv ERRA: empty input size!\n");
    return 0.;
  }
  return dnormv(x,sz);
}

double dcopynzidx(double* pin,double* pidx,double* pout,int sz) {
  int szout,i;
  szout=0;
  for(i=0;i<sz;i++) if(pidx[i]) pout[szout++]=pin[i];
  return (double) szout;
}

static double copynzidx (void* v) {
  double *x,*y,*z,ret; int sz;
  sz = vector_instance_px(v,&x);
  if(sz!=vector_arg_px(1,&y)) y=vector_newsize(vector_arg(1),sz);
  if(sz!=vector_arg_px(2,&z)) z=vector_newsize(vector_arg(2),sz);
  ret=dcopynzidx(x,y,z,sz);
  vector_resize(vector_arg(2),(int)ret);
  return ret;
}
ENDVERBATIM

: assumes v1,v2 are normalized
FUNCTION tsmul () {
  VERBATIM
  double* x, *y, *px, *py, sum = 0.;
  int xsz,ysz,i;
  xsz = vector_arg_px(1,&x);
  ysz = vector_arg_px(2,&y);
  px=x; 
  py=y;
  if(xsz!=ysz || xsz<1 || ysz<1){
    printf("tsmul ERRB: input vecs are invalid sizes: %d %d!\n",xsz,ysz);
    return -2.;
  }
  for(i=0;i<xsz;i++,px++,py++) sum += *px * *py;
  return sum / (double) xsz;
  ENDVERBATIM
}

PROCEDURE install () {
  if (INSTALLED==1) {
    printf("already installed tsa.mod")
  } else {
    INSTALLED=1
    VERBATIM
      //    install_vector_method("myspct", myspct);
    install_vector_method("inrange",inrange);
    install_vector_method("normv",normv);
    install_vector_method("copynzidx",copynzidx);
    ENDVERBATIM
  }
}
