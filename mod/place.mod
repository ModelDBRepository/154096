: $Id: place.mod,v 1.141 2011/07/06 15:24:40 samn Exp $
 
:* COMMENT
COMMENT
CODE modified from ~/nrniv/place/fenton/
SimulateTS/MakeGaussianField.c,SimulateTS.c
CorrelateSpikeTrains/CorrelateSpikeTrains.c,Kendall.c
ENDCOMMENT

NEURON {
  SUFFIX PLACE
  GLOBAL  INSTALLED, MAKE_GAUSSIAN, REFRAC
  GLOBAL  DEBUG, PEAK, THRESHOLD
  GLOBAL  samps_per_sec,scale_y,first_time_stamp,scale_x
  GLOBAL  MaxTime,AnalysisInterval,BinSize,IntervalStart,VERBOSE
}

PARAMETER {
  INSTALLED=0
  VERBOSE=0
  DEBUG=0
  MAKE_GAUSSIAN=1
  REFRAC=3
  PEAK = 5.0
  THRESHOLD= 1e-5
  samps_per_sec=0
  scale_y=0
  first_time_stamp=0
  scale_x=0
  MaxTime=1200          : seconds
  AnalysisInterval=1200 : seconds
  IntervalStart=0       : seconds
  BinSize=250
}

ASSIGNED { }

VERBATIM
#include "bpf.h"
#include "misc.h"
void KendallTauPairs(si1 *dx, si1 *dy, ui4 NumEvents, ui4 NumPairs, sf8 *tau, sf8 *z, sf8 *prob);

static float   **rate_map;
static double  storage[10];

#define MAP_YSZ		32
#define MAP_XSZ		32
#define YSZ		256
#define XSZ		256
#define	MAX_RATE	100.0
#define	MAX_SPKS_PER_SAMP	100
#define	TIMESTAMP_UNITS_PER_SEC	10000
#define	TRUE		1
#define	FALSE		0
#define         PARTS_PER_SAMP  2

//* from MakeGaussianField.c
static double  gaussaf (mu, sd, x)
        double mu, sd, x;
{ double exponent, num, denom, z;
 double value, max_value;

 z = (x - mu) / sd;
 exponent = -0.5 * z * z;

 max_value = (1.0 / (sd * (double)SQRT2PI));
 value = (exp(exponent) / (sd * (double)SQRT2PI));
 value /= max_value;
 return(value);
}

// MakeGaussianField(int x, int y, double sd, int xsz, int ysz, float peak, float
// threshold, float **map)
static void MakeGaussianField (x, y, sd, xsz, ysz, peak, threshold, map)
        int x, y, xsz, ysz;
        double  sd;
        float   peak, threshold, **map;
{ int   i, j;
  double        hypot(), d;

  for(i = 0; i < ysz; i++){
    for(j = 0; j < xsz; j++){
      d = hypot((double)(i - y), (double)(j - x));
      map[i][j] = peak * gaussaf(0.0, sd, d);
      if(map[i][j] < threshold) map[i][j] = 0.0;
      if (DEBUG==2 && map[i][j]!=0) printf("MAP: %d %d %g\n", i,j,map[i][j]);
    }
  }
  return;
}

//* From generate_spks.c
static int generate_spks(rate, max_rate, samps_per_sec)
        float rate, max_rate;
        double samps_per_sec;
{ register int i;
 register double prob;
 double  drand48();
 int     reduction_factor = 1, parts_per_samp = PARTS_PER_SAMP;
 int     n_spks;

 while(rate >= max_rate){
   rate /= 2.0;
   reduction_factor *= 2;
 }
 prob = (double)rate / (double)samps_per_sec / (double)parts_per_samp;

 n_spks = 0;
 for(i = parts_per_samp * reduction_factor; i--; )
   if(prob > (-1.0 * log(drand48())))
     n_spks++;
 return n_spks;
}

//* original SimulateTC.c
// why not just use fgets?
static char    *get_line(FILE* f,char* s,int iMaxLen)
{ int     i=0;
 while(fscanf(f,"%c",&s[i]) != EOF && i<iMaxLen-1){
   if(s[i] == '\n'){
     s[i+1] = '\0';
     return s;
   }
   i++;
 }
 return(NULL);
}

static void    instruct()
{
  printf("Call with a time_series format file in TS_DIR.\n"); 
  printf("Output to STS_DIR.\n"); 
  printf("Simulates a location-specific spike train using the track in the TS file. Positional rates are derived from a 2-D Gaussian rate map\n"); 
  printf("or a rate map stored as a file. See options.\n"); 
  printf("Options:\n"); 
  printf("\tF <filename> Don't use a Gaussian field for the simulation. Use the rate map stored in RATE_DIR/filename.\n"); 
  printf("\t\tRequires options X and Y\n");
  printf("\tP <cycles per minute> Modulate rate with a periodic function. Give the period (cycles per min).\n"); 
  printf("\tS <Seed the random number generator with your own integer seed. Default is a random seed.\n"); 
  printf("\tX <x resolution of rate map> The x dimension resolution of the rate map array. Necessary for option F. Default: %d\n", MAP_XSZ); 
  printf("\tY <y resolution of rate map> The y dimension resolution of the rate map array. Necessary for option F. Default: %d\n", MAP_YSZ); 
  hxe();
}

static int  ascending_int(unsigned int *a, unsigned int *b)
{ if (*a > *b)
    return(1);
  if (*a < *b)
    return(-1);
  return(0);
}
ENDVERBATIM

:* v1.glob(L_BOUND,R_BOUND,cut) -- group values together and give bounds (indices)
VERBATIM
static double glob (void* vv) {
  int i, j, nx, ny, nz;
  double *x, *y, *z, last, cut, diff;
  nx = vector_instance_px(vv, &x);
  ny = vector_arg_px(1, &y);
  nz = vector_arg_px(2, &z);
  if (ny!=nz || ny<2) {printf("Place:glob ERR0: ny!=nz %d %d\n",ny,nz); hxe();}
  cut=*getarg(3);
  for (last=x[0],i=1,j=0,y[0]=0; i<nx; i++) { 
    diff=x[i]-last;
    if (diff<0) {printf("Place:glob ERRA: not monotonic %g < gd\n",x[i],last); hxe();}
    if (diff>cut) {
      z[j]=(double)i-1;
      j++;
      if (j>=ny) { printf("Place:glob ERRB: OOR %d\n",ny); hxe(); }
      y[j]=(double)i;
    }
    last=x[i];
  }
  z[j]=nx-1; // index at end of vector
  vector_resize(vector_arg(1), j+1);
  vector_resize(vector_arg(2), j+1);
  return (double)j+1;
}
ENDVERBATIM

:* v1.cumul(out) 
VERBATIM
static double cumul (void* vv) {
  int i, j, nx, ny, nyn, fli;
  unsigned int* scr;
  double *x, *y, cut,fl;
  void *vvy;
  nx = vector_instance_px(vv, &x);
  ny = vector_arg_px(1, &y); vvy=vector_arg(1);
  cut=*getarg(2);
  if (nx<2) { printf("cumul:ERRA nx too small %d\n",nx); hxe();}
  vector_resize(vvy,ny=vector_buffer_size(vvy));
  scr=scrset(nx);
  for (i=0;i<nx;i++) scr[i]=i;
  nrn_mlh_gsort(x, scr, nx, cmpdfn);
  fl=floor(x[scr[nx-1]]); // this is the largest value
  if ((fli=(int)fl)>=ny-1) { printf("cumul:ERRB fli>=ny %d %d\n",fli+1,ny); hxe();}
  if (fl>=cut) {
    nyn=fli+2;
    for (i=cut;i<nyn-1;i++) y[i]=1;
    for (i=0;i<cut;i++) y[i]=0;
    y[fli+1]=0;
    for (i=nx-2;i>=0 && x[scr[i]]>=cut;i--) {
      if (x[scr[i]]<fl) fli=(int)(fl=floor(x[scr[i]]));
      for (j=cut;j<=fli;j++) y[j]++;
    }
  } else nyn=0;
  vector_resize(vvy, nyn);
  return (double)nyn;
}
ENDVERBATIM

:* main from SimulateTS.c
VERBATIM
static double simts (void* vv) {
  FILE    *fopen(), *ifp, *ofp, *ofq, *rfp;
  Object* ob;
  unsigned char xys[3];
  int    i,j,k,s,mt,fflag;
  int     op, AngularModulation;
  int     x, y, nmb_header_lines;
  double  center_x, center_y, sd, value;
  double  time;
  double   max_rate;
  char    line[256], keyword[64], *get_line(), rate_file[256];
  char    in_file[256], out_file[256], file_name[256];
  double   modulate, peak, threshold;
  double   CyclesPerMin;
  double   map_x, map_y, collapse_x, collapse_y;
  void    instruct(), MakeGaussianField();
  int     found, n_samps;
  double     sample_time, timestamp_units_per_samp;
  int     SampsPerCycle, SampsPerHalfCycle, samps;

  double *xx, *yy, jit, *loc[2]; 
  int nx,ny,maxsz,ii, q, sflag, nloc[2];
  void *vv1, *vloc[2];

  nx = vector_instance_px(vv, &xx);
  maxsz=vector_buffer_size(vv);
  vector_resize(vv, maxsz);
  sflag=0; jit=0.;
  ob =   *hoc_objgetarg(1);
  if (strncmp(hoc_object_name(ob),"List",4)==0) {
    fflag=0;
    // pick up list of 2 vectors
    for (i=0;i<2;i++) nloc[i]=list_vector_px2(ob, i, &loc[i], &vloc[i]);
    if (nloc[0]!=nloc[1]) {printf("simtst_PLACE ERRA: %d %d\n",nloc[0],nloc[1]); hxe();}
  } else {
    fflag=1;
    ifp=hoc_obj_file_arg(1);
  }
  if (ifarg(2)) {
    if (hoc_is_double_arg(2)) {
      jit= *getarg(2)*10.;  // convert from ms to 100 mus
    } else { // vector for storage
      ny = vector_arg_px(2, &yy);
      vv1=vector_arg(2); 
      if ((ii=vector_buffer_size(vv1))<maxsz) {
        fprintf(stderr, "2nd vec maxsz (%d) should be at least %d\n",ii,maxsz); hxe(); }
        vector_resize(vv1, maxsz);
        sflag=1; 
    }
  }

  q=found=0;
  AngularModulation = FALSE;
  max_rate = MAX_RATE;
  map_x = MAP_XSZ; map_y = MAP_YSZ;

  while (0) switch (op) {
    case    'F':    ;// Gaussian_flag -> replaced with MAKE_GAUSSIAN_PLACE
    // strcpy(rate_file, optarg);
    break;
    case    'P':    AngularModulation = TRUE;
    //    CyclesPerMin = atof(optarg);
    break;
    case    'S':    ; // Seed = atoi(optarg);
    break;
    case    'X':    map_x = MAP_XSZ; //atof(optarg);
    break;
    case    'Y':    map_y = MAP_YSZ; // atof(optarg);
    break;
    default:
    instruct();
  }

  if (rate_map==NULL) {
    if ((rate_map = (float **)calloc((size_t)YSZ, sizeof(float *))) == NULL) {
      fprintf(stderr, "Can't allocate\n"); hxe(); }
    for(i=0; i < YSZ; i++) if ((rate_map[i]=(float *)calloc(XSZ, sizeof(float))) == NULL){
        fprintf(stderr, "Can't allocate\n"); hxe(); }
  }

  if (MAKE_GAUSSIAN_PLACE==1) {
    // pick the center of the field
    center_x = center_y = 0.0;
    while(hypot(center_x - (double)(XSZ/2.0), center_y - (double)(YSZ/2.0)) > (XSZ/2.0)){
      center_x = (double)(XSZ * drand48());
      center_y = (double)(YSZ * drand48());
    }

    peak=PEAK; // AP/s
    threshold = THRESHOLD; // minimum rate AP/s 
    sd = 0.2 * XSZ/2.0;     // sd set to 20% of radius
    MakeGaussianField((int)center_x, (int)center_y, sd, (int)XSZ, (int)YSZ, peak, threshold, rate_map);
    collapse_x = collapse_y = 1.0;
    storage[0]=center_x;storage[1]=center_y;storage[2]=collapse_x;storage[3]=collapse_y;
    storage[4]=peak;storage[5]=threshold;storage[6]=sd;
  } else if (MAKE_GAUSSIAN==2) { // will use to read from vector
    for(y = 0; y < (int)map_y; y++)
      for(x = 0; x < (int)map_x; x++)
	fscanf(rfp,"%f", &(rate_map[y][x]));
    collapse_x = (double)XSZ / (double)map_x;
    collapse_y = (double)YSZ / (double)map_y;
  } else { // MAKE_GAUSSIAN==0: resuse existing map
    center_x=storage[0];center_y=storage[1];collapse_x=storage[2];collapse_y=storage[3];
    peak=storage[4];threshold=storage[5];sd=storage[6];
  }

  nmb_header_lines = found = 0;

  if (fflag) {
    fseek(ifp,0,SEEK_SET);
    fscanf(ifp, "%d\n", &nmb_header_lines);
    for (i = 1; i < nmb_header_lines ; i++) {
      get_line(ifp, line, 256);
      // fprintf(ofp, "%s", line);

      sscanf(line,"%s%lf", keyword, &value);
      if(!strcmp(keyword,"%SAMPLING_INTERVAL(samps/sec)")){
        samps_per_sec = value;
        found++;
      }else if(!strcmp(keyword, "%SCALE_Y(RatioTracktoMapPixels)")){
        scale_y = value;
        found++;
      }else if(!strcmp(keyword, "%FIRST_TIMESTAMP(100usUnits)")){
        first_time_stamp = value;
        found++;
      }else if(!strcmp(keyword, "%SCALE_X(RatioTracktoMapPixels)")){
        scale_x = value;
        found++;
      }
    }
    if (found != 4) {
      printf("PLACE simts() Only found %d/4. File is corrupt\n",found);
      hxe(); }
  }

  // deliberately rounding down
  timestamp_units_per_samp= (int)(TIMESTAMP_UNITS_PER_SEC/samps_per_sec);
  if(AngularModulation) {
    SampsPerCycle = (int)(samps_per_sec * 60.0 / CyclesPerMin);
    SampsPerHalfCycle = (int)(SampsPerCycle / 2.0);
  }

  n_samps=0; modulate = 1.0;
  while (1) {
    if (fflag) {
      if (fread(xys,sizeof(unsigned char), 3, ifp) != 3) break;
      for(i = 0; i < xys[2]; i++) fread(&mt,sizeof(int), 1, ifp); // throw times away
      x = xys[0]; y = xys[1];
    } else {
      if (n_samps>=nloc[0]) break;
      x = loc[0][n_samps]; y=loc[1][n_samps];
    }
    n_samps++;
    // determine rate modulation
    if (AngularModulation){
      samps = (n_samps % SampsPerCycle);
      if(samps > SampsPerHalfCycle)
	samps = SampsPerCycle - samps;
      modulate = (double)samps / (double)SampsPerHalfCycle;
    }
    // use the collapse factors in case the rate map is a different resolution than the track
    s = generate_spks(modulate*rate_map[(int)(y/collapse_y)][(int)(x/ collapse_x)],\
                      max_rate, samps_per_sec);
    // xys[2] = s;
    // write x y s -- CAN KILL (x loc, y loc, s # of spikes) ****************
    // fwrite(xys, sizeof(unsigned char), (size_t)3, ofp);
    if (sflag) {
      if (s>0) {
        xx[q]=(double)(n_samps + drand48())*timestamp_units_per_samp*0.1; // 0.1 converts to ms
        yy[q]=(double)s;
        q++;
      }
    } else {
      sample_time = n_samps * timestamp_units_per_samp;
      if (q>0 && s>0 && jit>0.) { 
        sample_time+=(jit*(1.-2.*drand48()));
        while(sample_time<xx[q-1]) { // don't allow regress in time
          sample_time=n_samps*timestamp_units_per_samp+(jit*(1.-2.*drand48())); }
      }
      for (i=0; i<s && q<maxsz; i++,q++) {
        if (i==0) {
          xx[q]=0.1*(sample_time+drand48()*timestamp_units_per_samp); 
          while (q>0 && xx[q]<xx[q-1]+REFRAC) xx[q]+=drand48()*timestamp_units_per_samp; 
        } else {
          xx[q]=xx[q-1]+REFRAC+0.1*drand48()*timestamp_units_per_samp;
        }
        if (DEBUG==1) printf("%d:%20.2lf:%20.2lf %d\n",q,sample_time,xx[q],s);
      }
    }
    if (q>maxsz) break;
    // write spike times -- CHANGE HERE ****************
    // fwrite(time, sizeof(int), (size_t)s, ofp); // this is the time in time
  }
  vector_resize(vv, q);
  if (sflag) vector_resize(vv1, q);
  return (double)q;
}
ENDVERBATIM

:* outvec.kendall(srcvecs,tmpvecs)
VERBATIM
static double kendall (void* vv) {
  ui4 i, j, k, m, n, p, q, maxsz, nx, av[VRRY], bv[VRRY], num, numb, indflag, match;
  ui4	MaxBins, NumPairs, *SpikeCounts;
  si1	**Diffs, dif;
  Object *ob, *ob2;
  double *x, *avo[VRRY], *bvo[VRRY], val[VRRY];
  double IntStart, IntervalStop, AnInt, MaxT;
  double tau, z, prob;
  void *vva[VRRY],*vvb[VRRY];
  void KendallTauPairs();
  nx = vector_instance_px(vv, &x);
  ob = *hoc_objgetarg(1);
  ob2 = *hoc_objgetarg(2);
  num = ivoc_list_count(ob);
  numb = ivoc_list_count(ob2);
  if (num>VRRY) hoc_execerror("kendall ****ERRA****: can only handle VRRY vectors", 0);
  if (numb!=4) hoc_execerror("kendall ****ERRB****: 2nd list count should be 4:i,j,tau,p", 0);
  for (i=0;i<num;i++) av[i]=list_vector_px2(ob, i, &avo[i], &vva[i]);
  for (i=0;i<numb;i++) { 
    bv[i]=list_vector_px3(ob2, i, &bvo[i], &vvb[i],1);
    if (bv[0]!=bv[i]) { printf("kendall ****ERRC2**** %d %d %d\n",i,bv[0],bv[i]);
      hoc_execerror("Vectors must all be same size: ", 0); }
  }
  Diffs = (si1 **)calloc(num, sizeof(si1 *));
  AnInt=AnalysisInterval*1e3; // convert from sec to ms
  MaxT=MaxTime*1e3; // convert from sec to ms
  MaxBins = (ui4) AnInt/BinSize;
  NumPairs = (ui4)MaxBins*(MaxBins-1); // not just neighboring bins

  SpikeCounts = (ui4 *)calloc(MaxBins, sizeof(ui4));
  for(i=0;i<num;i++) {
    Diffs[i] = (si1 *)calloc(NumPairs, sizeof(si1));
  }
  
  // do correlation over segments of the data between Start and Stop
  // MaxTime is the full duration of the recording
  IntStart = IntervalStart*1e3; // s to ms
  IntervalStop = IntStart + AnInt;
  while (IntervalStop <= MaxT) {
    // Initialize the Spike Count bins
    for (i=0;i<num;i++) {
      for(j=0;j<MaxBins;j++) SpikeCounts[j]=0;
      for (j=0;j<av[i];j++) {
        if (avo[i][j]<IntStart) continue;
        if (avo[i][j]>IntervalStop) break;
        SpikeCounts[(int)((avo[i][j]-IntStart)/BinSize)]++;
      }

      // PairwiseDiffSigns(SpikeCounts[i], Diffs[i], MaxBins, NumPairs);
      // SpikeCounts->data, Diffs->signs, MaxBins->NumEvents, NumPairs->NumPairs
      for (k=0,n=0;k<MaxBins;k++) {
        for (j=(k+1);j<MaxBins;j++,n++) {
          if (dif=(SpikeCounts[k]-SpikeCounts[j])) {
            Diffs[i][n]=(dif>0)?1:-1; 
          } else Diffs[i][n]=0;
        }
      }
    }
    // calculate Kendall's tau for each pair of time series
    for(i=0,k=0; i<num; i++) {
      for(j=i+1; j<num; j++,k++) {
        KendallTauPairs(Diffs[i], Diffs[j], MaxBins, NumPairs, &tau, &z, &prob);
        if (VERBOSE) printf("%g  %g  %d  %d  %0.6lf  %0.3lf  %0.3lf\n",\
               IntStart/1e3, IntervalStop/1e3, i, j, tau, z, prob);
        // i j tau prob -- save i, j, tau, prob in an NQS given by 2nd list arg
        if (k>=bv[0]){printf("PLACE kendall() ERRD out of room: %d %d\n",k,bv[0]); hxe();}
        bvo[0][k]=(double)i;bvo[1][k]=(double)j;bvo[2][k]=tau;bvo[3][k]=prob;
      }
    }
    IntStart+=AnInt;
    IntervalStop +=AnInt;
  }
  for (j=0;j<numb;j++) vector_resize(vvb[j], k);
  return (double)k;
}

static double kendal2 (void* vv) {
  ui4 i, n, j, k, MaxBins, NumPairs;
  double *x, *y, tau, z, prob;;
  si1	*Dx, *Dy; 
  si4 difx, dify;
  ui4 LeftNum = 0, RightNum = 0;
  si4 Numerator = 0;
  si1 dxdy;
  sf8 Variance;
  ui4 NumEvents;
  MaxBins = vector_instance_px(vv, &x);
  if ((i=vector_arg_px(1,&y)) != MaxBins ) {printf("kenall ERRA: %d %d\n",n,MaxBins); hxe();}
  NumPairs = (ui4)MaxBins*(MaxBins-1)/2; // not just neighboring bins
  if(VERBOSE>1) printf("NumPairs=%d, MaxBins=%d\n",NumPairs,MaxBins);
  Dx = (si1 *)calloc(NumPairs, sizeof(si1)); // sets memory to 0
  if(!Dx) { printf("kendal2 ERR0: out of mem!\n"); hxe(); }
  Dy = (si1 *)calloc(NumPairs, sizeof(si1));
  if(!Dy) { printf("kendal2 ERR1: out of mem!\n"); hxe(); }
  for (k=0,n=0;k<MaxBins;k++) for (j=(k+1);j<MaxBins;j++,n++) {
    if (difx=x[k]-x[j]) Dx[n]=(difx>0)?1:-1; 
    if (dify=y[k]-y[j]) Dy[n]=(dify>0)?1:-1; 
  }
  if(VERBOSE) printf("n=%d, NumPairs=%d\n",n,NumPairs);
  KendallTauPairs(Dx, Dy, MaxBins, NumPairs, &tau, &z, &prob);
  if(VERBOSE) printf("tau:%g z:%g p:%g\n",tau, z, prob);
  free(Dx); free(Dy);
  return tau;
}

// kend2() is like kendall but only looks at neighboring pairs
static double kend2 (void* vv) {
  ui4 i, n, j, k, pw, MaxBins, NumPairs;
  double *x, *y, tau, z, prob;;
  si1	*Dx, *Dy; 
  si4 difx, dify;
  MaxBins = vector_instance_px(vv, &x);
  if ((i=vector_arg_px(1,&y)) != MaxBins ) {printf("kend2 ERRA: %d %d\n",n,MaxBins); hxe();}
  if (ifarg(2)) pw=1; else pw=0; // pairwise flag
  NumPairs = (ui4)(MaxBins-1); // just neighboring bins
  Dx = (si1 *)calloc(NumPairs, sizeof(si1)); // sets memory to 0
  Dy = (si1 *)calloc(NumPairs, sizeof(si1));
  for (k=0,n=0;k<MaxBins-1;k++,n++) { j=(k+1);
    if (difx=x[k]-x[j]) Dx[n]=(difx>0)?1:-1; 
    if (dify=y[k]-y[j]) Dy[n]=(dify>0)?1:-1; 
  }
  KendallTauPairs(Dx, Dy, MaxBins, NumPairs, &tau, &z, &prob);
  if(VERBOSE>0) printf("tau:%g z:%g p:%g\n",tau, z, prob);
  free(Dx); free(Dy);
  return tau;
}

// Given data arrays x[1..n] and y[1..n], this program returns Kendall's correlation tau,
// its number of standard deviations from zero as z, and its two-tailed probability level as prob.
// Small values of prob indicate a significant correlation (tau positive) or anticorrelation (tau
// negative).
//
// the necessary formulas are:
// tau = numerator / denominator
// 	numerator = (concordant_pairs - discordant_pairs)
// 	denominator = sqrt(concordant_pairs + discordant_pairs+ extra-y) * sqrt(concordant_pairs + discordant_pairs+ extra-x)
// 		denominator = sqrt(LeftNum) * sqrt(RightNum)
//	where concordant pairs are those differences (xi - xj) that have the same sign as the differences (yi - yj)
//	extra-y are those non-zero differences (yi - yj) where (xi - xj) = 0
//	extra-x are those non-zero differences (xi - xj) where (yi - yj) = 0
// variance(tau) = (4N + 10) / 9N(N - 1)
// z = tau / sqrt(variance)
// Given arrays dx[1..n] and dy[1..n] which represent the signs (-1, 0, 1) of all
// pairs in the data sets x and y,
// this program returns Kendall's correlation tau. dx and dy are the result of
// calling PairwiseDiffSigns
// its number of standard deviations from zero as z, and its two-tailed probability level as prob.
// Small values of prob indicate a significant correlation (tau positive) or
// anticorrelation (tau negative).
void KendallTauPairs (si1 *dx, si1 *dy, ui4 NumEvents, ui4 NumPairs, sf8 *tau, sf8 *z, sf8 *prob)
{
  sf8 erfcc(sf8 x);
  ui4 i;
  ui4 LeftNum = 0, RightNum = 0;
  si4 Numerator=0;
  si1 dxdy;
  sf8 Variance;

  for (i=0; i < NumPairs; i++) {
    dxdy = dx[i] * dy[i];

    if (dxdy) { // Neither array has a tie.
      LeftNum++;
      RightNum++;
      dxdy > 0.0 ? Numerator++ : Numerator--;
    } else { // One or both arrays have ties.
      if(dy[i])
	LeftNum++;	//an extra y event
      else if(dx[i])
	RightNum++;	//an extra x event
    }
  }

  if(LeftNum && RightNum){
    *tau = (sf8) Numerator / (sqrt((sf8)LeftNum) * sqrt((sf8)RightNum));
    Variance = (sf8)(4.0 * NumEvents + 10.0) / (sf8)(9.0 * NumEvents * (NumEvents - 1.0));
    *z = (*tau) / sqrt(Variance);
    *prob = erfcc(fabs(*z)/1.4142136);
  }else{
    *tau = *z = 0.0;
    *prob = 1.0;
  }
  return;
}

// Returns the complementary error function erfc(x) with fractional error
//    everywhere less than 1.2 x 10^-7.
double erfcc(double x)
{
	double	mt,z,ans;
	z=fabs(x);
	mt=1.0/(1.0+0.5*z);
	ans=mt*exp(-z*z-1.26551223+mt*(1.00002368+mt*(0.37409196+mt*(0.09678418+ mt*(-0.18628806+mt*(0.27886807+mt*(-1.13520398+mt*(1.48851587+ mt*(-0.82215223+mt*0.17087277)))))))));
	return x >= 0.0 ? ans : 2.0-ans;
}
ENDVERBATIM

VERBATIM
static double mkgaussfield(void* vv) {
// NOT DEBUGGED
  double x,y,peak,threshold,sd,d, *map;
  int i,j,nx,maxsz;

  nx = vector_instance_px(vv, &map);
  maxsz=vector_buffer_size(vv);
  if (maxsz<XSZ*YSZ){printf("mkgasussfield ERR vector too small %d<%d",maxsz,XSZ*YSZ); hxe();}
  vector_resize(vv, XSZ*YSZ);

  // pick the center of the field
  x = y = 0.0;
  while(hypot(x - (double)(XSZ/2.0), y - (double)(YSZ/2.0)) > (XSZ/2.0)){
    x = (double)(XSZ * drand48());
    y = (double)(YSZ * drand48());
  }

  peak = 60.0; // AP/s
  threshold = 0.5; // minimum rate AP/s 
  sd = 0.2 * XSZ/2.0;     // sd set to 20% of radius

  for (i=0; i<YSZ; i++) { // YSZ rows
    for (j=0; j<XSZ; j++) { // XSZ cols
      d = hypot((double)(i-y), (double)(j-x));
      map[i*XSZ+j] = peak * gaussaf(0.0, sd, d);
      // if (map[i*XSZ+j]!=0) printf("BB: %d %d %f\n", i,j,map[i*XSZ+j]);
      if(map[i*XSZ+j] < threshold)
        map[i*XSZ+j] = 0.0;
    }
  }
}
ENDVERBATIM

VERBATIM
static double dumpratemap(void* vv) {
  int i,j,k,nx,maxsz;
  double *map;

  if (rate_map==NULL) {printf("rate_map not allocated\n"); hxe();}
  nx = vector_instance_px(vv, &map);
  maxsz=vector_buffer_size(vv);
  if (maxsz<XSZ*YSZ){printf("dumpratemap ERR vector too small %d<%d\n",maxsz,XSZ*YSZ); hxe();}
  vector_resize(vv, XSZ*YSZ);

  for (i=0; i<YSZ; i++) { // YSZ rows
    for (j=0; j<XSZ; j++) { // XSZ cols
      map[i*XSZ+j]=rate_map[i][j];
    }
  }
}
ENDVERBATIM

: from /u/billl/nrniv/place/fenton/read_time_series/read_time_series.c
: took out the BIG_ENDIAN byte reversal stuff

:* v.rd()
VERBATIM
static double rd (void* vv) {
  unsigned char x, y, s;
  unsigned int	i, j, n, time, beg, nx,ny, maxnx, locflag, linflag;
  unsigned int maxsz, k, nmb_header_lines, found;
  int err;
  double *xx, *xloc, *yloc, *xlin, value, errtol;
  void *vv1, *vv2, *vv4;
  char	line[4196], keyword[4196];
  FILE* fp, *hoc_obj_file_arg();
  fp = hoc_obj_file_arg(1);
  vector_instance_px(vv, &xx);
  if((maxsz=vector_buffer_size(vv))<1) {maxsz=100; xx=vector_newsize(vv,maxsz);} else vector_resize(vv, maxsz);
  err=-1;  maxnx=locflag=linflag=0;
  if (ifarg(3)) { 
    nx = vector_arg_px(2, &xloc); vv1=vector_arg(2); 
    maxnx=vector_buffer_size(vv1); 
    ny = vector_arg_px(3, &yloc); vv2=vector_arg(3); 
    k=vector_buffer_size(vv2); 
    if (k!=maxnx) {printf("rd_PLACE: ERRA not same size: %d %d\n",k,maxnx); hxe();}
    vector_resize(vv1, maxnx); vector_resize(vv2, maxnx);
    locflag=1; 
  } else if (ifarg(2) && hoc_is_double_arg(2)) {
    errtol=*getarg(2); err=0;
  } else if (ifarg(2)) { linflag=1;
    vv1=vector_arg(2); 
    xlin=vector_newsize(vv1,maxsz);
  }
  fseek(fp,0,SEEK_SET);
  found=0;
  fscanf(fp, "%d\n", &nmb_header_lines);
  for (i = 1; i < nmb_header_lines ; i++) {
    get_line(fp, line, 256);
    sscanf(line,"%s%lf", keyword, &value);
    if(!strcmp(keyword,"%SAMPLING_INTERVAL(samps/sec)")){
      samps_per_sec = value;
      found++;
    }else if(!strcmp(keyword, "%SCALE_Y(RatioTracktoMapPixels)")){
      scale_y = value;
      found++;
    }else if(!strcmp(keyword, "%FIRST_TIMESTAMP(100usUnits)")){
      first_time_stamp = value;
      found++;
    }else if(!strcmp(keyword, "%SCALE_X(RatioTracktoMapPixels)")){
      scale_x = value;
      found++;
    }
  }
  if (found != 4) {
    printf("PLACE rd() Only found %d/4. File is corrupt\n",found);
    hxe(); }

  for (n=0,j=0,k=0;fread(&x,sizeof(unsigned char), 1, fp) != 0; n++) {
    fread(&y,sizeof(unsigned char), 1, fp);
    fread(&s,sizeof(unsigned char), 1, fp);
    for (i=0; i<s; i++) {
      if(!fread(&time,sizeof(unsigned int), 1, fp))
        printf("rd_PLACE: WARNZ couldnt read spiketime %d of %d at byte %d\n",i+1,s,ftell(fp)); 
      if (k>=maxsz) {
        xx=vector_newsize(vv,maxsz*=2);
        if (linflag) xlin=vector_newsize(vv1,maxsz);
      }
      if (linflag) xlin[k]=(double)n+1.;
      xx[k++]=0.1*(double)(time); // ms
      if (err>=0 && 
          (xx[k-1]<((n-1)*1e3/samps_per_sec-errtol)||xx[k-1]>n*1e3/samps_per_sec+errtol)) {
        if (VERBOSE) {
   printf("rdERRB:%d:%g (%lf:%lf) ",n,xx[k-1],((n-1)*1e3/samps_per_sec),n*1e3/samps_per_sec);
        }
        err++; 
      }
    }
    if (locflag) {
      if (j>=maxnx) {
        maxnx*=2; 
        xloc=vector_newsize(vv1,maxnx); 
        yloc=vector_newsize(vv2,maxnx);
      }
      xloc[j]=(double)x; yloc[j]=(double)y; j++;
    }
  }
  vector_resize(vv, k);
  if (linflag) vector_resize(vv1,k);
  if (locflag) {vector_resize(vv1, j); vector_resize(vv2, j);}
  if (err>0) return (double)-err; else return (double)n*1e3/samps_per_sec;
}
ENDVERBATIM

: modified from rd()
:* spks.rdtts(file,[,xloc,yloc,tloc])
:  read spike times into spks vector; optionally read locations and times into 3 other vecs
VERBATIM
static double rdtts (void* vv) {
  unsigned char x, y, s;
  unsigned int	i, j, n, time, tt, beg, nx,ny, maxnx, locflag;
  unsigned int maxsz, k, nmb_header_lines, found;
  int err;
  double *xx, *tvec, *xloc, *yloc, value, errtol;
  void *vv1, *vv2, *vv3;
  char	line[4196], keyword[4196];
  FILE* fp, *hoc_obj_file_arg();
  fp = hoc_obj_file_arg(1);
  vector_instance_px(vv, &xx);
  if((maxsz=vector_buffer_size(vv))<1) {maxsz=100; xx=vector_newsize(vv,maxsz);} else vector_resize(vv, maxsz);
  err=-1;
  if (ifarg(3)) { 
    nx = vector_arg_px(2, &xloc); vv1=vector_arg(2); 
    maxnx=vector_buffer_size(vv1); 
    ny = vector_arg_px(3, &yloc); vv2=vector_arg(3); 
    ny = vector_arg_px(4, &tvec); vv3=vector_arg(4); 
    xloc=vector_newsize(vv1,maxnx); yloc=vector_newsize(vv2,maxnx); 
    tvec=vector_newsize(vv3,maxnx); locflag=1; 
  } else { 
    maxnx=locflag=0;
    if (ifarg(2)) errtol=*getarg(2); err=0;
  }
  fseek(fp,0,SEEK_SET);
  found=0;
  fscanf(fp, "%d\n", &nmb_header_lines);
  for (i = 1; i < nmb_header_lines ; i++) {
    get_line(fp, line, 256);
    sscanf(line,"%s%lf", keyword, &value);
    if(!strcmp(keyword,"%SAMPLING_INTERVAL(samps/sec)")){
      samps_per_sec = value;      found++;
    }else if(!strcmp(keyword, "%SCALE_Y(RatioTracktoMapPixels)")){
      scale_y = value;      found++;
    }else if(!strcmp(keyword, "%FIRST_TIMESTAMP(100usUnits)")){
      first_time_stamp = value;      found++;
    }else if(!strcmp(keyword, "%SCALE_X(RatioTracktoMapPixels)")){
      scale_x = value;      found++;
    }
  }
  if (found!=4){printf("PLACE rd() Only found %d/4. File is corrupt\n",found); hxe();}
  for (n=0,j=0,k=0;fread(&tt,sizeof(unsigned int), 1, fp) != 0; n++) {
    fread(&x,sizeof(unsigned char), 1, fp);
    fread(&y,sizeof(unsigned char), 1, fp);
    fread(&s,sizeof(unsigned char), 1, fp);
    for (i=0; i<s; i++) {
      if(!fread(&time,sizeof(unsigned int), 1, fp))
        printf("rd_PLACE: WARNZ couldnt read spiketime %d of %d at byte %d\n",i+1,s,ftell(fp)); 
      if (k>=maxsz) xx=vector_newsize(vv,maxsz*=2);
      xx[k++]=0.1*(double)(time); // ms
      if (err>=0 && 
          (xx[k-1]<((n-1)*1e3/samps_per_sec-errtol)||xx[k-1]>n*1e3/samps_per_sec+errtol)) {
        if (VERBOSE) {
   printf("rdERRB:%d:%g (%lf:%lf) ",n,xx[k-1],((n-1)*1e3/samps_per_sec),n*1e3/samps_per_sec);
        }
        err++; 
      }
    }
    if (locflag) {
      if (j>=maxnx) {maxnx*=2;  tvec=vector_newsize(vv3,maxnx); 
        xloc=vector_newsize(vv1,maxnx); yloc=vector_newsize(vv2,maxnx);
      }
      tvec[j]=0.1*(double)tt; xloc[j]=(double)x; yloc[j]=(double)y; j++;
    }
  }
  vector_resize(vv, k);
  if (locflag) {vector_resize(vv1, j); vector_resize(vv2, j); vector_resize(vv3, j);}
  if (err>0) return (double)-err; else return (double)n*1e3/samps_per_sec;
}
ENDVERBATIM

VERBATIM

//more bpf related stuff

si1	*sscanfLine(si1 *dp, si1 *s)
{
	int	i = 0;

	*s = (si1) 0x0;
	while (sscanf(dp, "%c", s + i) != EOF){
		dp++;
		if(*(s + i) == '\n'){
			*(s + i + 1) = '\0';
			return s;
		}
		++i;
	}
	*(s + i) = '\0';
	return s;
}


#define DEBUGPR(x) printf("HERE %d\n", x); fflush(stdout);
		
si4     ReadBPFHeader (ui1 *dp, ui4 *BPFRecordSizes, si2 *Gains, si4 *AcquisitionSystem, ui4 *EEGChList, ui4 *BPFRecordTypeNumbers )
{
	ui1	VerifyRecordFormat(ui1 *, ui4 *, ui1 *, si1 *, si1 *, ui1);
	ui1	*FindHeaderStart(), *FindHeaderEnd(), *FindSection(), *FindSectionEnd();
	void	GetRecordSizes(); 
	ui1	BPFRecordTypes[256];
	ui1	*HeaderStart, *HeaderEnd, *SectionStart, *SectionEnd;
	si1	KeyWordString[256], Encoding[256], **Parameters;
	si4	Type, NumberOfParameters, NumEEGChannels;
	si4	GetRecordTypesUsed(), FindKeyWord();

	// Look for header beginning
	HeaderStart = FindHeaderStart(dp);
	if(HeaderStart == NULL){
		fprintf(MyStdErr,"\nCan't find start of header (%s).\nWrong File Format.\n", HEADER_BEGIN_LABEL);
		return(-1);
	}
	dp = (HeaderStart + strlen(HEADER_BEGIN_LABEL));
	// Check to see if it's the BPF header or not
	sscanf((si1 *)dp,"%s", KeyWordString);
	if(strstr(KeyWordString, BPF_HEADER_TYPE) == NULL){
		fprintf(MyStdErr,"This is not a BPF file. It is a %s file.\n", KeyWordString);
		return(-1);
	}
        // Section: DATABASE_INFORMATION
        dp = HeaderStart;
        SectionStart = FindSection(dp, SECTION_DATABASE_INFORMATION);
        if(SectionStart == NULL){
                fprintf(MyStdErr, "\nCan't find Header section %s.\nWrong File Format.\n", SECTION_DATABASE_INFORMATION);
                return(-1);
        }

        dp = SectionStart;

        // Acquisition System
        *AcquisitionSystem = ACX;
        if(!FindKeyWord(SectionStart, ACQUISITION_SYSTEM_LABEL, KeyWordString, Encoding)){
                printf("AcX Acquisition System\n");
        }else{
                Type = GetType(KeyWordString);
                switch(Type){
                case -1:
                        fprintf(MyStdErr,"KeyWord.integer format not respected in: %s\nWrong File Format.\n", KeyWordString);
                        return(0);

                case 0: // ASCII encoded
                        NumberOfParameters = GetASCIIEncoding(Encoding, &Parameters);
                        if(NumberOfParameters == 1){
                                if(strcasestr(Parameters[0], "Axona")){
                                        *AcquisitionSystem = AXONA;
                                        if(VERBOSE) printf("AXONA Acquisition System\n");
                                }
                        }
                        break;

                default:
                        fprintf(MyStdErr,"%s must be ASCII encoded\nWrong File Format.\n", ACQUISITION_SYSTEM_LABEL);
                return(0);
                }
        }

	// Section: SETUP_INFORMATION
	dp = HeaderStart;
	SectionStart = FindSection(dp, SECTION_SETUP_INFORMATION);
	if(SectionStart == NULL){
		fprintf(MyStdErr, "\nCan't find Header section %s.\nWrong File Format.\n", SECTION_SETUP_INFORMATION);
		return(-1);
	}

        dp = SectionStart;
        NumEEGChannels = GetRecordTypesUsed(BPFRecordTypeNumbers, dp, EEGChList);

	if(!GetGainsList(Gains, dp)) return -1;

	SectionEnd = FindSectionEnd(SectionStart);
	if(SectionEnd == NULL)
		fprintf(MyStdErr,"\nCan't find end of Header Section (%s).\nWrong File Format.\n", SECTION_SETUP_INFORMATION);

	// Section: RECORD_FORMAT
	dp = HeaderStart;
	SectionStart = FindSection(dp, SECTION_RECORD_FORMAT_INFORMATION);
	if(SectionStart == NULL){ 
		fprintf(MyStdErr, "\nCan't find section %s. Wrong data format\n", SECTION_RECORD_FORMAT_INFORMATION);
		return (-1);
	}

	// EEG RECORD ID
	if(BPFRecordTypeNumbers[EEG_BPF_REC_TYPE]){
		BPFRecordSizes[EEG_BPF_REC_TYPE] = BPFRecordTypeNumbers[EEG_BPF_REC_TYPE] * EEG_BPF_REC_DATA_SIZE + EEG_BPF_REC_ID_SIZE;
		if(!VerifyRecordFormat(BPFRecordTypes, BPFRecordSizes, SectionStart, EEG_ID, EEG_RECORD_FORMAT, EEG_BPF_REC_TYPE))
			return(-1);
	}

	// SINGLE ELECTRODE RECORD ID
	if(BPFRecordTypeNumbers[SINGLE_BPF_REC_TYPE]){
		if(!VerifyRecordFormat(BPFRecordTypes, BPFRecordSizes, SectionStart, SINGLE_ID, SINGLE_RECORD_FORMAT, SINGLE_BPF_REC_TYPE))
			return(-1);
	}


	// STEREOTRODE RECORD ID
	if(BPFRecordTypeNumbers[STEREO_BPF_REC_TYPE]){
		if(!VerifyRecordFormat(BPFRecordTypes, BPFRecordSizes, SectionStart, STEREO_ID, STEREO_RECORD_FORMAT, STEREO_BPF_REC_TYPE))
			return(-1);
	}

	// TETROTRODE RECORD ID
	if(BPFRecordTypeNumbers[TETRODE_BPF_REC_TYPE]){
		if(!VerifyRecordFormat(BPFRecordTypes, BPFRecordSizes, SectionStart, TETRODE_ID, TETRODE_RECORD_FORMAT, TETRODE_BPF_REC_TYPE))
			return(-1);
	}

	// SYNC RECORD ID
	if(BPFRecordTypeNumbers[SYNC_BPF_REC_TYPE]){
		if(!VerifyRecordFormat(BPFRecordTypes, BPFRecordSizes, SectionStart, SYNC_ID, SYNC_RECORD_FORMAT, SYNC_BPF_REC_TYPE))
			return(-1);
	}

	// ROOM POSITION RECORD ID
	if(BPFRecordTypeNumbers[ROOM_POSITION_BPF_REC_TYPE]){
		if(!VerifyRecordFormat(BPFRecordTypes, BPFRecordSizes, SectionStart, ROOM_POSITION_ID, ROOM_POSITION_RECORD_FORMAT, ROOM_POSITION_BPF_REC_TYPE))
			return(-1);
	}

	// ARENA POSITION RECORD ID
	if(BPFRecordTypeNumbers[ARENA_POSITION_BPF_REC_TYPE]){
		if(!VerifyRecordFormat(BPFRecordTypes, BPFRecordSizes, SectionStart, ARENA_POSITION_ID, ARENA_POSITION_RECORD_FORMAT, ARENA_POSITION_BPF_REC_TYPE))
			return(-1);
	}

        // KEY EVENT RECORD ID
        if(BPFRecordTypeNumbers[KEY_EVENT_BPF_REC_TYPE]){
                VerifyRecordFormat(BPFRecordTypes, BPFRecordSizes, SectionStart, KEY_EVENT_ID, KEY_EVENT_RECORD_FORMAT, KEY_EVENT_BPF_REC_TYPE);
		// This record may not be defined
        }

        // INPUT EVENT RECORD ID
        if(BPFRecordTypeNumbers[INPUT_EVENT_BPF_REC_TYPE]){
                VerifyRecordFormat(BPFRecordTypes, BPFRecordSizes, SectionStart, INPUT_EVENT_ID, INPUT_EVENT_RECORD_FORMAT, INPUT_EVENT_BPF_REC_TYPE);
		// This record may not be defined
        }

        // OUTPUT EVENT RECORD ID
        if(BPFRecordTypeNumbers[OUTPUT_EVENT_BPF_REC_TYPE]){
                VerifyRecordFormat(BPFRecordTypes, BPFRecordSizes, SectionStart, OUTPUT_EVENT_ID, OUTPUT_EVENT_RECORD_FORMAT, OUTPUT_EVENT_BPF_REC_TYPE);
		// This record may not be defined
        }

	SectionEnd = FindSectionEnd(SectionStart);
	if(SectionEnd == NULL)
		fprintf(MyStdErr,"\nCan't find end of Header Section (%s).\nWrong File Format.\n", SECTION_RECORD_FORMAT_INFORMATION);


	HeaderEnd = FindHeaderEnd(SectionEnd);
	if(HeaderEnd == NULL){
		fprintf(MyStdErr, "\nCan't find end of header. Wrong data format\n");
		return (-1);
	}

	return(NumEEGChannels);
}

ui1     *FindHeaderStart(ui1 *mdata)
{
        si1 String[4196], Line[4196], *HeaderStart;
	ui1 *dp;

	dp = mdata;
        while(sscanf((si1 *)dp, "%s", String) == 1){
                HeaderStart = (si1 *)strstr(String, HEADER_BEGIN_LABEL);
		if(HeaderStart != NULL){
                        return(dp + (HeaderStart - String));
                }       
        	sscanfLine(dp, Line);
                dp += strlen(Line);
        }       

        fprintf(MyStdErr,"Can't find BPF Header begin label: %s\n", HEADER_BEGIN_LABEL);
	return(NULL);
}       

ui1	*FindHeaderEnd(ui1 *mdata)
{
	si1 String[4196], Line[4196];
	ui1 *dp;
	
	dp = mdata;
	
	while(sscanfLine((si1 *)dp, Line) != NULL){
		dp = (ui1 *)strstr((si1 *)dp, Line);
		dp += strlen(Line);
		sscanf(Line, "%s", String);
		if(strstr(String, HEADER_END_LABEL) != NULL){
			return(dp);
		}
	}
	
	fprintf(MyStdErr,"Can't find BPF Header end label: %s\n", HEADER_END_LABEL);
	return(NULL);
}

ui1	*FindSection(ui1 *dp, ui1 *Section)
{
	si1	Flag[4196];

	while (sscanf((si1 *)dp,"%s", Flag) != EOF){
		dp = (ui1 *)strstr((si1 *)dp, Flag);
		dp += strlen(Flag);
		if(Comment((si1 *)dp, Flag)) continue;

		if(!strcmp(Flag, SECTION_BEGIN_LABEL)){
			sscanf((si1 *)dp,"%s", Flag);
			dp = (ui1 *)strstr((si1 *)dp, Flag);
			dp += strlen(Flag);

			if(!strcmp(Flag, (si1 *)Section))
				return(dp);
		}
	}
	return NULL;
}

si4	FindKeyWord(si1 *dp, si1 *KeyWord, si1 *KeyWordString, si1 *Encoding)
{
	si4	SectionEnd(), n_chars;

	while (sscanf(dp,"%s", KeyWordString) != EOF){
		dp = (si1 *)strstr(dp, KeyWordString);
		dp += strlen(KeyWordString);
		if(Comment(dp, KeyWordString)) continue;
		if(SectionEnd(KeyWordString)) break;

		if(strstr(KeyWordString, KeyWord) != NULL){
  // search until the "." for a paren; 
  // Do this because some BPF files don't have a space between the Keyword and the encoding info
  // for example ListOfGains.0(2 2 2 2 etc
			n_chars = 0;
			while(*(dp-n_chars) != '.'){
				if(*(dp-n_chars) == '('){
					dp -= n_chars;
					break;
				}
				n_chars++;
			}
	
			(void)sscanfLine(dp, Encoding);
			return(1);
		}
	}
	return 0;
}

ui1	*FindSectionEnd(ui1 *dp)
{
	si1	Line[4196], String[4196];

	while (sscanfLine(dp, Line) != EOF){
		dp = (ui1 *)strstr((si1 *)dp, Line);
		dp+= strlen(Line);
		sscanf(Line,"%s",String);
		if(SectionEnd(String)){
			return(dp);
		}
	}
	return(NULL);
}

si4	Comment(ui1 *dp, si1 *String)
{
	si1	Line[4196];

	if(!strcmp(String, COMMENT)){
		sscanfLine(dp, Line);
		dp += strlen(Line);
		return(1);
	}else	
		return(0);
}

si4	SectionEnd(si1 *String)
{
	if(!strcmp(String, SECTION_END_LABEL))
		return(1);
	return(0);
}

si4	GetType(si1 *KeyWord)
{
	si1 *TypeString;
	si4	Type = -1;
	
	TypeString = (si1 *)strchr(KeyWord, TYPE_PREFIX);
	if(TypeString == NULL)
		return(-1);
	
	sscanf(TypeString+1,"%d", &Type);	// +1 to step past the '.'
	return(Type);
}

si4	GetBinaryEncoding(si1 *Encoding, si1 *P)
{ 
	static si1 String[4196], *TerminateString, *Line;

	//  The character pair (' must start the Binary Encoding
	if((Line = (si1 *)strstr(Encoding, BINARY_STARTING_MARK)) == NULL){
		fprintf(MyStdErr,"a (' must start binary encodings.\n");
		return (0);
	}
	Line+= strlen(BINARY_STARTING_MARK);;
	sscanf(Line,"%c",P);
	if(P == NULL){
		fprintf(MyStdErr,"Can't read Binary code.\n");
		return (0);
	}
	return(1);
}

si4	GetASCIIEncoding(si1 *Encoding, si1 ***P)
{
	si4 n; 
	static si1 String[4196], *TerminateString, *Line, **Parameters;

	
	Parameters = (si1 **)calloc(MAX_NUMBER_BPF_FORMAT_PARAMETERS, sizeof(si1 *));
	if(Parameters == NULL){
		fprintf(MyStdErr, "\nCan't allocate Parameters\n");
		return (0);
	}
	*P = Parameters;

	// a '(' must start the ASCII Encoding
	if((Line = (si1 *)strstr(Encoding, ASCII_STARTING_MARK )) == NULL){
		fprintf(MyStdErr,"a '(' must start ASCII encodings. Found: %s\n", Encoding);
		return (0);
	}
	Line+= strlen(ASCII_STARTING_MARK);;

	if(sscanf(Line, "%s", String) != 1){
		fprintf(MyStdErr,"No ASCII string to decode\n");
		return (0);
	}

	Parameters[0] = (si1 *)calloc(strlen(String + 1), sizeof(si1));
	strcpy(Parameters[0], String);

	if((TerminateString = (si1 *)strstr(Parameters[0], ASCII_TERMINATING_MARK )) != NULL){	// a ')' ends the encoding
		*TerminateString = '\0';
		return(1);
	}

	n = 1;

	Line += (strlen(String) + 1);
	while((*Line == ' ') || (*Line == '\t')){	// move to start of next string
		Line ++;
	}
	while(sscanf(Line, "%s", String) != EOF){
		if(!strcmp(String, ")"))
			return (n);
		Parameters[n] = (si1 *)calloc(strlen(String) + 1, sizeof(si1));
		strcpy(Parameters[n], String);
		if((TerminateString = (si1 *)strchr(Parameters[n], ')')) != NULL){ // a ')' indicates the end of the ASCII encoding
			*TerminateString = '\0' ;
			return(++n);
		}
		Line += (strlen(String) + 1);
		while((*Line == ' ') || (*Line == '\t')){	// move to start of next string
			Line ++;
		}
		n++;

		if(n == MAX_NUMBER_BPF_FORMAT_PARAMETERS){
			fprintf(MyStdErr,"Too many Parameters in Encoding of:\t");
			return(0);
		}
	}

	return (n);
}

ui1	VerifyRecordFormat(ui1 *RecordTypes, ui4 *RecordSizes, ui1 *SectionStart, si1 *id, si1 *RecordFormat, ui1 RecType){

	ui1	*dp;
	si1	KeyWordString[4196], Encoding[4196];
	si4	Type;
	
	// Get record type
        dp = SectionStart;
        if(!FindKeyWord((si1 *)dp, id, KeyWordString, Encoding)){
                fprintf(MyStdErr,"\nCan't find keyWord (%s).\nWrong File Format.\n", id);
                return 0;
        }
        Type = GetType(KeyWordString);   
        switch(Type){   
        case -1:        
                fprintf(MyStdErr,"KeyWord.integer format not respected in: %s\nWrong File Format.\n", KeyWordString);
                return 0;
                                
        case 0: // ASCII encoded
                fprintf(MyStdErr,"%s\nWrong File Format.\n", id);
                return 0;       

        case RECORD_ID_ENCODING: // 1 byte Binary encoded
                if(!GetBinaryEncoding(Encoding, (si1 *)RecordTypes+RecType)){
                        fprintf(MyStdErr,"\nCan't get RECORD ID from header.\n");
                        return 0;
                }
                if(RecordTypes[RecType] != RecType){
                        fprintf(MyStdErr,"\nRecord type (%c) does not match line in  the header.\n", RecType);
                        return(0);
                }
                break;

        default:
                fprintf(MyStdErr,"Record id (%s) must be Binary encoded\nWrong File Format.\n", Type, KeyWordString);
                return (0);
        }

	//Get record size
        dp = SectionStart;
        if(!FindKeyWord((si1 *)dp, RecordFormat, KeyWordString, Encoding)){
                fprintf(MyStdErr,"\nCan't find keyWord (%s).\nWrong File Format.\n", RecordFormat);
                return 0;
        }
	RecordSizes[RecType] = (ui4)GetType(KeyWordString); 
	/* if (RecordSizes[RecType] != (ui4)GetType(KeyWordString)){ 
               fprintf(MyStdErr,"\nRecord size (%d) does not match line in the header (%s).\n", RecordSizes[RecType], KeyWordString);
		return(0);
	}
	*/
	return(1);
}
si4     GetRecordTypesUsed(ui4 *RecordTypeNumbers, si1 *SectionStart, ui4 *BPFEEGChannels){

        ui4     GetNumberOfChannels(), GetIfPositions();
        si1     *dp;
        ui4     i, n;
        si4     NumEEGChannels;
        si4     GetEEGChannelList();

        for(i = 0; i < MAX_BPF_RECORD_TYPES; i++)
                RecordTypeNumbers[i] = 0;

        n = 0;
        dp = SectionStart;
        n += GetNumberOfChannels(RecordTypeNumbers, dp, NUMBER_OF_EEG_CHANNELS, EEG_BPF_REC_TYPE);
        if(n){  // There are EEG channels
                NumEEGChannels = GetEEGChannelList(BPFEEGChannels, dp);
        }else
                BPFEEGChannels = NULL;
        n += GetNumberOfChannels(RecordTypeNumbers, dp, NUMBER_OF_SINGLE_CHANNELS, SINGLE_BPF_REC_TYPE);
        n += GetNumberOfChannels(RecordTypeNumbers, dp, NUMBER_OF_STEREO_CHANNELS, STEREO_BPF_REC_TYPE);
        n += GetNumberOfChannels(RecordTypeNumbers, dp, NUMBER_OF_TETRODE_CHANNELS, TETRODE_BPF_REC_TYPE);
        if(!n){
                fprintf(MyStdErr,"\nNo indication in header that there are electrophysiological data.\nWrong File Format.\n");
                return (NumEEGChannels);
        }

        (void) GetNumberOfChannels(RecordTypeNumbers, dp, NUMBER_OF_SYNC_CHANNELS, SYNC_BPF_REC_TYPE);
        (void) GetIfPositions(RecordTypeNumbers, dp, ROOM_POSITIONS, ROOM_POSITION_BPF_REC_TYPE);
        (void) GetIfPositions(RecordTypeNumbers, dp, ARENA_POSITIONS, ARENA_POSITION_BPF_REC_TYPE);

        return (NumEEGChannels);
}

si4	GetEEGChannelList(ui4 *BPFEEGChannels, si1 *SectionStart){
	ui1	*dp;
	si1	KeyWordString[4196], Encoding[4196], **Parameters;
	si4	NumberOfParameters, Type, i;

	dp = (ui1 *)SectionStart;
        if(!FindKeyWord((si1 *)dp, LIST_OF_EEG_CHANNELS, KeyWordString, Encoding)){
                fprintf(MyStdErr,"\nCan't find keyWord (%s).\n", LIST_OF_EEG_CHANNELS);
			return -1;
        }
        Type = GetType(KeyWordString);
        switch(Type){
        case -1:
                fprintf(MyStdErr,"KeyWord.integer format not respected in: %s\nWrong File Format.\n", KeyWordString);
                return -1; 

        case 0: // ASCII encoded
                NumberOfParameters = GetASCIIEncoding(Encoding, &Parameters);
		if(NumberOfParameters){
			for(i=0; i < NumberOfParameters; i++){
                                BPFEEGChannels[i] = atoi(Parameters[i]) - 1; // The header count begins at 1 but the channel count begins at 0 for processing
			}
		}
                return NumberOfParameters;

        default:
                fprintf(MyStdErr,"%s must be ASCII encoded\nWrong File Format.\n", KeyWordString);
                return 0;
        }
}	

ui4	GetNumberOfChannels(BPFRecordTypeNumbers, SectionStart, KeyWord, RecType)
	ui4 *BPFRecordTypeNumbers;
	ui1 *SectionStart;
	si1* KeyWord;
	ui1 RecType;
{
	ui1	*dp;
	si1	KeyWordString[4196], Encoding[4196], **Parameters;
	si4	Type, NumberOfParameters;

        dp = SectionStart;
        if(!FindKeyWord((si1 *)dp, KeyWord, KeyWordString, Encoding)){
                fprintf(MyStdErr,"\nCan't find keyWord (%s).\n", KeyWord);
			return(-1);
        }

        Type = GetType(KeyWordString);
        switch(Type){
        case -1:
                fprintf(MyStdErr,"KeyWord.integer format not respected in: %s\nWrong File Format.\n", KeyWordString);
                return(0);   

        case 0: // ASCII encoded
                NumberOfParameters = GetASCIIEncoding(Encoding, &Parameters);
                if(NumberOfParameters != 1){    
                        fprintf(MyStdErr,"(%s).\nWrong File Format.\n", KeyWord);
                        return 0;
                }
                sscanf(Parameters[0],"%d", &(BPFRecordTypeNumbers[RecType]));
                return(1);

        default:
                fprintf(MyStdErr,"%s must be ASCII encoded\nWrong File Format.\n", KeyWord);
                return(0);
        }
}	

ui4	GetIfPositions(BPFRecordTypeNumbers, SectionStart, KeyWord, RecType)
	ui4 *BPFRecordTypeNumbers;
	ui1 *SectionStart;
	si1* KeyWord;
	ui1 RecType;
{

	ui1	*dp;
	si1	KeyWordString[4196], Encoding[4196], **Parameters;
	si4	Type, NumberOfParameters;

        dp = SectionStart;
        if(!FindKeyWord((si1 *)dp, KeyWord, KeyWordString, Encoding)){
            // fprintf(MyStdErr,"Can't find keyWord (%s). This is not critical.\n", KeyWord);
			BPFRecordTypeNumbers[RecType] = 0;
            return 1;
        }

        Type = GetType(KeyWordString);
        switch(Type){
        case -1:
                fprintf(MyStdErr,"KeyWord.integer format not respected in: %s\nWrong File Format.\n", KeyWordString);
                return(0);   

        case 0: // ASCII encoded
                NumberOfParameters = GetASCIIEncoding(Encoding, &Parameters);
                if(NumberOfParameters != 1){    
                        fprintf(MyStdErr,"(%s).\nWrong File Format.\n", KeyWord);
                        return 0;
                }
                sscanf(Parameters[0],"%d", &(BPFRecordTypeNumbers[RecType]));
                return(1);

        default:
                fprintf(MyStdErr,"%s must be ASCII encoded\nWrong File Format.\n", KeyWord);
                return(0);
        }
}	

int	GetGainsList(si2 *Gains, ui1 *SectionStart){
	ui1	*dp;
	si1	KeyWordString[4196], Encoding[4196], **Parameters;
	si4	NumberOfParameters, Type, i;

	dp = SectionStart;
        if(!FindKeyWord((si1 *)dp, LIST_OF_GAINS, KeyWordString, Encoding)){
                fprintf(MyStdErr,"\nCan't find keyWord (%s).\n", LIST_OF_GAINS);
			return 1;
        }

        Type = GetType(KeyWordString);
        switch(Type){
        case -1:
                fprintf(MyStdErr,"KeyWord.integer format not respected in: %s\nWrong File Format.\n", KeyWordString);
                return 1; 

        case 0: // ASCII encoded
                NumberOfParameters = GetASCIIEncoding(Encoding, &Parameters);
		if(NumberOfParameters){
                  if(NumberOfParameters>MAX_NUMBER_OF_BPF_CHANNELS){
                        fprintf(MyStdErr,"GetGainsList ERRA: NumberOfParameters=%d > MAX_NUMBER_OF_BPF_CHANNELS=%d!!!\n",NumberOfParameters,MAX_NUMBER_OF_BPF_CHANNELS);
                     return 0;
                  }
			// Gains = (ui2 *)calloc(MAX_NUMBER_OF_BPF_CHANNELS, sizeof(si2));
			for(i=0; i < NumberOfParameters; i++){
                		// sscanf(Parameters[i],"%d", &(Gains[i]));
                		Gains[i] = atoi(Parameters[i]);
			}
		}
                return 1;

        default:
                fprintf(MyStdErr,"%s must be ASCII encoded\nWrong File Format.\n", KeyWordString);
	        return 1;
        }
        return 1;
}


void	bpf_instruct(){

	printf("\n\nCall with a brain potential file (.bpf) format file in the pwd.\n"); 
	printf("Outputs the EEG waveform for each channel in the BPF file to a .wfm file in the pwd.\n"); 
	printf("The output voltage is in mV units.\n"); 
	printf("Options:\n"); 
	printf("-e# give the number of samples in the EEG record. Default is %d\n", NUMBER_OF_BPF_EEG_SAMPLES); 
	printf("-g# give the factor by which to multiply the voltages so as adjust the gain of the recorded data. The default is 1.0\n");
	printf("\t\tUse this option for example if there is external gain that the recording system does not register in the BPF file header.\n");
	printf("\t\tOne example is if an external amplifier was used and set to 1000x. In this case you will need to set option g to 0.001.\n"); 
	printf("\t\tAnother example is the wireless DT recording system which has a hidden gain of 2.08\n");
}

#define	PREAMP_GAIN	100.0
#define	DSP_GAIN	20.0

int	optind;
char	*optarg;

ENDVERBATIM

VERBATIM


inline ui4 reverse_ui4(ui1* b)
{       
  ui1     *f;     
  static  ui4     u;      

  f = (ui1 *) &u; 

//        b += 3;
//
//        *f++ = *b--;
//        *f++ = *b--;
//        *f++ = *b--;
//        *f = *b;

  *f++ = *b++;
  *f++ = *b++;
  *f++ = *b++;
  *f = *b;

  return(u); 
}       

inline ui2 reverse_ui2(ui1* b)
{       
  ui1     *f;     
  static  ui2     u;      

  f = (ui1 *) &u; 

//        b += 1;
//
//        *f++ = *b--;
//        *f = *b;

  *f++ = *b++;
  *f = *b;

  return(u); 
}       
inline si2 reverse_si2(ui1* b)
{       
  ui1     *f;     
  static  si2     u;      

  f = (ui1 *) &u; 

//        b += 1;
//
//        *f++ = *b--;
//        *f = *b;

  *f++ = *b++;
  *f = *b;

  return(u); 
}       

//* probenumvec.readbpfunits("filename.bpf",outveclist,channumvec[,waveformflag])
static double readbpfunits (void* vv) {
  Object* pList; void *vw;
  int iMaxCols, iSz, npvec, ncvec, vlen[VRRY], ix, ix2, waveformflag;
  char* fname; //input bpf file name
  double *pvec, *cvec, retval, cl, pr, *vvo[VRRY];
  FILE 	*fp;
  si1	d, type, dp[MAX_BPF_REC_SIZE], *hp;
  si4	ReadBPFHeader();
  ui4	data_offset, RecordSizes[256];
  ui1     dummy[5], key;
  si1     line[256], string[256];
  ui1	probe, clust, x, y;
  ui2	ang, reverse_ui2();
  si2	reverse_si2();
  ui4	i, j, time_stamp, reverse_ui4();
  si2     *ADCValue1, *ADCValue2, Gains[MAX_NUMBER_OF_BPF_CHANNELS];
  si4	channel, BytesPerChannel, TetrodeBPFRecSize, ChannelsPerTetrode, SamplesPerWaveform;
  sf8     ddt, volts1, volts2, energy;
  si4     AcquisitionSystem;
  ui4     EEGChList[MAX_NUMBER_OF_BPF_CHANNELS];
  ui4	  BPFRecordTypeNumbers[256];
  ChannelsPerTetrode=4; SamplesPerWaveform=NUMBER_OF_BPF_TETRODE_SAMPLES;

  npvec = vector_instance_px(vv,&pvec);
  fname = gargstr(1);//input bpf filename
  pList = *hoc_objgetarg(2); //output list of vectors
  ncvec = vector_arg_px(3, &cvec);
  waveformflag=(ifarg(4)?(int)*getarg(4):1);
  iMaxCols = ivoc_list_count(pList);
  if (iMaxCols>VRRY) {printf("readbpfunits; Can only handle VRRY==%d vecs\n",VRRY); hxe();}
  retval=0.;
  if (npvec!=ncvec) {
    printf("pv.readbpfunits(file,outvecl,cv) ERRA: need probe#s in pv, clust#s in cv %d %d\n",\
           npvec,ncvec); hxe();}
  if (waveformflag) {
    if (npvec*ChannelsPerTetrode!=iMaxCols) { 
      printf("readbpfunits ERRA0: waveform outveclist should be size %d * %d (#chans), not %d\n",\
           ChannelsPerTetrode,npvec,iMaxCols); hxe();}
  } else if (npvec!=iMaxCols) if (npvec*ChannelsPerTetrode!=iMaxCols) { 
      printf("readbpfunits ERRA0: outveclist should be size %d (#chans), not %d\n",\
             npvec,iMaxCols); hxe();
  }
  //set pointers to list vectors
  for(i=0;i<iMaxCols;i++) {
    vlen[i]=0;
    if (i==0) iSz=list_vector_px3(pList,i,&vvo[i],&vw); else {
      if ((j=list_vector_px(pList,i,&vvo[i]))!=iSz) {
        fprintf(stderr,"readbpfunits ERRF: different size vectors %d %d %d\n",i,j,iSz);hxe(); } 
    }
  }
  // set these here so that they can be determined by optional values
  RecordSizes[EEG_BPF_REC_TYPE] = EEG_BPF_REC_SIZE;
  RecordSizes[SINGLE_BPF_REC_TYPE] = SINGLE_BPF_REC_SIZE;
  RecordSizes[STEREO_BPF_REC_TYPE] = STEREO_BPF_REC_SIZE;
  RecordSizes[TETRODE_BPF_REC_TYPE] = TETRODE_BPF_REC_SIZE;
  RecordSizes[SYNC_BPF_REC_TYPE] = SYNC_BPF_REC_SIZE;
  RecordSizes[ROOM_POSITION_BPF_REC_TYPE] = ROOM_POSITION_BPF_REC_SIZE;
  RecordSizes[ARENA_POSITION_BPF_REC_TYPE] = ARENA_POSITION_BPF_REC_SIZE;
  RecordSizes[KEY_EVENT_BPF_REC_TYPE] = KEY_EVENT_BPF_REC_SIZE;
  RecordSizes[INPUT_EVENT_BPF_REC_TYPE] = INPUT_EVENT_BPF_REC_SIZE;
  RecordSizes[OUTPUT_EVENT_BPF_REC_TYPE] = OUTPUT_EVENT_BPF_REC_SIZE;

  fp = fopen(fname, "r");
  if(fp == NULL){printf("readbpfunits ERRZ: can't open %s for reading!\n", fname); goto dofree; }
  while(fgets(line, 256, fp) != NULL){ sscanf(line, "%s", string);
    if(!strcmp(string, HEADER_END_LABEL)) {data_offset=ftell(fp); break;}
  }
  hp = (si1 *)calloc(data_offset, sizeof(si1));
  if(hp == NULL){fprintf(stderr,"readbpfunits ERRM: calloc failed\n"); goto dofree;}
  rewind(fp);
  if(fread(hp,sizeof(ui1), data_offset, fp) != data_offset){
    fprintf(stderr,"readbpfunits ERRN: fread header failed\n"); goto dofree; }
  fflush(stdout);
  ReadBPFHeader(hp, RecordSizes, Gains, &AcquisitionSystem, EEGChList, BPFRecordTypeNumbers);
  BytesPerChannel = 2 * SamplesPerWaveform;
  ddt = 1.0 / (sf8)SamplesPerWaveform;
  while(1) {
    fseek(fp, data_offset,0);// do because in Windows fp doesn't increment properly after fread
    if (fread(((void *)dp),sizeof(ui1), (size_t)1, fp) == EOF) break;
    type = (si1)*dp;
    data_offset += RecordSizes[type];
    if(fread((void *)(dp+1),sizeof(ui1),(size_t)RecordSizes[type]-1,fp) != (RecordSizes[type]-1)){
      fprintf(stderr,"Couldn't read complete record of type %c\n", type);
      goto resize;
    }
    time_stamp = reverse_ui4(dp+1);
    // printf("%c\t%d", type, time_stamp);
    switch (type) {
      case EEG_BPF_REC_TYPE: // ignore these
      break;
      case SYNC_BPF_REC_TYPE:
      break;
      case INPUT_EVENT_BPF_REC_TYPE:
       break;
      case OUTPUT_EVENT_BPF_REC_TYPE:
       break;

      case KEY_EVENT_BPF_REC_TYPE:
      key = (si1) *(dp + BPF_KEY_EVENT_REC_OFFSET);
      // printf("\t%c\n",key);
      break;

      case ROOM_POSITION_BPF_REC_TYPE:
      case ARENA_POSITION_BPF_REC_TYPE:
      // x = *(dp + BPF_POS_REC_X_OFFSET);
      // y = *(dp + BPF_POS_REC_Y_OFFSET);
      // ang = reverse_ui2(dp + BPF_POS_REC_ANG_OFFSET);
      // printf("\t%d\t%d\t%d",x,y,ang);
      break;

      case SINGLE_BPF_REC_TYPE:
      case STEREO_BPF_REC_TYPE:
      case TETRODE_BPF_REC_TYPE:
      probe = *(dp + BPF_RECORD_PROBE_OFFSET);
      clust = *(dp + BPF_SPK_REC_CLUST_OFFSET);
      pr=(double)probe; cl=(double)clust;
      // if vec contains probe and in contains clust then save it to that vec
      for (i=0.;i<npvec;i++) if (pvec[i]==pr&&cvec[i]==cl) break;
      if ((ix=i)<npvec) { // add to the chosen vector
        if (!waveformflag) {
          if (vlen[ix]>=iSz) {  // lengthen all vectors
            iSz=((iSz==0)?1e4:iSz*2);
            for (i=0;i<iMaxCols;i++) vvo[i]=list_vector_resize(pList,i,iSz);
          }
          vvo[ix][vlen[ix]]=(double)time_stamp/10.0;
          vlen[ix]++;
        } else {
          ix*=ChannelsPerTetrode;
          if (vlen[ix]+SamplesPerWaveform+10>=iSz) {  // lengthen all vectors
            iSz=((iSz==0)?1e4:iSz*2)+SamplesPerWaveform+10;
            for (i=0;i<iMaxCols;i++) vvo[i]=list_vector_resize(pList,i,iSz);
          }
          for (channel=0; channel<ChannelsPerTetrode; channel++) {
            ix2=ix+channel;
            for(j=0; j<SamplesPerWaveform; j++) {
              ADCValue2 = (si2 *)(dp + 7 + (2 * j) + (channel * BytesPerChannel));
              *ADCValue2 = (si2)reverse_si2(ADCValue2);
              // convert ADC value to voltage -- REMOVED
              volts2 = 10.0 * (sf8)*ADCValue2 / (sf8)SHRT_MAX / (sf8)Gains[channel];
              vvo[ix2][vlen[ix2]]=volts2;
              vlen[ix2]++;
            }
          }
        }
      }
      // printf("\t%d\t%d\t%0.4lf",probe, clust, energy);
      break;

      default:
      printf("%c UNKNOWN RECORD TYPE", type); fflush(stdout);
      break;
    }
    // printf("\n"); fflush(stdout);
  }
resize:
  for (i=0;i<iMaxCols;i++) list_vector_resize(pList,i,vlen[i]);
  retval=1.0; //success
dofree:   //  free memory
  if(VERBOSE>1) printf("readbpfunits: freeing memory\n");
  if(fp) fclose(fp);
  if(hp) free(hp);
  return retval;
}

void SetBPFRecSizes(ui4 RecordSizes[256]) {
  // set these here so that they can be determined by optional values
  RecordSizes[EEG_BPF_REC_TYPE] = EEG_BPF_REC_SIZE;
  RecordSizes[SINGLE_BPF_REC_TYPE] = SINGLE_BPF_REC_SIZE;
  RecordSizes[STEREO_BPF_REC_TYPE] = STEREO_BPF_REC_SIZE;
  RecordSizes[TETRODE_BPF_REC_TYPE] = TETRODE_BPF_REC_SIZE;
  RecordSizes[SYNC_BPF_REC_TYPE] = SYNC_BPF_REC_SIZE;
  RecordSizes[ROOM_POSITION_BPF_REC_TYPE] = ROOM_POSITION_BPF_REC_SIZE;
  RecordSizes[ARENA_POSITION_BPF_REC_TYPE] = ARENA_POSITION_BPF_REC_SIZE;
  RecordSizes[KEY_EVENT_BPF_REC_TYPE] = KEY_EVENT_BPF_REC_SIZE;
  RecordSizes[INPUT_EVENT_BPF_REC_TYPE] = INPUT_EVENT_BPF_REC_SIZE;
  RecordSizes[OUTPUT_EVENT_BPF_REC_TYPE] = OUTPUT_EVENT_BPF_REC_SIZE;
}

//* nq.fcdv.rdbpfu("filename.bpf",nq.vl)
// nq cols should be eg "timestamp", "probenum", "clustnum", "waveform" with vdec("waveform")
static double rdbpfu (void* vv) {
  Object* pList; void *vw;
  int iMaxCols, iSz, novec, cnt, vpr;
  char* fname; //input bpf file name
  double *ovec, retval, *vvo[4];
  FILE 	*fp;
  si1	d, type, dp[MAX_BPF_REC_SIZE], *hp;
  si4	ReadBPFHeader();
  ui4	data_offset, RecordSizes[256], filebytes;
  ui1     dummy[5], key;
  si1     line[256], string[256];
  ui1	probe, clust, x, y;
  ui2	ang, reverse_ui2();
  si2	reverse_si2();
  ui4	i, j, time_stamp, reverse_ui4();
  si2     *ADCValue1, *ADCValue2, Gains[MAX_NUMBER_OF_BPF_CHANNELS], *wavep;
  si4	chan, BytesPerChannel, TetrodeBPFRecSize, ChannelsPerTetrode, SamplesPerWaveform;
  sf8     ddt, volts1, volts2, energy;
  si4     AcquisitionSystem;
  ui4     EEGChList[MAX_NUMBER_OF_BPF_CHANNELS];
  ui4	  BPFRecordTypeNumbers[256];
  ChannelsPerTetrode=4; SamplesPerWaveform=NUMBER_OF_BPF_TETRODE_SAMPLES;
  novec = vector_instance_px(vv,&ovec);
  fname = gargstr(1);//input bpf filename
  pList = *hoc_objgetarg(2); //output list of vectors
  iMaxCols = ivoc_list_count(pList);
  if (iMaxCols!=4) {printf("rdbpfu: NQS veclist should be size 4\n"); hxe();}
  retval=0.;
  //set pointers to list vectors
  for(i=0;i<iMaxCols;i++) {
    if (i==0) iSz=list_vector_px3(pList,i,&vvo[i],&vw); else {
      if ((j=list_vector_px3(pList,i,&vvo[i],&vw))!=iSz) {
        fprintf(stderr,"rdbpfu ERRF: different size vectors %d %d %d\n",i,j,iSz);hxe(); } 
    }
  }  
  SetBPFRecSizes(RecordSizes);// set these here so that they can be determined by optional values
  if((fp=fopen(fname, "r"))==NULL){printf("rdbpfu ERRZ: can't open %s for reading!\n", fname); goto dofree; }
  while(fgets(line, 256, fp) != NULL){ sscanf(line, "%s", string);
    if(!strcmp(string, HEADER_END_LABEL)) {data_offset=ftell(fp); break;}
  }
  hp = (si1 *)calloc(data_offset, sizeof(si1));
  if(hp == NULL){fprintf(stderr,"rdbpfu ERRM: calloc failed\n"); goto dofree;}
  fseek(fp,0,SEEK_END); // go to end of file
  filebytes = ftell(fp); if(VERBOSE>1) printf("filebytes=%d\n",filebytes); // get the # of bytes
  rewind(fp); // go back to beginning of file
  if(fread(hp,sizeof(ui1), data_offset, fp) != data_offset){ // read header into buffer
    fprintf(stderr,"rdbpfu ERRN: fread header failed\n"); goto dofree; }
  fflush(stdout);
  ReadBPFHeader(hp, RecordSizes, Gains, &AcquisitionSystem, EEGChList, BPFRecordTypeNumbers);//parse header
  BytesPerChannel = 2 * SamplesPerWaveform;
  ddt = 1.0/(sf8)SamplesPerWaveform;
  cnt=vpr=0;
  while (1) { //read structures from file
    fseek(fp, data_offset,0);// do because in Windows fp doesn't increment properly after fread (??)
    if (fread(((void *)dp),sizeof(ui1), (size_t)1, fp) == EOF) break;
    type=(si1)*dp;
    data_offset += RecordSizes[type]; // increment pointer
    if(data_offset >= filebytes) { // check for EOF here before reading next record, which could be past end of file
      if(VERBOSE>1) printf("data_offset=%d >= filebytes=%d\n",data_offset,filebytes);
      break; 
    }
    if(fread((void *)(dp+1),sizeof(ui1),(size_t)RecordSizes[type]-1,fp)!=(RecordSizes[type]-1)){
      fprintf(stderr,"rdbpfu: Couldn't read complete record of type %c, @ offset=%d\n", type,data_offset);
      goto resize;
    }
    switch (type) {
      case INPUT_EVENT_BPF_REC_TYPE: // ignore these
      case OUTPUT_EVENT_BPF_REC_TYPE:
      case EEG_BPF_REC_TYPE: 
      case SYNC_BPF_REC_TYPE:
      case ROOM_POSITION_BPF_REC_TYPE:
      case ARENA_POSITION_BPF_REC_TYPE:
      case KEY_EVENT_BPF_REC_TYPE:
        break;
      case SINGLE_BPF_REC_TYPE:
      case STEREO_BPF_REC_TYPE:
      case TETRODE_BPF_REC_TYPE:
        time_stamp = reverse_ui4(dp+1);
        probe = *(dp + BPF_RECORD_PROBE_OFFSET);
        clust = *(dp + BPF_SPK_REC_CLUST_OFFSET);
        // if vec contains probe and in contains clust then save it to that vec
        if (cnt>=iSz) {  // lengthen all vectors
          iSz=((iSz==0)?1e4:iSz*2);
          for (i=0;i<iMaxCols;i++) vvo[i]=list_vector_resize(pList,i,iSz);
        }
        vvo[0][cnt]=(double)time_stamp/10.0;
        vvo[1][cnt]=(double)probe;
        vvo[2][cnt]=(double)clust;
        vvo[3][cnt]=(double)vpr;
        cnt++;
        if (vpr+ChannelsPerTetrode*SamplesPerWaveform>=novec) {
          novec=(novec==0?10000:10*novec);
          ovec=vector_newsize(vv,novec);
        }
        wavep = (si2*) (dp + 7);//temporary si2 pointer to waveform
        for (chan=0;chan<ChannelsPerTetrode;chan++) for(j=0;j<SamplesPerWaveform;j++,wavep++) {
          ovec[vpr++] = (double) 10.0 * wavep[0] / SHRT_MAX; 
        }
        // printf("\t%d\t%d\t%0.4lf",probe, clust, energy);
        break;
      default:
        printf("%c UNKNOWN RECORD TYPE", type); fflush(stdout);
        break;
    }
  }
resize:
  if(VERBOSE>1) printf("vpr = %d, novec = %d\n",vpr,novec);
  for (i=0;i<iMaxCols;i++) list_vector_resize(pList,i,cnt);
  vector_resize(vv,vpr);
  retval=1.0; //success
dofree:   //  free memory
  if(VERBOSE>1) printf("rdbpfu: freeing memory\n");
  if(fp) fclose(fp);
  if(hp) free(hp);
  return retval;
}

// vector.readbpfeeg("filename.bpf",list_of_output_vectors [,records_per_sample,gain_adjust])
static double readbpfeeg(void* vv) {

  Object* pList = 0;
  int iMaxCols = 0;   

  double** vvo = 0;

  double retval = 0.0;

  FILE 	*fp = 0;
  si1	file[4196], *c;
  si1	type, dp[MAX_BPF_REC_SIZE], *hp = 0;
  si4	d, channel, EEGBPFRecSize, EEGChannels = NUMBER_OF_BPF_EEG_CHANNELS, SamplesPerRecord = NUMBER_OF_BPF_EEG_SAMPLES;
  si4	BytesPerChannel;
  ui4	data_offset, RecordSizes[256];
  si1	line[4096], string[4096];
  ui1	probe, clust, x, y;
  ui2	ang;
  ui4	i, time_stamp;
  si2	*ADCValue, Gains[MAX_NUMBER_OF_BPF_CHANNELS];
  sf4	volts, AdjustGain = 1.0;
  si4	AcquisitionSystem;
  ui4     EEGChList[MAX_NUMBER_OF_BPF_CHANNELS];
  ui4	BPFRecordTypeNumbers[256];
  char* fname; //input bpf file name
  double* pThisV=0; //stores eeg channel #s
  int iThisVSz=0 , iSz=0, iOffsetStart=0,iOffsetCur=0;
  iThisVSz = vector_instance_px(vv,&pThisV);
 
  // set these here so that they can be determined by optional values
  RecordSizes[EEG_BPF_REC_TYPE] = EEG_BPF_REC_SIZE;
  RecordSizes[SINGLE_BPF_REC_TYPE] = SINGLE_BPF_REC_SIZE;
  RecordSizes[STEREO_BPF_REC_TYPE] = STEREO_BPF_REC_SIZE;
  RecordSizes[TETRODE_BPF_REC_TYPE] = TETRODE_BPF_REC_SIZE;
  RecordSizes[SYNC_BPF_REC_TYPE] = SYNC_BPF_REC_SIZE;
  RecordSizes[ROOM_POSITION_BPF_REC_TYPE] = ROOM_POSITION_BPF_REC_SIZE;
  RecordSizes[ARENA_POSITION_BPF_REC_TYPE] = ARENA_POSITION_BPF_REC_SIZE;
  RecordSizes[INPUT_EVENT_BPF_REC_TYPE] = INPUT_EVENT_BPF_REC_SIZE;
  RecordSizes[OUTPUT_EVENT_BPF_REC_TYPE] = OUTPUT_EVENT_BPF_REC_SIZE;

  fname = gargstr(1);//input bpf filename
  pList = *hoc_objgetarg(2); //output list of vectors
  iMaxCols = ivoc_list_count(pList); //output list count (should be >= EEGChannels+1)

  if(ifarg(3)) SamplesPerRecord = *getarg(3); //optional argument for SamplesPerRecord
  if(ifarg(4)) AdjustGain = (sf4) *getarg(4); //optional argument to adjust gain

  fp = fopen(fname, "r");
  if(fp == NULL){
    printf("readbpfeeg ERRZ: can't open %s for reading!\n", fname);
    goto dofree;
  }	

  while(get_line(fp, line, 4096) != NULL){
    sscanf(line, "%s", string);
    if(!strcmp(string, HEADER_END_LABEL)){
      data_offset = ftell(fp);
      break;
    }
  }

  hp = (si1 *)calloc(data_offset, sizeof(si1));
  if(hp == NULL){
    fprintf(stderr,"readbpfeeg ERRM: calloc failed\n");
    goto dofree;
  }
  rewind(fp);
  if(fread(hp,sizeof(ui1), data_offset, fp) != data_offset){
    fprintf(stderr,"readbpfeeg ERRN: fread header failed\n");
    goto dofree;
  }

  fflush(stdout);

  //something wrong here
  EEGChannels = ReadBPFHeader(hp, RecordSizes, Gains, &AcquisitionSystem, EEGChList, BPFRecordTypeNumbers);

  if(EEGChannels < 1){
    printf("readbpfeeg ERRB: No EEG Channels.\n");
    goto dofree;
  } else if(VERBOSE) printf("readbpfeeg: EEGChannels=%d\n",EEGChannels);

  EEGBPFRecSize =  EEG_BPF_REC_INFO_SIZE + (SamplesPerRecord * 2 * EEGChannels);
		
  if(iMaxCols<EEGChannels+1){
    fprintf(stderr,"readbpfeeg ERRD: need EEGChannels+1=%d for eegdata+time, List.count=%d!!\n",EEGChannels+1,iMaxCols);
    goto dofree;
  } else if(VERBOSE) printf("readbpfeeg: iMaxCols=%d\n",iMaxCols);

  if(iThisVSz<EEGChannels){
    printf("readbpfeeg ERRG: This vector sz %d < EEGChannels %d!\n",iThisVSz,EEGChannels);
    goto dofree;
  } else if(VERBOSE) printf("readbpfeeg: iThisVSz=%d\n",iThisVSz);

  vvo = (double**)malloc(iMaxCols*sizeof(double*)); //alloc mem for pointers to list vectors
  if(!vvo){
    fprintf(stderr,"readbpfeeg ERRE: couldn't alloc mem for vvo!\n");
    goto dofree;
  }
 
  //set pointers to list vectors
  for(i=0;i<EEGChannels;i++){
    pThisV[i]=1+EEGChList[i]; //store channel #s
    if((iSz=list_vector_px(pList,i,&vvo[i]))<SamplesPerRecord){
      fprintf(stderr,"readbpfeeg ERRF: list vectors size %d < BPF SamplesPerRecord %d!\n",iSz,SamplesPerRecord);
      goto dofree;
    } else if(VERBOSE>1) printf("readbpfeeg: iSz %d = %d\n",i,iSz);
  } 
  if((iSz=list_vector_px(pList,i,&vvo[i]))<SamplesPerRecord)
  {  fprintf(stderr,"readbpfeeg ERRT: time Vector invalid!\n");
     goto dofree;
  }

  BytesPerChannel = 2 * SamplesPerRecord;
  if(VERBOSE>1) printf("readbpfeeg: BytesPerChannel=%d SamplesPerRecord=%d\n",BytesPerChannel,SamplesPerRecord);
  i = 0;

  clock_t startt,endt;
  startt=clock();
  if(AcquisitionSystem == ACX)
  { while(1)
    {  if(fread(((void *)dp),sizeof(ui1), (size_t)1, fp) == EOF)
         break;

       type = (si1)*dp;

       if(fread((void *)(dp+1),sizeof(ui1),(size_t)RecordSizes[type]-1,fp) != (RecordSizes[type]-1))
         break;
       
       if(type!=EEG_BPF_REC_TYPE) continue;

       for(channel=0; channel < EEGChannels; channel++)
       { iOffsetCur=iOffsetStart;
         for(i=0; i < SamplesPerRecord; i++, iOffsetCur++)
         {  ADCValue = (si2 *)(dp + 5 + (2 * i) + (channel * BytesPerChannel));
     	   *ADCValue = (si2)reverse_si2(ADCValue);
      	   // convert ADC value to voltage in mV
            vvo[channel][iOffsetCur]=AdjustGain*10.*(sf4)*ADCValue/(sf4)SHRT_MAX/(sf4)Gains[EEGChList[channel]]*1e3; 
         }
       }
       iOffsetStart += SamplesPerRecord; //increment offset
    }
  }
  else if(AcquisitionSystem == AXONA)
  { while(1)
    {  if(fread(((void *)dp),sizeof(ui1), (size_t)1, fp) == EOF)
         break;

       type = (si1)*dp;

       if(fread((void *)(dp+1),sizeof(ui1),(size_t)RecordSizes[type]-1,fp) != (RecordSizes[type]-1))
         break;

       if(type!=EEG_BPF_REC_TYPE) continue;

       for(channel=0; channel < EEGChannels; channel++)
       { iOffsetCur=iOffsetStart;
         for(i=0; i < SamplesPerRecord; i++,iOffsetCur++)
         {  ADCValue = (si2 *)(dp + 5 + (2 * i) + (channel * BytesPerChannel));
   	   *ADCValue = (si2)reverse_si2(ADCValue);
  	    // convert ADC value to voltage in mV
            vvo[channel][iOffsetCur]=1000.0*(AdjustGain*(sf4)*ADCValue/(sf4)SHRT_MAX);
         }
       }
       iOffsetStart += SamplesPerRecord; //increment offset
    }
  }
  int iT = EEGChannels, iLim = iOffsetStart; //create time index in minutes
  double tcur = 0.0, tinc = (1/2e3)/60.0;    //assumes sampling rate == 2Khz, which is default
  for(i=0;i<iLim;i++)
  { vvo[iT][i] = tcur;
    tcur += tinc;
  }
  endt=clock();
  printf("main loop time = %gms\n", 1000.0*((double)(endt-startt))/CLOCKS_PER_SEC);
  if(VERBOSE) printf("offset=%d\n",iOffsetStart);
  retval=1.0;//success
//free memory
dofree: 
  return retval; // ???????????????????? why return before free ?????
  if(VERBOSE>1) printf("readbpfeeg: freeing memory\n");
  if(fp) fclose(fp);
  if(vvo) free(vvo);
  if(hp) free(hp);
  return retval;
}

ENDVERBATIM

: bpfeegchlab(filepath,numchannels,str1,str2,...)
: return # of channel labels found (includes empty labels which are usually for SYNC channels)
FUNCTION bpfeegchlab () {
  VERBATIM
  double dRet;
  dRet=-1.0;
  FILE* fp = 0;
  char* fname = 0, **plabs = 0, line[4096],string[4096],*pch;
  fname = gargstr(1); //input bpf filename
  fp = fopen(fname,"r");
  if(!fp){
    printf("bpfeegchlab ERRA: couldn't open %s for reading!\n",fname);
    goto CHLCLEANUP;
  }
  int iChans = (int)*getarg(2), idx , jdx, iLine = 0, iMaxLine = 100 , iLabLen = strlen(CHANNEL_LABEL), iCH = 0;
  if(iChans < 1){
    printf("bpfeegchlab ERRB: invalid number of channels: %d!\n",iChans);
    goto CHLCLEANUP;
  }
  plabs = (char**)calloc(iChans,sizeof(char*));
  if(!plabs){
    printf("bpfeegchlab ERRC: out of memory!\n");
    goto CHLCLEANUP;
  }
  for(idx=0;idx<iChans;idx++){
    plabs[idx] = gargstr(idx+3);
    if(!plabs[idx]){
      printf("bpfeegchlab ERRD: couldn't get input string %d!\n",idx+3);
      goto CHLCLEANUP;
    }
    plabs[idx][0]=0; //set to empty string
  }
  while(iLine++<iMaxLine && iCH<iChans && get_line(fp,line,4096)!=NULL){
    if(!strncmp(line,CHANNEL_LABEL,iLabLen)){
      if(VERBOSE) printf("bpfeegchlab: %s\n",line);
      pch=strchr(line,'\"');
      if(!pch) continue;
      pch++;
      if(!*pch) continue;
      idx=0;
      while(*pch && *pch!='\"') plabs[iCH][idx++]=*pch++;
      plabs[iCH][idx]=0; //terminate string with 0
      if(VERBOSE) printf("bpfeegchlab: label %d = %s\n",iCH,plabs[iCH]);
      iCH++;
    } else if(VERBOSE>1) printf("bpfeegchlab: L%d=%s\n",iLine,line);
  }
  dRet=(double)iCH; //return # of channel labels found
CHLCLEANUP:
  if(fp) fclose(fp);
  if(plabs) free(plabs);
  return dRet;
  ENDVERBATIM
}

: get # of SYNC channels in a BPF file
FUNCTION bpfsyncchan () {
  VERBATIM
  FILE 	*fp = 0;
  si1	file[4196], *c;
  si1	type, dp[MAX_BPF_REC_SIZE], *hp = 0;
  si4	d, channel, EEGBPFRecSize, EEGChannels = NUMBER_OF_BPF_EEG_CHANNELS;
  ui4	data_offset, RecordSizes[256];
  si1	line[4096], string[4096];
  ui4	i;
  si2	Gains[MAX_NUMBER_OF_BPF_CHANNELS];
  sf4	AdjustGain = 1.0;
  si4	AcquisitionSystem;
  ui4   EEGChList[MAX_NUMBER_OF_BPF_CHANNELS];
  ui4	BPFRecordTypeNumbers[256];
  char* fname; //input bpf file name
 
  // set these here so that they can be determined by optional values
  RecordSizes[EEG_BPF_REC_TYPE] = EEG_BPF_REC_SIZE;
  RecordSizes[SINGLE_BPF_REC_TYPE] = SINGLE_BPF_REC_SIZE;
  RecordSizes[STEREO_BPF_REC_TYPE] = STEREO_BPF_REC_SIZE;
  RecordSizes[TETRODE_BPF_REC_TYPE] = TETRODE_BPF_REC_SIZE;
  RecordSizes[SYNC_BPF_REC_TYPE] = SYNC_BPF_REC_SIZE;
  RecordSizes[ROOM_POSITION_BPF_REC_TYPE] = ROOM_POSITION_BPF_REC_SIZE;
  RecordSizes[ARENA_POSITION_BPF_REC_TYPE] = ARENA_POSITION_BPF_REC_SIZE;
  RecordSizes[INPUT_EVENT_BPF_REC_TYPE] = INPUT_EVENT_BPF_REC_SIZE;
  RecordSizes[OUTPUT_EVENT_BPF_REC_TYPE] = OUTPUT_EVENT_BPF_REC_SIZE;

  fname = gargstr(1);//input bpf filename

  fp = fopen(fname, "r");
  if(fp == NULL){
    printf("bpfsyncchan ERRA: Can't open %s for reading\n", fname);
    goto dofree;
  }	

  while(get_line(fp, line, 4096) != NULL){
    sscanf(line, "%s", string);
    if(!strcmp(string, HEADER_END_LABEL)){
      data_offset = ftell(fp);
      break;
    }
  }

  hp = (si1 *)calloc(data_offset, sizeof(si1));
  if(hp == NULL){
    fprintf(stderr,"bpfsyncchan ERRB: calloc failed\n");
    goto dofree;
  }
  rewind(fp);
  if(fread(hp,sizeof(ui1), data_offset, fp) != data_offset){
    fprintf(stderr,"bpfsyncchan ERRC: fread header failed\n");
    goto dofree;
  }

  fflush(stdout);

  EEGChannels = ReadBPFHeader(hp, RecordSizes, Gains, &AcquisitionSystem, EEGChList, BPFRecordTypeNumbers); //something wrong here

 if(VERBOSE) printf("bpfsyncchan : SYNCChannels=%d\n",BPFRecordTypeNumbers[SYNC_BPF_REC_TYPE]);
dofree:
  if(VERBOSE>1) printf("bpfsyncchan: freeing memory\n");
  if(fp) fclose(fp);
  if(hp) free(hp);
  return (double) (BPFRecordTypeNumbers[SYNC_BPF_REC_TYPE] >= 0 ? BPFRecordTypeNumbers[SYNC_BPF_REC_TYPE] : 0);
  ENDVERBATIM
}

: get # of eeg channels in a BPF file
FUNCTION bpfeegchan(){
  VERBATIM
  FILE 	*fp = 0;
  si1	file[4196], *c;
  si1	type, dp[MAX_BPF_REC_SIZE], *hp = 0;
  si4	d, channel, EEGBPFRecSize, EEGChannels = NUMBER_OF_BPF_EEG_CHANNELS;
  ui4	data_offset, RecordSizes[256];
  si1	line[4096], string[4096];
  ui4	i;
  si2	Gains[MAX_NUMBER_OF_BPF_CHANNELS];
  sf4	AdjustGain = 1.0;
  si4	AcquisitionSystem;
  ui4   EEGChList[MAX_NUMBER_OF_BPF_CHANNELS];
  ui4	BPFRecordTypeNumbers[256];
  char* fname; //input bpf file name
 
  // set these here so that they can be determined by optional values
  RecordSizes[EEG_BPF_REC_TYPE] = EEG_BPF_REC_SIZE;
  RecordSizes[SINGLE_BPF_REC_TYPE] = SINGLE_BPF_REC_SIZE;
  RecordSizes[STEREO_BPF_REC_TYPE] = STEREO_BPF_REC_SIZE;
  RecordSizes[TETRODE_BPF_REC_TYPE] = TETRODE_BPF_REC_SIZE;
  RecordSizes[SYNC_BPF_REC_TYPE] = SYNC_BPF_REC_SIZE;
  RecordSizes[ROOM_POSITION_BPF_REC_TYPE] = ROOM_POSITION_BPF_REC_SIZE;
  RecordSizes[ARENA_POSITION_BPF_REC_TYPE] = ARENA_POSITION_BPF_REC_SIZE;
  RecordSizes[INPUT_EVENT_BPF_REC_TYPE] = INPUT_EVENT_BPF_REC_SIZE;
  RecordSizes[OUTPUT_EVENT_BPF_REC_TYPE] = OUTPUT_EVENT_BPF_REC_SIZE;

  fname = gargstr(1);//input bpf filename

  fp = fopen(fname, "r");
  if(fp == NULL){
    printf("bpfeegchan ERRA: Can't open %s for reading\n", fname);
    goto dofree;
  }	

  while(get_line(fp, line, 4096) != NULL){
    sscanf(line, "%s", string);
    if(!strcmp(string, HEADER_END_LABEL)){
      data_offset = ftell(fp);
      break;
    }
  }

  hp = (si1 *)calloc(data_offset, sizeof(si1));
  if(hp == NULL){
    fprintf(stderr,"bpfeegchan ERRB: calloc failed\n");
    goto dofree;
  }
  rewind(fp);
  if(fread(hp,sizeof(ui1), data_offset, fp) != data_offset){
    fprintf(stderr,"bpfeegchan ERRC: fread header failed\n");
    goto dofree;
  }

  fflush(stdout);

  EEGChannels = ReadBPFHeader(hp, RecordSizes, Gains, &AcquisitionSystem, EEGChList, BPFRecordTypeNumbers); //something wrong here

  if(EEGChannels < 1){
    printf("bpfeegchan ERRD: No EEG Channels.\n");
    goto dofree;
  } if(VERBOSE) printf("bpfeegchan : EEGChannels=%d\n",EEGChannels);
dofree:
  if(VERBOSE>1) printf("bpfeegchan: freeing memory\n");
  if(fp) fclose(fp);
  if(hp) free(hp);
  return (double)EEGChannels;
  ENDVERBATIM
}

: get total # of samples for a single EEG channel for the full duration of recording of a BPF file
FUNCTION bpfeegsamp(){
  VERBATIM
  FILE 	*fp = 0;
  si1	file[4196], *c;
  si1	type, dp[MAX_BPF_REC_SIZE], *hp = 0;
  si4	d, channel, EEGBPFRecSize, EEGChannels = NUMBER_OF_BPF_EEG_CHANNELS;
  //samples per channel for single unit of time - if 2000Hz recording, then SamplesPerRecord = 2000
  si4   SamplesPerRecord = NUMBER_OF_BPF_EEG_SAMPLES;
  si4	BytesPerChannel;
  ui4	data_offset, RecordSizes[256];
  si1	line[4096], string[4096];
  ui1	probe, clust, x, y;
  ui2	ang;
  ui4	i;
  si2	*ADCValue, Gains[MAX_NUMBER_OF_BPF_CHANNELS];
  sf4	volts, **EEGdata = 0, AdjustGain = 1.0;
  si4	AcquisitionSystem;
  ui4   EEGChList[MAX_NUMBER_OF_BPF_CHANNELS];
  ui4	BPFRecordTypeNumbers[256];
  char* fname; //input bpf file name
  int iTotalSamplesPerChannel=0; //samples per eeg channel for full BPF recording
 
  // set these here so that they can be determined by optional values
  RecordSizes[EEG_BPF_REC_TYPE] = EEG_BPF_REC_SIZE;
  RecordSizes[SINGLE_BPF_REC_TYPE] = SINGLE_BPF_REC_SIZE;
  RecordSizes[STEREO_BPF_REC_TYPE] = STEREO_BPF_REC_SIZE;
  RecordSizes[TETRODE_BPF_REC_TYPE] = TETRODE_BPF_REC_SIZE;
  RecordSizes[SYNC_BPF_REC_TYPE] = SYNC_BPF_REC_SIZE;
  RecordSizes[ROOM_POSITION_BPF_REC_TYPE] = ROOM_POSITION_BPF_REC_SIZE;
  RecordSizes[ARENA_POSITION_BPF_REC_TYPE] = ARENA_POSITION_BPF_REC_SIZE;
  RecordSizes[INPUT_EVENT_BPF_REC_TYPE] = INPUT_EVENT_BPF_REC_SIZE;
  RecordSizes[OUTPUT_EVENT_BPF_REC_TYPE] = OUTPUT_EVENT_BPF_REC_SIZE;

  fname = gargstr(1);//input bpf filename
  if(ifarg(2)) SamplesPerRecord = *getarg(2); //optional argument for SamplesPerRecord

  fp = fopen(fname, "r");
  if(fp == NULL)
  { printf("bpfeegsamp ERRA: Can't open %s for reading\n", fname);
    goto dofree;
  }	

  while(get_line(fp, line, 4096) != NULL)
  { sscanf(line, "%s", string);
    if(!strcmp(string, HEADER_END_LABEL))
    { data_offset = ftell(fp);
      break;
    }
  }

  hp = (si1 *)calloc(data_offset, sizeof(si1));
  if(hp == NULL)
  { fprintf(stderr,"bpfeegsamp ERRB: calloc failed\n");
    goto dofree;
  }
  rewind(fp);
  if(fread(hp,sizeof(ui1), data_offset, fp) != data_offset)
  { fprintf(stderr,"bpfeegsamp ERRC: fread header failed\n");
    goto dofree;
  }

  fflush(stdout);

  EEGChannels = ReadBPFHeader(hp, RecordSizes, Gains, &AcquisitionSystem, EEGChList, BPFRecordTypeNumbers); //something wrong here

  if(EEGChannels < 1)
  { printf("bpfeegsamp ERRD: No EEG Channels.\n");
    goto dofree;
  } if(VERBOSE) printf("bpfeegsamp: EEGChannels=%d\n",EEGChannels);

  EEGBPFRecSize =  EEG_BPF_REC_INFO_SIZE + (SamplesPerRecord * 2 * EEGChannels);
		 
  BytesPerChannel = 2 * SamplesPerRecord;
  i = 0;

  while(1)
  {
    if(fread(((void *)dp),sizeof(ui1), (size_t)1, fp) == EOF)
      break;

    type = (si1)*dp;

     data_offset += RecordSizes[type];
     if(fread((void *)(dp+1),sizeof(ui1),(size_t)RecordSizes[type]-1,fp) != (RecordSizes[type]-1))
       break;

     if(type == EEG_BPF_REC_TYPE)
       iTotalSamplesPerChannel+=SamplesPerRecord;
  }
  if(VERBOSE) printf("bpfeegsamp: TotalSamplesPerChannel=%d\n",iTotalSamplesPerChannel);
//free memory
dofree: 
  if(VERBOSE>1) printf("bpfeegsamp: freeing memory\n");
  if(fp) fclose(fp);
  if(hp) free(hp);
  return (double)iTotalSamplesPerChannel;
  ENDVERBATIM
}

VERBATIM
// * readbdat - read the extra/intra-cellular .dat files provided by buzsaki
// vector.readbdat("filename.dat",list_of_output_vectors,numchannels)
static double readbdat (void* vv) {
  FILE* fp;
  ListVec* pList;
  double** pLV, retval;
  int i,j,iChans,*pLen,iMinSz;
  short* pbuf;
  char* fname;
  retval=0.0; fp=0; pList=0; pbuf=0;
  fname=gargstr(1);
  fp=fopen(fname,"rb");
  if(!fp){
    printf("readbdat ERRA: couldn't open %s for reading!\n",fname);
    goto CLEANBD;
  }
  pList=AllocListVec(*hoc_objgetarg(2));
  if(!pList){
    printf("readbdat ERRB: couldn't initialize list vec arg 2!\n");
    goto CLEANBD;
  }
  pLV = pList->pv;  pLen = pList->plen;
  iChans = (int)*getarg(3);
  if(iChans < 1) {
    printf("readbdat ERRC: num channels must be > 0!\n",iChans);
    goto CLEANBD;
  }
  if(pList->isz < iChans) {
    printf("readbdat ERRD: ListVec arg 2 size %d < # of channels %d!\n",pList->isz,iChans);
    goto CLEANBD;
  }
  j=0;
  pbuf = (short*)malloc(sizeof(short)*iChans);
  iMinSz=pLen[0];
  for(i=0;i<iChans;i++) if(pLen[i]<iMinSz) iMinSz=pLen[i];
  while(1) {
    if(fread(pbuf,sizeof(short), iChans, fp)!=iChans) break;
    if(j < iMinSz) for(i=0;i<iChans;i++) pLV[i][j]=pbuf[i];
    j++;
  }
  if(j > iMinSz) printf("readbdat WARNA: could only read %d samples, file has %d samples\n",iMinSz,j);
  retval=j;
CLEANBD:
  if(pList) FreeListVec(&pList); if(fp) fclose(fp); if(pbuf) free(pbuf);
  return retval;
}
ENDVERBATIM

: get the probe number from a TS file-name -- convention is filename ends in
: P%dC%d , so the # after P is the probe number. probe,cluster are digits
FUNCTION GetProbeNum () {
  VERBATIM
  char* fname,buf[10]; int i,j,sz,ploc;
  fname=gargstr(1);
  sz=strlen(fname);
  for(i=sz-1;i>=0;i--) if(fname[i]=='P') break;
  if(i-1<=0 || fname[i-1] != '.') return -1;
  ploc=i; buf[0]=j=0;
  for(i=ploc+1;i<sz && fname[i]>='0' && fname[i]<='9' && j<10;i++) buf[j++]=fname[i];
  if(i>=sz || fname[i]!='C') return -1;
  if(j<10) buf[j]=0; else buf[9]=0;
  return atoi(buf);
  ENDVERBATIM
}

: get the cluster number from a TS file-name -- convention is filename ends in
: P%dC%d , so the # after C is the cluster number. probe,cluster are digits
FUNCTION GetClustNum () {
  VERBATIM
  char* fname,buf[10]; int i,j,sz,cloc;
  fname=gargstr(1);
  sz=strlen(fname);
  for(i=sz-2;i>=0;i--) if(fname[i]=='C') break;
  cloc=i; buf[0]=j=0;
  if(i+1>=sz) return -1;
  for(i=cloc+1;i<sz && fname[i]>='0' && fname[i]<='9' && j<10;i++) buf[j++]=fname[i];
  if(j<10) buf[j]=0; else buf[9]=0;
  return atoi(buf);
  ENDVERBATIM
}


:* PROCEDURE install()
PROCEDURE install () {
  if (INSTALLED==1) {
    printf("$Id: place.mod,v 1.141 2011/07/06 15:24:40 samn Exp $\n")
  } else {
  INSTALLED=1
VERBATIM
  install_vector_method("simts_PLACE", simts); // need to give _PLACE explicitly
  install_vector_method("rd_PLACE", rd);
  install_vector_method("rdtts_PLACE", rdtts);
  install_vector_method("kendall",kendal2);
  install_vector_method("kend2",kend2);
  install_vector_method("kendall_PLACE", kendall);
  install_vector_method("mkgaussfield_PLACE", mkgaussfield);
  install_vector_method("dumpratemap_PLACE", dumpratemap);
  install_vector_method("glob", glob); // not particularly placeish
  install_vector_method("cumul", cumul); // not particularly placeish
  install_vector_method("readbpfunits",readbpfunits);
  install_vector_method("rdbpfu",rdbpfu);
  install_vector_method("readbpfeeg",readbpfeeg);
  install_vector_method("readbdat",readbdat);
ENDVERBATIM
  }
}
