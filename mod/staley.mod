: $Id: staley.mod,v 1.75 2010/04/30 18:59:54 samn Exp $ 

NEURON {
  SUFFIX staley
  GLOBAL installed,verbose,samprate,mindist,minspikes,abovebth
  GLOBAL Lintdur,Bintdur,spkup
}

::* NEURON stuff 
PARAMETER {
  installed=0
  verbose=0
  samprate=2000    : sampling rate of original data, assumed to be 2KHz as in BPF default
  Lintdur=120      : little interval: duration in ms
  Bintdur=12       : big interval: duration in sec
  Sintdur=2        : smaller interval for checking spike counts for beg and end of sz
  mindist=0.1      : (sec) minimum duration between seizures (for concatenating them)
  minspikes=40     : minimum # of spikes in a 'seizure'
  abovebth=0       : not used currently
  endwting=1       : to weight avgnumspikes for criterion for end of sz
  flag=0           : choose set of criteria
  spkup=12         : duration of upswing to look for
  spklim=0.25      : used in upswing detection - larger means bigger spikes required
  useavgtots=0     : use average max_maxS-min_minS over all intervals for upvalue threshold

  usesharp=0       : use 2nd deriv to test spike
  sharpoff=4       : iff usesharp==1, uses x{tt+sharpoff} - 2*x{tt} + x{tt-sharpoff}
  sharpth=-500     : threshold for sharpness to be considered a spike - only used when usesharp=1

  incby1=0         : inc by 1 in inteval checks
}

VERBATIM

static double gsz(), gszspk();

#include "misc.h"
static  ListVec* pL;

typedef struct SEIZURE {
  int startIntervalIndex;
  int endIntervalIndex;
  int totalSpikesCount;
  int endSpikesCount;
  int ID;
  int startIndex;
  int endIndex;
} sSeizure;


//* utility functions and main struct
//note: not sure initialization of sSeizure is done correctly, need to see original class
// constructor
sSeizure* AllocSeizure (int intervalIndex,int spikesCount)
{
  sSeizure* p;
  if(!(p = calloc(1,sizeof(sSeizure)))){
    printf("AllocSeizure ERR: out of memory!\n"); hxe();
  }
  p->startIntervalIndex = p->endIntervalIndex = intervalIndex;
  p->endIntervalIndex+=1; // when looking at pairs
  p->totalSpikesCount = spikesCount;
  p->ID = -1; // invalid identifier
  return p;
}

typedef struct LSEIZURE {
  int bufsz;
  int count;
  sSeizure** pp;
} LSeizure;

int InitLSeizure (LSeizure* pl,int sz) {
  pl->bufsz = sz;
  pl->count = 0;
  if(sz==0) pl->pp=NULL; else pl->pp=calloc(sz,sizeof(sSeizure*));
  return pl->bufsz;
}

void FreeLSeizure(LSeizure* pl) {
  int i;
  //  for(i=0;i<pl->count;i++) free(pl->pp[i]);
  free(pl->pp);
}

int AddSeizure(LSeizure* pl, sSeizure* ps) {
  if(0) printf("pl=%p , pl->count=%d, pl->bufsz=%d\n",pl,pl->count,pl->bufsz);
  if( pl->count + 1 >= pl->bufsz ) {
    pl->bufsz *= 4;
    if(0)printf("pl=%p, realloc pl->count=%d, pl->bufsz=%d\n",pl,pl->count,pl->bufsz);
    if(! (pl->pp=realloc(pl->pp,sizeof(sSeizure*)*pl->bufsz)) ) {
      printf("AddSeizure: out of memory!\n"); hxe();
    }
  }
  pl->pp[pl->count++] = ps;
  return 1;
}

int RemoveSeizureAt(LSeizure* pl, int idx) {
  int i,j;
  if( idx < 0 || idx >= pl->count) {
    printf("RemoveSeizureAt: invalid index=%d, count=%d!\n",idx,pl->count); hxe();
  }
  for(i=idx+1;i<pl->count;i++) pl->pp[i-1]=pl->pp[i];
  pl->count--;
  return 1;
}

void printSeizure (sSeizure* p) {
  printf("ID:%d, startIntervalIndex:%d, endIntervalIndex:%d, totalSpikes:%d, endSpikes:%d, startIndex:%d, endIndex:%d\n",
	 p->ID,p->startIntervalIndex,p->endIntervalIndex,p->totalSpikesCount,p->endSpikesCount,p->startIndex,p->endIndex);
}

void printSeizures (LSeizure* lp) {
  int i;
  for(i=0;i<lp->count;i++) printSeizure(lp->pp[i]);
}

//* dgetseizures(double* channelData,int channelsLength) -- main routine
static LSeizure* dgetseizures (double* channelData,int channelsLength) {
  // based on Andy White & Kevin Staley routine for seizure detection
  // currently works for recording freq of 250Hz
  // This subroutine creates a measure of the density of the signal
  // We first find the maxima and minima for groups of 25 points
  //** declarations and allocations
  int intervalIndex, dataIndex, seizureIndex, iS, stopIndex, uplim;
  int intervalLength, channelIndex, i,j, intervalsCount, SgroupCount, true;
  int seizuresCount, spikesCount, totsLen, bufszStart, dbxi, ii, jj, LintCnt;
  double value, min_minS, max_maxS, HV, LV; //int value, min_minS, max_maxS, HV, LV;
  double diffc, totsum, diffthresh;
  int upflag, upcount, avgnumspikes, *spks;
  double upvalue, limit, *diffcor, *percentc, *tots, sharp;
  LSeizure tmpSeizures, *pSeizuresOut;  
  sSeizure *tmpSeizure, *currentSeizure, *nextSeizure; 
  tmpSeizure=currentSeizure=nextSeizure=0x0;
  diffcor=percentc=tots=0x0; spks=0x0;
  if(!(pSeizuresOut=calloc(1,sizeof(LSeizure)))){printf("getseizures ERR: out of memory!\n");hxe();}
  bufszStart=400;
  InitLSeizure(pSeizuresOut,bufszStart); InitLSeizure(&tmpSeizures,bufszStart);
  intervalLength = (int)Bintdur*samprate;
  intervalsCount = channelsLength / intervalLength; // this is the # of intervals
  uplim=(int)(0.5+spkup*samprate/1e3); // round up
  if(intervalsCount<1) {
    printf("getseizures ERR:invalid intervalsCount:%d %d\n",channelsLength,intervalLength);hxe();}
  SgroupCount = (int)Bintdur/Lintdur*1e3; // Bintdur/Lintdur; intervalLength/30; 
  LintCnt=(int)Lintdur*samprate/1000;
  if(intervalsCount<1 || SgroupCount<1) { 
    printf("Error: Data length too short, cannot process."); hxe(); }
  if((diffcor = (double*) calloc(intervalsCount,sizeof(double)))==0x0) {
    printf("getseizures ERRdfc: out of memory!\n"); hxe();  }
  if(!(percentc = (double*) calloc(intervalsCount,sizeof(double)))) { // not used currently
    printf("getseizures ERRpcc: out of memory!\n"); hxe();  }            // but still calculated
  if(!(tots = (double*) calloc(intervalsCount,sizeof(double)))) {
    printf("getseizures ERRtts: out of memory!\n"); hxe();  }
  if(!(spks=(unsigned int*) calloc((size_t)(jj=intervalsCount*Bintdur/Sintdur),sizeof(int)))) {
    printf("getseizures ERRspk: out of memory!\n"); hxe();  }
  for (ii=0;ii<jj;ii++) spks[ii]=0; // clear
  totsLen = intervalsCount;
  seizuresCount=0;
  // gsz() calculates the correlation measures
  diffthresh=\
    gsz(channelData,diffcor,percentc,tots, channelsLength, intervalLength,SgroupCount,LintCnt);
  //** 2nd full loop: spike counting; find spikes for cases with high diffcor
  // Now that we have a suspected (intervalLength) 3000 pts interval, 
  // we test how many spikes there are in the interval
  gszspk(channelData,spks,tots,channelsLength,intervalLength);
  if (verbose==13) {
    if (pL->isz<1||pL->plen[0]!=intervalsCount) {
      printf("For verbose 13 need 1 vec of %d each\n",intervalsCount); hxe(); } 
    for (ii=0;ii<intervalsCount;ii++) pL->pv[0][ii]=(double)spks[ii];
  }
  if(incby1) for (ii=0; ii<intervalsCount; ii++) { // whole trace; ii=intervalIndex
    if (flag==0) {
      true=(spks[ii]>minspikes && diffcor[ii]>diffthresh);
    } else if (flag==1) {true=(spks[ii]>minspikes);
    } else if (flag==2) {true=(diffcor[ii]>diffthresh);
    }
    if (true) {
      if(seizuresCount==0 || tmpSeizure->endIntervalIndex!=ii-1)  { // new one
        seizuresCount++;
        AddSeizure(&tmpSeizures,tmpSeizure=AllocSeizure(ii,spks[ii]));
        AddSeizure(pSeizuresOut,AllocSeizure(ii,spks[ii]));
      } else { // add to current sz
        tmpSeizure->endIntervalIndex = ii; //ending index upddate
        tmpSeizure->totalSpikesCount += spks[ii]; //total spikes count
        tmpSeizure->endSpikesCount = spks[ii]; //end spikes count
      }
    }
  } // full trace loop; next intervalIndex pair
  else  for (ii=0; ii+1<intervalsCount; ii+=2) { // whole trace; ii=intervalIndex
    if (flag==0) {
      true=(spks[ii]>minspikes && spks[ii+1]>minspikes && diffcor[ii]>diffthresh && diffcor[ii+1]>diffthresh);
    } else if (flag==1) {true=(spks[ii]>minspikes && spks[ii+1]>minspikes);
    } else if (flag==2) {true=(diffcor[ii]>diffthresh && diffcor[ii+1]>diffthresh);
    }
    if (true) {
      if(seizuresCount==0 || tmpSeizure->endIntervalIndex!=ii-1)  { // new one
        seizuresCount++;
        AddSeizure(&tmpSeizures,tmpSeizure=AllocSeizure(ii,spks[ii]+spks[ii+1]));
        AddSeizure(pSeizuresOut,AllocSeizure(ii,spks[ii]+spks[ii+1]));
      } else { // add to current sz
        tmpSeizure->endIntervalIndex = ii+1; //ending index upddate
        tmpSeizure->totalSpikesCount += spks[ii]+spks[ii+1]; //total spikes count
        tmpSeizure->endSpikesCount = spks[ii+1]; //end spikes count
      }
    }
  } // full trace loop; next intervalIndex pair
  if(verbose==-1) {
    printf("before tmpSeizures %d:\n",tmpSeizures.count);  printSeizures(&tmpSeizures); }
  //** find the end points of the seizure. to do this compute the number of spikes in 2 second
  // intervals.  There will be a dramatic drop at the end of a seizure
  for (seizureIndex=0; seizureIndex<seizuresCount; seizureIndex++) { 
    tmpSeizure = tmpSeizures.pp[seizureIndex];
    currentSeizure = pSeizuresOut->pp[seizureIndex];
    currentSeizure->ID = seizureIndex;
    currentSeizure->totalSpikesCount = tmpSeizure->totalSpikesCount;
    if(seizureIndex < seizuresCount-1) { //start of next seizure (to prevent overlap)
      //potential problem for seizure overlap:
      //endIntervalIndex marks the beginning of interval and stop condition looks at it
      //so seizures overlap within 4s; next seizure start dataIndex
      stopIndex = (tmpSeizures.pp[seizureIndex+1])->startIntervalIndex * intervalLength; 
    } else stopIndex = channelsLength;
    intervalIndex = tmpSeizure->endIntervalIndex + 1; //current seizure end point 
    // (next seizure doesn't start in next interval), otherwise would be the same seizure
    dataIndex = intervalIndex * intervalLength; 
    if (intervalIndex >= totsLen) intervalIndex = totsLen-1;
    limit = tots[intervalIndex]/50.; // 50 hardcoded; count many more spikes than orig
    avgnumspikes=tmpSeizure->totalSpikesCount/\
      (Bintdur*(tmpSeizure->endIntervalIndex - tmpSeizure->startIntervalIndex+1));
    if (verbose==-2) printf("avgnumspikes: %d\n",avgnumspikes);
    upcount=upvalue=upflag=0;
    for (j=1; j<minspikes; j++) { // j not used
      spikesCount = 0;
      for (i=0; i<samprate*Sintdur; i++, dataIndex++) { // 2s interval; i not used
        if(dataIndex>=stopIndex) { //start of next seizure reached -> stop
          dataIndex = stopIndex-1; break; }
        if (dataIndex+1 >= channelsLength) break;
          value=channelData[dataIndex+1]-channelData[dataIndex];
          if(value>0) { upflag = 1; // true;
            upcount++; //count of how many periods was increasing
            upvalue+=value; // total value of increase
          } else	{ 
            if (upflag && upvalue>limit && upcount>uplim) {
              if(usesharp && dataIndex-sharpoff>=0.0 && dataIndex+sharpoff<channelsLength) {
                sharp = channelData[dataIndex+(int)sharpoff]-2*channelData[dataIndex]+channelData[dataIndex-(int)sharpoff];
                if(sharp < sharpth) spikesCount++; // only sharp spikes count
              } else spikesCount++;
            }
            upflag=upcount=upvalue=0;
          }
      } //Loop
      if (spikesCount/Sintdur < endwting*avgnumspikes) break; // spikes/sec
    } //Loop
    currentSeizure->endIndex=dataIndex; //update seizure end dataIndex
    currentSeizure->endIntervalIndex = (int)dataIndex/intervalLength; 
  }

  //** search/update the start of the seizure index. 
  // do this by computing the number of significant spikes and
  //determining where that changes significantly. use 2 second intervals
  seizuresCount = pSeizuresOut->count; 
  for (seizureIndex = 0; seizureIndex<seizuresCount; seizureIndex++) { 
    tmpSeizure = tmpSeizures.pp[seizureIndex];
    currentSeizure = pSeizuresOut->pp[seizureIndex];
    if (seizureIndex>0) { //end of previous seizure (to prevent overlap)
      //!! ?? here is potential problem for seizure overlap !!!
      //endIntervalIndex marks the beginning of interval and stop condition looks at it
      //so seizures overlap within 4s !!!
      //stopIndex = ((seizureTmp^)tmpSeizures[seizureIndex-1])->endIntervalIndex * intervalLength; 
      stopIndex = (pSeizuresOut->pp[seizureIndex-1])->endIndex;
    }
    else stopIndex = 0;
    intervalIndex = tmpSeizure->startIntervalIndex; //current seizure start point
    if(intervalIndex >= totsLen) intervalIndex = totsLen-1;
    dataIndex = intervalIndex * intervalLength; 
    limit= tots[intervalIndex]/50; // *0.02 hardcoded 50
    avgnumspikes = tmpSeizure->totalSpikesCount/\
      (Bintdur*(tmpSeizure->endIntervalIndex-tmpSeizure->startIntervalIndex+1));
    upcount=0.0; upvalue=upflag=0; 
    for(j=1; j<minspikes; j++) { 
      spikesCount = 0;
      for(i=0; i<samprate*Sintdur; i++, dataIndex--) { //2s interval; hardcoded 2
	if(dataIndex<=stopIndex) { //end of previous seizure reached -> stop
	  dataIndex = stopIndex + 1; break; }
        if(dataIndex+1 >= channelsLength) break;
        value=channelData[dataIndex+1]-channelData[dataIndex];
        if(value>0) { 
          upflag=1; upcount++; upvalue+=value;
        } else { 
          if (upflag && upvalue>limit && upcount>uplim) {
            if(usesharp && dataIndex-sharpoff>=0.0 && dataIndex+sharpoff<channelsLength) {
              sharp = channelData[dataIndex+(int)sharpoff]-2*channelData[dataIndex]+channelData[dataIndex-(int)sharpoff];
              if(sharp < sharpth) spikesCount++; // only count sharp spikes
            } else spikesCount++;
          }
          upflag=upcount=upvalue=0;
        }
      } //Loop
      if (spikesCount/Sintdur<endwting*avgnumspikes) break;
    }
    currentSeizure->startIndex = dataIndex; // update seizure start dataIndex
    currentSeizure->startIntervalIndex = (int)dataIndex/intervalLength; 
  }
  //** connect overlapping seizures
  for(seizureIndex = seizuresCount-1; seizureIndex>0; seizureIndex--) { 
    currentSeizure = pSeizuresOut->pp[seizureIndex-1];
    nextSeizure =    pSeizuresOut->pp[seizureIndex]; 
    if(currentSeizure->endIndex >= (nextSeizure->startIndex - mindist)) {//overlap ??
      //connect overlapping seizures
      currentSeizure->endIndex = nextSeizure->endIndex; //copy end to start
      currentSeizure->totalSpikesCount += nextSeizure->totalSpikesCount;
      free(nextSeizure);     // comment by sam : do we need to delete this seizure here?
      RemoveSeizureAt(pSeizuresOut,seizureIndex);
    }
  }
  //** deallocations
  if(verbose==-1) {
    printf("after tmpSeizures %d:\n",tmpSeizures.count);  printSeizures(&tmpSeizures);
    printf("after pSeizuresOut %d:\n",pSeizuresOut->count); printSeizures(pSeizuresOut); }
 STALEY_DOFREE:
  if(diffcor) free(diffcor);
  if(percentc) free(percentc);
  if(tots) free(tots);
  for(i=0; i<tmpSeizures.count; i++) { //delete tmpSeizures
    if(tmpSeizures.pp[i]==0x0) printf("tmpSeizures.pp[%d]=0x0",i);
    free(tmpSeizures.pp[i]); 
  }
  FreeLSeizure(&tmpSeizures); 
  return pSeizuresOut;
}

//* gsz()
static double gsz (double* channelData,double* diffcor,double* percentc,double* tots,\
          int channelsLength,int intervalLength, int SgroupCount,int LintCnt) {
  int intervalIndex, dataIndex, iS, intervalsCount;
  int i,j, ii, SgroupCount4, dbxi;
  double value, min_minS, max_maxS, HV, LV; //int value, min_minS, max_maxS, HV, LV;
  double diffc, totsum, avg1, sdev1, avg2, sdev2;
  double upvalue, limit, *minS, *maxS;
  double bintdur,lintdur; // lower case are local versions
  //** redund variable definitions from calling routine
  intervalsCount = channelsLength / intervalLength; // this is the # of intervals
  SgroupCount4 = SgroupCount/4;  //and this is the # of sub-sub-groups?? hardcoded 4
  //minS=gcnew array<float>(SgroupCount+2);//minS=gcnew array<int>(SgroupCount+2);
  //100 => (averages over 30 values) * 100 = 3000 
  if((minS = (double*) calloc(SgroupCount+2,sizeof(double)))==0x0) {
    printf("getseizures ERR: out of memory!\n"); hxe();  }
  //maxS=gcnew array<float>(SgroupCount+2);//maxS=gcnew array<int>(SgroupCount+2);
  //102 to make correct sum at the end of 100 (maxS(i+1) maxS(i+2))
  if((maxS = (double*) calloc(SgroupCount+2,sizeof(double)))==0x0){
    printf("getseizures ERR: out of memory!\n"); hxe();  }

  //** first loop: search for seizures; for each period of bintdur detect correlations
  for (dbxi=0,intervalIndex=0; intervalIndex<intervalsCount; intervalIndex++) { //thru trace
    diffcor[intervalIndex]=percentc[intervalIndex]=tots[intervalIndex]=totsum=diffc=value=0;
    //find minS(i) and maxS(i)  also max(maxS) and min(minS) for all intervals			
    // (100+2)*30 == 3060; we need 2 more since maxS(i+1) maxS(i+2)
    for (dataIndex=intervalIndex*intervalLength, iS=0; iS<SgroupCount+2; iS++) { // Big Int.
      minS[iS]= 1e22; maxS[iS]= -1e22; 
      for(i=0; i<LintCnt; i++, dataIndex++) {  // find min/max on Little Interval
	if(dataIndex < channelsLength) value = channelData[dataIndex];
        if (minS[iS] > value)  minS[iS] = value;
        if (maxS[iS] < value)  maxS[iS] = value;
      }
    }
    if (verbose==11) {
      if (dbxi==0 && (pL->isz<2 || pL->plen[0]!=intervalsCount*(SgroupCount+2)\
          || pL->plen[1]!=intervalsCount*(SgroupCount+2))) {
        printf("For verbose 11 need 2 vecs of %d each\n",intervalsCount*(SgroupCount+2)); hxe();} 
      for (iS=0; iS<SgroupCount+2; iS++,dbxi++) {
        pL->pv[0][dbxi]=minS[iS]; pL->pv[1][dbxi]=maxS[iS];
      }
    }
    // calculate metrics for each group
    for (iS=0; iS<SgroupCount; iS++) {
      HV = MIN(maxS[iS], MAX(maxS[iS+1],maxS[iS+2])); // compare to how it ends
      LV = MAX(minS[iS], MIN(minS[iS+1],minS[iS+2]));
      //metric3 = sum of differences over 100 intervals Si 
      diffc += (HV - LV); //!!! metric3
      totsum += (maxS[iS]-minS[iS]); //sum (max(Si)-min(Si))
    }    
    //Now produce the metrics
    diffcor[intervalIndex] = diffc;	//metric3
    percentc[intervalIndex] = diffc/totsum; // almost metric4; should be SUM(a/b) not SUM(a)/SUM(b)
    // recompute totsum using quartile maxima
    for(totsum=0, iS=0, j=0; j<4; j++)  {  // smoothing over 4
      min_minS = 1e22;  max_maxS = -1e22; 
      for(i=0; i<SgroupCount4; i++, iS++) { 
        if (max_maxS < maxS[iS])  max_maxS = maxS[iS];
        if (min_minS > minS[iS])  min_minS = minS[iS];
      }
      totsum += (max_maxS - min_minS);
    }
    tots[intervalIndex] = totsum/4;
    if (verbose==12) {
      if (dbxi==0) {
        for (ii=0;ii<3;ii++) if (pL->isz<3||pL->plen[ii]!=intervalsCount) {
          printf("For verbose 12 need 3 vecs of %d each\n",intervalsCount); hxe(); } 
        printf("Verbose 12: diffcor,percentc,tots\n");
      }
      ii=intervalIndex;
      pL->pv[0][dbxi]=diffcor[ii]; pL->pv[1][dbxi]=percentc[ii]; pL->pv[2][dbxi]=tots[ii];
      dbxi++;
    }
  } // whole trace; next intervalIndex
  //** Calculate the standard deviation of the high and low values
  // This is a measure of the correlation of the items - for a seizure the value should be low -
  // it will be higher for random processes.
  //  Review the results for this slice of time
  //  calculate mean and std-dev for percentc and diffcor
  avg1=avg2=sdev1=sdev2=0.0;
  for(intervalIndex=0; intervalIndex<intervalsCount; intervalIndex++)  { 
    avg1 += diffcor[intervalIndex];  // avg1 for diffcor 
    avg2 += percentc[intervalIndex]; // avg2 for percentc
    sdev1 += diffcor[intervalIndex]*diffcor[intervalIndex];   //for standard-deviations
    sdev2 += percentc[intervalIndex]*percentc[intervalIndex];
  }
  avg1 = avg1/intervalsCount; avg2 = avg2/intervalsCount; 
  sdev1 = sdev1/intervalsCount - avg1*avg1; // standard-deviations
  if(sdev1>0.) sdev1=sqrt(sdev1); else sdev1=avg1;
  sdev2 = sdev2/intervalsCount - avg2*avg2;
  if(sdev2>0.) sdev2=sqrt(sdev2); else sdev2=avg2;
  if (verbose>1) printf("diffcor: %g (%g,%g), percentc: %g (%g)\n",\
                        avg1,sdev1,avg1+abovebth*sdev1,avg2,sdev2);
  GSZ_DOFREE:
  if(minS) free(minS);
  if(maxS) free(maxS);
  return avg1+abovebth*sdev1;
}

#ifdef MYSPUD

int dspud (double* src, int nsrc, int lc) {
  int i, k, m, n, nqsz, nsrc, jj[UDSL], f[UDSL], lc, dsz[UDSL], nqmax, thsz, lc2, done, dbn;
  double *src, *tvec, *th, *dest[UDSL], *nq[UDNQ], *tmp, *dbx, lt, thdist;
  Object *ob, *ob2;
  void *vvd[UDSL], *vvth, *vnq[UDNQ];
  //** read in vectors and verify sizes, etc
  //nsrc = vector_instance_px(vv, &src); // trace to analyze
  thsz = vector_arg_px(1, &th);        // vector of thresholds to check
  ob =  *hoc_objgetarg(2);             // storage for values for each threshold
  ob2 = *hoc_objgetarg(3);             // list of NQS vectors for returning values
  tmp = (double *)ecalloc(nsrc, sizeof(double));  // tmp is size of trace
  lc =  ivoc_list_count(ob);
  lc2 = ivoc_list_count(ob2);
  if (lc>UDSL) {printf("updown ERRF mismatch: max slice list:%d %d\n",UDSL,lc); hxf(tmp);}
  if (lc2!=UDNQ){printf("updown ERRB mismatch: NQS sz is %d (%d in list)\n",UDNQ,lc2);hxf(tmp);}
  if (nsrc<lc) {printf("updown ERRC mismatch: %d %d\n",lc,nsrc); hxf(tmp);} // ??
  if (lc!=thsz) {printf("updown ERRA mismatch: %d %d\n",lc,thsz); hxf(tmp);}
  if (!ismono1(th,thsz,-1)) {printf("updown ERRD: not mono dec %g %d\n",th[0],thsz); hxf(tmp);}
  // thdist=(th[thsz-2]-th[thsz-1])/2; // NOT BEING USED: the smallest spike we will accept
  for (k=0;k <lc;k++)  dsz[k] =list_vector_px3(ob , k, &dest[k], &vvd[k]);
  for (k=0;k<lc2;k++) {
    i=list_vector_px3(ob2, k, &nq[k],   &vnq[k]);
    if (k==0) nqmax=i; else if (i!=nqmax) { // all NQ vecs same size
      printf("updown ERRE mismatch: %d %d %d\n",k,i,nqmax); hxf(tmp); }
  }
  //** store crossing points and midpoints in dest[k]
  // dest vectors dest[k] will store crossing points and midpoints at each th[k] slice location
  // as triplets: up/max/down
  for (k=0; k<lc; k++) {   // iterate thru thresholds
    jj[k]=f[k]=0; // jj[k] is ind into dest[k]; f[k] is flag for threshold  crossings
    for (i=0;i<nsrc && src[i]>th[k];i++) {} // start somewhere below this thresh th[k]
    for (; i<nsrc; i++) { // iterate through trace
      if (src[i]>th[k]) { 
        if (f[k]==0) { // ? passing thresh 
          if (jj[k]>=dsz[k]){printf("(%d,%d,%d) :: ",k,jj[k],dsz[k]);
            hoc_execerror("Dest vec too small in updown ", 0); }
          dest[k][jj[k]++] = (i-1) + (th[k]-src[i-1])/(src[i]-src[i-1]); // interpolate
          f[k]=1; 
          tmp[k]=-1e9; dest[k][jj[k]]=-1.; // flag in tmp says that a thresh found here
        }
        if (f[k]==1 && src[i]>tmp[k]) { // use tmp[] even more temporarily
          tmp[k]=src[i]; // pick out max
          dest[k][jj[k]] = (double)i; // location of this peak
        }
      } else {          // below thresh 
        if (f[k]==1) {  // just passed going down 
          jj[k]++;      // triplet will be indices of cross-up/peak/cross-down
          dest[k][jj[k]++] = (i-1) + (src[i-1]-th[k])/(src[i-1]-src[i]);
          f[k]=0; 
        }
      }
    }
  }
  //** truncate dest vectors to multiples of 3:
  for (k=0;k<lc;k++) vector_resize(vvd[k],(int)(floor((double)jj[k]/3.)*3.));
  for (i=0; i<nsrc; i++) tmp[i]=0.; // clear temp space
  //** go through all the slices to find identical peaks and save widths and locations
  // tmp[] uses triplets centered around a location corresponding to a max loc in the
  // original vector; the widest flanks for each are then on either side of this loc
  for (k=0;k<lc;k++) { // need to go from top to bottom to widen flanks
    for (i=1;i<jj[k];i+=3) { // through centers (peaks)
      m=(int)dest[k][i]; // hash: place center at location
      if (tmp[m-2]<0 || tmp[m-1]<0 || tmp[m+1]<0 || tmp[m+2]<0) continue; // ignore; too crowded
      tmp[m]--;  // count how many slices have found this peak (use negative)
      tmp[m-1]=dest[k][i-1]; tmp[m+1]=dest[k][i+1]; // flanks
    }
  }
  //** 1st (of 2) loops through tmp[] -- pick up flanks
  // step through tmp[] looking for negatives which indicate the slice count and pick up 
  // flanks from these
  // nq=new NQS("LOC","PEAK","WIDTH","BASE","HEIGHT","START","SLICES","SHARP","INDEX","FILE")
  for (i=0,k=0; i<nsrc; i++) if (tmp[i]<0.) { // tmp holds neg of count of slices
    if (k>=nqmax) { printf("updown ERRG OOR in NQ db: %d %d\n",k,nqmax); hxf(tmp); }
    LOC[k]=(double)i;  // approx location of the peak of the spike
    WIDTH[k]=tmp[i+1]; // location of right side -- temp storage
    START[k]=tmp[i-1]; // start of spike (left side)
    SLICES[k]=-tmp[i];  // # of slices
    k++;
  }
  nqsz=k;   // k ends up as size of NQS db
  if (DEBUG_UPDOWN && ifarg(4)) { dbn=vector_arg_px(4, &dbx); // DEBUG -- save tmp vector
    if (dbn<nsrc) printf("updown ERRH: Insufficient room in debug vec (%d<%d)\n",dbn,nsrc); 
    else for (i=0;i<nsrc;i++) dbx[i]=tmp[i]; 
  }
  //** adjust flanks to handle nested bumps
  // 3 ways to handle spike nested in a spike or elongated base:
  // NB always using same slice for both L and R flanks; NOV_UPDOWN flag: (no-overlap)
  //   0. nested spike(s) share flanks determined by shared base
  //   1. nested spike(s) have individual bases, 1st and last use flanks from base
  //   2. nested spike(s) have individual bases, base flanks listed separately w/out peak
  // here use 
  // search nq vecs to compare flanks to neighboring centers
  // if flanks overlap the centers on LT or RT side,
  // correct them by going back to original slice loc info (in dest[])
  //*** look at left side -- is this flank to left of center of another bump?
  if (NOV_UPDOWN) for (i=0;i<nqsz;i++) { // iterate through NQS db
    if ((i-1)>0 && START[i] < LOC[i-1]) { // flank is to left of prior center
      if (DEBUG_UPDOWN) printf("LT problem %d %g %g<%g\n",i,LOC[i],START[i],LOC[i-1]);
      for (m=lc-1,done=0;m>=0 && !done;m--) { // m:go from bottom (widest) to top
        for (n=1;n<jj[m] && !done;n+=3) {     // n:through centers
          // pick out lowest slice with this peak LOC whose flank is to RT of prior peak
          if (floor(dest[m][n])==LOC[i] && dest[m][n-1]>LOC[i-1]) {
            // ??[i]=START[i]; // temp storage for L end of this overlap
            // replace both left and right flanks at this level -- #1 above
            START[i]=dest[m][n-1]; WIDTH[i]=dest[m][n+1]; done=1; 
          }
        }
      }
    }
    //*** now look at RT side
    if ((i+1)<nqsz && WIDTH[i]>LOC[i+1]) {
      if (DEBUG_UPDOWN) printf("RT problem %d %g %g>%g\n",i,LOC[i],WIDTH[i],LOC[i+1]);
      for (m=lc-1,done=0;m>=0 && !done;m--) { // m: go from bottom to top
        for (n=1;n<jj[m] && !done;n+=3) {     // n: through centers
          // pick out lowest slice with this peak LOC whose flank is to LT of next peak
          if (floor(dest[m][n])==LOC[i] && dest[m][n+1]<LOC[i+1]) {
            // ??[i]=WIDTH[i]; // end of overlap
            START[i]=dest[m][n-1]; WIDTH[i]=dest[m][n+1]; done=1;
          }
        }
      }        
    }
  }

  //make sure left and right sides of bump occur at local minima
  //shouldn't creeping be before NOV_UPDOWN=1 overlap check???
  //creeping can result only in equal borders btwn two bumps
  //on one side, so it should be ok here...
  if(CREEP_UPDOWN) for(i=0,k=0;i<nsrc;i++) if(tmp[i]<0.){

    //move left side to local minima
    int idx = (int)START[k];
    while(idx >= 1 && src[idx] >= src[idx-1]) idx--;
    START[k] = idx;

    //move right side to local minima
    idx = (int)WIDTH[k];
    while(idx < nsrc-1 && src[idx] >= src[idx+1]) idx++;
    WIDTH[k] = idx;

    k++;
  }

  //** 2nd loop through tmp[] used to fill in the rest of NQS
  // needed to split into 2 loops so that could check for overlaps and correct those
  // before filling in the rest of nq
  for (i=0,k=0; i<nsrc; i++) if (tmp[i]<0.) { // tmp holds neg of count of slices
    // calculate a base voltage lt as interpolated value on left side
    lt=src[(int)floor(START[k])]+(START[k]-floor(START[k]))*\
      (src[(int)floor(START[k]+1.)]-src[(int)floor(START[k])]);
    BASE[k]=lt;         // base voltage
    PEAK[k]=src[i];     // peak voltage
    WIDTH[k] = WIDTH[k] - START[k]; // width = RT_flank-LT_flank
    HEIGHT[k]=PEAK[k]-BASE[k]; // redund measure -- can eliminate
    // measure of sharpness diff of 1st derivs btwn peak and SHM_UPDOWN dist from peak
    // to get 2nd deriv would be normalized by 2*SHM_UPDOWN*tstep
    // ??could take an ave. or max first deriv for certain distance on either side
    SHARP[k]=(src[i]-src[i-(int)SHM_UPDOWN])-(src[i+(int)SHM_UPDOWN]-src[i]);
    INDEX[k]=(double)k;
    k++;
  }
  
  int iNumBumps = k;

  //count # of other bumps nested within each bump
  if(!NOV_UPDOWN){
    for(i=0; i<iNumBumps; i++){
      NESTED[i] = 0;
      int j = 0;
      for(;j<iNumBumps;j++){
        if(i!=j && LOC[j] >= START[i] && LOC[j] <= START[i]+WIDTH[i]){
          NESTED[i]+=1.0;
        }
      }
    }
  } else for(i=0;i<iNumBumps;i++) NESTED[i]=0.0;

  //** finish up
  for (i=0;i<lc2;i++) vector_resize(vnq[i], nqsz);
  if (k!=nqsz) { printf("updown ERRI INT ERR: %d %d\n",k,nqsz); hxf(tmp); }
  free(tmp);
  return jj[0];
}

#endif

static double gszspk (double* channelData, int* spks,\
                      double* tots, int channelsLength, int intervalLength) {
  int intervalIndex, dataIndex, intervalsCount, foundspk;
  int i,j, ii, upflag, upcount, spikesCount, uplim, dbgSpikes,didpr,cnt;
  double value, upvalue, limit, sum, sharp; // lower case are local versions
  //** redund variable definitions from calling routine
  intervalsCount = channelsLength / intervalLength; // this is the # of intervals
  uplim=(int)(0.5+spkup*samprate/1e3);
  dbgSpikes=didpr=0;

  if(useavgtots) {
    sum=0.0;
    for(intervalIndex=0; intervalIndex<intervalsCount; intervalIndex++) sum += tots[intervalIndex];
    sum /= (double) intervalsCount;
    limit= sum*spklim;
  }

  for (intervalIndex=0; intervalIndex<intervalsCount; intervalIndex++) { // whole trace
    dataIndex = intervalIndex * intervalLength;
    if(!useavgtots) limit= tots[intervalIndex]*spklim; // 0.25 of average difference max_maxS-min_minS
    if(verbose>=15 && !didpr){ printf("limit = %g\n",limit); didpr=1; }
    upcount=upvalue=spikesCount=upflag=0;
    for (j=0; j<intervalLength; j++, dataIndex++)  {  // Bintdur ~12sec
      if (dataIndex > channelsLength) continue; // bounds-check added by Sam
      value = channelData[dataIndex+1]-channelData[dataIndex];
      if(value>0) { // increasing value
        upflag=1; upcount++; upvalue+=value;
      } else { // not increasing value
        foundspk=0; sharp=0.0;
        if(upflag && upvalue>limit && upcount>uplim) {
          foundspk=1;
          if(usesharp && dataIndex-sharpoff>=0.0 && dataIndex+sharpoff<channelsLength) {
            sharp = channelData[dataIndex+(int)sharpoff]-2*channelData[dataIndex]+channelData[dataIndex-(int)sharpoff];
            if(sharp > sharpth) foundspk=0; // only count sharp spikes
          }
        }
        if (foundspk) { // found a spike
          spikesCount++;  // increase the spike count
          if(verbose>=14) { // save spike x,y ?
            if(pL->isz < 3 || pL->plen[0]<channelsLength || pL->plen[1]<channelsLength || pL->plen[2]<1) {
              printf("need at least 2 vectors of size %d for verbose==14!\n",channelsLength); hxe();
            } else { 
              if(usesharp && dataIndex-sharpoff>=0.0 && dataIndex+sharpoff<channelsLength) {
                sharp = channelData[dataIndex+(int)sharpoff]-2*channelData[dataIndex]+channelData[dataIndex-(int)sharpoff];
                if(verbose>=15) printf("spike found (x,y,upcount,upvalue,sharp)=(%d,%g,%d,%g,%g)\n",dataIndex,channelData[dataIndex],upcount,upvalue,sharp);
              } else {
                sharp=0.0;
                if(verbose>=15) printf("spike found (x,y,upcount,upvalue)=(%d,%g,%d,%g)\n",dataIndex,channelData[dataIndex],upcount,upvalue);
              }              
              pL->pv[0][dbgSpikes]=dataIndex; pL->pv[1][dbgSpikes++]=channelData[dataIndex];}
          }
        }
        upflag=upcount=upvalue=0; // not increasing so reset counts
      }
    } // Bintdur loop; next j,dataIndex
    spks[intervalIndex]=spikesCount;
  }
  if(verbose>=14) pL->pv[2][0]=dbgSpikes;
  return 0.;
}

// veceeg.getseizures(totalSpikesCount,startIndex,endIndex)
// the 3 input args are Vectors to store results
static double getseizures (void* vv) {
  int n, cnt,i; LSeizure* pSeizures;
  double *p,*totalSpikesCount,*startIndex,*endIndex;
  n = vector_instance_px(vv,&p);
  if(verbose>10) {
    if (!ifarg(4)) { printf("Use veclist for dbx with verbose>10\n");hxe();
    } else pL=AllocListVec(*hoc_objgetarg(4));
  }
  if(!(pSeizures=dgetseizures(p,n))) return 0.0;
  cnt=pSeizures->count;
  totalSpikesCount=vector_newsize(vector_arg(1),cnt);
  startIndex=vector_newsize(vector_arg(2),cnt);
  endIndex=vector_newsize(vector_arg(3),cnt);

  for(i=0;i<cnt;i++) {
    totalSpikesCount[i] =   (double)pSeizures->pp[i]->totalSpikesCount;
    startIndex[i] =         (double)pSeizures->pp[i]->startIndex;
    endIndex[i] =           (double)pSeizures->pp[i]->endIndex;
  } 
  for(i=0;i<pSeizures->count;i++) free(pSeizures->pp[i]);
  FreeLSeizure(pSeizures);
  if (pL) FreeListVec(&pL);
  return (double)cnt;
}
ENDVERBATIM

PROCEDURE install () {
  VERBATIM
  if(!installed) {
    install_vector_method("getseizures",getseizures);
  } else printf("%s\n","$Id: staley.mod,v 1.75 2010/04/30 18:59:54 samn Exp $");
  ENDVERBATIM
  installed=1
}
