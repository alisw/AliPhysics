//$Id$

// Author: Anders Vestbo <mailto:vestbo$fi.uib.no>
//*-- Copyright &copy ASV

#include <stream.h>
#include <string.h>

#include "AliL3ModelTrack.h"
#include "AliL3Defs.h"

ClassImp(AliL3ModelTrack)

AliL3ModelTrack::AliL3ModelTrack()
{
  fNClusters = 0;
  fClusters = 0;
  fOverlap = -1;
  fPad=0;
  fTime=0;
  fClusterCharge=0;
}


AliL3ModelTrack::~AliL3ModelTrack()
{
  if(fClusters)
    delete [] fClusters;
  if(fPad)
    delete [] fPad;
  if(fTime)
    delete [] fTime;
}

void AliL3ModelTrack::Init(Int_t slice,Int_t patch)
{
  fNClusters = 0;
  Int_t nrows = NumRows[patch];
  fClusters = new ClusterComp[nrows];
  memset((void*)fClusters,0,nrows*sizeof(ClusterComp));
  
  fPad = new Float_t[NRowsSlice];
  fTime = new Float_t[NRowsSlice];
  
  fClusterCharge = 100;
  
}


void AliL3ModelTrack::SetCluster(Float_t fpad,Float_t ftime,Float_t charge,Float_t sigmaY2,Float_t sigmaZ2)
{
  ClusterComp *cl = &fClusters[fNClusters];
  if(!charge)
    cl->fEmpty = kTRUE;
  else
    {
      cl->fEmpty = kFALSE;
      cl->fDTime = ftime - GetTimeHit(fNClusters);
      cl->fDPad = fpad - GetPadHit(fNClusters);
      cl->fDCharge = charge - fClusterCharge;
    }
  cout<<"DPad "<<fpad<<" dtime "<<ftime<<" charge "<<charge<<" sigmaY2 "<<sigmaY2<<" sigmaZ2 "<<sigmaZ2<<endl;
  fNClusters++;
}

