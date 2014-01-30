/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *    
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  * 
 **************************************************************************/

/* $Id:  $ */

//_________________________________________________________________________//
//_________________________________________________________________________//

#include "AliESDTOFcluster.h"

ClassImp(AliESDTOFcluster)

AliESDTOFcluster::AliESDTOFcluster() :
  AliVTOFcluster(),
  fNTOFhits(1),
  fClusterIndex(new TArrayI(1)),
  fTOFchannel(new TArrayI(1)),
  fTime(new TArrayF(1)),
  fTimeRaw(new TArrayF(1)),
  fTOT(new TArrayF(1)),
  fTOFlabel(new TArrayI(3)),
  fDeltaBC(new TArrayI(1)),
  fL0L1Latency(new TArrayI(1)),
  fStatus(0),
  fZ(0),
  fPhi(0),
  fR(0),
  fNmatchableTracks(0),
  fTrackIndex(new TArrayI(1)),
  fDx(new TArrayF(1)),
  fDy(new TArrayF(1)),
  fDz(new TArrayF(1)),
  fTrackLength(new TArrayF(1)),
  fIntegratedTimes(new TArrayD(9))
{
  //
  // default ctor
  //
}

AliESDTOFcluster::AliESDTOFcluster(Int_t clusterIndex,Int_t tofChannel,Float_t tofTime,Float_t timeRaw,Float_t tofTot,Int_t label[3],Int_t deltaBC,Int_t l0l1Latency,
				   Bool_t status,Float_t zClu,Float_t phiClu,Float_t rClu,
				   Int_t trackIndex,Float_t dX,Float_t dY,Float_t dZ,Float_t length,Double_t expTimes[9]) :
  AliVTOFcluster(),
  fNTOFhits(1),
  fClusterIndex(new TArrayI(1)),
  fTOFchannel(new TArrayI(1)),
  fTime(new TArrayF(1)),
  fTimeRaw(new TArrayF(1)),
  fTOT(new TArrayF(1)),
  fTOFlabel(new TArrayI(3)),
  fDeltaBC(new TArrayI(1)),
  fL0L1Latency(new TArrayI(1)),
  fStatus(status),
  fZ(zClu),
  fPhi(phiClu),
  fR(rClu),
  fNmatchableTracks(1),
  fTrackIndex(new TArrayI(1)),
  fDx(new TArrayF(1)),
  fDy(new TArrayF(1)),
  fDz(new TArrayF(1)),
  fTrackLength(new TArrayF(1)),
  fIntegratedTimes(new TArrayD(9))
{
  //
  // Constructor of AliESDTOFcluster object
  //
  fClusterIndex->AddAt(clusterIndex,0);
  fTOFchannel->AddAt(tofChannel,0);
  fTime->AddAt(tofTime,0);
  fTimeRaw->AddAt(timeRaw,0);
  fTOT->AddAt(tofTot,0);
  fTOFlabel->AddAt(label[0],0);
  fTOFlabel->AddAt(label[1],1);
  fTOFlabel->AddAt(label[2],2);
  fDeltaBC->AddAt(deltaBC,0);
  fL0L1Latency->AddAt(l0l1Latency,0);
  fTrackIndex->AddAt(trackIndex,0);
  fDx->AddAt(dX,0);
  fDy->AddAt(dY,0);
  fDz->AddAt(dZ,0);
  fTrackLength->AddAt(length,0);
  for (Int_t ii=0; ii<9; ii++) fIntegratedTimes->AddAt(expTimes[ii],ii);

}

AliESDTOFcluster::AliESDTOFcluster(Int_t clusterIndex,Int_t tofChannel,Float_t tofTime,Float_t timeRaw,Float_t tofTot,Int_t label[3],Int_t deltaBC,Int_t l0l1Latency,
				   Bool_t status,Float_t zClu,Float_t phiClu,Float_t rClu) :
  AliVTOFcluster(),
  fNTOFhits(1),
  fClusterIndex(new TArrayI(1)),
  fTOFchannel(new TArrayI(1)),
  fTime(new TArrayF(1)),
  fTimeRaw(new TArrayF(1)),
  fTOT(new TArrayF(1)),
  fTOFlabel(new TArrayI(3)),
  fDeltaBC(new TArrayI(1)),
  fL0L1Latency(new TArrayI(1)),
  fStatus(status),
  fZ(zClu),
  fPhi(phiClu),
  fR(rClu),
  fNmatchableTracks(0),
  fTrackIndex(new TArrayI(1)),
  fDx(new TArrayF(1)),
  fDy(new TArrayF(1)),
  fDz(new TArrayF(1)),
  fTrackLength(new TArrayF(1)),
  fIntegratedTimes(new TArrayD(9))
{
  //
  // Constructor of AliESDTOFcluster object
  //
  fClusterIndex->AddAt(clusterIndex,0);
  fTOFchannel->AddAt(tofChannel,0);
  fTime->AddAt(tofTime,0);
  fTimeRaw->AddAt(timeRaw,0);
  fTOT->AddAt(tofTot,0);
  fTOFlabel->AddAt(label[0],0);
  fTOFlabel->AddAt(label[1],1);
  fTOFlabel->AddAt(label[2],2);
  fDeltaBC->AddAt(deltaBC,0);
  fL0L1Latency->AddAt(l0l1Latency,0);
}

AliESDTOFcluster::AliESDTOFcluster(const AliESDTOFcluster & source) :
  AliVTOFcluster(source),
  fNTOFhits(source.fNTOFhits),
  fClusterIndex(new TArrayI(source.fNTOFhits)),
  fTOFchannel(new TArrayI(source.fNTOFhits)),
  fTime(new TArrayF(source.fNTOFhits)),
  fTimeRaw(new TArrayF(source.fNTOFhits)),
  fTOT(new TArrayF(source.fNTOFhits)),
  fTOFlabel(new TArrayI(3*source.fNTOFhits)),
  fDeltaBC(new TArrayI(source.fNTOFhits)),
  fL0L1Latency(new TArrayI(source.fNTOFhits)),
  fStatus(source.fStatus),
  fZ(source.fZ),
  fPhi(source.fPhi),
  fR(source.fR),
  fNmatchableTracks(source.fNmatchableTracks),
  fTrackIndex(NULL),
  fDx(NULL),
  fDy(NULL),
  fDz(NULL),
  fTrackLength(NULL),
  fIntegratedTimes(NULL)
{
  // 
  // copy ctor for AliESDTOFcluster object
  //

  printf("Copy constructor (matchable tracks = %i)\n",fNmatchableTracks);

  for(Int_t i=0;i < source.fNTOFhits;i++){
    fClusterIndex->AddAt((*source.fClusterIndex)[i],i);
    fTOFchannel->AddAt((*source.fTOFchannel)[i],i);
    fTime->AddAt((*source.fTime)[i],i);
    fTimeRaw->AddAt((*source.fTimeRaw)[i],i);
    fTOT->AddAt((*source.fTOT)[i],i);
    fTOFlabel->AddAt((*source.fTOFlabel)[0+3*i],0+3*i);
    fTOFlabel->AddAt((*source.fTOFlabel)[1+3*i],1+3*i);
    fTOFlabel->AddAt((*source.fTOFlabel)[2+3*i],2+3*i);
    fDeltaBC->AddAt((*source.fDeltaBC)[i],i);
    fL0L1Latency->AddAt((*source.fL0L1Latency)[i],i);
  }

  if (fNmatchableTracks>0) {
    fTrackIndex = new TArrayI(fNmatchableTracks);
    fDx = new TArrayF(fNmatchableTracks);
    fDy = new TArrayF(fNmatchableTracks);
    fDz = new TArrayF(fNmatchableTracks);
    fTrackLength = new TArrayF(fNmatchableTracks);
    fIntegratedTimes = new TArrayD(fNmatchableTracks*9);
    for(Int_t i=0;i<fNmatchableTracks;i++) {
      (*fTrackIndex)[i]=source.fTrackIndex->At(i);
      (*fDx)[i]=source.fDx->At(i);
      (*fDy)[i]=source.fDy->At(i);
      (*fDz)[i]=source.fDz->At(i);
      (*fTrackLength)[i]=source.fTrackLength->At(i);
      for(Int_t j=0;j<9;j++)
	(*fIntegratedTimes)[i*9+j]=source.fIntegratedTimes->At(i*9+j);
    }
  }
  else{
    fTrackIndex = new TArrayI(1);
    fDx = new TArrayF(1);
    fDy = new TArrayF(1);
    fDz = new TArrayF(1);
    fTrackLength = new TArrayF(1);
    fIntegratedTimes = new TArrayD(9);
  }

  printf("END -> Copy constructor\n");

}

AliESDTOFcluster & AliESDTOFcluster::operator=(const AliESDTOFcluster & source)
{
  // 
  // copy ctor for AliESDTOFcluster object
  //
  if (this == &source) return *this;
  AliVTOFcluster::operator=(source);

  fNTOFhits = source.fNTOFhits;

  if(fClusterIndex) delete fClusterIndex;
  if(fTOFchannel->GetArray()) delete fTOFchannel;
  if(fTime) delete fTime;
  if(fTimeRaw) delete fTimeRaw;
  if(fTOT) delete fTOT;
  if(fTOFlabel) delete fTOFlabel;
  if(fDeltaBC) delete fDeltaBC;
  if(fL0L1Latency) delete fL0L1Latency;

  if(fNTOFhits){
    fClusterIndex = new TArrayI(source.fNTOFhits);
    fTOFchannel = new TArrayI(source.fNTOFhits);
    fTime = new TArrayF(source.fNTOFhits);
    fTimeRaw = new TArrayF(source.fNTOFhits);
    fTOT = new TArrayF(source.fNTOFhits);
    fTOFlabel = new TArrayI(3*source.fNTOFhits);
    fDeltaBC = new TArrayI(source.fNTOFhits);
    fL0L1Latency = new TArrayI(source.fNTOFhits);
      
    
    for(Int_t i=0;i < source.fNTOFhits;i++){
      fClusterIndex->AddAt(source.fClusterIndex->At(i),i);
      fTOFchannel->AddAt(source.fTOFchannel->At(i),i);
      fTime->AddAt(source.fTime->At(i),i);
      fTimeRaw->AddAt(source.fTimeRaw->At(i),i);
      fTOT->AddAt(source.fTOT->At(i),i);
      fTOFlabel->AddAt(source.fTOFlabel->At(0+3*i),0+3*i);
      fTOFlabel->AddAt(source.fTOFlabel->At(1+3*i),1+3*i);
      fTOFlabel->AddAt(source.fTOFlabel->At(2+3*i),2+3*i);
      fDeltaBC->AddAt(source.fDeltaBC->At(i),i);
      fL0L1Latency->AddAt(source.fL0L1Latency->At(i),i);
    }
  }
  else{
    fClusterIndex = new TArrayI(1);
    fTOFchannel = new TArrayI(1);
    fTime = new TArrayF(1);
    fTimeRaw = new TArrayF(1);
    fTOT = new TArrayF(1);
    fTOFlabel = new TArrayI(3);
    fDeltaBC = new TArrayI(1);
    fL0L1Latency = new TArrayI(1);    
  }

  fStatus=source.fStatus;
  fZ=source.fZ;
  fPhi=source.fPhi;
  fR=source.fR;

  if(fTrackIndex) delete fTrackIndex;
  if(fDx) delete fDx;
  if(fDy) delete fDy;
  if(fDz) delete fDz;
  if(fTrackLength) delete fTrackLength;
  if(fIntegratedTimes) delete fIntegratedTimes;

  fNmatchableTracks=source.fNmatchableTracks;

  if (fNmatchableTracks>0) {
    fTrackIndex = new TArrayI(fNmatchableTracks);
    fDx = new TArrayF(fNmatchableTracks);
    fDy = new TArrayF(fNmatchableTracks);
    fDz = new TArrayF(fNmatchableTracks);
    fTrackLength = new TArrayF(fNmatchableTracks);
    fIntegratedTimes = new TArrayD(fNmatchableTracks*9);
    for(Int_t i=0;i<fNmatchableTracks;i++) {
      (*fTrackIndex)[i]=source.fTrackIndex->At(i);
      (*fDx)[i]=source.fDx->At(i);
      (*fDy)[i]=source.fDy->At(i);
      (*fDz)[i]=source.fDz->At(i);
      (*fTrackLength)[i]=source.fTrackLength->At(i);
      for(Int_t j=0;j<9;j++)
	(*fIntegratedTimes)[i*9+j]=source.fIntegratedTimes->At(i*9+j);
    }
  } else {
    fTrackIndex= new TArrayI(1);
    fDx = new TArrayF(1);
    fDy = new TArrayF(1);
    fDz = new TArrayF(1);
    fTrackLength = new TArrayF(1);
    fIntegratedTimes = new TArrayD(9);
  }

  return *this;
}

Int_t AliESDTOFcluster::Update(Int_t trackIndex,Float_t dX,Float_t dY,Float_t dZ,Float_t length,Double_t expTimes[9])
{

  if (fNmatchableTracks==0) {
    fNmatchableTracks++;
    fTrackIndex = new TArrayI(fNmatchableTracks);
    fDx = new TArrayF(fNmatchableTracks);
    fDy = new TArrayF(fNmatchableTracks);
    fDz = new TArrayF(fNmatchableTracks);
    fTrackLength = new TArrayF(fNmatchableTracks);
    fIntegratedTimes = new TArrayD(fNmatchableTracks*9);
    fTrackIndex->AddAt(trackIndex,fNmatchableTracks-1);
    fDx->AddAt(dX,fNmatchableTracks-1);
    fDy->AddAt(dY,fNmatchableTracks-1);
    fDz->AddAt(dZ,fNmatchableTracks-1);
    fTrackLength->AddAt(length,fNmatchableTracks-1);
    for (Int_t ii=0; ii<9; ii++) fIntegratedTimes->AddAt(expTimes[ii],(fNmatchableTracks-1)*9+ii);

  } else {
    Bool_t updatedOneAlreadyStored = kFALSE;
    for (Int_t ii=0; ii<fNmatchableTracks; ii++) {
      if (trackIndex==fTrackIndex->At(ii)) {
	// evitare di mettere piu' trackPoints della stessa traccia,
	// scegliere solo quello piu' vicino al piano centrale del pad
	updatedOneAlreadyStored=kTRUE;
// 	if (TMath::Abs(dX)<TMath::Abs(fDx->At(ii))) {
// 	  fDx->AddAt(dX,ii);
// 	  fDy->AddAt(dY,ii);
// 	  fDz->AddAt(dZ,ii);
// 	  fTrackLength->AddAt(length,ii);
// 	  for (Int_t jj=0; jj<9; jj++) fIntegratedTimes->AddAt(expTimes[jj],ii*9+jj);
// 	}
	return 1;
      }
    }
    if (!updatedOneAlreadyStored) {
      fNmatchableTracks++;
      fTrackIndex->Set(fNmatchableTracks);
      fDx->Set(fNmatchableTracks);
      fDy->Set(fNmatchableTracks);
      fDz->Set(fNmatchableTracks);
      fTrackLength->Set(fNmatchableTracks);
      fIntegratedTimes->Set(fNmatchableTracks*9);
      fTrackIndex->AddAt(trackIndex,fNmatchableTracks-1);
      fDx->AddAt(dX,fNmatchableTracks-1);
      fDy->AddAt(dY,fNmatchableTracks-1);
      fDz->AddAt(dZ,fNmatchableTracks-1);
      fTrackLength->AddAt(length,fNmatchableTracks-1);
      for (Int_t jj=0; jj<9; jj++) fIntegratedTimes->AddAt(expTimes[jj],(fNmatchableTracks-1)*9+jj);
    }
  }
  return 0;
}

AliESDTOFcluster::~AliESDTOFcluster()
{
  //
  // dtr
  //
  if(fClusterIndex) delete fClusterIndex;
  if(fTOFchannel) delete fTOFchannel;
  if(fTime) delete fTime;
  if(fTimeRaw) delete fTimeRaw;
  if(fTOT) delete fTOT;
  if(fTOFlabel) delete fTOFlabel;
  if(fDeltaBC) delete fDeltaBC;
  if(fL0L1Latency) delete fL0L1Latency;

  if (fTrackIndex) delete fTrackIndex;
  if (fDx) delete fDx;
  if (fDy) delete fDy;
  if (fDz) delete fDz;
  if (fTrackLength) delete fTrackLength;
  if (fIntegratedTimes) delete fIntegratedTimes;

}

void AliESDTOFcluster::AddTOFhit(Int_t clusterIndex,Int_t tofChannel,Float_t tofTime,Float_t timeRaw,Float_t tofTot,Int_t label[3],Int_t deltaBC,Int_t l0l1Latency, Bool_t status,Float_t zClu,Float_t phiClu,Float_t rClu){

  TArrayI *ClusterIndexOld = fClusterIndex;
  TArrayI *TOFchannelOld = fTOFchannel;
  TArrayF *TimeOld = fTime;
  TArrayF *TimeRawOld = fTimeRaw;
  TArrayF *TOTOld = fTOT;
  TArrayI *TOFlabelOld = fTOFlabel;
  TArrayI *DeltaBCOld = fDeltaBC;
  TArrayI *L0L1LatencyOld = fL0L1Latency;
  
  fClusterIndex = new TArrayI(fNTOFhits+1);
  fTOFchannel = new TArrayI(fNTOFhits+1);
  fTime = new TArrayF(fNTOFhits+1);
  fTimeRaw = new TArrayF(fNTOFhits+1);
  fTOT = new TArrayF(fNTOFhits+1);
  fTOFlabel = new TArrayI(3*fNTOFhits+3);
  fDeltaBC = new TArrayI(fNTOFhits+1);
  fL0L1Latency = new TArrayI(fNTOFhits+1);

  for(Int_t i=0;i < fNTOFhits;i++){
    fClusterIndex->AddAt(ClusterIndexOld->At(i),i);
    fTOFchannel->AddAt(TOFchannelOld->At(i),i);
    fTime->AddAt(TimeOld->At(i),i);
    fTimeRaw->AddAt(TimeRawOld->At(i),i);
    fTOT->AddAt(TOTOld->At(i),i);
    fTOFlabel->AddAt(TOFlabelOld->At(0+3*i),0+3*i);
    fTOFlabel->AddAt(TOFlabelOld->At(1+3*i),1+3*i);
    fTOFlabel->AddAt(TOFlabelOld->At(2+3*i),2+3*i);
    fDeltaBC->AddAt(DeltaBCOld->At(i),i);
    fL0L1Latency->AddAt(L0L1LatencyOld->At(i),i);
  }
  
  if(ClusterIndexOld) delete ClusterIndexOld;
  if(TOFchannelOld) delete TOFchannelOld;
  if(TimeOld) delete TimeOld;
  if(TimeRawOld) delete TimeRawOld;
  if(TOTOld) delete TOTOld;
  if(TOFlabelOld) delete TOFlabelOld;
  if(DeltaBCOld) delete DeltaBCOld;
  if(L0L1LatencyOld) delete L0L1LatencyOld;

  fClusterIndex->AddAt(clusterIndex,fNTOFhits);
  fTOFchannel->AddAt(tofChannel,fNTOFhits);
  fTime->AddAt(tofTime,fNTOFhits);
  fTimeRaw->AddAt(timeRaw,0);
  fTOT->AddAt(tofTot,fNTOFhits);
  fTOFlabel->AddAt(label[0],3*fNTOFhits+0);
  fTOFlabel->AddAt(label[1],3*fNTOFhits+1);
  fTOFlabel->AddAt(label[2],3*fNTOFhits+2);
  fDeltaBC->AddAt(deltaBC,fNTOFhits);
  fL0L1Latency->AddAt(l0l1Latency,fNTOFhits);

  if(status){
    fZ = (fZ*fNTOFhits + zClu)/(fNTOFhits+1);
    fPhi = (fPhi*fNTOFhits + phiClu)/(fNTOFhits+1);
    fR = (fR*fNTOFhits + rClu)/(fNTOFhits+1);
  }

  fNTOFhits++;
}
