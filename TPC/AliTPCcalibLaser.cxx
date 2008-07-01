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

/*
  laser track clasification;
  TCut cutT("cutT","abs(Tr.fP[3])<0.06");
  TCut cutPt("cutPt","abs(Tr.fP[4])<0.1");
  TCut cutN("cutN","fTPCncls>100");
  
  TCut cutA = cutT+cutPt+cutN;

  TCut cutFi("cutZB","");
  TCut cutFi("cutFi","abs((180*atan2(x1,x0)/pi-20)%40)<5");
  TChain * chain = tool.MakeChain("chain.txt","Track",0,1000000)
*/



#include "TLinearFitter.h"
#include "AliTPCcalibLaser.h"
#include "AliExternalTrackParam.h"
#include "AliESDEvent.h"
#include "AliESDfriend.h"
#include "AliESDtrack.h"
#include "AliTPCTracklet.h"
#include "TH1D.h"
#include "TVectorD.h"
#include "TTreeStream.h"
#include "TFile.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "AliTPCclusterMI.h"
#include "AliTPCseed.h"
#include "AliTracker.h"
#include "TClonesArray.h"


#include "TTreeStream.h"
#include <iostream>
#include <sstream>
#include "AliTPCLaserTrack.h"

using namespace std;

ClassImp(AliTPCcalibLaser)

AliTPCcalibLaser::AliTPCcalibLaser():
  AliTPCcalibBase(),
  fESD(0),
  fESDfriend(),
  fTracksMirror(1000),
  fTracksEsd(1000),
  fTracksEsdParam(1000),
  fTracksTPC(1000),
  fRun(0)
{
  //
  // Constructor
  //
  fTracksEsdParam.SetOwner(kTRUE);
}

AliTPCcalibLaser::AliTPCcalibLaser(const Text_t *name, const Text_t *title):
  AliTPCcalibBase(),
  fESD(0),
  fESDfriend(0),
  fTracksMirror(1000),
  fTracksEsd(1000),
  fTracksEsdParam(1000),
  fTracksTPC(1000),
  fRun(0)
{
  SetName(name);
  SetTitle(title);
  //
  // Constructor
  //
  fTracksEsdParam.SetOwner(kTRUE);  
}

AliTPCcalibLaser::~AliTPCcalibLaser() {
  //
  // destructor
  //
}



void AliTPCcalibLaser::Process(AliESDEvent * event) {
  //
  // 
  // Loop over tracks and call  Process function
  //
  fESD = event;
  if (!fESD) {
    return;
  }
  fESDfriend=static_cast<AliESDfriend*>(fESD->FindListObject("AliESDfriend"));
  if (!fESDfriend) {
    return;
  }
  fTracksTPC.Clear();
  fTracksEsd.Clear();
  fTracksEsdParam.Delete();
  //
  Int_t n=fESD->GetNumberOfTracks();
  Int_t run = fESD->GetRunNumber();
  fRun = run;
  for (Int_t i=0;i<n;++i) {
    AliESDfriendTrack *friendTrack=fESDfriend->GetTrack(i);
    AliESDtrack *track=fESD->GetTrack(i);
    TObject *calibObject=0;
    AliTPCseed *seed=0;
    for (Int_t j=0;(calibObject=friendTrack->GetCalibObject(j));++j)
      if ((seed=dynamic_cast<AliTPCseed*>(calibObject)))
	break;
    if (track&&seed) FindMirror(track,seed);
    //
  }
  
  FitDriftV();
  
  //
  for (Int_t id=0; id<1000; id++){
    //
    //
    if (!fTracksEsdParam.At(id)) continue;
    DumpLaser(id);
    RefitLaser(id);    
  }
}

Int_t  AliTPCcalibLaser::FindMirror(AliESDtrack *track, AliTPCseed *seed){
  //
  // Find corresponding mirror
  // add the corresponding tracks
  //
  Float_t kRadius0 = 252;
  Float_t kRadius  = 254.25;
  if (!track->GetOuterParam()) return -1;
  AliExternalTrackParam param(*(track->GetOuterParam()));
  AliTracker::PropagateTrackTo(&param,kRadius0,0.10566,3,kTRUE);
  AliTracker::PropagateTrackTo(&param,kRadius,0.10566,0.1,kTRUE);
  AliTPCLaserTrack ltr;
  AliTPCLaserTrack *ltrp=0x0;
  //
  Int_t id = AliTPCLaserTrack::IdentifyTrack(&param);
  if (id!=-1 && (AliTPCLaserTrack::GetTracks()->UncheckedAt(id))) 
    ltrp=(AliTPCLaserTrack*)AliTPCLaserTrack::GetTracks()->UncheckedAt(id);
  else 
    ltrp=&ltr;
  
  if (id>=0){
    if (!fTracksMirror.At(id)) fTracksMirror.AddAt(ltrp,id);
    fTracksEsdParam.AddAt(param.Clone(),id);
    fTracksEsd.AddAt(track,id);
    fTracksTPC.AddAt(seed,id);
  }
  return id;
}



void AliTPCcalibLaser::DumpLaser(Int_t id) {
  //
  //  Dump Laser info to the tree
  //
  AliESDtrack   *track    = (AliESDtrack*)fTracksEsd.At(id);
  AliExternalTrackParam *param=(AliExternalTrackParam*)fTracksEsdParam.At(id);
  AliTPCLaserTrack *ltrp = ( AliTPCLaserTrack*)fTracksMirror.At(id);
  //
  // Fast laser ID
  //
  Double_t xyz[3];
  Double_t pxyz[3];
  param->GetXYZ(xyz);
  param->GetPxPyPz(pxyz);
  if (fStreamLevel>0){
    TTreeSRedirector *cstream = GetDebugStreamer();
    if (cstream){
      (*cstream)<<"Track"<<
	"run="<<fRun<<
	"id="<<id<<
	//
        "LTr.="<<ltrp<<
	"Esd.="<<track<<
	"Tr.="<<param<<
	"x0="<<xyz[0]<<
	"x1="<<xyz[1]<<
	"x2="<<xyz[2]<<
	"px0="<<pxyz[0]<<
	"px1="<<pxyz[1]<<
	"px2="<<pxyz[2]<<
	"\n";
    }
  }
}


void AliTPCcalibLaser::RefitLaser(Int_t id){
  //
  // Refit the track store residuals
  //

  AliTPCseed *track    = (AliTPCseed*)fTracksTPC.At(id);
  AliExternalTrackParam *param=(AliExternalTrackParam*)fTracksEsdParam.At(id);
  AliTPCLaserTrack *ltrp = (AliTPCLaserTrack*)fTracksMirror.At(id);
			     
  //
  static TLinearFitter fy2(3,"hyp2");
  static TLinearFitter fz2(3,"hyp2");
  static TLinearFitter fy1(2,"hyp1");
  static TLinearFitter fz1(2,"hyp1");
  static TVectorD vecy2,vecz2,vecy1,vecz1;

  const Int_t kMinClusters=20;
  Int_t nclusters[72]; 
  //
  for (Int_t i=0;i<72;++i) nclusters[i]=0;

  for (Int_t i=0;i<160;++i) {    
    AliTPCclusterMI *c=track->GetClusterPointer(i);
    if (c) nclusters[c->GetDetector()]++;
  }
   
  for (Int_t isec=0; isec<72;isec++){
    if (nclusters[isec]<kMinClusters) continue;
    fy2.ClearPoints();
    fz2.ClearPoints();
    fy1.ClearPoints();
    fz1.ClearPoints();
    //
    for (Int_t irow=0;irow<160;++irow) {      
      AliTPCclusterMI *c=track->GetClusterPointer(irow);
      //if (c && RejectCluster(c)) continue;
      if (c&&c->GetDetector()==isec) {
	Double_t xd = c->GetX()-120;;
	Double_t x[2]={xd,xd*xd};
	fy2.AddPoint(x,c->GetY());
	fz2.AddPoint(x,c->GetZ());
	//
	fy1.AddPoint(x,c->GetY());
	fz1.AddPoint(x,c->GetZ());
      }
    }
    fy2.Eval();
    fz2.Eval();
    fy1.Eval();
    fz1.Eval();
    fy1.GetParameters(vecy1);
    fy2.GetParameters(vecy2);
    fz1.GetParameters(vecz1);
    fz2.GetParameters(vecz2);
    
    if (fStreamLevel>0){
      TTreeSRedirector *cstream = GetDebugStreamer();
      if (cstream){
	Float_t dedx = track->GetdEdx();
	(*cstream)<<"Tracklet"<<
	  "LTr.="<<ltrp<<
	  "Tr.="<<param<<
	  "sec="<<isec<<
	  "ncl="<<nclusters[isec]<<
	  "dedx="<<dedx<<
	  "dedx="<<dedx<<
	  "vy1.="<<&vecy1<<
	  "vy2.="<<&vecy2<<
	  "vz1.="<<&vecz1<<
	  "vz2.="<<&vecz2<<
	  "\n";
      }
    }
  }
  //
  //
  //
  //   for (Int_t irow=0;irow<160;++irow) {      
  //       AliTPCclusterMI *c=track->GetClusterPointer(irow);
  //       if (c && RejectCluster(c)) continue;
  //       if (c&&c->GetDetector()==isec) {
  // 	Double_t x[2]={c->GetX(),c->GetX()*c->GetX()};
  // 	fy2.AddPoint(&x,c->GetY());
  // 	fz2.AddPoint(&x,c->GetZ());
  // 	//
  // 	fy1.AddPoint(&x,c->GetY());
  // 	fz1.AddPoint(&x,c->GetZ());
  //       }
  //     }    
  
}


void AliTPCcalibLaser::Analyze(){
  //
  //
  //
}




