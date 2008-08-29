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


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//   Component for redoing the reconstruction from the clusters and tracks
//      
//   The new calibration data used 
//
//   In reality it overwrites the content of the ESD 
//    

/*

  gSystem->Load("libANALYSIS");
  gSystem->Load("libTPCcalib");
  //
  gSystem->AddIncludePath("-I$ALICE_ROOT/TPC/macros");
  gROOT->LoadMacro("$ALICE_ROOT/TPC/macros/AliXRDPROOFtoolkit.cxx+")
  AliXRDPROOFtoolkit tool; 
  TChain * chain = tool.MakeChain("esd.txt","esdTree",0,5000);
  chain->Lookup();



*/



//  marian.ivanov@cern.ch
// 
#include "AliTPCcalibCalib.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTreeStream.h"
#include "AliLog.h"
#include "TTimeStamp.h"
#include "AliESDEvent.h"
#include "AliESDfriend.h"
#include "AliESDtrack.h"
#include "AliTracker.h"

#include "AliTPCcalibDB.h"
#include "AliTPCTransform.h"
#include "AliTPCclusterMI.h"
#include "AliTPCseed.h"

ClassImp(AliTPCcalibCalib)

AliTPCcalibCalib::AliTPCcalibCalib():
    AliTPCcalibBase()
{
  //
  // Constructor
  //
}


AliTPCcalibCalib::AliTPCcalibCalib(const Text_t *name, const Text_t *title) 
  :AliTPCcalibBase()
{  
  SetName(name);
  SetTitle(title);
}


AliTPCcalibCalib::AliTPCcalibCalib(const AliTPCcalibCalib&calib):
  AliTPCcalibBase(calib)
{
  //
  // copy constructor
  //
}

AliTPCcalibCalib &AliTPCcalibCalib::operator=(const AliTPCcalibCalib&calib){
  //
  //
  //
  ((AliTPCcalibBase *)this)->operator=(calib);
  return *this;
}


AliTPCcalibCalib::~AliTPCcalibCalib() {
  //
  // destructor
  //
}


void     AliTPCcalibCalib::Process(AliESDEvent *event){
  //
  // 
  //
  if (!event) {
    return;
  }  
  AliESDfriend *ESDfriend=static_cast<AliESDfriend*>(event->FindListObject("AliESDfriend"));
  if (!ESDfriend) {
   return;
  }
  
  if (GetDebugLevel()>20) printf("Hallo world: Im here\n");
  Int_t ntracks=event->GetNumberOfTracks();   
  //
  //
  //

  for (Int_t i=0;i<ntracks;++i) {
    AliESDtrack *track = event->GetTrack(i);  
    const AliExternalTrackParam * trackIn = track->GetInnerParam();
    const AliExternalTrackParam * trackOut = track->GetOuterParam();
    if (!trackIn) continue;
    if (!trackOut) continue;
   
    AliESDfriendTrack *friendTrack = ESDfriend->GetTrack(i);
    TObject *calibObject;
    AliTPCseed *seed = 0;
    for (Int_t l=0;(calibObject=friendTrack->GetCalibObject(l));++l) {
      if ((seed=dynamic_cast<AliTPCseed*>(calibObject))) break;
    }
    if (!seed) continue;
    RefitTrack(track, seed);
  }
  return;
}

Bool_t  AliTPCcalibCalib::RefitTrack(AliESDtrack * track, AliTPCseed *seed){
  //
  // Refit track
  //

  //
  // First apply calibration
  //

  AliTPCTransform *transform = AliTPCcalibDB::Instance()->GetTransform();
  for (Int_t irow=0;irow<159;irow++) {
    AliTPCclusterMI *cluster=seed->GetClusterPointer(irow);
    if (!cluster) continue; 
    AliTPCclusterMI cl0(*cluster);
    Double_t x[3]={cluster->GetRow(),cluster->GetPad(),cluster->GetTimeBin()};
    Int_t i[1]={cluster->GetDetector()};
    transform->Transform(x,i,0,1);
    //
    cluster->SetX(x[0]);
    cluster->SetY(x[1]);
    cluster->SetZ(x[2]);
    if (fStreamLevel>2){
      TTreeSRedirector *cstream = GetDebugStreamer();
      if (cstream){
	(*cstream)<<"Clusters"<<
	  "cl0.="<<&cl0<<
	  "cl.="<<cluster<<
	  "\n";
      }
    }
  }
  Int_t ncl = seed->GetNumberOfClusters();  
  // 
  // And now do refit
  //
  AliExternalTrackParam * trackInOld  = (AliExternalTrackParam*)track->GetTPCInnerParam();
  AliExternalTrackParam * trackOutOld = (AliExternalTrackParam*)track->GetOuterParam();

  AliExternalTrackParam trackIn  = *trackOutOld;
  AliExternalTrackParam trackOut = *trackInOld;
  trackIn.ResetCovariance(10.);
  trackOut.ResetCovariance(10.);
  Double_t xyz[3];
  Int_t nclIn=0,nclOut=0;
  //
  // Refit out
  //
  for (Int_t irow=0; irow<160; irow++){
    AliTPCclusterMI *cl=seed->GetClusterPointer(irow);
    if (!cl) continue;
    if (cl->GetX()<80) continue;
    Int_t sector = cl->GetDetector();
    Float_t dalpha = TMath::DegToRad()*(sector%18*20.+10.)-trackOut.GetAlpha();

    if (TMath::Abs(dalpha)>0.01)
      trackOut.Rotate(TMath::DegToRad()*(sector%18*20.+10.));
    if (RejectCluster(cl,&trackOut)) continue;
    Double_t r[3]={cl->GetX(),cl->GetY(),cl->GetZ()};
    Double_t cov[3]={0.01,0.,0.01}; //TODO: correct error parametrisation
    trackOut.GetXYZ(xyz);
    Double_t bz = AliTracker::GetBz(xyz);
    if (trackOut.PropagateTo(r[0],bz)) nclOut++;
    trackOut.Update(&r[1],cov);    
    
  }
  //
  // Refit in
  //

  for (Int_t irow=159; irow>0; irow--){
    AliTPCclusterMI *cl=seed->GetClusterPointer(irow);
    if (!cl) continue;
    if (cl->GetX()<80) continue;
    Int_t sector = cl->GetDetector();
    Float_t dalpha = TMath::DegToRad()*(sector%18*20.+10.)-trackIn.GetAlpha();
    if (TMath::Abs(dalpha)>0.01)
      trackIn.Rotate(TMath::DegToRad()*(sector%18*20.+10.));
    if (RejectCluster(cl,&trackIn)) continue;
    Double_t r[3]={cl->GetX(),cl->GetY(),cl->GetZ()};
    Double_t cov[3]={0.01,0.,0.01}; //TODO: correct error parametrisation
    trackOut.GetXYZ(xyz);
    Double_t bz = AliTracker::GetBz(xyz);

    if (trackIn.PropagateTo(r[0],bz)) nclIn++;
    trackIn.Update(&r[1],cov);    
  }
  trackIn.Rotate(trackInOld->GetAlpha());
  trackOut.Rotate(trackOutOld->GetAlpha());
  //
  trackInOld->GetXYZ(xyz);
  Double_t bz = AliTracker::GetBz(xyz);
  trackIn.PropagateTo(trackInOld->GetX(),bz);
  //
  trackOutOld->GetXYZ(xyz);
  bz = AliTracker::GetBz(xyz);  
  trackOut.PropagateTo(trackOutOld->GetX(),bz);
  

  if (fStreamLevel>0){
    TTreeSRedirector *cstream = GetDebugStreamer();
    if (cstream){
      (*cstream)<<"Tracks"<<
	"nclIn="<<nclIn<<
	"nclOut="<<nclOut<<
	"ncl="<<ncl<<
	"TrIn0.="<<trackInOld<<
	"TrOut0.="<<trackOutOld<<
	"TrIn1.="<<&trackIn<<
	"TrOut1.="<<&trackOut<<
	"\n";
    }
  }
  //
  // And now rewrite ESDtrack
  //
 
  (*trackInOld)  = trackIn;
  (*trackOutOld) = trackOut;
  AliExternalTrackParam *t = &trackIn;
  track->Set(t->GetX(),t->GetAlpha(),t->GetParameter(),t->GetCovariance());
  return kTRUE;
}



Bool_t AliTPCcalibCalib::RejectCluster(AliTPCclusterMI* cl, AliExternalTrackParam * param){
  //
  // check the acceptance of cluster
  // Cut on edge effects
  //
  Bool_t isReject = kFALSE;
  Float_t edgeY = cl->GetX()*TMath::Tan(TMath::Pi()/18);
  Float_t dist  = edgeY - TMath::Abs(cl->GetY());
  if (param)  dist  = TMath::Abs(edgeY - TMath::Abs(param->GetY()));
  if (dist<3) isReject=kTRUE;
  if (cl->GetType()<0) isReject=kTRUE;
  return isReject;
}


