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
  //
  // FUNCTIONALITY:
  //
  // 1. The laser track is associated with the mirror
  //    see function FindMirror
  //
  // 2. The laser track is accepted for the analysis under certain condition
  //    (see function Accpet laser)
  // 
  // 3. The drift velocity and jitter is calculated event by event
  //    (see function drift velocity) 
  //
  //
  //
  gSystem->Load("libANALYSIS");
  gSystem->Load("libTPCcalib");
  TFile fcalib("CalibObjects.root");
  TObjArray * array = (TObjArray*)fcalib.Get("TPCCalib");
  AliTPCcalibLaser * laser = ( AliTPCcalibLaser *)array->FindObject("laserTPC");
  laser->DumpMeanInfo(-0.4)
  TFile fmean("laserMean.root")

  //  laser track clasification;

  TCut cutT("cutT","abs(Tr.fP[3])<0.06");
  TCut cutPt("cutPt","abs(Tr.fP[4])<0.1");
  TCut cutN("cutN","fTPCncls>70");
  TCut cutP("cutP","abs(atan2(x1,x0)-atan2(lx1,lx0))<0.03")
  TCut cutA = cutT+cutPt+cutP;

  TFile f("laserTPCDebug.root");
  TTree * treeT = (TTree*)f.Get("Track");


  // LASER scan 
  gSystem->AddIncludePath("-I$ALICE_ROOT/TPC/macros");
  gROOT->LoadMacro("$ALICE_ROOT/TPC/macros/AliXRDPROOFtoolkit.cxx+")
  AliXRDPROOFtoolkit tool; 
  TChain * chain = tool.MakeChain("laserScan.txt","Mean",0,10200);
  chain->Lookup();


  treeT->Draw("(atan2(x1,x0)-atan2(lx1,lx0))*250.:fBundle","fSide==1&&fRod==0"+cutA,"prof") 

  gSystem->Load("libSTAT.so")
  TStatToolkit toolkit;
  Double_t chi2;
  TVectorD fitParam;
  TMatrixD covMatrix;
  Int_t npoints;
  TString *strq0 = toolkit.FitPlane(treeT,"Tr.fP[1]-LTr.fP[1]","lx1++lx2", "fSide==1"+cutA, chi2,npoints,fitParam,covMatrix);

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
  fTracksMirror(336),
  fTracksEsd(336),
  fTracksEsdParam(336),
  fTracksTPC(336),
  fDeltaZ(336),          // array of histograms of delta z for each track
  fDeltaPhi(336),          // array of histograms of delta z for each track
  fDeltaPhiP(336),          // array of histograms of delta z for each track
  fSignals(336),           // array of dedx signals
  fFitAside(new TVectorD(3)),        // drift fit - A side
  fFitCside(new TVectorD(3)),        // drift fit - C- side
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
  fTracksMirror(336),
  fTracksEsd(336),
  fTracksEsdParam(336),
  fTracksTPC(336),
  fDeltaZ(336),          // array of histograms of delta z for each track
  fDeltaPhi(336),          // array of histograms of delta z for each track
  fDeltaPhiP(336),          // array of histograms of delta z for each track
  fSignals(336),           // array of dedx signals
  fFitAside(new TVectorD(3)),        // drift fit - A side
  fFitCside(new TVectorD(3)),        // drift fit - C- side
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
  MakeDistHisto();
  //
  for (Int_t id=0; id<336; id++){
    //
    //
    if (!fTracksEsdParam.At(id)) continue;
    DumpLaser(id);
    RefitLaser(id);    
    
  }
}

void AliTPCcalibLaser::MakeDistHisto(){
  //
  //
  //
  for (Int_t id=0; id<336; id++){
    //
    //
    if (!fTracksEsdParam.At(id)) continue;  
    if (!AcceptLaser(id)) continue;
    //
    //
    TH1F * hisdz = (TH1F*)fDeltaZ.At(id);
    TH1F * hisdphi = (TH1F*)fDeltaPhi.At(id);
    TH1F * hisdphiP = (TH1F*)fDeltaPhiP.At(id);
    TH1F * hisSignal = (TH1F*)fSignals.At(id);

    if (!hisdz){      
      hisdz = new TH1F(Form("hisdz%d",id),Form("hisdz%d",id),500,-10,10);
      fDeltaZ.AddAt(hisdz,id);
      //
      hisdphi = new TH1F(Form("hisdphi%d",id),Form("hisdphi%d",id),500,-1,1);
      fDeltaPhi.AddAt(hisdphi,id);
      //
      hisdphiP = new TH1F(Form("hisdphiP%d",id),Form("hisdphiP%d",id),500,-0.01,0.01);
      fDeltaPhiP.AddAt(hisdphiP,id);
      hisSignal = new TH1F(Form("hisSignal%d",id),Form("hisSignal%d",id),500,0,1000);
      fSignals.AddAt(hisSignal,id);
    }

    AliExternalTrackParam *param=(AliExternalTrackParam*)fTracksEsdParam.At(id);
    AliTPCLaserTrack *ltrp = ( AliTPCLaserTrack*)fTracksMirror.At(id);
    AliESDtrack   *track    = (AliESDtrack*)fTracksEsd.At(id);
    Double_t xyz[3];
    Double_t pxyz[3];
    Double_t lxyz[3];
    Double_t lpxyz[3];
    param->GetXYZ(xyz);
    param->GetPxPyPz(pxyz);
    ltrp->GetXYZ(lxyz);
    ltrp->GetPxPyPz(lpxyz);
    //
    Float_t dz   = param->GetZ()-ltrp->GetZ();
    Float_t dphi = (TMath::ATan2(xyz[1],xyz[0])- TMath::ATan2(lxyz[1],lxyz[0]))*254.;
    Float_t dphiP = param->GetParameter()[2]-ltrp->GetParameter()[2];
    hisdz->Fill(dz);
    hisdphi->Fill(dphi);
    hisdphiP->Fill(dphiP); 
    hisSignal->Fill(track->GetTPCsignal());
  }
}

void AliTPCcalibLaser::FitDriftV(){
  //
  // Fit drift velocity - linear approximation in the z and global y
  //
  static TLinearFitter fdriftA(3,"hyp2");
  static TLinearFitter fdriftC(3,"hyp2");
  fdriftA.ClearPoints();
  fdriftC.ClearPoints();
  //
  for (Int_t id=0; id<336; id++){
    if (!fTracksEsdParam.At(id)) continue;  
    if (!AcceptLaser(id)) continue;
    AliExternalTrackParam *param=(AliExternalTrackParam*)fTracksEsdParam.At(id);
    AliTPCLaserTrack *ltrp = ( AliTPCLaserTrack*)fTracksMirror.At(id);
    Double_t xyz[3];
    Double_t pxyz[3];
    Double_t lxyz[3];
    Double_t lpxyz[3];
    param->GetXYZ(xyz);
    param->GetPxPyPz(pxyz);
    ltrp->GetXYZ(lxyz);
    ltrp->GetPxPyPz(lpxyz);
    Double_t xxx[2] = {lxyz[2],lxyz[1]};
    if (ltrp->GetSide()==0){
      fdriftA.AddPoint(xxx,xyz[2],1);
    }else{
      fdriftC.AddPoint(xxx,xyz[2],1);
    }
  }
  Float_t chi2A = 0;
  Float_t chi2C = 0;
  Int_t npointsA=0;
  Int_t npointsC=0;
  //
  if (fdriftA.GetNpoints()>10){
    fdriftA.Eval();
    fdriftA.EvalRobust(0.8);
    fdriftA.GetParameters(*fFitAside);
    npointsA= fdriftA.GetNpoints();
    chi2A = fdriftA.GetChisquare()/fdriftA.GetNpoints();
  }
  if (fdriftC.GetNpoints()>10){
    fdriftC.Eval();
    fdriftC.EvalRobust(0.8);
    fdriftC.GetParameters(*fFitCside);
    npointsC= fdriftC.GetNpoints();
    chi2C = fdriftC.GetChisquare()/fdriftC.GetNpoints();
  }
  
  if (fStreamLevel>0){
    TTreeSRedirector *cstream = GetDebugStreamer();
    Int_t time = fESD->GetTimeStamp();
    if (cstream){
      (*cstream)<<"driftv"<<
	"driftA.="<<fFitAside<<
	"driftC.="<<fFitCside<<
	"chi2A="<<chi2A<<
	"chi2C="<<chi2C<<
	"nA="<<npointsA<<
	"nC="<<npointsC<<
	"time="<<time<<
	"\n";
    }
  }
  //
}


Bool_t  AliTPCcalibLaser::AcceptLaser(Int_t id){
  //
  //
  //
  /*
  TCut cutT("cutT","abs(Tr.fP[3])<0.06");
  TCut cutPt("cutPt","abs(Tr.fP[4])<0.1");
  TCut cutN("cutN","fTPCncls>70");
  TCut cutP("cutP","abs(atan2(x1,x0)-atan2(lx1,lx0))<0.03")
  TCut cutA = cutT+cutPt+cutP;
  */
  AliExternalTrackParam *param =(AliExternalTrackParam*)fTracksEsdParam.At(id);
  AliTPCLaserTrack *ltrp  = ( AliTPCLaserTrack*)fTracksMirror.At(id);
  AliESDtrack   *track    = (AliESDtrack*)fTracksEsd.At(id);

  if (TMath::Abs(param->GetParameter()[4])>0.03) return kFALSE;
  if (TMath::Abs(param->GetParameter()[3])>0.06) return kFALSE;  
  if (TMath::Abs(param->GetParameter()[2]-ltrp->GetParameter()[2])>0.06) return kFALSE;
  if (TMath::Abs(param->GetParameter()[1]-ltrp->GetParameter()[1])>10) return kFALSE;
  //
  // dedx cut
  //
  if (TMath::Abs(track->GetTPCsignal())<20) return kFALSE;
  if (TMath::Abs(track->GetTPCsignal())>800) return kFALSE;
  //
  return kTRUE;  
}

Int_t  AliTPCcalibLaser::FindMirror(AliESDtrack *track, AliTPCseed *seed){
  //
  // Find corresponding mirror
  // add the corresponding tracks
  //
  Float_t kRadius0 = 252;
  Float_t kRadius  = 253.4;
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
    //
    //
    Float_t radius=TMath::Abs(ltrp->GetX());
    AliTracker::PropagateTrackTo(&param,radius,0.10566,0.01,kTRUE);
    //
    if (!fTracksMirror.At(id)) fTracksMirror.AddAt(ltrp,id);
    fTracksEsdParam.AddAt(param.Clone(),id);
    fTracksEsd.AddAt(track,id);
    fTracksTPC.AddAt(seed,id);
    //
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
  Double_t lxyz[3];
  Double_t lpxyz[3];
  param->GetXYZ(xyz);
  param->GetPxPyPz(pxyz);
  ltrp->GetXYZ(lxyz);
  ltrp->GetPxPyPz(lpxyz);

  if (fStreamLevel>0){
    TTreeSRedirector *cstream = GetDebugStreamer();
    Int_t time = fESD->GetTimeStamp();
    Bool_t accept = AcceptLaser(id);
    if (cstream){
      (*cstream)<<"Track"<<
	"run="<<fRun<<
	"id="<<id<<
	"accept="<<accept<<
	"driftA.="<<fFitAside<<
	"driftC.="<<fFitCside<<
	"time="<<time<<
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
	//
	"lx0="<<lxyz[0]<<
	"lx1="<<lxyz[1]<<
	"lx2="<<lxyz[2]<<
	"lpx0="<<lpxyz[0]<<
	"lpx1="<<lpxyz[1]<<
	"lpx2="<<lpxyz[2]<<
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


void AliTPCcalibLaser::DumpMeanInfo(Float_t bfield){
  //
  //
  //
  AliTPCcalibLaser *laser = this;
  TTreeSRedirector *pcstream = new TTreeSRedirector("laserMean.root");
  TF1 fg("fg","gaus");
  //
  //
  for (Int_t id=0; id<336; id++){
    TH1F * hisphi  = (TH1F*)laser->fDeltaPhi.At(id);
    TH1F * hisphiP = (TH1F*)laser->fDeltaPhiP.At(id);
    TH1F * hisZ    = (TH1F*)laser->fDeltaZ.At(id);
    //    TH1F * hisS    = laser->fDeltaZ->At(id);
    if (!hisphi) continue;;
    Double_t entries = hisphi->GetEntries();
    if (entries<30) continue;
    //
    AliTPCLaserTrack *ltrp = (AliTPCLaserTrack*)fTracksMirror.At(id);
    if (!ltrp) {
     AliTPCLaserTrack::LoadTracks();
      ltrp =(AliTPCLaserTrack*)AliTPCLaserTrack::GetTracks()->UncheckedAt(id);
    }
    Float_t meanphi = hisphi->GetMean();
    Float_t rmsphi = hisphi->GetRMS();
    Float_t meanphiP = hisphiP->GetMean();
    Float_t rmsphiP = hisphiP->GetRMS();
    Float_t meanZ = hisZ->GetMean();
    Float_t rmsZ = hisZ->GetRMS();
    hisphi->Fit(&fg,"","",hisphi->GetMean()-4*hisphi->GetRMS(),hisphi->GetMean()+4*hisphi->GetRMS());
    Double_t gphi1 = fg.GetParameter(1); 
    Double_t gphi2 = fg.GetParameter(2); 
    hisphiP->Fit(&fg,"","",hisphiP->GetMean()-4*hisphiP->GetRMS(),hisphiP->GetMean()+4*hisphiP->GetRMS());
    Double_t gphiP1 = fg.GetParameter(1); 
    Double_t gphiP2 = fg.GetParameter(2); 
    hisZ->Fit(&fg,"","",hisZ->GetMean()-4*hisZ->GetRMS(),hisZ->GetMean()+4*hisZ->GetRMS());
    Double_t gz1 = fg.GetParameter(1); 
    Double_t gz2 = fg.GetParameter(2); 

    //
    Double_t lxyz[3];
    Double_t lpxyz[3];
    ltrp->GetXYZ(lxyz);
    ltrp->GetPxPyPz(lpxyz);
    //
    (*pcstream)<<"Mean"<<
      "entries="<<entries<<      // number of entries
      "bz="<<bfield<<            // bfield
      "LTr.="<<ltrp<<             // refernece track
      "lx0="<<lxyz[0]<<          // reference x
      "lx1="<<lxyz[1]<<          // reference y
      "lx2="<<lxyz[2]<<          // refernece z      
      "lpx0="<<lpxyz[0]<<          // reference x
      "lpx1="<<lpxyz[1]<<          // reference y
      "lpx2="<<lpxyz[2]<<          // refernece z            
      "mphi="<<meanphi<<         //
      "rmsphi="<<rmsphi<<        //
      "gphi1="<<gphi1<<
      "gphi2="<<gphi2<<
      "mphiP="<<meanphiP<<       //
      "rmsphiP="<<rmsphiP<<      //
      "gphiP1="<<gphiP1<<
      "gphiP2="<<gphiP2<<
      "meanZ="<<meanZ<<
      "rmsZ="<<rmsZ<<
      "gz1="<<gz1<<
      "gz2="<<gz2<<

      "\n";
  }
  delete pcstream;
}

void AliTPCcalibLaser::Analyze(){
  //
  //
  //
}




