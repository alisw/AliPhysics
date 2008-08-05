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
  // To make laser scan the user interaction neccessary
  //
  .x ~/UliStyle.C
  gSystem->Load("libANALYSIS");
  gSystem->Load("libTPCcalib");
  TFile fcalib("CalibObjects.root");
  TObjArray * array = (TObjArray*)fcalib.Get("TPCCalib");
  AliTPCcalibLaser * laser = ( AliTPCcalibLaser *)array->FindObject("laserTPC");
  laser->DumpMeanInfo(-0.4)
  TFile fmean("laserMean.root")
  //
  //  laser track clasification;
  //
  TCut cutT("cutT","abs(Tr.fP[3])<0.06");
  TCut cutPt("cutPt","abs(Tr.fP[4])<0.1");
  TCut cutN("cutN","fTPCncls>70");
  TCut cutP("cutP","abs(atan2(x1,x0)-atan2(lx1,lx0))<0.03")
  TCut cutA = cutT+cutPt+cutP;
  TFile f("laserTPCDebug.root");
  TTree * treeT = (TTree*)f.Get("Track");
  //
  //
  // Analyze  LASER scan 
  //
  gSystem->AddIncludePath("-I$ALICE_ROOT/TPC/macros");
  gROOT->LoadMacro("$ALICE_ROOT/TPC/macros/AliXRDPROOFtoolkit.cxx+")
  AliXRDPROOFtoolkit tool;
  TChain * chain = tool.MakeChain("laserScan.txt","Mean",0,10200);
  chain->Lookup();
  AliTPCcalibLaser::DumpScanInfo(chain)
  TFile fscan("laserScan.root")
  TTree * treeT = (TTree*)fscan.Get("Mean")
  //
  // Analyze laser 
  //
  gSystem->AddIncludePath("-I$ALICE_ROOT/TPC/macros");
  gROOT->LoadMacro("$ALICE_ROOT/TPC/macros/AliXRDPROOFtoolkit.cxx+")
  AliXRDPROOFtoolkit tool;
  TChain * chain = tool.MakeChain("laser.txt","Residuals",0,10200);
  chain->Lookup();

*/



#include "TLinearFitter.h"
#include "AliTPCcalibLaser.h"
#include "AliExternalTrackParam.h"
#include "AliESDEvent.h"
#include "AliESDfriend.h"
#include "AliESDtrack.h"
#include "AliTPCTracklet.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TVectorD.h"
#include "TTreeStream.h"
#include "TFile.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "AliTPCclusterMI.h"
#include "AliTPCseed.h"
#include "AliTracker.h"
#include "AliLog.h"
#include "TClonesArray.h"
#include "TPad.h"


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
  fDeltaZ(336),
  fDeltaPhi(336),
  fDeltaPhiP(336),
  fSignals(336),
  fDeltaYres(336),
  fDeltaZres(336),  
  fPol2Par2InY(336),
  fDiffPar1InY(336),
  fPol2Par2OutY(336),
  fDiffPar1OutY(336),
  fPol2Par2InZ(336),
  fDiffPar1InZ(336),
  fPol2Par2OutZ(336),
  fDiffPar1OutZ(336),
  fFitAside(new TVectorD(3)),
  fFitCside(new TVectorD(3)),      
  fEdgeXcuts(5),    
  fEdgeYcuts(5),    
  fNClCuts(5),      
  fNcuts(0),        
  fRun(0),
  fEvent(0)
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
  fDeltaYres(336),
  fDeltaZres(336),  
  fPol2Par2InY(336),
  fDiffPar1InY(336),
  fPol2Par2OutY(336),
  fDiffPar1OutY(336),
  fPol2Par2InZ(336),
  fDiffPar1InZ(336),
  fPol2Par2OutZ(336),
  fDiffPar1OutZ(336),
  fFitAside(new TVectorD(3)),        // drift fit - A side
  fFitCside(new TVectorD(3)),        // drift fit - C- side
  fEdgeXcuts(5),       // cuts in local x direction; used in the refit of the laser tracks
  fEdgeYcuts(5),       // cuts in local y direction; used in the refit of the laser tracks
  fNClCuts(5),         // cuts on the number of clusters per tracklet; used in the refit of the laser tracks
  fNcuts(0),           // number of cuts
  fRun(0),
  fEvent(0)
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
  AliDebug(4,Form("Event number in current file: %d",event->GetEventNumberInFile()));
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
//    RefitLaser(id);
    RefitLaserJW(id);

  }
//  fEvent++;
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
      hisdz = new TH1F(Form("hisdz%d",id),Form("hisdz%d",id),1000,-10,10);
      hisdz->SetDirectory(0);
      fDeltaZ.AddAt(hisdz,id);
      //
      hisdphi = new TH1F(Form("hisdphi%d",id),Form("hisdphi%d",id),1000,-1,1);
      hisdphi->SetDirectory(0);
      fDeltaPhi.AddAt(hisdphi,id);
      //
      hisdphiP = new TH1F(Form("hisdphiP%d",id),Form("hisdphiP%d",id),1000,-0.01,0.01);
      hisdphiP->SetDirectory(0);
      fDeltaPhiP.AddAt(hisdphiP,id);
      hisSignal = new TH1F(Form("hisSignal%d",id),Form("hisSignal%d",id),1000,0,1000);
      hisSignal->SetDirectory(0);
      fSignals.AddAt(hisSignal,id);
    }

    AliExternalTrackParam *param=(AliExternalTrackParam*)fTracksEsdParam.At(id);
    AliTPCLaserTrack *ltrp = ( AliTPCLaserTrack*)fTracksMirror.At(id);
    AliESDtrack   *track    = (AliESDtrack*)fTracksEsd.At(id);
    if (!param) return;
    if (!ltrp) return;
    if (!track) return;
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
    if (hisdz) hisdz->Fill(dz);
    if (hisdphi) hisdphi->Fill(dphi);
    if (hisdphiP) hisdphiP->Fill(dphiP);
    if (hisSignal) hisSignal->Fill(track->GetTPCsignal());
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
  //  AliTPCLaserTrack *ltrpjw=0x0;
  //
  Int_t id   = AliTPCLaserTrack::IdentifyTrack(&param);
 // Int_t idjw = AliTPCLaserTrack::IdentifyTrackJW(&param);
  //AliDebug(4,Form("Identified Track: %03d (%03d)",id,idjw));

  if (id!=-1 && (AliTPCLaserTrack::GetTracks()->UncheckedAt(id)))
    ltrp=(AliTPCLaserTrack*)AliTPCLaserTrack::GetTracks()->UncheckedAt(id);
  else
    ltrp=&ltr;
  /*
    if (idjw!=-1 && (AliTPCLaserTrack::GetTracks()->UncheckedAt(idjw)))
    ltrpjw=(AliTPCLaserTrack*)AliTPCLaserTrack::GetTracks()->UncheckedAt(idjw);
  else
    ltrpjw=&ltr;


    if (fStreamLevel>0){
    TTreeSRedirector *cstream = GetDebugStreamer();
    if (cstream){
	(*cstream)<<"idcmp"<<
	    "id=" << id <<
	    "idjw=" << idjw <<
	    "tr.="  << ltrp <<
	    "trjw.="<< ltrpjw <<
	    "seed.="<<seed<<
            "event="<<fEvent <<
	    "\n";
    }
  }

      */


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

void AliTPCcalibLaser::RefitLaserJW(Int_t id){
  //
  // Refit the track with different tracklet models:
  // 1. Per ROC using the kalman filter, different edge cuts
  // 2. Per ROC linear in y and z
  // 3. Per ROC quadratic in y and z
  // 4. Per track offset for each sector, linear for each sector, common quadratic
  // store x, y, z information for all models and the cluster to calculate the residuals
  //
  AliTPCseed *track      = (AliTPCseed*)fTracksTPC.At(id);
  AliExternalTrackParam *extparam=(AliExternalTrackParam*)fTracksEsdParam.At(id);
  AliTPCLaserTrack *ltrp = (AliTPCLaserTrack*)fTracksMirror.At(id);

  AliTPCclusterMI dummyCl;

  //two tracklets
  Int_t kMaxTracklets=2;
  //=============================================//
  // Linear Fitters for the Different Approaches //
  //=============================================//
  //linear fit model in y and z; inner - outer sector
  static TLinearFitter fy1I(2,"hyp1");
  static TLinearFitter fy1O(2,"hyp1");
  static TLinearFitter fz1I(2,"hyp1");
  static TLinearFitter fz1O(2,"hyp1");
  //quadratic fit model in y and z; inner - sector
  static TLinearFitter fy2I(3,"hyp2");
  static TLinearFitter fy2O(3,"hyp2");
  static TLinearFitter fz2I(3,"hyp2");
  static TLinearFitter fz2O(3,"hyp2");
  //common quadratic fit for IROC and OROC in y and z
  static TLinearFitter fy4(5,"hyp4");
  static TLinearFitter fz4(5,"hyp4");


  //set standard cuts
  if ( fNcuts==0 ){
      fNcuts=1;
      fEdgeXcuts[0]=4;
      fEdgeYcuts[0]=3;
      fNClCuts[0]=20;
  }
  //=============================//
  // Loop over all Tracklet Cuts //
  //=============================//
  for (Int_t icut=0; icut<fNcuts; icut++){
      AliDebug(4,Form("Processing cut %d for track with ID %d",icut,id));
      //cut parameters
      Double_t edgeCutX = fEdgeXcuts[icut];
      Double_t edgeCutY = fEdgeYcuts[icut];
      Int_t    nclCut   = fNClCuts[icut];
      //=========================//
      // Parameters to calculate //
      //=========================//
      //fit parameter inner, outer tracklet and combined fit
      TVectorD vecy1resInner(2),vecz1resInner(2); //pol1 fit parameters inner
      TVectorD vecy2resInner(3),vecz2resInner(3); //pol2 fit parameters inner
      //
      TVectorD vecy1resOuter(2),vecz1resOuter(2); //pol1 fit parameters outer
      TVectorD vecy2resOuter(3),vecz2resOuter(3); //pol2 fit parameters outer
      TVectorD vecy4res(5),vecz4res(5);
      // cluster and track positions for each row - used for residuals
      TVectorD vecX(159);        // x is the same for all (row center)
      TVectorD vecYkalman(159);  // y from kalman fit
      TVectorD vecZkalman(159);  // z from kalman fit
      TVectorD vecY1(159);       // y from pol1 fit per ROC
      TVectorD vecZ1(159);       // z from pol1 fit per ROC
      TVectorD vecY2(159);       // y from pol2 fit per ROC
      TVectorD vecZ2(159);       // z from pol2 fit per ROC
      TVectorD vecY4(159);       // y from sector fit
      TVectorD vecZ4(159);       // z from sector fit
      TVectorD vecClY(159);      // y cluster position
      TVectorD vecClZ(159);      // z cluster position
      TVectorD vecSec(159);      // sector for each row
      //chi2 of fits
      Double_t chi2I1z=-1;       // chi2 of pol1 fit in z (inner)
      Double_t chi2I1y=-1;       // chi2 of pol1 fit in y (inner)
      Double_t chi2O1z=-1;       // chi2 of pol1 fit in z (outer)
      Double_t chi2O1y=-1;       // chi2 of pol1 fit in y (outer)
      Double_t chi2I2z=-1;       // chi2 of pol2 fit in z (inner)
      Double_t chi2I2y=-1;       // chi2 of pol2 fit in y (inner)
      Double_t chi2O2z=-1;       // chi2 of pol2 fit in z (outer)
      Double_t chi2O2y=-1;       // chi2 of pol2 fit in y (outer)
      Double_t chi2IOz=-1;       // chi2 of hyp4 fit in z (inner+outer)
      Double_t chi2IOy=-1;       // chi2 of hyp4 fit in y (inner+outer)
      //more
      Int_t innerSector = -1;    // number of inner sector
      Int_t outerSector = -1;    // number of outer sector
      Int_t nclI=0;              // number of clusters (inner)
      Int_t nclO=0;              // number of clusters (outer)
      //=================================================//
      // Perform the Kalman Fit using the Tracklet Class //
      //=================================================//
      AliTPCTracklet::SetEdgeCut(edgeCutX,edgeCutY);
      TObjArray tracklets=
	  AliTPCTracklet::CreateTracklets(track,AliTPCTracklet::kKalman,
					  kFALSE,nclCut,kMaxTracklets);
      // tracklet pointers
      AliTPCTracklet *trInner = (AliTPCTracklet*)tracklets.At(0);
      AliTPCTracklet *trOuter = (AliTPCTracklet*)tracklets.At(1);
      AliTPCTracklet *tr=0x0;
      AliTPCTracklet dummy;
      //continue if we didn't find a tracklet
      if ( !trInner && !trOuter ) continue;
      //================================================================================//
      // Swap Inner and Outer if necessary (inner sector is definde by smaller local x) //
      //================================================================================//
      if ( trInner && trOuter ){
	  if ( !trInner->GetInner() || !trOuter->GetInner() ) continue;
	  if ( trInner->GetInner()->GetX() > trOuter->GetInner()->GetX() ){
	      tr = trInner;
	      trInner=trOuter;
	      trOuter=tr;
	      AliDebug(5,Form("Swapped Sectors: %02d (%f) <-> %02d (%f)", trOuter->GetSector(), trOuter->GetInner()->GetX(), trInner->GetSector(), trInner->GetInner()->GetX()));
	  }
      } else {
	  if ( trInner ){
              if ( !trInner->GetInner() ) continue;
	      trOuter=&dummy;
	      if ( trInner->GetSector()>35 ){
		  trOuter=trInner;
                  trInner=&dummy;
	      }
	  } else { //trOuter
              if ( !trOuter->GetInner() ) continue;
              trInner=&dummy;
	      if ( trOuter->GetSector()<36 ){
		  trInner=trOuter;
		  trOuter=&dummy;
	      }
	  }
      }
      innerSector = trInner->GetSector();
      if ( innerSector>=0 ) AliDebug(5,Form("Found inner Sector %02d at X %.2f", innerSector, trInner->GetInner()->GetX()));
      outerSector = trOuter->GetSector();
      if ( outerSector>=0 ) AliDebug(5,Form("Found outer Sector %02d at X %.2f", outerSector, trOuter->GetInner()->GetX()));

      // array of clusters
      TClonesArray arrCl("AliTPCclusterMI",159);
      arrCl.ExpandCreateFast(159);
      //=======================================//
      // fill fitters with cluster information //
      //=======================================//
      AliDebug(3,"Filing Cluster Information");
      for (Int_t irow=158;irow>-1;--irow) {
	  AliTPCclusterMI *c=track->GetClusterPointer(irow);
	  AliTPCclusterMI &cl = (AliTPCclusterMI&) (*arrCl[irow]);
	  cl=dummyCl;
          //
          vecSec[irow]=-1;
	  if (!c) continue;
          //
	  Int_t roc = static_cast<Int_t>(c->GetDetector());
	  if ( roc!=innerSector && roc!=outerSector ) continue;
          vecSec[irow]=roc;
	  //store clusters in clones array
	  cl=*c;
	  //cluster position
	  vecX[irow]   = c->GetX();
	  vecClY[irow] = c->GetY();
	  vecClZ[irow] = c->GetZ();
          //
          Double_t x = vecX[irow]-133.4; //reference is between IROC and OROC
          Double_t y = vecClY[irow];
	  Double_t z = vecClZ[irow];
          //
	  Double_t x2[2]={x,x*x};   //linear and parabolic parameters
	  Double_t x4[4]={0,0,0,0}; //hyp4 parameters
	  if ( roc == innerSector ) {
	      x4[0]=1; //offset inner - outer sector
	      x4[1]=x; //slope parameter inner sector
	  } else {
	      x4[2]=x; //slope parameter outer sector
	  }
	  x4[3]=x*x;   //common parabolic shape
          //
	  if ( roc==innerSector ){
	      fy1I.AddPoint(x2,y);
	      fz1I.AddPoint(x2,z);
	      fy2I.AddPoint(x2,y);
	      fz2I.AddPoint(x2,z);
              ++nclI;
	  }
	  if ( roc==outerSector ){
	      fy1O.AddPoint(x2,y);
	      fz1O.AddPoint(x2,z);
	      fy2O.AddPoint(x2,y);
	      fz2O.AddPoint(x2,z);
              ++nclO;
	  }
	  fy4.AddPoint(x4,y);
	  fz4.AddPoint(x4,z);
      }
      //======================================//
      // Evaluate and retrieve fit parameters //
      //======================================//
      AliDebug(5,Form("Evaluating fits with inner (outer) Sec: %02d (%02d)",innerSector,outerSector));
      //inner sector
      if (  (innerSector>-1) && (fy1I.GetNpoints()>0) ){
	  fy1I.Eval();
	  fz1I.Eval();
	  fy2I.Eval();
	  fz2I.Eval();
	  fy1I.GetParameters(vecy1resInner);
	  fz1I.GetParameters(vecz1resInner);
	  fy2I.GetParameters(vecy2resInner);
	  fz2I.GetParameters(vecz2resInner);
	  chi2I1y=fy1I.GetChisquare()/(fy1I.GetNpoints()-2);
	  chi2I1z=fz1I.GetChisquare()/(fz1I.GetNpoints()-2);
	  chi2I2y=fy2I.GetChisquare()/(fy2I.GetNpoints()-3);
	  chi2I2z=fz2I.GetChisquare()/(fz2I.GetNpoints()-3);
      }
      //outer sector
      if (  (outerSector>-1) && (fy1O.GetNpoints()>0) ){
	  fy1O.Eval();
	  fz1O.Eval();
	  fy2O.Eval();
	  fz2O.Eval();
	  fy1O.GetParameters(vecy1resOuter);
	  fz1O.GetParameters(vecz1resOuter);
	  fy2O.GetParameters(vecy2resOuter);
	  fz2O.GetParameters(vecz2resOuter);
	  chi2O1y=fy1O.GetChisquare()/(fy1O.GetNpoints()-2);
	  chi2O1z=fz1O.GetChisquare()/(fz1O.GetNpoints()-2);
	  chi2O2y=fy2O.GetChisquare()/(fy2O.GetNpoints()-3);
	  chi2O2z=fz2O.GetChisquare()/(fz2O.GetNpoints()-3);
      }
      //combined hyp4 fit
      if ( innerSector>0 && outerSector>0 ){
	  if (fy4.GetNpoints()>0) {
	      fy4.Eval();
	      fy4.GetParameters(vecy4res);
	      chi2IOy=fy4.GetChisquare()/(fy4.GetNpoints()-5);
	  }
	  if (fz4.GetNpoints()>0) {
	      fz4.Eval();
	      fz4.GetParameters(vecz4res);
	      chi2IOz=fz4.GetChisquare()/(fz4.GetNpoints()-5);
	  }
      }
      //clear points
      fy4.ClearPoints();  fz4.ClearPoints();
      fy1I.ClearPoints(); fy1O.ClearPoints();
      fz1I.ClearPoints(); fz1O.ClearPoints();
      fy2I.ClearPoints(); fy2O.ClearPoints();
      fz2I.ClearPoints(); fz2O.ClearPoints();
      //==============================//
      // calculate tracklet positions //
      //==============================//
      AliDebug(4,"Calculate tracklet positions");
      for (Int_t irow=158;irow>-1;--irow) {
	  if ( vecSec[irow]==-1 ) continue;  //no cluster info
          if ( vecSec[irow]!=innerSector && vecSec[irow]!=outerSector ) continue;
	  tr=&dummy;
	  Double_t x=vecX[irow];
	  Double_t xref=x-133.4;
          //
          Double_t yoffInner=0;
          Double_t zoffInner=0;
          Double_t yslopeInner=0;
          Double_t yslopeOuter=0;
          Double_t zslopeInner=0;
	  Double_t zslopeOuter=0;
          //positions of hyperplane fits
	  if ( vecSec[irow] == outerSector ) {
	      tr=trOuter;
              vecY1[irow]=vecy1resOuter[0]+vecy1resOuter[1]*xref;
              vecZ1[irow]=vecz1resOuter[0]+vecz1resOuter[1]*xref;
              vecY2[irow]=vecy2resOuter[0]+vecy2resOuter[1]*xref+vecy2resOuter[2]*xref*xref;
	      vecZ2[irow]=vecz2resOuter[0]+vecz2resOuter[1]*xref+vecz2resOuter[2]*xref*xref;
              yslopeOuter=vecy4res[3];
	      zslopeOuter=vecz4res[3];
	  } else {
	      tr=trInner;
              vecY1[irow]=vecy1resInner[0]+vecy1resInner[1]*xref;
              vecZ1[irow]=vecz1resInner[0]+vecz1resInner[1]*xref;
              vecY2[irow]=vecy2resInner[0]+vecy2resInner[1]*xref+vecy2resInner[2]*xref*xref;
              vecZ2[irow]=vecz2resInner[0]+vecz2resInner[1]*xref+vecz2resInner[2]*xref*xref;
              yoffInner=vecy4res[1];
	      zoffInner=vecz4res[1];
              yslopeInner=vecy4res[2];
	      zslopeInner=vecz4res[2];
	  }
	  vecY4[irow]=vecy4res[0]+yoffInner+yslopeInner*xref+yslopeOuter*xref+vecy4res[4]*xref*xref;
	  vecZ4[irow]=vecz4res[0]+zoffInner+zslopeInner*xref+zslopeOuter*xref+vecz4res[4]*xref*xref;
          //positions of kalman fits
	  Double_t gxyz[3],xyz[3];
	  AliExternalTrackParam *param = 0x0;
          //
	  param=tr->GetInner();
	  if (param){
	      param->GetXYZ(gxyz);
	      Float_t bz = AliTracker::GetBz(gxyz);
	      param->GetYAt(x, bz, xyz[1]);
	      param->GetZAt(x, bz, xyz[2]);
	      vecYkalman[irow]=xyz[1];
	      vecZkalman[irow]=xyz[2];
	  }
      }
      //=====================================================================//
      // write results from the different tracklet fits with debug streamers //
      //=====================================================================//
      if (fStreamLevel>4){
	  TTreeSRedirector *cstream = GetDebugStreamer();
	  if (cstream){
	      Float_t dedx = track->GetdEdx();
	      (*cstream)<<"FitModels"<<
		  "cutNr="      << icut <<
                  "edgeCutX="   << edgeCutX <<
		  "edgeCutY="   << edgeCutY <<
		  "nclCut="     << nclCut <<
                  "innerSector="<< innerSector <<
		  "outerSector="<< outerSector <<
                  "dEdx="       << dedx <<
		  "LTr.="       << ltrp <<
		  "Tr.="        << extparam <<
                  "yPol1In.="   << &vecy1resInner <<
                  "zPol1In.="   << &vecz1resInner <<
                  "yPol2In.="   << &vecy2resInner <<
                  "zPol2In.="   << &vecz2resInner <<
                  "yPol1Out.="  << &vecy1resOuter <<
                  "zPol1Out.="  << &vecz1resOuter <<
                  "yPol2Out.="  << &vecy2resOuter <<
                  "zPol2Out.="  << &vecz2resOuter <<
		  "yInOut.="    << &vecy4res <<
		  "zInOut.="    << &vecz4res <<
                  "chi2y1In="   << chi2I1y <<
                  "chi2z1In="   << chi2I1z <<
                  "chi2y1Out="  << chi2O1y <<
                  "chi2z1Out="  << chi2O1z <<
                  "chi2y2In="   << chi2I2y <<
                  "chi2z2In="   << chi2I2z <<
                  "chi2y2Out="  << chi2O2y <<
                  "chi2z2Out="  << chi2O2z <<
                  "chi2yInOut=" << chi2IOy <<
                  "chi2zInOut=" << chi2IOz <<
		  "trletIn.="   << trInner <<
		  "trletOut.="  << trOuter <<
		  "nclI="       << nclI <<
                  "nclO="       << nclO <<
		  "\n";
	  }
      }

      // wirte residuals information
      if (fStreamLevel>5){
	  TTreeSRedirector *cstream = GetDebugStreamer();
	  if (cstream){
	      Float_t dedx = track->GetdEdx();
	      (*cstream)<<"Residuals"<<
		  "cutNr="      << icut <<
                  "edgeCutX="   << edgeCutX <<
		  "edgeCutY="   << edgeCutY <<
		  "nclCut="     << nclCut   <<
		  "LTr.="       << ltrp <<
		  "Tr.="        << extparam<<
		  "dEdx="       << dedx <<
		  "Cl.="        << &arrCl <<
		  "TrX.="       << &vecX <<
		  "TrYpol1.="   << &vecY1 <<
		  "TrZpol1.="   << &vecZ1 <<
		  "TrYpol2.="   << &vecY2 <<
		  "TrZpol2.="   << &vecZ2 <<
		  "TrYInOut.="  << &vecY4 <<
		  "TrZInOut.="  << &vecZ4 <<
		  "ClY.="       << &vecClY <<
		  "ClZ.="       << &vecClZ <<
                  "sec.="       << &vecSec <<
		  "nclI="       << nclI <<
                  "nclO="       << nclO <<
		  "yInOut.="    << &vecy4res <<
		  "zInOut.="    << &vecz4res <<
                  "chi2y1In="   << chi2I1y <<
                  "chi2z1In="   << chi2I1z <<
                  "chi2y1Out="  << chi2O1y <<
                  "chi2z1Out="  << chi2O1z <<
                  "chi2y2In="   << chi2I2y <<
                  "chi2z2In="   << chi2I2z <<
                  "chi2y2Out="  << chi2O2y <<
                  "chi2z2Out="  << chi2O2z <<
                  "chi2yInOut=" << chi2IOy <<
                  "chi2zInOut=" << chi2IOz <<
		  "\n";

	  }
      }
      //==========================//
      // Fill Residual Histograms //
      //==========================//
      TProfile *profy = (TProfile*)fDeltaYres.UncheckedAt(id);
      TProfile *profz = (TProfile*)fDeltaZres.UncheckedAt(id);
      if (!profy){
	profy=new TProfile(Form("pry%03d",id),Form("Y Residuals for Laser Beam %03d",id),115,80,250);
	profy->SetDirectory(0);
	fDeltaYres.AddAt(profy,id);
      }
      if (!profz){
	profz=new TProfile(Form("prz%03d",id),Form("Z Residuals for Laser Beam %03d",id),115,80,250);
	profz->SetDirectory(0);
	fDeltaZres.AddAt(profz,id);
      }
      for (Int_t irow=158;irow>-1;--irow) {
	if (vecSec[irow]==-1)continue; //no cluster info
	Double_t x    = vecX[irow];
	Double_t ycl  = vecClY[irow];
	Double_t yfit = vecY1[irow];
	Double_t zcl  = vecClZ[irow];
	Double_t zfit = vecZ1[irow];
	if (profy) 
	  if (profy->GetEntries()<1000000) 
	    profy->Fill(x,yfit-ycl);
	if (profz) 
	  if (profz->GetEntries()<1000000) 
	    profz->Fill(x,zfit-zcl);
      }
      //===============================//
      // Fill Fit Parameter Histograms //
      //===============================//
      TH1F *pol2InnerY =(TH1F*)fPol2Par2InY.UncheckedAt(id);
      TH1F *diff1InnerY=(TH1F*)fDiffPar1InY.UncheckedAt(id);
      TH1F *pol2OuterY =(TH1F*)fPol2Par2OutY.UncheckedAt(id);
      TH1F *diff1OuterY=(TH1F*)fDiffPar1OutY.UncheckedAt(id);
      TH1F *pol2InnerZ =(TH1F*)fPol2Par2InZ.UncheckedAt(id);
      TH1F *diff1InnerZ=(TH1F*)fDiffPar1InZ.UncheckedAt(id);
      TH1F *pol2OuterZ =(TH1F*)fPol2Par2OutZ.UncheckedAt(id);
      TH1F *diff1OuterZ=(TH1F*)fDiffPar1OutZ.UncheckedAt(id);
      //create histograms if the do not already exist
      if (!pol2InnerY){
	  pol2InnerY =new TH1F(Form("pol2par2inY%03d",id),
			       Form("2nd derivative from pol2 fit in Y for Laser beam %03d (inner sector)",id),
			       500,-.005,.005);
	  pol2InnerY->SetDirectory(0);
	  pol2OuterY =new TH1F(Form("pol2par2outY%03d",id),
			       Form("2nd derivative from pol2 fit in Y for Laser beam %03d (outer sector)",id),
			       500,0.01,.01);
	  pol2OuterY->SetDirectory(0);

	  diff1InnerY=new TH1F(Form("diff1inY%03d",id),
			       Form("diff of 1st derivative from pol1 and pol2 in Y fit for Laser beam %03d (inner sector)",id),
			       500,-.5,.5);
	  diff1InnerY->SetDirectory(0);

	  diff1OuterY=new TH1F(Form("diff1outY%03d",id),
			       Form("diff of 1st derivative from pol1 and pol2 in Yfit for Laser beam %03d (outer sector)",id),
			       500,-1,1);
	  diff1OuterY->SetDirectory(0);

	  pol2InnerZ =new TH1F(Form("pol2par2inZ%03d",id),
			       Form("2nd derivative from pol2 fit in Z for Laser beam %03d (inner sector)",id),
			       500,-.002,.002);
	  pol2InnerZ->SetDirectory(0);

	  pol2OuterZ =new TH1F(Form("pol2par2outZ%03d",id),
			       Form("2nd derivative from pol2 fit in Z for Laser beam %03d (outer sector)",id),
			       500,-.005,.005);
	  pol2OuterZ->SetDirectory(0);
	  diff1InnerZ=new TH1F(Form("diff1inZ%03d",id),
			       Form("diff of 1st derivative from pol1 and pol2 in Z fit for Laser beam %03d (inner sector)",id),
			       500,-.02,.02);
	  diff1InnerZ->SetDirectory(0);
	  diff1OuterZ=new TH1F(Form("diff1outZ%03d",id),
			       Form("diff of 1st derivative from pol1 and pol2 in Zfit for Laser beam %03d (outer sector)",id),
			       500,-.03,.03);
	  diff1OuterZ->SetDirectory(0);
          //add
	  fPol2Par2InY.AddAt(pol2InnerY,id);
	  fDiffPar1InY.AddAt(diff1InnerY,id);
	  fPol2Par2OutY.AddAt(pol2OuterY,id);
	  fDiffPar1OutY.AddAt(diff1OuterY,id);
	  fPol2Par2InZ.AddAt(pol2InnerZ,id);
	  fDiffPar1InZ.AddAt(diff1InnerZ,id);
	  fPol2Par2OutZ.AddAt(pol2OuterZ,id);
	  fDiffPar1OutZ.AddAt(diff1OuterZ,id);
      }
      //fill histograms
      pol2InnerY ->Fill(vecy2resInner[2]);
      pol2OuterY ->Fill(vecy2resOuter[2]);
      diff1InnerY->Fill(vecy2resInner[1]-vecy1resInner[1]);
      diff1OuterY->Fill(vecy2resOuter[1]-vecy1resOuter[1]);
      pol2InnerZ ->Fill(vecz2resInner[2]);
      pol2OuterZ ->Fill(vecz2resOuter[2]);
      diff1InnerZ->Fill(vecz2resInner[1]-vecz1resInner[1]);
      diff1OuterZ->Fill(vecz2resOuter[1]-vecz1resOuter[1]);



  }
}


void AliTPCcalibLaser::RefitLaser(Int_t id){
  //
  // Refit the track store residuals
  //

  AliTPCseed *track    = (AliTPCseed*)fTracksTPC.At(id);
  AliExternalTrackParam *param=(AliExternalTrackParam*)fTracksEsdParam.At(id);
  AliTPCLaserTrack *ltrp = (AliTPCLaserTrack*)fTracksMirror.At(id);

  //linear fit model in y and z per sector
  static TLinearFitter fy1(2,"hyp1");
  static TLinearFitter fz1(2,"hyp1");
  //quadratic fit model in y and z per sector
  static TLinearFitter fy2(3,"hyp2");
  static TLinearFitter fz2(3,"hyp2");
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
      Double_t xd = c->GetX()-133.4; // reference x is beteen iroc and oroc
      if (c&&c->GetDetector()==isec) {
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


void AliTPCcalibLaser::DumpMeanInfo(Float_t bfield,Int_t minEntries){
  //
  //  Dump information about laser beams
  //  isOK variable indicates usability of the beam  
  //  Beam is not usable if:
  //     a.  No entries in range (krmsCut0)
  //     b.  Big sperad          (krmscut1)
  //     c.  RMSto fit sigma bigger then (kmultiCut)
  //     d.  Too big angular spread 
  //  

  const Float_t krmsCut0=0.001;
  const Float_t krmsCut1=0.16;
  const Float_t kmultiCut=2;
  const Float_t kcutP0=0.002;
  //
  AliTPCcalibLaser *laser = this;
  TTreeSRedirector *pcstream = new TTreeSRedirector("laserMean.root");
  TF1 fg("fg","gaus");
  //
  //
  for (Int_t id=0; id<336; id++){
    Bool_t isOK=kTRUE;
    TH1F * hisphi  = (TH1F*)laser->fDeltaPhi.At(id);
    TH1F * hisphiP = (TH1F*)laser->fDeltaPhiP.At(id);
    TH1F * hisZ    = (TH1F*)laser->fDeltaZ.At(id);
    TH1F * hisS    = (TH1F*)laser->fSignals.At(id);
    if (!hisphi) continue;;
    Double_t entries = hisphi->GetEntries();
    if (entries<minEntries) continue;
    //
    AliTPCLaserTrack *ltrp = (AliTPCLaserTrack*)fTracksMirror.At(id);
    if (!ltrp) {
     AliTPCLaserTrack::LoadTracks();
      ltrp =(AliTPCLaserTrack*)AliTPCLaserTrack::GetTracks()->UncheckedAt(id);
    }
    Float_t meanphi = hisphi->GetMean();
    Float_t rmsphi = hisphi->GetRMS();
    //
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
    Float_t meanS=hisS->GetMean();
    //
    Double_t lxyz[3];
    Double_t lpxyz[3];
    ltrp->GetXYZ(lxyz);
    ltrp->GetPxPyPz(lpxyz);

    if (rmsphi<krmsCut0) isOK=kFALSE; // empty in range - not entries inside
    if (rmsphi>krmsCut1) isOK=kFALSE; // empty in range - not entries inside
    if (rmsphi>krmsCut0) if (gphi2/rmsphi>kmultiCut) isOK=kFALSE;   // multi peak structure
    if (gphiP2>kcutP0) isOK=kFALSE;
    //
    (*pcstream)<<"Mean"<<
      "isOK="<<isOK<<
      "entries="<<entries<<      // number of entries
      "bz="<<bfield<<            // bfield
      "LTr.="<<ltrp<<             // refernece track
      //
      "lx0="<<lxyz[0]<<          // reference x
      "lx1="<<lxyz[1]<<          // reference y
      "lx2="<<lxyz[2]<<          // refernece z
      "lpx0="<<lpxyz[0]<<          // reference x
      "lpx1="<<lpxyz[1]<<          // reference y
      "lpx2="<<lpxyz[2]<<          // refernece z
      //
      "msig="<<meanS<<
      //
      "mphi="<<meanphi<<         //
      "rmsphi="<<rmsphi<<        //
      "gphi1="<<gphi1<<
      "gphi2="<<gphi2<<
      //
      "mphiP="<<meanphiP<<       //
      "rmsphiP="<<rmsphiP<<      //
      "gphiP1="<<gphiP1<<
      "gphiP2="<<gphiP2<<
      //
      "meanZ="<<meanZ<<
      "rmsZ="<<rmsZ<<
      "gz1="<<gz1<<
      "gz2="<<gz2<<

      "\n";
  }
  delete pcstream;
}



void AliTPCcalibLaser::DumpScanInfo(TTree * chain){
  //
  //
  //
  TTreeSRedirector *pcstream = new TTreeSRedirector("laserScan.root");
  TFile * f = pcstream->GetFile();
  f->mkdir("dirphi");
  f->mkdir("dirphiP");
  f->mkdir("dirZ");
  TF1 fp("p1","pol1");
  //
  //
  char    cut[1000];
  char grname[1000];
  char grnamefull[1000];

  Double_t mphi[100];
  Double_t mphiP[100];
  Double_t smphi[100];
  Double_t smphiP[100];
  Double_t mZ[100];
  Double_t smZ[100];
  Double_t bz[100];
  Double_t sbz[100];
  // fit parameters
  Double_t pphi[3];
  Double_t pphiP[3];
  Double_t pmZ[3];
  //
  for (Int_t id=0; id<336; id++){
    // id =205;
    sprintf(cut,"isOK&&fId==%d",id);
    Int_t entries = chain->Draw("bz",cut,"goff");
    if (entries<3) continue;
    AliTPCLaserTrack *ltrp = 0;;
    if (!AliTPCLaserTrack::GetTracks()) AliTPCLaserTrack::LoadTracks();
    ltrp =(AliTPCLaserTrack*)AliTPCLaserTrack::GetTracks()->UncheckedAt(id);
    Double_t lxyz[3];
    Double_t lpxyz[3];
    ltrp->GetXYZ(lxyz);
    ltrp->GetPxPyPz(lpxyz);

    chain->Draw("bz",cut,"goff");
    memcpy(bz, chain->GetV1(), entries*sizeof(Double_t));
    chain->Draw("0.01*abs(bz)+0.02",cut,"goff");
    memcpy(sbz, chain->GetV1(), entries*sizeof(Double_t));
    //
    chain->Draw("gphi1",cut,"goff");
    memcpy(mphi, chain->GetV1(), entries*sizeof(Double_t));
    chain->Draw("0.05*abs(mphi)+gphi2",cut,"goff");
    memcpy(smphi, chain->GetV1(), entries*sizeof(Double_t));
    //
    chain->Draw("gphiP1",cut,"goff");
    memcpy(mphiP, chain->GetV1(), entries*sizeof(Double_t));
    chain->Draw("0.05*abs(mphiP)+gphiP2",cut,"goff");
    memcpy(smphiP, chain->GetV1(), entries*sizeof(Double_t));
    //
    chain->Draw("gz1",cut,"goff");
    memcpy(mZ, chain->GetV1(), entries*sizeof(Double_t));
    chain->Draw("0.01*abs(meanZ)+gz2",cut,"goff");
    memcpy(smZ, chain->GetV1(), entries*sizeof(Double_t));
    //
    //
    sprintf(grnamefull,"Side_%d_Bundle_%d_Rod_%d_Beam_%d",
	    ltrp->GetSide(),  ltrp->GetBundle(), ltrp->GetRod(), ltrp->GetBeam());
    // store data  
    // phi
    f->cd("dirphi");
    TGraphErrors *grphi = new TGraphErrors(entries,bz,mphi,sbz,smphi);
    grphi->Draw("a*");
    grphi->Fit(&fp);
    pphi[0] = fp.GetParameter(0);                          // offset
    pphi[1] = fp.GetParameter(1);                          // slope
    pphi[2] = TMath::Sqrt(fp.GetChisquare()/(entries-2.));  // normalized chi2
    sprintf(grname,"phi_id%d",id);
    grphi->SetName(grname);  grphi->SetTitle(grnamefull);
    grphi->GetXaxis()->SetTitle("b_{z} (T)");
    grphi->GetYaxis()->SetTitle("#Delta r#phi (cm)");
    grphi->SetMaximum(1.2);
    grphi->SetMinimum(-1.2);
    grphi->Draw("a*");

    grphi->Write();
    gPad->SaveAs(Form("pic/phi/phi_%s.gif",grnamefull));
    // phiP
    f->cd("dirphiP)");
    TGraphErrors *grphiP = new TGraphErrors(entries,bz,mphiP,sbz,smphiP);
    grphiP->Draw("a*");
    grphiP->Fit(&fp);
    pphiP[0] = fp.GetParameter(0);                          // offset
    pphiP[1] = fp.GetParameter(1);                          // slope
    pphiP[2] = TMath::Sqrt(fp.GetChisquare()/(entries-2.));  // normalized chi2
    sprintf(grname,"phiP_id%d",id);
    grphiP->SetName(grname);  grphiP->SetTitle(grnamefull);
    grphiP->GetXaxis()->SetTitle("b_{z} (T)");
    grphiP->GetYaxis()->SetTitle("#Delta #phi (rad)");
    grphiP->SetMaximum(pphiP[0]+0.005);
    grphiP->SetMinimum(pphiP[0]-0.005);

    gPad->SaveAs(Form("pic/phiP/phiP_%s.gif",grnamefull));
    grphiP->Write();
    //
    //Z
    f->cd("dirZ");
    TGraphErrors *grmZ = new TGraphErrors(entries,bz,mZ,sbz,smZ);
    grmZ->Draw("a*");
    grmZ->Fit(&fp);
    pmZ[0] = fp.GetParameter(0);                          // offset
    pmZ[1] = fp.GetParameter(1);                          // slope
    pmZ[2] = TMath::Sqrt(fp.GetChisquare()/(entries-2.));  // normalized chi2
    sprintf(grname,"mZ_id%d",id);
    grmZ->SetName(grname);  grmZ->SetTitle(grnamefull);
    grmZ->GetXaxis()->SetTitle("b_{z} (T)");
    grmZ->GetYaxis()->SetTitle("#Delta z (cm)");

    gPad->SaveAs(Form("pic/z/z_%s.gif",grnamefull));
    grmZ->Write();
    

    for (Int_t ientry=0; ientry<entries; ientry++){
      (*pcstream)<<"Mean"<<
	"id="<<id<<
	"LTr.="<<ltrp<<
	"entries="<<entries<<
	"bz="<<bz[ientry]<<
	"lx0="<<lxyz[0]<<          // reference x
	"lx1="<<lxyz[1]<<          // reference y
	"lx2="<<lxyz[2]<<          // refernece z      
	"lpx0="<<lpxyz[0]<<          // reference x
	"lpx1="<<lpxyz[1]<<          // reference y
	"lpx2="<<lpxyz[2]<<          // refernece z            
	//values
	"gphi1="<<mphi[ientry]<< // mean - from gaus fit
	"pphi0="<<pphi[0]<<   // offset
	"pphi1="<<pphi[1]<<   // mean
	"pphi2="<<pphi[2]<<   // norm chi2
	//
	"gphiP1="<<mphiP[ientry]<< // mean - from gaus fit
	"pphiP0="<<pphiP[0]<< // offset
	"pphiP1="<<pphiP[1]<< // mean
	"pphiP2="<<pphiP[2]<< // norm chi2
	//
	"gz1="<<mZ[ientry]<<
	"pmZ0="<<pmZ[0]<<     // offset
	"pmZ1="<<pmZ[1]<<     // mean
	"pmZ2="<<pmZ[2]<<     // norm chi2
	"\n";
    }
  }
  
  delete pcstream;
  
}


void AliTPCcalibLaser::Analyze(){
  //
  //
  //
}


Long64_t AliTPCcalibLaser::Merge(TCollection *li) {

  TIterator* iter = li->MakeIterator();
  AliTPCcalibLaser* cal = 0;

  while ((cal = (AliTPCcalibLaser*)iter->Next())) {
    if (!cal->InheritsFrom(AliTPCcalibLaser::Class())) {
      Error("Merge","Attempt to add object of class %s to a %s", cal->ClassName(), this->ClassName());
      return -1;
    }
    //
   //  fHistNTracks->Add(cal->fHistNTracks);
//     fClusters->Add(cal->fClusters);
//     fModules->Add(cal->fModules);
//     fHistPt->Add(cal->fHistPt);
//     fPtResolution->Add(cal->fPtResolution);
//     fDeDx->Add(cal->fDeDx);


    TH1F *h=0x0;
    TH1F *hm=0x0;
    TProfile *hp=0x0;
    TProfile *hpm=0x0;

    for (Int_t id=0; id<336; id++){
      // merge fDeltaZ histograms
      hm = (TH1F*)cal->fDeltaZ.At(id);
      h  = (TH1F*)fDeltaZ.At(id);
      if (!h) {
	h=new TH1F(Form("hisdz%d",id),Form("hisdz%d",id),1000,-10,10);
	h->SetDirectory(0);
	fDeltaZ.AddAt(h,id);
      }
      if (hm) h->Add(hm);
      // merge fDeltaPhi histograms
      hm = (TH1F*)cal->fDeltaPhi.At(id);
      h  = (TH1F*)fDeltaPhi.At(id);
      if (!h) {
	h= new TH1F(Form("hisdphi%d",id),Form("hisdphi%d",id),1000,-1,1);
	h->SetDirectory(0);
	fDeltaPhi.AddAt(h,id);
      }
      if (hm) h->Add(hm);
      // merge fDeltaPhiP histograms
      hm = (TH1F*)cal->fDeltaPhiP.At(id);
      h  = (TH1F*)fDeltaPhiP.At(id);
      if (!h) {
	h=new TH1F(Form("hisdphiP%d",id),Form("hisdphiP%d",id),1000,-0.01,0.01);
	h->SetDirectory(0);
	fDeltaPhiP.AddAt(h,id);
      }
      if (hm) h->Add(hm);
      // merge fSignals histograms
      hm = (TH1F*)cal->fSignals.At(id);
      h  = (TH1F*)fSignals.At(id);
      if (!h) {
	h=new TH1F(Form("hisSignal%d",id),Form("hisSignal%d",id),1000,0,1000);
	h->SetDirectory(0);
	fSignals.AddAt(h,id);
      }
      if (hm) h->Add(hm);
      //
      //
      // merge ProfileY histograms
      hpm = (TProfile*)cal->fDeltaYres.At(id);
      hp  = (TProfile*)fDeltaYres.At(id);
      if (!hp) {
	hp=new TProfile(Form("pry%03d",id),Form("Y Residuals for Laser Beam %03d",id),115,80,250);
	hp->SetDirectory(0);
	fDeltaYres.AddAt(hp,id);
      }
      if (hpm) hp->Add(hpm);
      //
      hpm = (TProfile*)cal->fDeltaZres.At(id);
      hp  = (TProfile*)fDeltaZres.At(id);
      if (!hp) {
	hp=new TProfile(Form("prz%03d",id),Form("Z Residuals for Laser Beam %03d",id),115,80,250);
	hp->SetDirectory(0);
	fDeltaZres.AddAt(hp,id);
      }
      if (hpm) hp->Add(hpm);
      //
      //
      //merge fit param histograms
      //local hists
      TH1F *pol2InnerY =(TH1F*)fPol2Par2InY.UncheckedAt(id);
      TH1F *diff1InnerY=(TH1F*)fDiffPar1InY.UncheckedAt(id);
      TH1F *pol2OuterY =(TH1F*)fPol2Par2OutY.UncheckedAt(id);
      TH1F *diff1OuterY=(TH1F*)fDiffPar1OutY.UncheckedAt(id);
      TH1F *pol2InnerZ =(TH1F*)fPol2Par2InZ.UncheckedAt(id);
      TH1F *diff1InnerZ=(TH1F*)fDiffPar1InZ.UncheckedAt(id);
      TH1F *pol2OuterZ =(TH1F*)fPol2Par2OutZ.UncheckedAt(id);
      TH1F *diff1OuterZ=(TH1F*)fDiffPar1OutZ.UncheckedAt(id);
      //hists to merge
      TH1F *pol2InnerYm =(TH1F*)cal->fPol2Par2InY.UncheckedAt(id);
      TH1F *diff1InnerYm=(TH1F*)cal->fDiffPar1InY.UncheckedAt(id);
      TH1F *pol2OuterYm =(TH1F*)cal->fPol2Par2OutY.UncheckedAt(id);
      TH1F *diff1OuterYm=(TH1F*)cal->fDiffPar1OutY.UncheckedAt(id);
      TH1F *pol2InnerZm =(TH1F*)cal->fPol2Par2InZ.UncheckedAt(id);
      TH1F *diff1InnerZm=(TH1F*)cal->fDiffPar1InZ.UncheckedAt(id);
      TH1F *pol2OuterZm =(TH1F*)cal->fPol2Par2OutZ.UncheckedAt(id);
      TH1F *diff1OuterZm=(TH1F*)cal->fDiffPar1OutZ.UncheckedAt(id);
      //create histos if they do not exist
      if (!pol2InnerY){
          pol2InnerY =new TH1F(Form("pol2par2inY%03d",id),
			       Form("2nd derivative from pol2 fit in Y for Laser beam %03d (inner sector)",id),
			       500,-.005,.005);
	  pol2InnerY->SetDirectory(0);

	  pol2OuterY =new TH1F(Form("pol2par2outY%03d",id),
			       Form("2nd derivative from pol2 fit in Y for Laser beam %03d (outer sector)",id),
			       500,0.01,.01);
	  pol2OuterY->SetDirectory(0);
	  diff1InnerY=new TH1F(Form("diff1inY%03d",id),
			       Form("diff of 1st derivative from pol1 and pol2 in Y fit for Laser beam %03d (inner sector)",id),
			       500,-.5,.5);
	  diff1OuterY=new TH1F(Form("diff1outY%03d",id),
			       Form("diff of 1st derivative from pol1 and pol2 in Yfit for Laser beam %03d (outer sector)",id),
			       500,-1,1);
	  diff1InnerY->SetDirectory(0);


	  pol2InnerZ =new TH1F(Form("pol2par2inZ%03d",id),
			       Form("2nd derivative from pol2 fit in Z for Laser beam %03d (inner sector)",id),
			       500,-.002,.002);
	  pol2InnerZ->SetDirectory(0);

	  pol2OuterZ =new TH1F(Form("pol2par2outZ%03d",id),
			       Form("2nd derivative from pol2 fit in Z for Laser beam %03d (outer sector)",id),
			       500,-.005,.005);
	  pol2OuterZ->SetDirectory(0);
	  diff1InnerZ=new TH1F(Form("diff1inZ%03d",id),
			       Form("diff of 1st derivative from pol1 and pol2 in Z fit for Laser beam %03d (inner sector)",id),
			       500,-.02,.02);
	  diff1InnerZ->SetDirectory(0);
	  diff1OuterZ=new TH1F(Form("diff1outZ%03d",id),
			       Form("diff of 1st derivative from pol1 and pol2 in Zfit for Laser beam %03d (outer sector)",id),
			       500,-.03,.03);
	  diff1OuterZ->SetDirectory(0);
      }
      pol2InnerYm =(TH1F*)cal->fPol2Par2InY.UncheckedAt(id);
      diff1InnerYm=(TH1F*)cal->fDiffPar1InY.UncheckedAt(id);
      pol2OuterYm =(TH1F*)cal->fPol2Par2OutY.UncheckedAt(id);
      diff1OuterYm=(TH1F*)cal->fDiffPar1OutY.UncheckedAt(id);
      pol2InnerZm =(TH1F*)cal->fPol2Par2InZ.UncheckedAt(id);
      diff1InnerZm=(TH1F*)cal->fDiffPar1InZ.UncheckedAt(id);
      pol2OuterZm =(TH1F*)cal->fPol2Par2OutZ.UncheckedAt(id);
      diff1OuterZm=(TH1F*)cal->fDiffPar1OutZ.UncheckedAt(id);
    }
  }
  return 0;
}




/*
 gSystem->Load("libSTAT.so")
 TStatToolkit toolkit;
 Double_t chi2;
 TVectorD fitParam;
 TMatrixD covMatrix;
 Int_t npoints;
 TCut cutA("entries>2&&pphi2<3&&abs(gphiP1-pphiP0)<0.003&&abs(gz1)<6");


TString fstring="";
//
fstring+="(abs(LTr.fP[1]/250)^3-1)*bz++";                               //1
fstring+="(abs(LTr.fP[1]/250)^3-1)*bz*LTr.fP[2]++";                     //2
fstring+="(abs(LTr.fP[1]/250)^1-1)*bz++";                               //3
fstring+="(abs(LTr.fP[1]/250)-1)*bz*LTr.fP[2]++";                       //4 
//
fstring+="(abs(LTr.fP[1]/250)^3-1)*bz*sin(atan2(lx1,lx0))++"            //5
fstring+="(abs(LTr.fP[1]/250)^3-1)*bz*sin(atan2(lx1,lx0))*LTr.fP[2]++"  //6
fstring+="(abs(LTr.fP[1]/250)-1)*bz*sin(atan2(lx1,lx0))++"              //7
fstring+="(abs(LTr.fP[1]/250)-1)*bz*sin(atan2(lx1,lx0))*LTr.fP[2]++"    //8
//   
fstring+="(abs(LTr.fP[1]/250)^3-1)*bz*cos(atan2(lx1,lx0))++"            //9
fstring+="(abs(LTr.fP[1]/250)^3-1)*bz*cos(atan2(lx1,lx0))*LTr.fP[2]++"  //10
fstring+="(abs(LTr.fP[1]/250)-1)*bz*cos(atan2(lx1,lx0))++"              //11
fstring+="(abs(LTr.fP[1]/250)-1)*bz*cos(atan2(lx1,lx0))*LTr.fP[2]++"    //12




 TString *strq0 = toolkit.FitPlane(treeT,"gphi1-pphi0",fstring->Data(), "fSide==1"+cutA, chi2,npoints,fitParam,covMatrix);

 treeT->SetAlias("fit",strq0->Data());
 

 TString *strqP = toolkit.FitPlane(treeT,"1000*(gphiP1-pphiP0)",fstring->Data(), "fSide==1"+cutA, chi2,npoints,fitParam,covMatrix);

 treeT->SetAlias("fitP",strqP->Data());


 TString *strqDrift = toolkit.FitPlane(treeT,"gz1","LTr.fP[1]++(1-2*(fSide==1))++lx1", cutA, chi2,npoints,fitParam,covMatrix);
 treeT->SetAlias("fitD",strqDrift->Data());


treeT->Draw("fit:LTr.fP[1]","abs(bz+0.4)<0.05"+cutA,""); 
{
for (Int_t i=0; i<6;i++){
treeT->SetLineColor(i+2);
treeT->SetMarkerSize(1);
treeT->SetMarkerStyle(22+i);
treeT->SetMarkerColor(i+2);

treeT->Draw("fit:LTr.fP[1]",Form("abs(bz+0.4)<0.05&fRod==%d",i)+cutA,"same"); 
}
} 
 */
 


/*
  TTree * tree = (TTree*)f.Get("FitModels");

  TEventList listLFit0("listLFit0","listLFit0");
  TEventList listLFit1("listLFit1","listLFit1");
  tree->Draw(">>listLFit0","seed.fdEdx<200&&seed.fdEdx>40");
  tree->SetEventList(&listLFit0);
  



  gSystem->Load("libSTAT.so")
  TStatToolkit toolkit;
  Double_t chi2;
  TVectorD fitParam;
  TMatrixD covMatrix;
  Int_t npoints;

  chain->SetAlias("dp","((Cl.fPad-int(Cl.fPad))*pi)");
  chain->SetAlias("dt","((Cl.fTimeBin-int(Cl.fTimeBin))*pi)");


  TString fstring="";  
  fstring+="cos(dp)++";
  fstring+="sin(dp)++";
  fstring+="cos(dt)++";
  fstring+="sin(dt)++";
  
  TString *str = toolkit.FitPlane(chain,"Cl.fZ-TrZInOut.fElements",fstring->Data(), "Cl.fDetector>35", chi2,npoints,fitParam,covMatrix,-1,0,200);


*/
