/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercialf purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

/* $Id: AliTRDresolution.cxx 27496 2008-07-22 08:35:45Z cblume $ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  TRD tracking resolution                                               //
//
// The class performs resolution and residual studies 
// of the TRD tracks for the following quantities :
//   - spatial position (y, [z])
//   - angular (phi) tracklet
//   - momentum at the track level
// 
// The class has to be used for regular detector performance checks using the official macros:
//   - $ALICE_ROOT/TRD/qaRec/run.C
//   - $ALICE_ROOT/TRD/qaRec/makeResults.C
// 
// For stand alone usage please refer to the following example: 
// {  
//   gSystem->Load("libANALYSIS.so");
//   gSystem->Load("libTRDqaRec.so");
//   AliTRDresolution *res = new AliTRDresolution();
//   //res->SetMCdata();
//   //res->SetVerbose();
//   //res->SetVisual();
//   res->Load();
//   if(!res->PostProcess()) return;
//   res->GetRefFigure(0);
// }  
//
//  Authors:                                                              //
//    Alexandru Bercuci <A.Bercuci@gsi.de>                                //
//    Markus Fasel <M.Fasel@gsi.de>                                       //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include <TSystem.h>

#include <TROOT.h>
#include <TObjArray.h>
#include <TH3.h>
#include <TH2.h>
#include <TH1.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TGaxis.h>
#include <TBox.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TMath.h>
#include <TMatrixT.h>
#include <TVectorT.h>
#include <TTreeStream.h>
#include <TGeoManager.h>
#include <TDatabasePDG.h>

#include "AliPID.h"
#include "AliLog.h"
#include "AliESDtrack.h"

#include "AliTRDresolution.h"
#include "AliTRDgeometry.h"
#include "AliTRDpadPlane.h"
#include "AliTRDcluster.h"
#include "AliTRDseedV1.h"
#include "AliTRDtrackV1.h"
#include "AliTRDReconstructor.h"
#include "AliTRDrecoParam.h"
#include "AliTRDpidUtil.h"

#include "info/AliTRDclusterInfo.h"

ClassImp(AliTRDresolution)

UChar_t const AliTRDresolution::fgNproj[kNviews] = {
  2, 2, 5, 5, 5,
  2, 5, 11, 11, 11
};
Char_t const * AliTRDresolution::fgPerformanceName[kNviews] = {
     "Charge"
    ,"Cluster2Track"
    ,"Tracklet2Track"
    ,"Tracklet2TRDin" 
    ,"Tracklet2TRDout" 
    ,"Cluster2MC"
    ,"Tracklet2MC"
    ,"TRDin2MC"
    ,"TRDout2MC"
    ,"TRD2MC"
};
// Configure segmentation for y resolution/residuals
Int_t const AliTRDresolution::fgkNresYsegm[3] = {
  AliTRDgeometry::kNsector
 ,AliTRDgeometry::kNsector*AliTRDgeometry::kNstack
 ,AliTRDgeometry::kNdet
};
Char_t const *AliTRDresolution::fgkResYsegmName[3] = {
  "Sector", "Stack", "Detector"};


//________________________________________________________
AliTRDresolution::AliTRDresolution()
  :AliTRDrecoTask()
  ,fStatus(0)
  ,fSegmentLevel(0)
  ,fIdxPlot(0)
  ,fIdxFrame(0)
  ,fPtThreshold(1.)
  ,fReconstructor(NULL)
  ,fGeo(NULL)
  ,fDBPDG(NULL)
  ,fGraphS(NULL)
  ,fGraphM(NULL)
  ,fCl(NULL)
  ,fMCcl(NULL)
/*  ,fTrklt(NULL)
  ,fMCtrklt(NULL)*/
{
  //
  // Default constructor
  //
  SetNameTitle("TRDresolution", "TRD spatial and momentum resolution");
  SetSegmentationLevel();
}

//________________________________________________________
AliTRDresolution::AliTRDresolution(char* name)
  :AliTRDrecoTask(name, "TRD spatial and momentum resolution")
  ,fStatus(0)
  ,fSegmentLevel(0)
  ,fIdxPlot(0)
  ,fIdxFrame(0)
  ,fPtThreshold(1.)
  ,fReconstructor(NULL)
  ,fGeo(NULL)
  ,fDBPDG(NULL)
  ,fGraphS(NULL)
  ,fGraphM(NULL)
  ,fCl(NULL)
  ,fMCcl(NULL)
/*  ,fTrklt(NULL)
  ,fMCtrklt(NULL)*/
{
  //
  // Default constructor
  //

  fReconstructor = new AliTRDReconstructor();
  fReconstructor->SetRecoParam(AliTRDrecoParam::GetLowFluxParam());
  fGeo = new AliTRDgeometry();

  InitFunctorList();
  SetSegmentationLevel();

  DefineOutput(kClToTrk, TObjArray::Class()); // cluster2track
  DefineOutput(kClToMC, TObjArray::Class()); // cluster2mc
/*  DefineOutput(kTrkltToTrk, TObjArray::Class()); // tracklet2track
  DefineOutput(kTrkltToMC, TObjArray::Class()); // tracklet2mc*/
}

//________________________________________________________
AliTRDresolution::~AliTRDresolution()
{
  //
  // Destructor
  //

  if(fGraphS){fGraphS->Delete(); delete fGraphS;}
  if(fGraphM){fGraphM->Delete(); delete fGraphM;}
  delete fGeo;
  delete fReconstructor;
  if(gGeoManager) delete gGeoManager;
  if(fCl){fCl->Delete(); delete fCl;}
  if(fMCcl){fMCcl->Delete(); delete fMCcl;}
/*  if(fTrklt){fTrklt->Delete(); delete fTrklt;}
  if(fMCtrklt){fMCtrklt->Delete(); delete fMCtrklt;}*/
}


//________________________________________________________
void AliTRDresolution::UserCreateOutputObjects()
{
  // spatial resolution

  if(!fReconstructor){
    fReconstructor = new AliTRDReconstructor();
    fReconstructor->SetRecoParam(AliTRDrecoParam::GetLowFluxParam());
  }
  if(!fGeo) fGeo = new AliTRDgeometry();

  if(!HasFunctorList()) InitFunctorList();
  fContainer = Histos();
  InitExchangeContainers();
}

//________________________________________________________
void AliTRDresolution::InitExchangeContainers()
{
  fCl = new TObjArray();
  fCl->SetOwner(kTRUE);
  fMCcl = new TObjArray();
  fMCcl->SetOwner(kTRUE);
/*  fTrklt = new TObjArray();
  fTrklt->SetOwner(kTRUE);
  fMCtrklt = new TObjArray();
  fMCtrklt->SetOwner(kTRUE);*/
}

//________________________________________________________
void AliTRDresolution::UserExec(Option_t *opt)
{
  //
  // Execution part
  //

  fCl->Delete();
  fMCcl->Delete();
/*  fTrklt->Delete();
  fMCtrklt->Delete();*/
  AliTRDrecoTask::UserExec(opt);
  PostData(kClToTrk, fCl);
  PostData(kClToMC, fMCcl);
/*  PostData(kTrkltToTrk, fTrklt);
  PostData(kTrkltToMC, fMCtrklt);*/
}

//________________________________________________________
Bool_t AliTRDresolution::Pulls(Double_t dyz[2], Double_t cov[3], Double_t tilt)
{
// Helper function to calculate pulls in the yz plane 
// using proper tilt rotation
// Uses functionality defined by AliTRDseedV1.

  Double_t t2(tilt*tilt);

  // rotate along pad
  Double_t cc[3];
  cc[0] = cov[0] - 2.*tilt*cov[1] + t2*cov[2]; 
  cc[1] = cov[1]*(1.-t2) + tilt*(cov[0] - cov[2]);
  cc[2] = t2*cov[0] + 2.*tilt*cov[1] + cov[2];
  // do sqrt
  Double_t sqr[3]={0., 0., 0.}; 
  if(AliTRDseedV1::GetCovSqrt(cc, sqr)) return kFALSE;
  Double_t invsqr[3]={0., 0., 0.}; 
  if(AliTRDseedV1::GetCovInv(sqr, invsqr)<1.e-40) return kFALSE;
  Double_t tmp(dyz[0]);
  dyz[0] = invsqr[0]*tmp + invsqr[1]*dyz[1];
  dyz[1] = invsqr[1]*tmp + invsqr[2]*dyz[1];
  return kTRUE;
}

//________________________________________________________
TH1* AliTRDresolution::PlotCharge(const AliTRDtrackV1 *track)
{
  //
  // Plots the charge distribution
  //

  if(track) fkTrack = track;
  if(!fkTrack){
    AliDebug(4, "No Track defined.");
    return NULL;
  }
  TObjArray *arr = NULL;
  if(!fContainer || !(arr = ((TObjArray*)fContainer->At(kCharge)))){
    AliWarning("No output container defined.");
    return NULL;
  }
  TH3S* h = NULL;

  AliTRDseedV1 *fTracklet = NULL;  
  AliTRDcluster *c = NULL;
  for(Int_t ily=0; ily<AliTRDgeometry::kNlayer; ily++){
    if(!(fTracklet = fkTrack->GetTracklet(ily))) continue;
    if(!fTracklet->IsOK()) continue;
    Float_t x0 = fTracklet->GetX0();
    Float_t dq, dl;
    for(Int_t itb=AliTRDseedV1::kNtb; itb--;){
      if(!(c = fTracklet->GetClusters(itb))){ 
        if(!(c = fTracklet->GetClusters(AliTRDseedV1::kNtb+itb))) continue;
      }
      dq = fTracklet->GetdQdl(itb, &dl);
      dl /= 0.15; // dl/dl0, dl0 = 1.5 mm for nominal vd
      (h = (TH3S*)arr->At(0))->Fill(dl, x0-c->GetX(), dq);
    }

//     if(!HasMCdata()) continue;
//     UChar_t s;
//     Float_t pt0, y0, z0, dydx0, dzdx0;
//     if(!fMC->GetDirections(x0, y0, z0, dydx0, dzdx0, pt0, s)) continue;

  }
  return h;
}


//________________________________________________________
TH1* AliTRDresolution::PlotCluster(const AliTRDtrackV1 *track)
{
  //
  // Plot the cluster distributions
  //

  if(track) fkTrack = track;
  if(!fkTrack){
    AliDebug(4, "No Track defined.");
    return NULL;
  }
  TObjArray *arr = NULL;
  if(!fContainer || !(arr = ((TObjArray*)fContainer->At(kCluster)))){
    AliWarning("No output container defined.");
    return NULL;
  }
  ULong_t status = fkESD ? fkESD->GetStatus():0;

  Int_t sgm[3];
  Double_t covR[7], cov[3], dy[2], dz[2];
  Float_t pt, x0, y0, z0, dydx, dzdx;
  AliTRDseedV1 *fTracklet(NULL);  TObjArray *clInfoArr(NULL);
  for(Int_t ily=0; ily<AliTRDgeometry::kNlayer; ily++){
    if(!(fTracklet = fkTrack->GetTracklet(ily))) continue;
    if(!fTracklet->IsOK()) continue;
    x0 = fTracklet->GetX0();
    pt = fTracklet->GetPt();
    sgm[2] = fTracklet->GetDetector();
    sgm[0] = AliTRDgeometry::GetSector(sgm[2]);
    sgm[1] = sgm[0] * AliTRDgeometry::kNstack + AliTRDgeometry::GetStack(sgm[2]);
    
    // retrive the track angle with the chamber
    y0   = fTracklet->GetYref(0);
    z0   = fTracklet->GetZref(0);
    dydx = fTracklet->GetYref(1);
    dzdx = fTracklet->GetZref(1);
    fTracklet->GetCovRef(covR);
    Double_t tilt(fTracklet->GetTilt())
            ,t2(tilt*tilt)
            ,corr(1./(1. + t2))
            ,cost(TMath::Sqrt(corr));
    AliTRDcluster *c = NULL;
    fTracklet->ResetClusterIter(kFALSE);
    while((c = fTracklet->PrevCluster())){
      Float_t xc = c->GetX();
      Float_t yc = c->GetY();
      Float_t zc = c->GetZ();
      Float_t dx = x0 - xc; 
      Float_t yt = y0 - dx*dydx;
      Float_t zt = z0 - dx*dzdx; 
      dy[0] = yc-yt; dz[0]= zc-zt;

      // rotate along pad
      dy[1] = cost*(dy[0] - dz[0]*tilt);
      dz[1] = cost*(dz[0] + dy[0]*tilt);
      if(pt>fPtThreshold && c->IsInChamber()) ((TH3S*)arr->At(0))->Fill(dydx, dy[1], sgm[fSegmentLevel]);

      // tilt rotation of covariance for clusters
      Double_t sy2(c->GetSigmaY2()), sz2(c->GetSigmaZ2());
      cov[0] = (sy2+t2*sz2)*corr;
      cov[1] = tilt*(sz2 - sy2)*corr;
      cov[2] = (t2*sy2+sz2)*corr;
      // sum with track covariance
      cov[0]+=covR[0]; cov[1]+=covR[1]; cov[2]+=covR[2];
      Double_t dyz[2]= {dy[1], dz[1]};
      Pulls(dyz, cov, tilt);
      ((TH3S*)arr->At(1))->Fill(sgm[fSegmentLevel], dyz[0], dyz[1]);
  
      // Get z-position with respect to anode wire
      Int_t istk = fGeo->GetStack(c->GetDetector());
      AliTRDpadPlane *pp = fGeo->GetPadPlane(ily, istk);
      Float_t row0 = pp->GetRow0();
      Float_t d  =  row0 - zt + pp->GetAnodeWireOffset();
      d -= ((Int_t)(2 * d)) / 2.0;
      if (d > 0.25) d  = 0.5 - d;

      AliTRDclusterInfo *clInfo = new AliTRDclusterInfo;
      clInfo->SetCluster(c);
      Float_t covF[] = {cov[0], cov[1], cov[2]};
      clInfo->SetGlobalPosition(yt, zt, dydx, dzdx, covF);
      clInfo->SetResolution(dy[1]);
      clInfo->SetAnisochronity(d);
      clInfo->SetDriftLength(dx);
      clInfo->SetTilt(tilt);
      if(fCl) fCl->Add(clInfo);
      else AliDebug(1, "Cl exchange container missing. Activate by calling \"InitExchangeContainers()\"");

      if(DebugLevel()>=1){
        if(!clInfoArr){ 
          clInfoArr=new TObjArray(AliTRDseedV1::kNclusters);
          clInfoArr->SetOwner(kFALSE);
        }
        clInfoArr->Add(clInfo);
      }
    }
    if(DebugLevel()>=1 && clInfoArr){
      (*DebugStream()) << "cluster"
        <<"status="  << status
        <<"clInfo.=" << clInfoArr
        << "\n";
      clInfoArr->Clear();
    }
  }
  if(clInfoArr) delete clInfoArr;
  return (TH3S*)arr->At(0);
}


//________________________________________________________
TH1* AliTRDresolution::PlotTracklet(const AliTRDtrackV1 *track)
{
// Plot normalized residuals for tracklets to track. 
// 
// We start from the result that if X=N(|m|, |Cov|)
// BEGIN_LATEX
// (Cov^{-1})^{1/2}X = N((Cov^{-1})^{1/2}*|m|, |1|)
// END_LATEX
// in our case X=(y_trklt - y_trk z_trklt - z_trk) and |Cov| = |Cov_trklt| + |Cov_trk| at the radial 
// reference position. 
  if(track) fkTrack = track;
  if(!fkTrack){
    AliDebug(4, "No Track defined.");
    return NULL;
  }
  TObjArray *arr = NULL;
  if(!fContainer || !(arr = (TObjArray*)fContainer->At(kTrack ))){
    AliWarning("No output container defined.");
    return NULL;
  }

  Int_t sgm[3];
  Double_t cov[3], covR[7]/*, sqr[3], inv[3]*/;
  Double_t pt, phi, tht, x, dx, dy[2], dz[2];
  AliTRDseedV1 *fTracklet(NULL);  
  for(Int_t il(0); il<AliTRDgeometry::kNlayer; il++){
    if(!(fTracklet = fkTrack->GetTracklet(il))) continue;
    if(!fTracklet->IsOK()) continue;
    sgm[2] = fTracklet->GetDetector();
    sgm[0] = AliTRDgeometry::GetSector(sgm[2]);
    sgm[1] = sgm[0] * AliTRDgeometry::kNstack + AliTRDgeometry::GetStack(sgm[2]);
    x   = fTracklet->GetX();
    dx  = fTracklet->GetX0() - x;
    pt  = fTracklet->GetPt();
    phi = fTracklet->GetYref(1);
    tht = fTracklet->GetZref(1);
    // compute dy and dz
    dy[0]= fTracklet->GetYref(0)-dx*fTracklet->GetYref(1) - fTracklet->GetY();
    dz[0]= fTracklet->GetZref(0)-dx*fTracklet->GetZref(1) - fTracklet->GetZ();
    Double_t tilt(fTracklet->GetTilt())
            ,t2(tilt*tilt)
            ,corr(1./(1. + t2))
            ,cost(TMath::Sqrt(corr));
    Bool_t rc(fTracklet->IsRowCross());

    // calculate residuals using tilt rotation
    dy[1]= cost*(dy[0] - dz[0]*tilt);
    dz[1]= cost*(dz[0] + dy[0]*tilt);
    ((TH3S*)arr->At(0))->Fill(phi, dy[1], sgm[fSegmentLevel]+rc*fgkNresYsegm[fSegmentLevel]);
    ((TH3S*)arr->At(2))->Fill(tht, dz[1], rc);

    // compute covariance matrix
    fTracklet->GetCovAt(x, cov);
    fTracklet->GetCovRef(covR);
    cov[0] += covR[0]; cov[1] += covR[1]; cov[2] += covR[2]; 
    Double_t dyz[2]= {dy[1], dz[1]};
    Pulls(dyz, cov, tilt);
    ((TH3S*)arr->At(1))->Fill(sgm[fSegmentLevel], dyz[0], dyz[1]);
    ((TH3S*)arr->At(3))->Fill(tht, dyz[1], rc);

    Double_t dphi((phi-fTracklet->GetYfit(1))/(1-phi*fTracklet->GetYfit(1)));
    Double_t dtht((tht-fTracklet->GetZfit(1))/(1-tht*fTracklet->GetZfit(1)));
    ((TH2I*)arr->At(4))->Fill(phi, TMath::ATan(dphi));

    if(DebugLevel()>=1){
      UChar_t err(fTracklet->GetErrorMsg());
      (*DebugStream()) << "tracklet"
        <<"pt="  << pt
        <<"phi=" << phi
        <<"tht=" << tht
        <<"det=" << sgm[2]
        <<"dy0="  << dy[0]
        <<"dz0="  << dz[0]
        <<"dy="  << dy[1]
        <<"dz="  << dz[1]
        <<"dphi="<< dphi
        <<"dtht="<< dtht
        <<"dyp=" << dyz[0]
        <<"dzp=" << dyz[1]
        <<"rc="  << rc
        <<"err=" << err
        << "\n";
    }
  }


  return (TH2I*)arr->At(0);
}


//________________________________________________________
TH1* AliTRDresolution::PlotTrackIn(const AliTRDtrackV1 *track)
{
// Store resolution/pulls of Kalman before updating with the TRD information 
// at the radial position of the first tracklet. The following points are used 
// for comparison  
//  - the (y,z,snp) of the first TRD tracklet
//  - the (y, z, snp, tgl, pt) of the MC track reference
// 
// Additionally the momentum resolution/pulls are calculated for usage in the 
// PID calculation. 

  if(track) fkTrack = track;
  if(!fkTrack){
    AliDebug(4, "No Track defined.");
    return NULL;
  }
  TObjArray *arr = NULL;
  if(!fContainer || !(arr = (TObjArray*)fContainer->At(kTrackIn))){
    AliWarning("No output container defined.");
    return NULL;
  }
  AliExternalTrackParam *tin = NULL;
  if(!(tin = fkTrack->GetTrackIn())){
    AliWarning("Track did not entered TRD fiducial volume.");
    return NULL;
  }
  TH1 *h = NULL;
  
  Double_t x = tin->GetX();
  AliTRDseedV1 *fTracklet = NULL;  
  for(Int_t ily=0; ily<AliTRDgeometry::kNlayer; ily++){
    if(!(fTracklet = fkTrack->GetTracklet(ily))) continue;
    break;
  }
  if(!fTracklet || TMath::Abs(x-fTracklet->GetX())>1.e-3){
    AliWarning("Tracklet did not match Track.");
    return NULL;
  }
  Int_t sgm[3];
  sgm[2] = fTracklet->GetDetector();
  sgm[0] = AliTRDgeometry::GetSector(sgm[2]);
  sgm[1] = sgm[0] * AliTRDgeometry::kNstack + AliTRDgeometry::GetStack(sgm[2]);
  Double_t tilt(fTracklet->GetTilt())
          ,t2(tilt*tilt)
          ,corr(1./(1. + t2))
          ,cost(TMath::Sqrt(corr));
  Bool_t rc(fTracklet->IsRowCross());

  const Int_t kNPAR(5);
  Double_t parR[kNPAR]; memcpy(parR, tin->GetParameter(), kNPAR*sizeof(Double_t));
  Double_t covR[3*kNPAR]; memcpy(covR, tin->GetCovariance(), 3*kNPAR*sizeof(Double_t));
  Double_t cov[3]; fTracklet->GetCovAt(x, cov);

  // define sum covariances
  TMatrixDSym COV(kNPAR); TVectorD PAR(kNPAR);
  Double_t *pc = &covR[0], *pp = &parR[0];
  for(Int_t ir=0; ir<kNPAR; ir++, pp++){
    PAR(ir) = (*pp);
    for(Int_t ic = 0; ic<=ir; ic++,pc++){ 
      COV(ir,ic) = (*pc); COV(ic,ir) = (*pc);
    }
  }
  PAR[4] = TMath::Abs(PAR[4]); // remove sign of pt !!
  //COV.Print(); PAR.Print();

  //TODO Double_t dydx =  TMath::Sqrt(1.-parR[2]*parR[2])/parR[2]; 
  Double_t dy[2]={parR[0] - fTracklet->GetY(), 0.}
          ,dz[2]={parR[1] - fTracklet->GetZ(), 0.}
          ,dphi(TMath::ASin(PAR[2])-TMath::ATan(fTracklet->GetYfit(1)));
  // calculate residuals using tilt rotation
  dy[1] = cost*(dy[0] - dz[0]*tilt);
  dz[1] = cost*(dz[0] + dy[0]*tilt);

  if(1./PAR[4]>fPtThreshold) ((TH3S*)arr->At(0))->Fill(fTracklet->GetYref(1), dy[1], sgm[fSegmentLevel]+rc*fgkNresYsegm[fSegmentLevel]);
  ((TH3S*)arr->At(2))->Fill(fTracklet->GetZref(1), dz[1], rc);
  ((TH2I*)arr->At(4))->Fill(fTracklet->GetYref(1), dphi);

  Double_t dyz[2] = {dy[1], dz[1]};
  Double_t cc[3] = {COV(0,0)+cov[0], COV(0,1)+cov[1], COV(1,1)+cov[2]};
  Pulls(dyz, cc, tilt);
  ((TH3S*)arr->At(1))->Fill(sgm[fSegmentLevel], dyz[0], dyz[1]);
  ((TH3S*)arr->At(3))->Fill(fTracklet->GetZref(1), dyz[1], rc);



  // register reference histo for mini-task
  h = (TH2I*)arr->At(0);

  if(DebugLevel()>=2){
    (*DebugStream()) << "trackIn"
      << "x="       << x
      << "P="       << &PAR
      << "C="       << &COV
      << "\n";

    Double_t y = fTracklet->GetY(); 
    Double_t z = fTracklet->GetZ(); 
    (*DebugStream()) << "trackletIn"
      << "y="       << y
      << "z="       << z
      << "Vy="      << cov[0]
      << "Cyz="     << cov[1]
      << "Vz="      << cov[2]
      << "\n";
  }


  if(!HasMCdata()) return h;
  UChar_t s;
  Float_t dx, pt0, x0=fTracklet->GetX0(), y0, z0, dydx0, dzdx0;
  if(!fkMC->GetDirections(x0, y0, z0, dydx0, dzdx0, pt0, s)) return h;
  Int_t pdg = fkMC->GetPDG(),
        sIdx(AliTRDpidUtil::Pdg2Pid(TMath::Abs(pdg))+1), // species index
        sign(0);
  if(!fDBPDG) fDBPDG=TDatabasePDG::Instance();
  TParticlePDG *ppdg(fDBPDG->GetParticle(pdg));
  if(ppdg) sign = ppdg->Charge() > 0. ? 1 : -1;

  // translate to reference radial position
  dx = x0 - x; y0 -= dx*dydx0; z0 -= dx*dzdx0;
  Float_t norm = 1./TMath::Sqrt(1.+dydx0*dydx0); // 1/sqrt(1+tg^2(phi))
  //Fill MC info
  TVectorD PARMC(kNPAR);
  PARMC[0]=y0; PARMC[1]=z0;
  PARMC[2]=dydx0*norm; PARMC[3]=dzdx0*norm;
  PARMC[4]=1./pt0;

//   TMatrixDSymEigen eigen(COV);
//   TVectorD evals = eigen.GetEigenValues();
//   TMatrixDSym evalsm(kNPAR);
//   for(Int_t ir=0; ir<kNPAR; ir++) for(Int_t ic=0; ic<kNPAR; ic++) evalsm(ir,ic) = (ir==ic ? evals(ir): 0.);
//   TMatrixD evecs = eigen.GetEigenVectors();
//   TMatrixD sqrcov(evecs, TMatrixD::kMult, TMatrixD(evalsm, TMatrixD::kMult, evecs.T()));
  
  // fill histos
  if(!(arr = (TObjArray*)fContainer->At(kMCtrackIn))) {
    AliWarning("No MC container defined.");
    return h;
  }

  // y resolution/pulls
  if(pt0>fPtThreshold) ((TH3S*)arr->At(0))->Fill(dydx0, PARMC[0]-PAR[0], sgm[fSegmentLevel]);
  ((TH3S*)arr->At(1))->Fill(sgm[fSegmentLevel], (PARMC[0]-PAR[0])/TMath::Sqrt(COV(0,0)), (PARMC[1]-PAR[1])/TMath::Sqrt(COV(1,1)));
  // z resolution/pulls
  ((TH3S*)arr->At(2))->Fill(dzdx0, PARMC[1]-PAR[1], 0);
  ((TH3S*)arr->At(3))->Fill(dzdx0, (PARMC[1]-PAR[1])/TMath::Sqrt(COV(1,1)), 0);
  // phi resolution/snp pulls
  ((TH2I*)arr->At(4))->Fill(dydx0, TMath::ASin(PARMC[2])-TMath::ASin(PAR[2]));
  ((TH2I*)arr->At(5))->Fill(dydx0, (PARMC[2]-PAR[2])/TMath::Sqrt(COV(2,2)));
  // theta resolution/tgl pulls
  ((TH2I*)arr->At(6))->Fill(dzdx0, TMath::ATan((PARMC[3]-PAR[3])/(1-PARMC[3]*PAR[3])));
  ((TH2I*)arr->At(7))->Fill(dzdx0, (PARMC[3]-PAR[3])/TMath::Sqrt(COV(3,3)));
  // pt resolution\\1/pt pulls\\p resolution/pull
  ((TH3S*)arr->At(8))->Fill(pt0, PARMC[4]/PAR[4]-1., sign*sIdx);
  ((TH3S*)arr->At(9))->Fill(PARMC[4], (PARMC[4]-PAR[4])/TMath::Sqrt(COV(4,4)), sign*sIdx);

  Double_t p0 = TMath::Sqrt(1.+ PARMC[3]*PARMC[3])*pt0, p;
  p = TMath::Sqrt(1.+ PAR[3]*PAR[3])/PAR[4];
  ((TH3S*)arr->At(10))->Fill(p0, p/p0-1., sign*sIdx);
//   Float_t sp = 
//     p*p*PAR[4]*PAR[4]*COV(4,4)
//    +2.*PAR[3]*COV(3,4)/PAR[4]
//    +PAR[3]*PAR[3]*COV(3,3)/p/p/PAR[4]/PAR[4]/PAR[4]/PAR[4];
//   if(sp>0.) ((TH3S*)arr->At(11))->Fill(p0, (p0-p)/TMath::Sqrt(sp), sign*sIdx);

  // fill debug for MC 
  if(DebugLevel()>=3){
    (*DebugStream()) << "trackInMC"
      << "P="   << &PARMC
      << "\n";
  }
  return h;
}

//________________________________________________________
TH1* AliTRDresolution::PlotTrackOut(const AliTRDtrackV1 *track)
{
// Store resolution/pulls of Kalman after last update with the TRD information 
// at the radial position of the first tracklet. The following points are used 
// for comparison  
//  - the (y,z,snp) of the first TRD tracklet
//  - the (y, z, snp, tgl, pt) of the MC track reference
// 
// Additionally the momentum resolution/pulls are calculated for usage in the 
// PID calculation. 

  if(track) fkTrack = track;
  if(!fkTrack){
    AliDebug(4, "No Track defined.");
    return NULL;
  }
  TObjArray *arr = NULL;
  if(!fContainer || !(arr = (TObjArray*)fContainer->At(kTrackOut))){
    AliWarning("No output container defined.");
    return NULL;
  }
  AliExternalTrackParam *tout = NULL;
  if(!(tout = fkTrack->GetTrackOut())){
    AliDebug(2, "Track did not exit TRD.");
    return NULL;
  }
  TH1 *h(NULL);
  
  Double_t x = tout->GetX();
  AliTRDseedV1 *fTracklet(NULL);  
  for(Int_t ily=0; ily<AliTRDgeometry::kNlayer; ily++){
    if(!(fTracklet = fkTrack->GetTracklet(ily))) continue;
    break;
  }
  if(!fTracklet || TMath::Abs(x-fTracklet->GetX())>1.e-3){
    AliWarning("Tracklet did not match Track position.");
    return NULL;
  }
  Int_t sgm[3];
  sgm[2] = fTracklet->GetDetector();
  sgm[0] = AliTRDgeometry::GetSector(sgm[2]);
  sgm[1] = sgm[0] * AliTRDgeometry::kNstack + AliTRDgeometry::GetStack(sgm[2]);
  Double_t tilt(fTracklet->GetTilt())
          ,t2(tilt*tilt)
          ,corr(1./(1. + t2))
          ,cost(TMath::Sqrt(corr));
  Bool_t rc(fTracklet->IsRowCross());

  const Int_t kNPAR(5);
  Double_t parR[kNPAR]; memcpy(parR, tout->GetParameter(), kNPAR*sizeof(Double_t));
  Double_t covR[3*kNPAR]; memcpy(covR, tout->GetCovariance(), 3*kNPAR*sizeof(Double_t));
  Double_t cov[3]; fTracklet->GetCovAt(x, cov);

  // define sum covariances
  TMatrixDSym COV(kNPAR); TVectorD PAR(kNPAR);
  Double_t *pc = &covR[0], *pp = &parR[0];
  for(Int_t ir=0; ir<kNPAR; ir++, pp++){
    PAR(ir) = (*pp);
    for(Int_t ic = 0; ic<=ir; ic++,pc++){ 
      COV(ir,ic) = (*pc); COV(ic,ir) = (*pc);
    }
  }
  PAR[4] = TMath::Abs(PAR[4]); // remove sign of pt !!
  //COV.Print(); PAR.Print();

  //TODO Double_t dydx =  TMath::Sqrt(1.-parR[2]*parR[2])/parR[2]; 
  Double_t dy[3]={parR[0] - fTracklet->GetY(), 0., 0.}
          ,dz[3]={parR[1] - fTracklet->GetZ(), 0., 0.}
          ,dphi(TMath::ASin(PAR[2])-TMath::ATan(fTracklet->GetYfit(1)));
  // calculate residuals using tilt rotation
  dy[1] = cost*(dy[0] - dz[0]*tilt);
  dz[1] = cost*(dz[0] + dy[0]*tilt);

  if(1./PAR[4]>fPtThreshold) ((TH3S*)arr->At(0))->Fill(fTracklet->GetYref(1), 1.e2*dy[1], sgm[fSegmentLevel]+rc*fgkNresYsegm[fSegmentLevel]); // scale to fit general residual range !!!
  ((TH3S*)arr->At(2))->Fill(fTracklet->GetZref(1), dz[1], rc);
  ((TH2I*)arr->At(4))->Fill(fTracklet->GetYref(1), dphi);

  Double_t dyz[2] = {dy[1], dz[1]};
  Double_t cc[3] = {COV(0,0)+cov[0], COV(0,1)+cov[1], COV(1,1)+cov[2]};
  Pulls(dyz, cc, tilt);
  ((TH3S*)arr->At(1))->Fill(sgm[fSegmentLevel], dyz[0], dyz[1]);
  ((TH3S*)arr->At(3))->Fill(fTracklet->GetZref(1), dyz[1], rc);

  // register reference histo for mini-task
  h = (TH2I*)arr->At(0);

  if(DebugLevel()>=2){
    (*DebugStream()) << "trackOut"
      << "x="       << x
      << "P="       << &PAR
      << "C="       << &COV
      << "\n";

    Double_t y = fTracklet->GetY(); 
    Double_t z = fTracklet->GetZ(); 
    (*DebugStream()) << "trackletOut"
      << "y="       << y
      << "z="       << z
      << "Vy="      << cov[0]
      << "Cyz="     << cov[1]
      << "Vz="      << cov[2]
      << "\n";
  }


  if(!HasMCdata()) return h;
  UChar_t s;
  Float_t dx, pt0, x0=fTracklet->GetX0(), y0, z0, dydx0, dzdx0;
  if(!fkMC->GetDirections(x0, y0, z0, dydx0, dzdx0, pt0, s)) return h;
  Int_t pdg = fkMC->GetPDG(),
        sIdx(AliTRDpidUtil::Pdg2Pid(TMath::Abs(pdg))+1), // species index
        sign(0);
  if(!fDBPDG) fDBPDG=TDatabasePDG::Instance();
  TParticlePDG *ppdg(fDBPDG->GetParticle(pdg));
  if(ppdg) sign = ppdg->Charge() > 0. ? 1 : -1;

  // translate to reference radial position
  dx = x0 - x; y0 -= dx*dydx0; z0 -= dx*dzdx0;
  Float_t norm = 1./TMath::Sqrt(1.+dydx0*dydx0); // 1/sqrt(1+tg^2(phi))
  //Fill MC info
  TVectorD PARMC(kNPAR);
  PARMC[0]=y0; PARMC[1]=z0;
  PARMC[2]=dydx0*norm; PARMC[3]=dzdx0*norm;
  PARMC[4]=1./pt0;

//   TMatrixDSymEigen eigen(COV);
//   TVectorD evals = eigen.GetEigenValues();
//   TMatrixDSym evalsm(kNPAR);
//   for(Int_t ir=0; ir<kNPAR; ir++) for(Int_t ic=0; ic<kNPAR; ic++) evalsm(ir,ic) = (ir==ic ? evals(ir): 0.);
//   TMatrixD evecs = eigen.GetEigenVectors();
//   TMatrixD sqrcov(evecs, TMatrixD::kMult, TMatrixD(evalsm, TMatrixD::kMult, evecs.T()));
  
  // fill histos
  if(!(arr = (TObjArray*)fContainer->At(kMCtrackOut))){
    AliWarning("No MC container defined.");
    return h;
  }
  // y resolution/pulls
  if(pt0>fPtThreshold) ((TH3S*)arr->At(0))->Fill(dydx0, PARMC[0]-PAR[0], sgm[fSegmentLevel]);
  ((TH3S*)arr->At(1))->Fill(sgm[fSegmentLevel], (PARMC[0]-PAR[0])/TMath::Sqrt(COV(0,0)), (PARMC[1]-PAR[1])/TMath::Sqrt(COV(1,1)));
  // z resolution/pulls
  ((TH3S*)arr->At(2))->Fill(dzdx0, PARMC[1]-PAR[1], 0);
  ((TH3S*)arr->At(3))->Fill(dzdx0, (PARMC[1]-PAR[1])/TMath::Sqrt(COV(1,1)), 0);
  // phi resolution/snp pulls
  ((TH2I*)arr->At(4))->Fill(dydx0, TMath::ASin(PARMC[2])-TMath::ASin(PAR[2]));
  ((TH2I*)arr->At(5))->Fill(dydx0, (PARMC[2]-PAR[2])/TMath::Sqrt(COV(2,2)));
  // theta resolution/tgl pulls
  ((TH2I*)arr->At(6))->Fill(dzdx0, TMath::ATan((PARMC[3]-PAR[3])/(1-PARMC[3]*PAR[3])));
  ((TH2I*)arr->At(7))->Fill(dzdx0, (PARMC[3]-PAR[3])/TMath::Sqrt(COV(3,3)));
  // pt resolution\\1/pt pulls\\p resolution/pull
  ((TH3S*)arr->At(8))->Fill(pt0, PARMC[4]/PAR[4]-1., sign*sIdx);
  ((TH3S*)arr->At(9))->Fill(PARMC[4], (PARMC[4]-PAR[4])/TMath::Sqrt(COV(4,4)), sign*sIdx);

  Double_t p0 = TMath::Sqrt(1.+ PARMC[3]*PARMC[3])*pt0, p;
  p = TMath::Sqrt(1.+ PAR[3]*PAR[3])/PAR[4];
  ((TH3S*)arr->At(10))->Fill(p0, p/p0-1., sign*sIdx);
//   Float_t sp = 
//     p*p*PAR[4]*PAR[4]*COV(4,4)
//    +2.*PAR[3]*COV(3,4)/PAR[4]
//    +PAR[3]*PAR[3]*COV(3,3)/p/p/PAR[4]/PAR[4]/PAR[4]/PAR[4];
//   if(sp>0.) ((TH3S*)arr->At(11))->Fill(p0, (p0-p)/TMath::Sqrt(sp), sign*sIdx);

  // fill debug for MC 
  if(DebugLevel()>=3){
    (*DebugStream()) << "trackOutMC"
      << "P="   << &PARMC
      << "\n";
  }
  return h;
}

//________________________________________________________
TH1* AliTRDresolution::PlotMC(const AliTRDtrackV1 *track)
{
  //
  // Plot MC distributions
  //

  if(!HasMCdata()){ 
    AliDebug(2, "No MC defined. Results will not be available.");
    return NULL;
  }
  if(track) fkTrack = track;
  if(!fkTrack){
    AliDebug(4, "No Track defined.");
    return NULL;
  }
  if(!fContainer){
    AliWarning("No output container defined.");
    return NULL;
  }
  // retriev track characteristics
  Int_t pdg = fkMC->GetPDG(),
        sIdx(AliTRDpidUtil::Pdg2Pid(TMath::Abs(pdg))+1), // species index
        sign(0),
        sgm[3],
        label(fkMC->GetLabel());
  if(!fDBPDG) fDBPDG=TDatabasePDG::Instance();
  TParticlePDG *ppdg(fDBPDG->GetParticle(pdg));
  if(ppdg) sign = ppdg->Charge() > 0. ? 1 : -1;

  TObjArray *arr(NULL);TH1 *h(NULL);
  UChar_t s;
  Double_t xAnode, x, y, z, pt, dydx, dzdx, dzdl;
  Float_t pt0, x0, y0, z0, dx, dy, dz, dydx0, dzdx0;
  Double_t covR[7]/*, cov[3]*/;

  if(DebugLevel()>=3){
    TVectorD dX(12), dY(12), dZ(12), Pt(12), dPt(12), cCOV(12*15);
    fkMC->PropagateKalman(&dX, &dY, &dZ, &Pt, &dPt, &cCOV);
    (*DebugStream()) << "MCkalman"
      << "pdg=" << pdg
      << "dx="  << &dX
      << "dy="  << &dY
      << "dz="  << &dZ
      << "pt="  << &Pt
      << "dpt=" << &dPt
      << "cov=" << &cCOV
      << "\n";
  }

  AliTRDReconstructor rec;
  AliTRDseedV1 *fTracklet(NULL); TObjArray *clInfoArr(NULL);
  for(Int_t ily=0; ily<AliTRDgeometry::kNlayer; ily++){
    if(!(fTracklet = fkTrack->GetTracklet(ily)))/* ||
       !fTracklet->IsOK())*/ continue;

    sgm[2] = fTracklet->GetDetector();
    sgm[0] = AliTRDgeometry::GetSector(sgm[2]);
    sgm[1] = sgm[0] * AliTRDgeometry::kNstack + AliTRDgeometry::GetStack(sgm[2]);
    Double_t tilt(fTracklet->GetTilt())
            ,t2(tilt*tilt)
            ,corr(1./(1. + t2))
            ,cost(TMath::Sqrt(corr));
    x0  = fTracklet->GetX0();
    //radial shift with respect to the MC reference (radial position of the pad plane)
    x= fTracklet->GetX();
    Bool_t rc(fTracklet->IsRowCross());
    if(!fkMC->GetDirections(x0, y0, z0, dydx0, dzdx0, pt0, s)) continue;
    xAnode  = fTracklet->GetX0();

    // MC track position at reference radial position
    dx  = x0 - x;
    if(DebugLevel()>=4){
      (*DebugStream()) << "MC"
        << "det="     << sgm[2]
        << "pdg="     << pdg
        << "sgn="     << sign
        << "pt="      << pt0
        << "x="       << x0
        << "y="       << y0
        << "z="       << z0
        << "dydx="    << dydx0
        << "dzdx="    << dzdx0
        << "\n";
    }
    Float_t ymc = y0 - dx*dydx0;
    Float_t zmc = z0 - dx*dzdx0;
    //p = pt0*TMath::Sqrt(1.+dzdx0*dzdx0); // pt -> p

    // Kalman position at reference radial position
    dx = xAnode - x;
    dydx = fTracklet->GetYref(1);
    dzdx = fTracklet->GetZref(1);
    dzdl = fTracklet->GetTgl();
    y  = fTracklet->GetYref(0) - dx*dydx;
    dy = y - ymc;
    z  = fTracklet->GetZref(0) - dx*dzdx;
    dz = z - zmc;
    pt = TMath::Abs(fTracklet->GetPt());
    fTracklet->GetCovRef(covR);

    arr = (TObjArray*)((TObjArray*)fContainer->At(kMCtrack))->At(ily);
    // y resolution/pulls
    if(pt0>fPtThreshold) ((TH3S*)arr->At(0))->Fill(dydx0, dy, sgm[fSegmentLevel]);
    ((TH3S*)arr->At(1))->Fill(sgm[fSegmentLevel], dy/TMath::Sqrt(covR[0]), dz/TMath::Sqrt(covR[2]));
    // z resolution/pulls
    ((TH3S*)arr->At(2))->Fill(dzdx0, dz, 0);
    ((TH3S*)arr->At(3))->Fill(dzdx0, dz/TMath::Sqrt(covR[2]), 0);
    // phi resolution/ snp pulls
    Double_t dtgp = (dydx - dydx0)/(1.- dydx*dydx0);
    ((TH2I*)arr->At(4))->Fill(dydx0, TMath::ATan(dtgp));
    Double_t dsnp = dydx/TMath::Sqrt(1.+dydx*dydx) - dydx0/TMath::Sqrt(1.+dydx0*dydx0);
    ((TH2I*)arr->At(5))->Fill(dydx0, dsnp/TMath::Sqrt(covR[3]));
    // theta resolution/ tgl pulls
    Double_t dzdl0 = dzdx0/TMath::Sqrt(1.+dydx0*dydx0),
              dtgl = (dzdl - dzdl0)/(1.- dzdl*dzdl0);
    ((TH2I*)arr->At(6))->Fill(dzdl0, 
    TMath::ATan(dtgl));
    ((TH2I*)arr->At(7))->Fill(dzdl0, (dzdl - dzdl0)/TMath::Sqrt(covR[4]));
    // pt resolution  \\ 1/pt pulls \\ p resolution for PID
    Double_t p0 = TMath::Sqrt(1.+ dzdl0*dzdl0)*pt0,
             p  = TMath::Sqrt(1.+ dzdl*dzdl)*pt;
    ((TH3S*)((TObjArray*)arr->At(8)))->Fill(pt0, pt/pt0-1., sign*sIdx);
    ((TH3S*)((TObjArray*)arr->At(9)))->Fill(1./pt0, (1./pt-1./pt0)/TMath::Sqrt(covR[6]), sign*sIdx);
    ((TH3S*)((TObjArray*)arr->At(10)))->Fill(p0, p/p0-1., sign*sIdx);

    // Fill Debug stream for Kalman track
    if(DebugLevel()>=4){
      (*DebugStream()) << "MCtrack"
        << "pt="      << pt
        << "x="       << x
        << "y="       << y
        << "z="       << z
        << "dydx="    << dydx
        << "dzdx="    << dzdx
        << "s2y="     << covR[0]
        << "s2z="     << covR[2]
        << "\n";
    }

    // recalculate tracklet based on the MC info
    AliTRDseedV1 tt(*fTracklet);
    tt.SetZref(0, z0 - (x0-xAnode)*dzdx0);
    tt.SetZref(1, dzdx0); 
    tt.SetReconstructor(&rec);
    tt.Fit(1);
    x= tt.GetX();y= tt.GetY();z= tt.GetZ();
    dydx = tt.GetYfit(1);
    dx = x0 - x;
    ymc = y0 - dx*dydx0;
    zmc = z0 - dx*dzdx0;
    dy = y-ymc;
    dz = z-zmc;
    Float_t dphi  = (dydx - dydx0);
    dphi /= (1.- dydx*dydx0);

    // add tracklet residuals for y and dydx
    arr = (TObjArray*)fContainer->At(kMCtracklet);

    if(pt0>fPtThreshold) ((TH3S*)arr->At(0))->Fill(dydx0, dy, sgm[fSegmentLevel]);
    if(tt.GetS2Y()>0. && tt.GetS2Z()>0.) ((TH3S*)arr->At(1))->Fill(sgm[fSegmentLevel], dy/TMath::Sqrt(tt.GetS2Y()), dz/TMath::Sqrt(tt.GetS2Z()));
    ((TH3S*)arr->At(2))->Fill(dzdl0, dz, rc);
    if(tt.GetS2Z()>0.) ((TH3S*)arr->At(3))->Fill(dzdl0, dz/TMath::Sqrt(tt.GetS2Z()), rc);
    ((TH2I*)arr->At(4))->Fill(dydx0, TMath::ATan(dphi));
  
    // Fill Debug stream for tracklet
    if(DebugLevel()>=4){
      Float_t s2y = tt.GetS2Y();
      Float_t s2z = tt.GetS2Z();
      (*DebugStream()) << "MCtracklet"
        << "rc="    << rc
        << "x="     << x
        << "y="     << y
        << "z="     << z
        << "dydx="  << dydx
        << "s2y="   << s2y
        << "s2z="   << s2z
        << "\n";
    }

    AliTRDpadPlane *pp = fGeo->GetPadPlane(ily, AliTRDgeometry::GetStack(sgm[2]));
    Float_t zr0 = pp->GetRow0() + pp->GetAnodeWireOffset();
    //Double_t exb = AliTRDCommonParam::Instance()->GetOmegaTau(1.5);

    arr = (TObjArray*)fContainer->At(kMCcluster);
    AliTRDcluster *c = NULL;
    tt.ResetClusterIter(kFALSE);
    while((c = tt.PrevCluster())){
      Float_t  q = TMath::Abs(c->GetQ());
      x = c->GetX(); y = c->GetY();z = c->GetZ();
      dx = x0 - x; 
      ymc= y0 - dx*dydx0;
      zmc= z0 - dx*dzdx0;
      dy = cost*(y - ymc - tilt*(z-zmc));
      dz = cost*(z - zmc + tilt*(y-ymc));
      
      // Fill Histograms
      if(q>20. && q<250. && pt0>fPtThreshold && c->IsInChamber()){ 
        ((TH3S*)arr->At(0))->Fill(dydx0, dy, sgm[fSegmentLevel]);
        ((TH3S*)arr->At(1))->Fill(sgm[fSegmentLevel], dy/TMath::Sqrt(c->GetSigmaY2()), dz/TMath::Sqrt(c->GetSigmaZ2()));
      }

      // Fill calibration container
      Float_t d = zr0 - zmc;
      d -= ((Int_t)(2 * d)) / 2.0;
      if (d > 0.25) d  = 0.5 - d;
      AliTRDclusterInfo *clInfo = new AliTRDclusterInfo;
      clInfo->SetCluster(c);
      clInfo->SetMC(pdg, label);
      clInfo->SetGlobalPosition(ymc, zmc, dydx0, dzdx0);
      clInfo->SetResolution(dy);
      clInfo->SetAnisochronity(d);
      clInfo->SetDriftLength(dx);
      clInfo->SetTilt(tilt);
      if(fMCcl) fMCcl->Add(clInfo);
      else AliDebug(1, "MCcl exchange container missing. Activate by calling \"InitExchangeContainers()\"");
      if(DebugLevel()>=5){ 
        if(!clInfoArr){ 
          clInfoArr=new TObjArray(AliTRDseedV1::kNclusters);
          clInfoArr->SetOwner(kFALSE);
        }
        clInfoArr->Add(clInfo);
      }
    }
    // Fill Debug Tree
    if(DebugLevel()>=5 && clInfoArr){
      (*DebugStream()) << "MCcluster"
        <<"clInfo.=" << clInfoArr
        << "\n";
      clInfoArr->Clear();
    }
  }
  if(clInfoArr) delete clInfoArr;
  return h;
}


//________________________________________________________
Bool_t AliTRDresolution::GetRefFigure(Int_t ifig)
{
  //
  // Get the reference figures
  //

  Float_t xy[4] = {0., 0., 0., 0.};
  if(!gPad){
    AliWarning("Please provide a canvas to draw results.");
    return kFALSE;
  }
  Int_t selection[100], n(0), selStart(0); // 
  Int_t ly0(0), dly(5);
  //Int_t ly0(1), dly(2); // used for SA
  TList *l(NULL); TVirtualPad *pad(NULL); 
  TGraphErrors *g(NULL);TGraphAsymmErrors *ga(NULL);
  switch(ifig){
  case 0: // charge resolution
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    ((TVirtualPad*)l->At(0))->cd();
    ga=((TGraphAsymmErrors*)((TObjArray*)fGraphM->At(kCharge))->At(0));
    if(ga->GetN()) ga->Draw("apl");
    ((TVirtualPad*)l->At(1))->cd();
    g = ((TGraphErrors*)((TObjArray*)fGraphS->At(kCharge))->At(0));
    if(g->GetN()) g->Draw("apl");
    break;
  case 1: // cluster2track residuals
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    xy[0] = -.3; xy[1] = -100.; xy[2] = .3; xy[3] = 1000.;
    pad = (TVirtualPad*)l->At(0); pad->cd();
    pad->SetMargin(0.125, 0.015, 0.1, 0.015);
    selStart=0; for(n=0; n<fgkNresYsegm[fSegmentLevel]/3; n++) selection[n]=selStart+n;
    if(!GetGraphArray(xy, kCluster, 0, 1, n, selection)) break;
    pad=(TVirtualPad*)l->At(1); pad->cd();
    pad->SetMargin(0.125, 0.015, 0.1, 0.015);
    selStart=fgkNresYsegm[fSegmentLevel]/3; for(n=0; n<fgkNresYsegm[fSegmentLevel]/3; n++) selection[n]=selStart+n;
    if(!GetGraphArray(xy, kCluster, 0, 1, n, selection)) break;
    return kTRUE;
  case 2: // cluster2track residuals
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    xy[0] = -.3; xy[1] = -100.; xy[2] = .3; xy[3] = 1000.;
    pad = (TVirtualPad*)l->At(0); pad->cd();
    pad->SetMargin(0.125, 0.015, 0.1, 0.015);
    selStart=2*fgkNresYsegm[fSegmentLevel]/3; for(n=0; n<fgkNresYsegm[fSegmentLevel]/3; n++) selection[n]=selStart+n;
    if(!GetGraphArray(xy, kCluster, 0, 1, n, selection)) break;
    xy[0] = -.5; xy[1] = -0.5; xy[2] = fgkNresYsegm[fSegmentLevel]-0.5; xy[3] = 2.5;
    pad=(TVirtualPad*)l->At(1); pad->cd();
    pad->SetMargin(0.125, 0.015, 0.1, 0.015);
    if(!GetGraphArray(xy, kCluster, 1, 1)) break;
    return kTRUE;
  case 3: // kTrack y
    gPad->Divide(3, 2, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    xy[0] = -.3; xy[1] = -20.; xy[2] = .3; xy[3] = 100.;
    ((TVirtualPad*)l->At(0))->cd();
    selStart=0; for(n=0; n<fgkNresYsegm[fSegmentLevel]/3; n++) selection[n]=selStart+n;
    if(!GetGraphArray(xy, kTrack, 0, 1, n, selection)) break;

    ((TVirtualPad*)l->At(1))->cd();
    selStart=fgkNresYsegm[fSegmentLevel]/3; for(n=0; n<fgkNresYsegm[fSegmentLevel]/3; n++) selection[n]=selStart+n;
    if(!GetGraphArray(xy, kTrack, 0, 1, n, selection)) break;

    ((TVirtualPad*)l->At(2))->cd();
    selStart=2*fgkNresYsegm[fSegmentLevel]/3; for(n=0; n<fgkNresYsegm[fSegmentLevel]/3; n++) selection[n]=selStart+n;
    if(!GetGraphArray(xy, kTrack, 0, 1, n, selection)) break;

    ((TVirtualPad*)l->At(3))->cd();
    selStart=fgkNresYsegm[fSegmentLevel]; for(n=0; n<fgkNresYsegm[fSegmentLevel]/3; n++) selection[n]=selStart+n;
    if(!GetGraphArray(xy, kTrack, 0, 1, n, selection, "[RC]")) break;

    ((TVirtualPad*)l->At(4))->cd();
    selStart=fgkNresYsegm[fSegmentLevel]/3+fgkNresYsegm[fSegmentLevel]; for(n=0; n<fgkNresYsegm[fSegmentLevel]/3; n++) selection[n]=selStart+n;
    if(!GetGraphArray(xy, kTrack, 0, 1, n, selection, "[RC]")) break;

    ((TVirtualPad*)l->At(5))->cd();
    selStart=2*fgkNresYsegm[fSegmentLevel]/3+fgkNresYsegm[fSegmentLevel]; for(n=0; n<fgkNresYsegm[fSegmentLevel]/3; n++) selection[n]=selStart+n;
    if(!GetGraphArray(xy, kTrack, 0, 1, n, selection, "[RC]")) break;
    return kTRUE;
  case 4: // kTrack z
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 

    xy[0] = -1.; xy[1] = -150.; xy[2] = 1.; xy[3] = 1000.;
    ((TVirtualPad*)l->At(0))->cd();
    selection[0]=1;
    if(!GetGraphArray(xy, kTrack, 2, 1, 1, selection)) break;

    xy[0] = -1.; xy[1] = -1500.; xy[2] = 1.; xy[3] = 10000.;
    ((TVirtualPad*)l->At(1))->cd();
    selection[0]=0;
    if(!GetGraphArray(xy, kTrack, 2, 1, 1, selection)) break;

    return kTRUE;
  case 5: // kTrack  pulls
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 

    xy[0] = -.5; xy[1] = -0.5; xy[2] = fgkNresYsegm[fSegmentLevel]-.5; xy[3] = 2.5;
    ((TVirtualPad*)l->At(0))->cd();
    if(!GetGraphArray(xy, kTrack, 1, 1)) break;

    xy[0] = -1.; xy[1] = -0.5; xy[2] = 1.; xy[3] = 2.5;
    ((TVirtualPad*)l->At(1))->cd();
    if(!GetGraphArray(xy, kTrack, 3, 1)) break;
    return kTRUE;
  case 6: // kTrack  phi
    xy[0] = -.3; xy[1] = -5.; xy[2] = .3; xy[3] = 50.;
    if(GetGraph(&xy[0], kTrack , 4)) return kTRUE;
    break;
  case 7: // kTrackIn y
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    xy[0] = -.3; xy[1] = -1500.; xy[2] = .3; xy[3] = 5000.;
    pad = ((TVirtualPad*)l->At(0)); pad->cd();
    pad->SetMargin(0.1, 0.1, 0.1, 0.01);
    selStart=0; for(n=0; n<fgkNresYsegm[fSegmentLevel]/3; n++) selection[n]=selStart+n;
    if(!GetGraphArray(xy, kTrackIn, 0, 1, n, selection)) break;
    pad=((TVirtualPad*)l->At(1)); pad->cd();
    pad->SetMargin(0.1, 0.1, 0.1, 0.01);
    selStart=fgkNresYsegm[fSegmentLevel]/3; for(n=0; n<fgkNresYsegm[fSegmentLevel]/3; n++) selection[n]=selStart+n;
    if(!GetGraphArray(xy, kTrackIn, 0, 1, n, selection)) break;
    return kTRUE;
  case 8: // kTrackIn y
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    xy[0] = -.3; xy[1] = -1500.; xy[2] = .3; xy[3] = 5000.;
    pad = ((TVirtualPad*)l->At(0)); pad->cd();
    pad->SetMargin(0.1, 0.1, 0.1, 0.01);
    selStart=2*fgkNresYsegm[fSegmentLevel]/3; for(n=0; n<fgkNresYsegm[fSegmentLevel]/3; n++) selection[n]=selStart+n;
    if(!GetGraphArray(xy, kTrackIn, 0, 1, n, selection)) break;
    xy[0] = -.5; xy[1] = -0.5; xy[2] = fgkNresYsegm[fSegmentLevel]-.5; xy[3] = 2.5;
    pad=((TVirtualPad*)l->At(1)); pad->cd();
    pad->SetMargin(0.1, 0.1, 0.1, 0.01);
    if(!GetGraphArray(xy, kTrackIn, 1, 1)) break;
    return kTRUE;
  case 9: // kTrackIn z
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    xy[0] = -1.; xy[1] = -1000.; xy[2] = 1.; xy[3] = 4000.;
    pad = ((TVirtualPad*)l->At(0)); pad->cd();
    pad->SetMargin(0.1, 0.1, 0.1, 0.01);
    selection[0]=1;
    if(!GetGraphArray(xy, kTrackIn, 2, 1, 1, selection)) break;
    xy[0] = -1.; xy[1] = -0.5; xy[2] = 1.; xy[3] = 2.5;
    pad = ((TVirtualPad*)l->At(1)); pad->cd();
    pad->SetMargin(0.1, 0.1, 0.1, 0.01);
    if(!GetGraphArray(xy, kTrackIn, 3, 1)) break;
    return kTRUE;
  case 10: // kTrackIn phi
    xy[0] = -.3; xy[1] = -5.; xy[2] = .3; xy[3] = 50.;
    if(GetGraph(&xy[0], kTrackIn, 4)) return kTRUE;
    break;
  case 11: // kTrackOut y
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    xy[0] = -.3; xy[1] = -50.; xy[2] = .3; xy[3] = 150.;
    pad = ((TVirtualPad*)l->At(0)); pad->cd();
    pad->SetMargin(0.1, 0.1, 0.1, 0.01);
    selStart=0; for(n=0; n<fgkNresYsegm[fSegmentLevel]/3; n++) selection[n]=selStart+n;
    if(!GetGraphArray(xy, kTrackOut, 0, 1, n, selection)) break;
    pad=((TVirtualPad*)l->At(1)); pad->cd();
    pad->SetMargin(0.1, 0.1, 0.1, 0.01);
    selStart=fgkNresYsegm[fSegmentLevel]/3; for(n=0; n<fgkNresYsegm[fSegmentLevel]/3; n++) selection[n]=selStart+n;
    if(!GetGraphArray(xy, kTrackOut, 0, 1, n, selection)) break;
    return kTRUE;
  case 12: // kTrackOut y
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    xy[0] = -.3; xy[1] = -50.; xy[2] = .3; xy[3] = 150.;
    pad = ((TVirtualPad*)l->At(0)); pad->cd();
    pad->SetMargin(0.1, 0.1, 0.1, 0.01);
    selStart=2*fgkNresYsegm[fSegmentLevel]/3; for(n=0; n<fgkNresYsegm[fSegmentLevel]/3; n++) selection[n]=selStart+n;
    if(!GetGraphArray(xy, kTrackOut, 0, 1, n, selection)) break;
    xy[0] = -.5; xy[1] = -0.5; xy[2] = fgkNresYsegm[fSegmentLevel]-.5; xy[3] = 2.5;
    pad=((TVirtualPad*)l->At(1)); pad->cd();
    pad->SetMargin(0.1, 0.1, 0.1, 0.01);
    if(!GetGraphArray(xy, kTrackOut, 1, 1)) break;
    return kTRUE;
  case 13: // kTrackOut z
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    xy[0] = -1.; xy[1] = -1000.; xy[2] = 1.; xy[3] = 4000.;
    pad = ((TVirtualPad*)l->At(0)); pad->cd();
    pad->SetMargin(0.1, 0.1, 0.1, 0.01);
    if(!GetGraphArray(xy, kTrackOut, 2, 1)) break;
    xy[0] = -1.; xy[1] = -0.5; xy[2] = 1.; xy[3] = 2.5;
    pad = ((TVirtualPad*)l->At(1)); pad->cd();
    pad->SetMargin(0.1, 0.1, 0.1, 0.01);
    if(!GetGraphArray(xy, kTrackOut, 3, 1)) break;
    return kTRUE;
  case 14: // kTrackOut phi
    xy[0] = -.3; xy[1] = -5.; xy[2] = .3; xy[3] = 50.;
    if(GetGraph(&xy[0], kTrackOut, 4)) return kTRUE;
    break;
  case 15: // kMCcluster
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    xy[0]=-.3; xy[1]=-50.; xy[2]=.3; xy[3]=650.;
    ((TVirtualPad*)l->At(0))->cd();
    selStart=0; for(n=0; n<fgkNresYsegm[fSegmentLevel]/3; n++) selection[n]=selStart+n;
    if(!GetGraphArray(xy, kMCcluster, 0, 1, n, selection)) break;
    ((TVirtualPad*)l->At(1))->cd();
    selStart=fgkNresYsegm[fSegmentLevel]/3; for(n=0; n<fgkNresYsegm[fSegmentLevel]/3; n++) selection[n]=selStart+n;
    if(!GetGraphArray(xy, kMCcluster, 0, 1, n, selection)) break;
    return kTRUE;
  case 16: // kMCcluster
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    xy[0]=-.3; xy[1]=-50.; xy[2]=.3; xy[3]=650.;
    ((TVirtualPad*)l->At(0))->cd();
    selStart=2*fgkNresYsegm[fSegmentLevel]/3; for(n=0; n<fgkNresYsegm[fSegmentLevel]/3; n++) selection[n]=selStart+n;
    if(!GetGraphArray(xy, kMCcluster, 0, 1, n, selection)) break;
    ((TVirtualPad*)l->At(1))->cd();
    xy[0]=-.5; xy[1]=-0.5; xy[2]=fgkNresYsegm[fSegmentLevel]-.5; xy[3]=2.5;
    if(!GetGraphArray(xy, kMCcluster, 1, 1)) break;
    return kTRUE;
  case 17: //kMCtracklet [y]
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    xy[0]=-.3; xy[1]=-50.; xy[2]=.3; xy[3] =500.;
    ((TVirtualPad*)l->At(0))->cd();
    selStart=0; for(n=0; n<fgkNresYsegm[fSegmentLevel]/3; n++) selection[n]=selStart+n;
    if(!GetGraphArray(xy, kMCtracklet, 0, 1, n, selection)) break;
    ((TVirtualPad*)l->At(1))->cd();
    selStart=fgkNresYsegm[fSegmentLevel]/3; for(n=0; n<fgkNresYsegm[fSegmentLevel]/3; n++) selection[n]=selStart+n;
    if(!GetGraphArray(xy, kMCtracklet, 0, 1, n, selection)) break;
    return kTRUE;
  case 18: //kMCtracklet [y]
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    xy[0]=-.3; xy[1]=-50.; xy[2]=.3; xy[3] =500.;
    ((TVirtualPad*)l->At(0))->cd();
    selStart=2*fgkNresYsegm[fSegmentLevel]/3; for(n=0; n<fgkNresYsegm[fSegmentLevel]/3; n++) selection[n]=selStart+n;
    if(!GetGraphArray(xy, kMCtracklet, 0, 1, n, selection)) break;
    ((TVirtualPad*)l->At(1))->cd();
    xy[0]=-.5; xy[1]=-0.5; xy[2]=fgkNresYsegm[fSegmentLevel]-.5; xy[3]=2.5;
    if(!GetGraphArray(xy, kMCtracklet, 1, 1)) break;
    return kTRUE;
  case 19: //kMCtracklet [z]
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    xy[0]=-1.; xy[1]=-100.; xy[2]=1.; xy[3] =2500.;
    ((TVirtualPad*)l->At(0))->cd();
    if(!GetGraphArray(xy, kMCtracklet, 2)) break;
    xy[0] = -1.; xy[1] = -0.5; xy[2] = 1.; xy[3] = 2.5;
    ((TVirtualPad*)l->At(1))->cd();
    if(!GetGraphArray(xy, kMCtracklet, 3)) break;
    return kTRUE;
  case 20: //kMCtracklet [phi]
    xy[0]=-.3; xy[1]=-3.; xy[2]=.3; xy[3] =25.;
    if(!GetGraph(&xy[0], kMCtracklet, 4)) break;
    return kTRUE;
  case 21: //kMCtrack [y] ly [0]
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    xy[0]=-.2; xy[1]=-50.; xy[2]=.2; xy[3] =400.;
    ((TVirtualPad*)l->At(0))->cd();
    selStart=Int_t(fgkNresYsegm[fSegmentLevel]*0.); for(n=0; n<fgkNresYsegm[fSegmentLevel]/2; n++) selection[n]=selStart+n;
    if(!GetGraphArray(xy, kMCtrack, 0, 1, n, selection, "Layer1")) break;
    ((TVirtualPad*)l->At(1))->cd();
    selStart=Int_t(fgkNresYsegm[fSegmentLevel]*0.5); for(n=0; n<fgkNresYsegm[fSegmentLevel]/2; n++) selection[n]=selStart+n;
    if(!GetGraphArray(xy, kMCtrack, 0, 1, n, selection, "Layer1")) break;
    return kTRUE;
  case 22: //kMCtrack [y] ly [1]
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    xy[0]=-.2; xy[1]=-50.; xy[2]=.2; xy[3] =400.;
    ((TVirtualPad*)l->At(0))->cd();
    selStart=Int_t(fgkNresYsegm[fSegmentLevel]*1.); for(n=0; n<fgkNresYsegm[fSegmentLevel]/2; n++) selection[n]=selStart+n;
    if(!GetGraphArray(xy, kMCtrack, 0, 1, n, selection, "Layer2")) break;
    ((TVirtualPad*)l->At(1))->cd();
    selStart=Int_t(fgkNresYsegm[fSegmentLevel]*1.5); for(n=0; n<fgkNresYsegm[fSegmentLevel]/2; n++) selection[n]=selStart+n;
    if(!GetGraphArray(xy, kMCtrack, 0, 1, n, selection, "Layer2")) break;
    return kTRUE;
  case 23: //kMCtrack [y] ly [2]
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    xy[0]=-.2; xy[1]=-50.; xy[2]=.2; xy[3] =400.;
    ((TVirtualPad*)l->At(0))->cd();
    selStart=Int_t(fgkNresYsegm[fSegmentLevel]*2.); for(n=0; n<fgkNresYsegm[fSegmentLevel]/2; n++) selection[n]=selStart+n;
    if(!GetGraphArray(xy, kMCtrack, 0, 1, n, selection, "Layer3")) break;
    ((TVirtualPad*)l->At(1))->cd();
    selStart=Int_t(fgkNresYsegm[fSegmentLevel]*2.5); for(n=0; n<fgkNresYsegm[fSegmentLevel]/2; n++) selection[n]=selStart+n;
    if(!GetGraphArray(xy, kMCtrack, 0, 1, n, selection, "Layer3")) break;
    return kTRUE;
  case 24: //kMCtrack [y] ly [3]
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    xy[0]=-.2; xy[1]=-50.; xy[2]=.2; xy[3] =400.;
    ((TVirtualPad*)l->At(0))->cd();
    selStart=Int_t(fgkNresYsegm[fSegmentLevel]*3.); for(n=0; n<fgkNresYsegm[fSegmentLevel]/2; n++) selection[n]=selStart+n;
    if(!GetGraphArray(xy, kMCtrack, 0, 1, n, selection, "Layer4")) break;
    ((TVirtualPad*)l->At(1))->cd();
    selStart=Int_t(fgkNresYsegm[fSegmentLevel]*3.5); for(n=0; n<fgkNresYsegm[fSegmentLevel]/2; n++) selection[n]=selStart+n;
    if(!GetGraphArray(xy, kMCtrack, 0, 1, n, selection, "Layer4")) break;
    return kTRUE;
  case 25: //kMCtrack [y] ly [4]
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    xy[0]=-.2; xy[1]=-50.; xy[2]=.2; xy[3] =400.;
    ((TVirtualPad*)l->At(0))->cd();
    selStart=Int_t(fgkNresYsegm[fSegmentLevel]*4.); for(n=0; n<fgkNresYsegm[fSegmentLevel]/2; n++) selection[n]=selStart+n;
    if(!GetGraphArray(xy, kMCtrack, 0, 1, n, selection, "Layer5")) break;
    ((TVirtualPad*)l->At(1))->cd();
    selStart=Int_t(fgkNresYsegm[fSegmentLevel]*4.5); for(n=0; n<fgkNresYsegm[fSegmentLevel]/2; n++) selection[n]=selStart+n;
    if(!GetGraphArray(xy, kMCtrack, 0, 1, n, selection, "Layer5")) break;
    return kTRUE;
  case 26: //kMCtrack [y] ly [5]
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    xy[0]=-.2; xy[1]=-50.; xy[2]=.2; xy[3] =400.;
    ((TVirtualPad*)l->At(0))->cd();
    selStart=Int_t(fgkNresYsegm[fSegmentLevel]*5.); for(n=0; n<fgkNresYsegm[fSegmentLevel]/2; n++) selection[n]=selStart+n;
    if(!GetGraphArray(xy, kMCtrack, 0, 1, n, selection, "Layer6")) break;
    ((TVirtualPad*)l->At(1))->cd();
    selStart=Int_t(fgkNresYsegm[fSegmentLevel]*5.5); for(n=0; n<fgkNresYsegm[fSegmentLevel]/2; n++) selection[n]=selStart+n;
    if(!GetGraphArray(xy, kMCtrack, 0, 1, n, selection, "Layer6")) break;
    return kTRUE;
  case 27: //kMCtrack [y pulls] 
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    xy[0] = -.5; xy[1] = -0.5; xy[2] = fgkNresYsegm[fSegmentLevel]-.5; xy[3] = 5.5;
    ((TVirtualPad*)l->At(0))->cd();
    selStart=0; for(n=0; n<6; n++) selection[n]=selStart+n;
    if(!GetGraphArray(xy, kMCtrack, 1, 1, n, selection)) break;
    ((TVirtualPad*)l->At(1))->cd();
    selStart=6; for(n=0; n<6; n++) selection[n]=selStart+n;
    if(!GetGraphArray(xy, kMCtrack, 1, 1, n, selection)) break;
    return kTRUE;
  case 28: //kMCtrack [z]
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    xy[0]=-1.; xy[1]=-1500.; xy[2]=1.; xy[3] =6000.;
    ((TVirtualPad*)l->At(0))->cd();
    if(!GetGraphArray(xy, kMCtrack, 2)) break;
    xy[0] = -1.; xy[1] = -1.5; xy[2] = 1.; xy[3] = 5.;
    ((TVirtualPad*)l->At(1))->cd();
    if(!GetGraphArray(xy, kMCtrack, 3)) break;
    return kTRUE;
  case 29: //kMCtrack [phi/snp]
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    xy[0]=-.2; xy[1]=-0.5; xy[2]=.2; xy[3] =10.;
    ((TVirtualPad*)l->At(0))->cd();
    if(!GetGraphArray(xy, kMCtrack, 4)) break;
    xy[0] = -.2; xy[1] = -1.5; xy[2] = .2; xy[3] = 5.;
    ((TVirtualPad*)l->At(1))->cd();
    if(!GetGraphArray(xy, kMCtrack, 5)) break;
    return kTRUE;
  case 30: //kMCtrack [theta/tgl]
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    xy[0]=-1.; xy[1]=-0.5; xy[2]=1.; xy[3] =5.;
    ((TVirtualPad*)l->At(0))->cd();
    if(!GetGraphArray(xy, kMCtrack, 6)) break;
    xy[0] = -.2; xy[1] = -0.5; xy[2] = .2; xy[3] = 2.5;
    ((TVirtualPad*)l->At(1))->cd();
    if(!GetGraphArray(xy, kMCtrack, 7)) break;
    return kTRUE;
  case 31: //kMCtrack [pt]
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    pad = (TVirtualPad*)l->At(0); pad->cd();
    pad->SetMargin(0.125, 0.015, 0.1, 0.015);
    // pi selection
    n=0;
    for(Int_t il(ly0); il<AliTRDgeometry::kNlayer; il+=dly){
      selection[n++] = il*11 + 2; // pi-
      selection[n++] = il*11 + 8; // pi+
    }
    xy[0] = 0.2; xy[1] = -.7; xy[2] = 7.; xy[3] = 4.;
    //xy[0] = 0.2; xy[1] = -1.; xy[2] = 7.; xy[3] = 10.; // SA
    if(!GetGraphArray(xy, kMCtrack, 8, kTRUE, n, selection, "#pi#pm")) break;
    pad->Modified(); pad->Update(); pad->SetLogx();
    pad = (TVirtualPad*)l->At(1); pad->cd();
    pad->SetMargin(0.125, 0.015, 0.1, 0.015);
    // mu selection
    n=0;    
    for(Int_t il(ly0); il<AliTRDgeometry::kNlayer; il+=dly){
      selection[n++] = il*11 + 3; // mu-
      selection[n++] = il*11 + 7; // mu+
    }
    if(!GetGraphArray(xy, kMCtrack, 8, kTRUE, n, selection, "#mu#pm")) break;
    pad->Modified(); pad->Update(); pad->SetLogx();
    return kTRUE;
  case 32: //kMCtrack [pt]
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    pad = (TVirtualPad*)l->At(0); pad->cd();
    pad->SetMargin(0.125, 0.015, 0.1, 0.015);
    // p selection
    n=0;
    for(Int_t il(ly0); il<AliTRDgeometry::kNlayer; il+=dly){
      selection[n++] = il*11 + 0; // p bar
      selection[n++] = il*11 + 10; // p
    }
    xy[0] = 0.2; xy[1] = -.7; xy[2] = 7.; xy[3] = 8.;
    //xy[0] = 0.2; xy[1] = -1.; xy[2] = 7.; xy[3] = 10.; // SA
    if(!GetGraphArray(xy, kMCtrack, 8, kTRUE, n, selection, "p&p bar")) break;
    pad->Modified(); pad->Update(); pad->SetLogx();
    pad = (TVirtualPad*)l->At(1); pad->cd();
    pad->SetMargin(0.125, 0.015, 0.1, 0.015);
    // e selection
    n=0;    
    for(Int_t il(ly0); il<AliTRDgeometry::kNlayer; il+=dly){
      selection[n++] = il*11 + 4; // e-
      selection[n++] = il*11 + 6; // e+
    }
    xy[0] = 0.2; xy[1] = -1.5; xy[2] = 7.; xy[3] = 12.;
    //xy[0] = 0.2; xy[1] = -1.5; xy[2] = 7.; xy[3] = 14.; // SA
    if(!GetGraphArray(xy, kMCtrack, 8, kTRUE, n, selection, "e#pm")) break;
    pad->Modified(); pad->Update(); pad->SetLogx();
    return kTRUE;
  case 33: //kMCtrack [1/pt] pulls
    xy[0] = 0.; xy[1] = -1.; xy[2] = 2.; xy[3] = 3.5;
    //xy[0] = 0.; xy[1] = -1.; xy[2] = 2.; xy[3] = 4.5; // SA
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    pad = (TVirtualPad*)l->At(0); pad->cd();
    pad->SetMargin(0.125, 0.015, 0.1, 0.015);
    // pi selection
    n=0;
    for(Int_t il(ly0); il<AliTRDgeometry::kNlayer; il+=dly){
      selection[n++] = il*11 + 2; // pi-
      selection[n++] = il*11 + 8; // pi+
    }
    if(!GetGraphArray(xy, kMCtrack, 9, kTRUE, n, selection, "#pi#pm")) break;
    pad = (TVirtualPad*)l->At(1); pad->cd();
    pad->SetMargin(0.125, 0.015, 0.1, 0.015);
    // mu selection
    n=0;    
    for(Int_t il(ly0); il<AliTRDgeometry::kNlayer; il+=dly){
      selection[n++] = il*11 + 3; // mu-
      selection[n++] = il*11 + 7; // mu+
    }
    if(!GetGraphArray(xy, kMCtrack, 9, kTRUE, n, selection, "#mu#pm")) break;
    return kTRUE;
  case 34: //kMCtrack [1/pt] pulls
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    pad = (TVirtualPad*)l->At(0); pad->cd();
    pad->SetMargin(0.125, 0.015, 0.1, 0.015);
    // p selection
    n=0;
    for(Int_t il(ly0); il<AliTRDgeometry::kNlayer; il+=dly){
      selection[n++] = il*11 + 0; // p bar
      selection[n++] = il*11 + 10; // p
    }
    xy[0] = 0.; xy[1] = -1.; xy[2] = 2.; xy[3] = 3.5;
    //xy[0] = 0.; xy[1] = -1.; xy[2] = 2.; xy[3] = 6.; // SA
    if(!GetGraphArray(xy, kMCtrack, 9, kTRUE, n, selection, "p & p bar")) break;
    pad = (TVirtualPad*)l->At(1); pad->cd();
    pad->SetMargin(0.125, 0.015, 0.1, 0.015);
    // e selection
    n=0;    
    for(Int_t il(ly0); il<AliTRDgeometry::kNlayer; il+=dly){
      selection[n++] = il*11 + 4; // e-
      selection[n++] = il*11 + 6; // e+
    }
    xy[0] = 0.; xy[1] = -2.; xy[2] = 2.; xy[3] = 4.5;
    if(!GetGraphArray(xy, kMCtrack, 9, kTRUE, n, selection, "e#pm")) break;
    return kTRUE;
  case 35: //kMCtrack [p]
    xy[0] = 0.2; xy[1] = -.7; xy[2] = 7.; xy[3] = 4.;
    //xy[0] = 0.2; xy[1] = -1.5; xy[2] = 7.; xy[3] = 10.;
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    pad = (TVirtualPad*)l->At(0); pad->cd();
    pad->SetMargin(0.125, 0.015, 0.1, 0.015);
    // pi selection
    n=0;
    for(Int_t il(ly0); il<AliTRDgeometry::kNlayer; il+=dly){
      selection[n++] = il*11 + 2; // pi-
      selection[n++] = il*11 + 8; // pi+
    }
    if(!GetGraphArray(xy, kMCtrack, 10, kTRUE, n, selection, "#pi#pm")) break;
    pad->Modified(); pad->Update(); pad->SetLogx();
    pad = (TVirtualPad*)l->At(1); pad->cd();
    pad->SetMargin(0.125, 0.015, 0.1, 0.015);
    // mu selection
    n=0;    
    for(Int_t il(ly0); il<AliTRDgeometry::kNlayer; il+=dly){
      selection[n++] = il*11 + 3; // mu-
      selection[n++] = il*11 + 7; // mu+
    }
    if(!GetGraphArray(xy, kMCtrack, 10, kTRUE, n, selection, "#mu#pm")) break;
    pad->Modified(); pad->Update(); pad->SetLogx();
    return kTRUE;
  case 36: //kMCtrack [p]
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    pad = (TVirtualPad*)l->At(0); pad->cd();
    pad->SetMargin(0.125, 0.015, 0.1, 0.015);
    // p selection
    n=0;
    for(Int_t il(ly0); il<AliTRDgeometry::kNlayer; il+=dly){
      selection[n++] = il*11 + 0; // p bar
      selection[n++] = il*11 + 10; // p
    }
    xy[0] = 0.2; xy[1] = -.7; xy[2] = 7.; xy[3] = 8.;
    //xy[0] = 0.2; xy[1] = -1.5; xy[2] = 7.; xy[3] = 12.; // SA
    if(!GetGraphArray(xy, kMCtrack, 10, kTRUE, n, selection, "p & p bar")) break;
    pad->Modified(); pad->Update(); pad->SetLogx();
    pad = (TVirtualPad*)l->At(1); pad->cd();
    pad->SetMargin(0.125, 0.015, 0.1, 0.015);
    // e selection
    n=0;    
    for(Int_t il(ly0); il<AliTRDgeometry::kNlayer; il+=dly){
      selection[n++] = il*11 + 4; // e-
      selection[n++] = il*11 + 6; // e+
    }
    xy[0] = 0.2; xy[1] = -1.5; xy[2] = 7.; xy[3] = 12.;
    //xy[0] = 0.2; xy[1] = -1.5; xy[2] = 7.; xy[3] = 14.; // SA
    if(!GetGraphArray(xy, kMCtrack, 10, kTRUE, n, selection, "e#pm")) break;
    pad->Modified(); pad->Update(); pad->SetLogx();
    return kTRUE;
  case 37: // kMCtrackIn [y]
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    xy[0]=-.25; xy[1]=-1000.; xy[2]=.25; xy[3] =3000.;
    ((TVirtualPad*)l->At(0))->cd();
    selStart=0; for(n=0; n<fgkNresYsegm[fSegmentLevel]/3; n++) selection[n]=selStart+n;
    if(!GetGraphArray(xy, kMCtrackIn, 0, 1, n, selection)) break;
    ((TVirtualPad*)l->At(1))->cd();
    selStart=fgkNresYsegm[fSegmentLevel]/3; for(n=0; n<fgkNresYsegm[fSegmentLevel]/3; n++) selection[n]=selStart+n;
    if(!GetGraphArray(&xy[0], kMCtrackIn, 0, 1, n, selection)) break;
    return kTRUE;
  case 38: // kMCtrackIn [y]
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    xy[0]=-.25; xy[1]=-1000.; xy[2]=.25; xy[3] =3000.;
    ((TVirtualPad*)l->At(0))->cd();
    selStart=2*fgkNresYsegm[fSegmentLevel]/3; for(n=0; n<fgkNresYsegm[fSegmentLevel]/3; n++) selection[n]=selStart+n;
    if(!GetGraphArray(xy, kMCtrackIn, 0, 1, n, selection)) break;
    xy[0] = -.5; xy[1] = -0.5; xy[2] = fgkNresYsegm[fSegmentLevel]-.5; xy[3] = 2.5;
    ((TVirtualPad*)l->At(1))->cd();
    if(!GetGraphArray(xy, kMCtrackIn, 1, 1)) break;
    return kTRUE;
  case 39: // kMCtrackIn [z]
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    xy[0]=-1.; xy[1]=-500.; xy[2]=1.; xy[3] =800.;
    ((TVirtualPad*)l->At(0))->cd();
    if(!GetGraphArray(xy, kMCtrackIn, 2, 1)) break;
    xy[0] = -1.; xy[1] = -0.5; xy[2] = 1.; xy[3] = 2.5;
    ((TVirtualPad*)l->At(1))->cd();
    if(!GetGraphArray(xy, kMCtrackIn, 3, 1)) break;
    return kTRUE;
  case 40: // kMCtrackIn [phi|snp]
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    xy[0]=-.25; xy[1]=-0.5; xy[2]=.25; xy[3] =2.5;
    ((TVirtualPad*)l->At(0))->cd();
    if(!GetGraph(&xy[0], kMCtrackIn, 4)) break;
    xy[0] = -.25; xy[1] = -0.5; xy[2] = .25; xy[3] = 1.5;
    ((TVirtualPad*)l->At(1))->cd();
    if(!GetGraph(&xy[0], kMCtrackIn, 5)) break;
    return kTRUE;
  case 41: // kMCtrackIn [theta|tgl]
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    xy[0]=-1.; xy[1]=-1.; xy[2]=1.; xy[3] =4.;
    ((TVirtualPad*)l->At(0))->cd();
    if(!GetGraph(&xy[0], kMCtrackIn, 6)) break;
    xy[0] = -1.; xy[1] = -0.5; xy[2] = 1.; xy[3] = 1.5;
    ((TVirtualPad*)l->At(1))->cd();
    if(!GetGraph(&xy[0], kMCtrackIn, 7)) break;
    return kTRUE;
  case 42: // kMCtrackIn [pt]
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    xy[0] = 0.2; xy[1] = -.8; xy[2] = 7.; xy[3] = 6.;
    //xy[0] = 0.2; xy[1] = -1.5; xy[2] = 7.; xy[3] = 10.; // SA
    pad=(TVirtualPad*)l->At(0); pad->cd(); pad->SetLogx();
    pad->SetMargin(0.125, 0.015, 0.1, 0.015);
    n=0; selection[n++]=2; selection[n++]=3; selection[n++]=7; selection[n++]=8;
    if(!GetGraphArray(xy, kMCtrackIn, 8, 1, n, selection)) break;
    pad = (TVirtualPad*)l->At(1); pad->cd(); pad->SetLogx();
    pad->SetMargin(0.125, 0.015, 0.1, 0.015);
    n=0; selection[n++]=0; selection[n++]=4; selection[n++]=6; selection[n++]=10;
    if(!GetGraphArray(xy, kMCtrackIn, 8, 1, n, selection)) break;
    return kTRUE;
  case 43: //kMCtrackIn [1/pt] pulls
    xy[0] = 0.; xy[1] = -1.; xy[2] = 2.; xy[3] = 3.5;
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    pad = (TVirtualPad*)l->At(0); pad->cd();
    pad->SetMargin(0.125, 0.015, 0.1, 0.015);
    n=0; selection[n++]=2; selection[n++]=3; selection[n++]=7; selection[n++]=8;
    if(!GetGraphArray(xy, kMCtrackIn, 9, 1, n, selection)) break;
    pad = (TVirtualPad*)l->At(1); pad->cd();
    pad->SetMargin(0.125, 0.015, 0.1, 0.015);
    n=0; selection[n++]=0; selection[n++]=4; selection[n++]=6; selection[n++]=10;
    if(!GetGraphArray(xy, kMCtrackIn, 9, 1, n, selection)) break;
    return kTRUE;
  case 44: // kMCtrackIn [p]
    xy[0] = 0.2; xy[1] = -.8; xy[2] = 7.; xy[3] = 6.;
    //xy[0] = 0.2; xy[1] = -1.5; xy[2] = 7.; xy[3] = 10.;
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    pad = ((TVirtualPad*)l->At(0));pad->cd();pad->SetLogx();
    pad->SetMargin(0.125, 0.015, 0.1, 0.015);
    n=0; selection[n++]=2; selection[n++]=3; selection[n++]=7; selection[n++]=8;
    if(!GetGraphArray(xy, kMCtrackIn, 10, 1, n, selection)) break;
    pad = ((TVirtualPad*)l->At(1)); pad->cd();pad->SetLogx();
    pad->SetMargin(0.125, 0.015, 0.1, 0.015);
    n=0; selection[n++]=0; selection[n++]=4; selection[n++]=6; selection[n++]=10;
    if(!GetGraphArray(xy, kMCtrackIn, 10, 1,  n, selection)) break;
    return kTRUE;
  case 45: // kMCtrackOut [y]
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    xy[0]=-.3; xy[1]=-50.; xy[2]=.3; xy[3] =400.;
    ((TVirtualPad*)l->At(0))->cd();
    selStart=0; for(n=0; n<fgkNresYsegm[fSegmentLevel]/3; n++) selection[n]=selStart+n;
    if(!GetGraphArray(xy, kMCtrackOut, 0, 1, n, selection)) break;
    ((TVirtualPad*)l->At(1))->cd();
    selStart=fgkNresYsegm[fSegmentLevel]/3; for(n=0; n<fgkNresYsegm[fSegmentLevel]/3; n++) selection[n]=selStart+n;
    if(!GetGraphArray(&xy[0], kMCtrackOut, 0, 1, n, selection)) break;
    return kTRUE;
  case 46: // kMCtrackOut [y]
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    xy[0]=-.3; xy[1]=-50.; xy[2]=.3; xy[3] =400.;
    ((TVirtualPad*)l->At(0))->cd();
    selStart=2*fgkNresYsegm[fSegmentLevel]/3; for(n=0; n<fgkNresYsegm[fSegmentLevel]/3; n++) selection[n]=selStart+n;
    if(!GetGraphArray(xy, kMCtrackOut, 0, 1, n, selection)) break;
    xy[0] = -.5; xy[1] = -0.5; xy[2] = fgkNresYsegm[fSegmentLevel]-.5; xy[3] = 2.5;
    ((TVirtualPad*)l->At(1))->cd();
    if(!GetGraphArray(xy, kMCtrackOut, 1, 1)) break;
    return kTRUE;
  case 47: // kMCtrackOut [z]
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    xy[0]=-1.; xy[1]=-500.; xy[2]=1.; xy[3] =1500.;
    ((TVirtualPad*)l->At(0))->cd();
    if(!GetGraphArray(xy, kMCtrackOut, 2, 1)) break;
    xy[0] = -1.; xy[1] = -0.5; xy[2] = 1.; xy[3] = 2.5;
    ((TVirtualPad*)l->At(1))->cd();
    if(!GetGraphArray(xy, kMCtrackOut, 3, 1)) break;
    return kTRUE;
  case 48: // kMCtrackOut [phi|snp]
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    xy[0]=-.25; xy[1]=-0.5; xy[2]=.25; xy[3] =2.5;
    ((TVirtualPad*)l->At(0))->cd();
    if(!GetGraph(&xy[0], kMCtrackOut, 4)) break;
    xy[0] = -.25; xy[1] = -0.5; xy[2] = .25; xy[3] = 1.5;
    ((TVirtualPad*)l->At(1))->cd();
    if(!GetGraph(&xy[0], kMCtrackOut, 5)) break;
    return kTRUE;
  case 49: // kMCtrackOut [theta|tgl]
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    xy[0]=-1.; xy[1]=-1.; xy[2]=1.; xy[3] =4.;
    ((TVirtualPad*)l->At(0))->cd();
    if(!GetGraph(&xy[0], kMCtrackOut, 6)) break;
    xy[0] = -1.; xy[1] = -0.5; xy[2] = 1.; xy[3] = 15.;
    ((TVirtualPad*)l->At(1))->cd();
    if(!GetGraph(&xy[0], kMCtrackOut, 7)) break;
    return kTRUE;
  case 50: // kMCtrackOut [pt]
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    xy[0] = 0.2; xy[1] = -.8; xy[2] = 7.; xy[3] = 6.;
    pad=(TVirtualPad*)l->At(0); pad->cd(); pad->SetLogx();
    pad->SetMargin(0.125, 0.015, 0.1, 0.015);
    n=0; selection[n++]=2; selection[n++]=3; selection[n++]=7; selection[n++]=8;
    if(!GetGraphArray(xy, kMCtrackOut, 8, 1, n, selection)) break;
    pad = (TVirtualPad*)l->At(1); pad->cd();pad->SetLogx();
    pad->SetMargin(0.125, 0.015, 0.1, 0.015);
    n=0; selection[n++]=0; selection[n++]=4; selection[n++]=6; selection[n++]=10;
    if(!GetGraphArray(xy, kMCtrackOut, 8, 1, n, selection)) break;
    return kTRUE;
  case 51: //kMCtrackOut [1/pt] pulls
    xy[0] = 0.; xy[1] = -1.; xy[2] = 2.; xy[3] = 3.5;
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    pad = (TVirtualPad*)l->At(0); pad->cd();
    pad->SetMargin(0.125, 0.015, 0.1, 0.015);
    n=0; selection[n++]=2; selection[n++]=3; selection[n++]=7; selection[n++]=8;
    if(!GetGraphArray(xy, kMCtrackOut, 9, 1, n, selection)) break;
    pad = (TVirtualPad*)l->At(1); pad->cd();
    pad->SetMargin(0.125, 0.015, 0.1, 0.015);
    n=0; selection[n++]=0; selection[n++]=4; selection[n++]=6; selection[n++]=10;
    if(!GetGraphArray(xy, kMCtrackOut, 9, 1, n, selection)) break;
    return kTRUE;
  case 52: // kMCtrackOut [p]
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    xy[0] = 0.2; xy[1] = -.8; xy[2] = 7.; xy[3] = 6.;
    pad = ((TVirtualPad*)l->At(0));pad->cd();pad->SetLogx();
    pad->SetMargin(0.125, 0.015, 0.1, 0.015);
    n=0; selection[n++]=2; selection[n++]=3; selection[n++]=7; selection[n++]=8;
    if(!GetGraphArray(xy, kMCtrackOut, 10, 1, n, selection)) break;
    pad = ((TVirtualPad*)l->At(1)); pad->cd();pad->SetLogx();
    pad->SetMargin(0.125, 0.015, 0.1, 0.015);
    n=0; selection[n++]=0; selection[n++]=4; selection[n++]=6; selection[n++]=10;
    if(!GetGraphArray(xy, kMCtrackOut, 10, 1, n, selection)) break;
    return kTRUE;
  }
  AliWarning(Form("Reference plot [%d] missing result", ifig));
  return kFALSE;
}

Char_t const *fgParticle[11]={
  " p bar", " K -", " #pi -", " #mu -", " e -",
  " No PID",
  " e +", " #mu +", " #pi +", " K +", " p",
};
const Color_t fgColorS[11]={
kOrange, kOrange-3, kMagenta+1, kViolet, kRed,
kGray,
kRed, kViolet, kMagenta+1, kOrange-3, kOrange
};
const Color_t fgColorM[11]={
kCyan-5, kAzure-4, kBlue-7, kBlue+2, kViolet+10,
kBlack,
kViolet+10, kBlue+2, kBlue-7, kAzure-4, kCyan-5
};
const Marker_t fgMarker[11]={
30, 30, 26, 25, 24,
28,
20, 21, 22, 29, 29
};
//________________________________________________________
Bool_t AliTRDresolution::PostProcess()
{
  //fContainer = dynamic_cast<TObjArray*>(GetOutputData(0));
  if (!fContainer) {
    AliError("ERROR: list not available");
    return kFALSE;
  }
  TGraph *gm= NULL, *gs= NULL;
  if(!fGraphS && !fGraphM){ 
    TObjArray *aM(NULL), *aS(NULL);
    Int_t n = fContainer->GetEntriesFast();
    fGraphS = new TObjArray(n); fGraphS->SetOwner();
    fGraphM = new TObjArray(n); fGraphM->SetOwner();
    for(Int_t ig(0), nc(0); ig<n; ig++){
      fGraphM->AddAt(aM = new TObjArray(fgNproj[ig]), ig);
      fGraphS->AddAt(aS = new TObjArray(fgNproj[ig]), ig);
      
      for(Int_t ic=0; ic<fgNproj[ig]; ic++, nc++){
        AliDebug(2, Form("building G[%d] P[%d] N[%d]", ig, ic, fNcomp[nc]));
        if(fNcomp[nc]>1){
          TObjArray *agS(NULL), *agM(NULL);
          aS->AddAt(agS = new TObjArray(fNcomp[nc]), ic); 
          aM->AddAt(agM = new TObjArray(fNcomp[nc]), ic); 
          for(Int_t is=fNcomp[nc]; is--;){
            agS->AddAt(gs = new TGraphErrors(), is);
            Int_t is0(is%11), il0(is/11);
            gs->SetMarkerStyle(fgMarker[is0]);
            gs->SetMarkerColor(fgColorS[is0]);
            gs->SetLineColor(fgColorS[is0]);
            gs->SetLineStyle(il0);gs->SetLineWidth(2);
            gs->SetName(Form("s_%d_%02d_%02d", ig, ic, is));

            agM->AddAt(gm = new TGraphErrors(), is);
            gm->SetMarkerStyle(fgMarker[is0]);
            gm->SetMarkerColor(fgColorM[is0]);
            gm->SetLineColor(fgColorM[is0]);
            gm->SetLineStyle(il0);gm->SetLineWidth(2);
            gm->SetName(Form("m_%d_%02d_%02d", ig, ic, is));
            // this is important for labels in the legend 
            if(ic==0) {
              gs->SetTitle(Form("%s %02d", fgkResYsegmName[fSegmentLevel], is%fgkNresYsegm[fSegmentLevel]));
              gm->SetTitle(Form("%s %02d", fgkResYsegmName[fSegmentLevel], is%fgkNresYsegm[fSegmentLevel]));
            } else if(ic==1) {
              gs->SetTitle(Form("%s Ly[%d]", is%2 ?"z":"y", is/2));
              gm->SetTitle(Form("%s Ly[%d]", is%2?"z":"y", is/2));
            } else if(ic==2||ic==3) {
              gs->SetTitle(Form("%s Ly[%d]", is%2 ?"RC":"no RC", is/2));
              gm->SetTitle(Form("%s Ly[%d]", is%2?"RC":"no RC", is/2));
            } else if(ic<=7) {
              gs->SetTitle(Form("Layer[%d]", is%AliTRDgeometry::kNlayer));
              gm->SetTitle(Form("Layer[%d]", is%AliTRDgeometry::kNlayer));
            } else {
              gs->SetTitle(Form("%s @ ly[%d]", fgParticle[is0], il0));
              gm->SetTitle(Form("%s @ ly[%d]", fgParticle[is0], il0));
            }
          }
        } else {
          aS->AddAt(gs = new TGraphErrors(), ic);
          gs->SetMarkerStyle(23);
          gs->SetMarkerColor(kRed);
          gs->SetLineColor(kRed);
          gs->SetNameTitle(Form("s_%d_%02d", ig, ic), "sigma");
  
          aM->AddAt(gm = ig ? (TGraph*)new TGraphErrors() : (TGraph*)new TGraphAsymmErrors(), ic);
          gm->SetLineColor(kBlack);
          gm->SetMarkerStyle(7);
          gm->SetMarkerColor(kBlack);
          gm->SetNameTitle(Form("m_%d_%02d", ig, ic), "mean");
        }
      }
    }
  }

/*  printf("\n\n\n"); fGraphS->ls();
  printf("\n\n\n"); fGraphM->ls();*/
  

  // DEFINE MODELS
  // simple gauss
  TF1 fg("fGauss", "gaus", -.5, .5);  
  // Landau for charge resolution
  TF1 fch("fClCh", "landau", 0., 1000.);  
  // Landau for e+- pt resolution
  TF1 fpt("fPt", "landau", -0.1, 0.2);  

  //PROCESS EXPERIMENTAL DISTRIBUTIONS
  // Charge resolution
  //Process3DL(kCharge, 0, &fl); 
  // Clusters residuals
  Process3D(kCluster, 0, &fg, 1.e4); 
  Process3Dlinked(kCluster, 1, &fg); 
  fNRefFigures = 3;
  // Tracklet residual/pulls
  Process3D(kTrack , 0, &fg, 1.e4);
  Process3Dlinked(kTrack , 1, &fg); 
  Process3D(kTrack , 2, &fg, 1.e4); 
  Process3D(kTrack , 3, &fg); 
  Process2D(kTrack , 4, &fg, 1.e3);
  fNRefFigures = 7;
  // TRDin residual/pulls
  Process3D(kTrackIn, 0, &fg, 1.e4); 
  Process3Dlinked(kTrackIn, 1, &fg); 
  Process3D(kTrackIn, 2, &fg, 1.e4); 
  Process3D(kTrackIn, 3, &fg); 
  Process2D(kTrackIn, 4, &fg, 1.e3); 
  fNRefFigures = 11;
  // TRDout residual/pulls
  Process3D(kTrackOut, 0, &fg, 1.e3); // scale to fit - see PlotTrackOut
  Process3Dlinked(kTrackOut, 1, &fg); 
  Process3D(kTrackOut, 2, &fg, 1.e4); 
  Process3D(kTrackOut, 3, &fg); 
  Process2D(kTrackOut, 4, &fg, 1.e3); 
  fNRefFigures = 15;

  if(!HasMCdata()) return kTRUE;


  //PROCESS MC RESIDUAL DISTRIBUTIONS

  // CLUSTER Y RESOLUTION/PULLS
  Process3D(kMCcluster, 0, &fg, 1.e4);
  Process3Dlinked(kMCcluster, 1, &fg, 1.);
  fNRefFigures = 17;

  // TRACKLET RESOLUTION/PULLS
  Process3D(kMCtracklet, 0, &fg, 1.e4); // y
  Process3Dlinked(kMCtracklet, 1, &fg, 1.);   // y pulls
  Process3D(kMCtracklet, 2, &fg, 1.e4); // z
  Process3D(kMCtracklet, 3, &fg, 1.);   // z pulls
  Process2D(kMCtracklet, 4, &fg, 1.e3); // phi
  fNRefFigures = 21;

  // TRACK RESOLUTION/PULLS
  Process3Darray(kMCtrack, 0, &fg, 1.e4);   // y
  Process3DlinkedArray(kMCtrack, 1, &fg);   // y PULL
  Process3Darray(kMCtrack, 2, &fg, 1.e4);   // z
  Process3Darray(kMCtrack, 3, &fg);         // z PULL
  Process2Darray(kMCtrack, 4, &fg, 1.e3);   // phi
  Process2Darray(kMCtrack, 5, &fg);         // snp PULL
  Process2Darray(kMCtrack, 6, &fg, 1.e3);   // theta
  Process2Darray(kMCtrack, 7, &fg);         // tgl PULL
  Process3Darray(kMCtrack, 8, &fg, 1.e2);   // pt resolution
  Process3Darray(kMCtrack, 9, &fg);         // 1/pt pulls
  Process3Darray(kMCtrack, 10, &fg, 1.e2);  // p resolution
  fNRefFigures+=16;

  // TRACK TRDin RESOLUTION/PULLS
  Process3D(kMCtrackIn, 0, &fg, 1.e4);// y resolution
  Process3Dlinked(kMCtrackIn, 1, &fg);      // y pulls
  Process3D(kMCtrackIn, 2, &fg, 1.e4);// z resolution
  Process3D(kMCtrackIn, 3, &fg);      // z pulls
  Process2D(kMCtrackIn, 4, &fg, 1.e3);// phi resolution
  Process2D(kMCtrackIn, 5, &fg);      // snp pulls
  Process2D(kMCtrackIn, 6, &fg, 1.e3);// theta resolution
  Process2D(kMCtrackIn, 7, &fg);      // tgl pulls
  Process3D(kMCtrackIn, 8, &fg, 1.e2);// pt resolution
  Process3D(kMCtrackIn, 9, &fg);      // 1/pt pulls
  Process3D(kMCtrackIn, 10, &fg, 1.e2);// p resolution
  fNRefFigures+=8;

  // TRACK TRDout RESOLUTION/PULLS
  Process3D(kMCtrackOut, 0, &fg, 1.e4);// y resolution
  Process3Dlinked(kMCtrackOut, 1, &fg);      // y pulls
  Process3D(kMCtrackOut, 2, &fg, 1.e4);// z resolution
  Process3D(kMCtrackOut, 3, &fg);      // z pulls
  Process2D(kMCtrackOut, 4, &fg, 1.e3);// phi resolution
  Process2D(kMCtrackOut, 5, &fg);      // snp pulls
  Process2D(kMCtrackOut, 6, &fg, 1.e3);// theta resolution
  Process2D(kMCtrackOut, 7, &fg);      // tgl pulls
  Process3D(kMCtrackOut, 8, &fg, 1.e2);// pt resolution
  Process3D(kMCtrackOut, 9, &fg);      // 1/pt pulls
  Process3D(kMCtrackOut, 10, &fg, 1.e2);// p resolution
  fNRefFigures+=8;

  return kTRUE;
}


//________________________________________________________
void AliTRDresolution::Terminate(Option_t *opt)
{
  AliTRDrecoTask::Terminate(opt);
  if(HasPostProcess()) PostProcess();
}

//________________________________________________________
void AliTRDresolution::AdjustF1(TH1 *h, TF1 *f)
{
// Helper function to avoid duplication of code
// Make first guesses on the fit parameters

  // find the intial parameters of the fit !! (thanks George)
  Int_t nbinsy = Int_t(.5*h->GetNbinsX());
  Double_t sum = 0.;
  for(Int_t jbin=nbinsy-4; jbin<=nbinsy+4; jbin++) sum+=h->GetBinContent(jbin); sum/=9.;
  f->SetParLimits(0, 0., 3.*sum);
  f->SetParameter(0, .9*sum);
  Double_t rms = h->GetRMS();
  f->SetParLimits(1, -rms, rms);
  f->SetParameter(1, h->GetMean());

  f->SetParLimits(2, 0., 2.*rms);
  f->SetParameter(2, rms);
  if(f->GetNpar() <= 4) return;

  f->SetParLimits(3, 0., sum);
  f->SetParameter(3, .1*sum);

  f->SetParLimits(4, -.3, .3);
  f->SetParameter(4, 0.);

  f->SetParLimits(5, 0., 1.e2);
  f->SetParameter(5, 2.e-1);
}

//________________________________________________________
TObjArray* AliTRDresolution::BuildMonitorContainerCluster(const char* name, Bool_t expand)
{
// Build performance histograms for AliTRDcluster.vs TRD track or MC
//  - y reziduals/pulls

  TObjArray *arr = new TObjArray(2);
  arr->SetName(name); arr->SetOwner();
  TH1 *h(NULL); char hname[100], htitle[300];

  // tracklet resolution/pull in y direction
  sprintf(hname, "%s_%s_Y", GetNameId(), name);
  sprintf(htitle, "Y res for \"%s\" @ %s;tg(#phi);#Delta y [cm];%s", GetNameId(), name, fgkResYsegmName[fSegmentLevel]);
  if(!(h = (TH3S*)gROOT->FindObject(hname))){
    Int_t nybins=fgkNresYsegm[fSegmentLevel];
    if(expand) nybins*=2;
    h = new TH3S(hname, htitle, 
                 48, -.48, .48, 60, -.15, .15, nybins, -0.5, nybins-0.5);
  } else h->Reset();
  arr->AddAt(h, 0);
  sprintf(hname, "%s_%s_YZpull", GetNameId(), name);
  sprintf(htitle, "YZ pull for \"%s\" @ %s;%s;#Delta y  / #sigma_{y};#Delta z  / #sigma_{z}", GetNameId(), name, fgkResYsegmName[fSegmentLevel]);
  if(!(h = (TH3S*)gROOT->FindObject(hname))){
    h = new TH3S(hname, htitle, fgkNresYsegm[fSegmentLevel], -0.5, fgkNresYsegm[fSegmentLevel]-0.5, 100, -4.5, 4.5, 100, -4.5, 4.5);
  } else h->Reset();
  arr->AddAt(h, 1);

  return arr;
}

//________________________________________________________
TObjArray* AliTRDresolution::BuildMonitorContainerTracklet(const char* name, Bool_t expand)
{
// Build performance histograms for AliExternalTrackParam.vs TRD tracklet
//  - y reziduals/pulls
//  - z reziduals/pulls
//  - phi reziduals
  TObjArray *arr = BuildMonitorContainerCluster(name, expand); 
  arr->Expand(5);
  TH1 *h(NULL); char hname[100], htitle[300];

  // tracklet resolution/pull in z direction
  sprintf(hname, "%s_%s_Z", GetNameId(), name);
  sprintf(htitle, "Z res for \"%s\" @ %s;tg(#theta);#Delta z [cm];row cross", GetNameId(), name);
  if(!(h = (TH3S*)gROOT->FindObject(hname))){
    h = new TH3S(hname, htitle, 50, -1., 1., 100, -1.5, 1.5, 2, -0.5, 1.5);
  } else h->Reset();
  arr->AddAt(h, 2);
  sprintf(hname, "%s_%s_Zpull", GetNameId(), name);
  sprintf(htitle, "Z pull for \"%s\" @ %s;tg(#theta);#Delta z  / #sigma_{z};row cross", GetNameId(), name);
  if(!(h = (TH3S*)gROOT->FindObject(hname))){
    h = new TH3S(hname, htitle, 50, -1., 1., 100, -5.5, 5.5, 2, -0.5, 1.5);
    h->GetZaxis()->SetBinLabel(1, "no RC");
    h->GetZaxis()->SetBinLabel(2, "RC");
  } else h->Reset();
  arr->AddAt(h, 3);

  // tracklet to track phi resolution
  sprintf(hname, "%s_%s_PHI", GetNameId(), name);
  sprintf(htitle, "#Phi res for \"%s\" @ %s;tg(#phi);#Delta #phi [rad];entries", GetNameId(), name);
  if(!(h = (TH2I*)gROOT->FindObject(hname))){
    h = new TH2I(hname, htitle, 21, -.33, .33, 100, -.5, .5);
  } else h->Reset();
  arr->AddAt(h, 4);

  return arr;
}

//________________________________________________________
TObjArray* AliTRDresolution::BuildMonitorContainerTrack(const char* name)
{
// Build performance histograms for AliExternalTrackParam.vs MC
//  - y resolution/pulls
//  - z resolution/pulls
//  - phi resolution, snp pulls
//  - theta resolution, tgl pulls
//  - pt resolution, 1/pt pulls, p resolution

  TObjArray *arr = BuildMonitorContainerTracklet(name); 
  arr->Expand(11);
  TH1 *h(NULL); char hname[100], htitle[300];
  TAxis *ax(NULL);

  // snp pulls
  sprintf(hname, "%s_%s_SNPpull", GetNameId(), name);
  sprintf(htitle, "SNP pull for \"%s\" @ %s;tg(#phi);#Delta snp  / #sigma_{snp};entries", GetNameId(), name);
  if(!(h = (TH2I*)gROOT->FindObject(hname))){
    h = new TH2I(hname, htitle, 60, -.3, .3, 100, -4.5, 4.5);
  } else h->Reset();
  arr->AddAt(h, 5);

  // theta resolution
  sprintf(hname, "%s_%s_THT", GetNameId(), name);
  sprintf(htitle, "#Theta res for \"%s\" @ %s;tg(#theta);#Delta #theta [rad];entries", GetNameId(), name);
  if(!(h = (TH2I*)gROOT->FindObject(hname))){
    h = new TH2I(hname, htitle, 100, -1., 1., 100, -5e-3, 5e-3);
  } else h->Reset();
  arr->AddAt(h, 6);
  // tgl pulls
  sprintf(hname, "%s_%s_TGLpull", GetNameId(), name);
  sprintf(htitle, "TGL pull for \"%s\" @ %s;tg(#theta);#Delta tgl  / #sigma_{tgl};entries", GetNameId(), name);
  if(!(h = (TH2I*)gROOT->FindObject(hname))){
    h = new TH2I(hname, htitle, 100, -1., 1., 100, -4.5, 4.5);
  } else h->Reset();
  arr->AddAt(h, 7);
  
  const Int_t kNpt(14);
  const Int_t kNdpt(150); 
  const Int_t kNspc = 2*AliPID::kSPECIES+1;
  Float_t Pt=0.1, DPt=-.1, Spc=-5.5;
  Float_t binsPt[kNpt+1], binsSpc[kNspc+1], binsDPt[kNdpt+1];
  for(Int_t i=0;i<kNpt+1; i++,Pt=TMath::Exp(i*.15)-1.) binsPt[i]=Pt;
  for(Int_t i=0; i<kNspc+1; i++,Spc+=1.) binsSpc[i]=Spc;
  for(Int_t i=0; i<kNdpt+1; i++,DPt+=2.e-3) binsDPt[i]=DPt;

  // Pt resolution
  sprintf(hname, "%s_%s_Pt", GetNameId(), name);
  sprintf(htitle, "P_{t} res for \"%s\" @ %s;p_{t} [GeV/c];#Delta p_{t}/p_{t}^{MC};SPECIES", GetNameId(), name);
  if(!(h = (TH3S*)gROOT->FindObject(hname))){
    h = new TH3S(hname, htitle, 
                 kNpt, binsPt, kNdpt, binsDPt, kNspc, binsSpc);
    ax = h->GetZaxis();
    for(Int_t ib(1); ib<=ax->GetNbins(); ib++) ax->SetBinLabel(ib, fgParticle[ib-1]);
  } else h->Reset();
  arr->AddAt(h, 8);
  // 1/Pt pulls
  sprintf(hname, "%s_%s_1Pt", GetNameId(), name);
  sprintf(htitle, "1/P_{t} pull for \"%s\" @ %s;1/p_{t}^{MC} [c/GeV];#Delta(1/p_{t})/#sigma(1/p_{t});SPECIES", GetNameId(), name);
  if(!(h = (TH3S*)gROOT->FindObject(hname))){
    h = new TH3S(hname, htitle, 
                 kNpt, 0., 2., 100, -4., 4., kNspc, -5.5, 5.5);
    ax = h->GetZaxis();
    for(Int_t ib(1); ib<=ax->GetNbins(); ib++) ax->SetBinLabel(ib, fgParticle[ib-1]);
  } else h->Reset();
  arr->AddAt(h, 9);
  // P resolution
  sprintf(hname, "%s_%s_P", GetNameId(), name);
  sprintf(htitle, "P res for \"%s\" @ %s;p [GeV/c];#Delta p/p^{MC};SPECIES", GetNameId(), name);
  if(!(h = (TH3S*)gROOT->FindObject(hname))){
    h = new TH3S(hname, htitle, 
                 kNpt, binsPt, kNdpt, binsDPt, kNspc, binsSpc);
    ax = h->GetZaxis();
    for(Int_t ib(1); ib<=ax->GetNbins(); ib++) ax->SetBinLabel(ib, fgParticle[ib-1]);
  } else h->Reset();
  arr->AddAt(h, 10);

  return arr;
}


//________________________________________________________
TObjArray* AliTRDresolution::Histos()
{
  //
  // Define histograms
  //

  if(fContainer) return fContainer;

  fContainer  = new TObjArray(kNviews);
  //fContainer->SetOwner(kTRUE);
  TH1 *h(NULL);
  TObjArray *arr(NULL);

  // binnings for plots containing momentum or pt
  const Int_t kNpt(14), kNphi(48), kNdy(60);
  Float_t Phi=-.48, Dy=-.3, Pt=0.1;
  Float_t binsPhi[kNphi+1], binsDy[kNdy+1], binsPt[kNpt+1];
  for(Int_t i=0; i<kNphi+1; i++,Phi+=.02) binsPhi[i]=Phi;
  for(Int_t i=0; i<kNdy+1; i++,Dy+=.01) binsDy[i]=Dy;
  for(Int_t i=0;i<kNpt+1; i++,Pt=TMath::Exp(i*.15)-1.) binsPt[i]=Pt;

  // cluster to track residuals/pulls
  fContainer->AddAt(arr = new TObjArray(2), kCharge);
  arr->SetName("Charge");
  if(!(h = (TH3S*)gROOT->FindObject("hCharge"))){
    h = new TH3S("hCharge", "Charge Resolution", 20, 1., 2., 24, 0., 3.6, 100, 0., 500.);
    h->GetXaxis()->SetTitle("dx/dx_{0}");
    h->GetYaxis()->SetTitle("x_{d} [cm]");
    h->GetZaxis()->SetTitle("dq/dx [ADC/cm]");
  } else h->Reset();
  arr->AddAt(h, 0);

  // cluster to track residuals/pulls
  fContainer->AddAt(BuildMonitorContainerCluster("Cl"), kCluster);
  // tracklet to TRD track
  fContainer->AddAt(BuildMonitorContainerTracklet("Trk", kTRUE), kTrack);
  // tracklet to TRDin
  fContainer->AddAt(BuildMonitorContainerTracklet("TrkIN", kTRUE), kTrackIn);
  // tracklet to TRDout
  fContainer->AddAt(BuildMonitorContainerTracklet("TrkOUT"), kTrackOut);


  // Resolution histos
  if(!HasMCdata()) return fContainer;

  // cluster resolution 
  fContainer->AddAt(BuildMonitorContainerCluster("MCcl"),  kMCcluster);

  // tracklet resolution
  fContainer->AddAt(BuildMonitorContainerTracklet("MCtracklet"), kMCtracklet);

  // track resolution
  fContainer->AddAt(arr = new TObjArray(AliTRDgeometry::kNlayer), kMCtrack);
  arr->SetName("MCtrk");
  for(Int_t il(0); il<AliTRDgeometry::kNlayer; il++) arr->AddAt(BuildMonitorContainerTrack(Form("MCtrk_Ly%d", il)), il);

  // TRDin TRACK RESOLUTION
  fContainer->AddAt(BuildMonitorContainerTrack("MCtrkIN"), kMCtrackIn);

  // TRDout TRACK RESOLUTION
  fContainer->AddAt(BuildMonitorContainerTrack("MCtrkOUT"), kMCtrackOut);

  return fContainer;
}

//________________________________________________________
Bool_t AliTRDresolution::Load(const Char_t *file, const Char_t *dir)
{
// Custom load function. Used to guess the segmentation level of the data.

  if(!AliTRDrecoTask::Load(file, dir)) return kFALSE;

  // look for cluster residual plot - always available
  TH3S* h3((TH3S*)((TObjArray*)fContainer->At(kClToTrk))->At(0));
  Int_t segmentation(h3->GetNbinsZ()/2);
  if(segmentation==fgkNresYsegm[0]){ // default segmentation. Nothing to do
    return kTRUE;
  } else if(segmentation==fgkNresYsegm[1]){ // stack segmentation.
    SetSegmentationLevel(1);
  } else if(segmentation==fgkNresYsegm[2]){ // detector segmentation.
    SetSegmentationLevel(2);
  } else {
    AliError(Form("Unknown segmentation [%d].", h3->GetNbinsZ()));
    return kFALSE;
  }

  AliDebug(2, Form("Segmentation set to level \"%s\"", fgkResYsegmName[fSegmentLevel]));
  return kTRUE;
}


//________________________________________________________
Bool_t AliTRDresolution::Process(TH2 * const h2, TF1 *f, Float_t k, TGraphErrors **g)
{
  //
  // Do the processing
  //

  Char_t pn[10]; sprintf(pn, "p%03d", fIdxPlot);
  Int_t n = 0;
  if((n=g[0]->GetN())) for(;n--;) g[0]->RemovePoint(n);
  if((n=g[1]->GetN())) for(;n--;) g[1]->RemovePoint(n);
  if(Int_t(h2->GetEntries())){ 
    AliDebug(4, Form("%s: g[%s %s]", pn, g[0]->GetName(), g[0]->GetTitle()));
  } else {
    AliDebug(2, Form("%s: g[%s %s]: Missing entries.", pn, g[0]->GetName(), g[0]->GetTitle()));
    fIdxPlot++;
    return kTRUE;
  }

  const Int_t kINTEGRAL=1;
  for(Int_t ibin = 0; ibin < Int_t(h2->GetNbinsX()/kINTEGRAL); ibin++){
    Int_t abin(ibin*kINTEGRAL+1),bbin(abin+kINTEGRAL-1),mbin(abin+Int_t(kINTEGRAL/2));
    Double_t x = h2->GetXaxis()->GetBinCenter(mbin);
    TH1D *h = h2->ProjectionY(pn, abin, bbin);
    if((n=(Int_t)h->GetEntries())<100){ 
      AliDebug(4, Form("  x[%f] range[%d %d] stat[%d] low statistics !", x, abin, bbin, n));
      continue;
    }
    h->Fit(f, "QN");
    Int_t ip = g[0]->GetN();
    AliDebug(4, Form("  x_%d[%f] range[%d %d] stat[%d] M[%f] Sgm[%f]", ip, x, abin, bbin, n, f->GetParameter(1), f->GetParameter(2)));
    g[0]->SetPoint(ip, x, k*f->GetParameter(1));
    g[0]->SetPointError(ip, 0., k*f->GetParError(1));
    g[1]->SetPoint(ip, x, k*f->GetParameter(2));
    g[1]->SetPointError(ip, 0., k*f->GetParError(2));
/*  
    g[0]->SetPoint(ip, x, k*h->GetMean());
    g[0]->SetPointError(ip, 0., k*h->GetMeanError());
    g[1]->SetPoint(ip, x, k*h->GetRMS());
    g[1]->SetPointError(ip, 0., k*h->GetRMSError());*/
  }
  fIdxPlot++;
  return kTRUE;
}

//________________________________________________________
Bool_t AliTRDresolution::Process2D(ETRDresolutionPlot plot, Int_t idx, TF1 *f, Float_t k, Int_t gidx)
{
  //
  // Do the processing
  //

  if(!fContainer || !fGraphS || !fGraphM) return kFALSE;

  // retrive containers
  TH2I *h2(NULL);
  if(idx<0){
    if(!(h2= (TH2I*)(fContainer->At(plot)))) return kFALSE; 
  } else{
    TObjArray *a0(NULL);
    if(!(a0=(TObjArray*)(fContainer->At(plot)))) return kFALSE;
    if(!(h2=(TH2I*)a0->At(idx))) return kFALSE;
  }
  if(Int_t(h2->GetEntries())){ 
    AliDebug(2, Form("p[%d] idx[%d] : h[%s] %s", plot, idx, h2->GetName(), h2->GetTitle()));
  } else {
    AliDebug(2, Form("p[%d] idx[%d] : Missing entries.", plot, idx));
    return kFALSE;
  }

  TGraphErrors *g[2];
  if(gidx<0) gidx=idx;
  if(!(g[0] = gidx<0 ? (TGraphErrors*)fGraphM->At(plot) : (TGraphErrors*)((TObjArray*)(fGraphM->At(plot)))->At(gidx))) return kFALSE;

  if(!(g[1] = gidx<0 ? (TGraphErrors*)fGraphS->At(plot) : (TGraphErrors*)((TObjArray*)(fGraphS->At(plot)))->At(gidx))) return kFALSE;

  return Process(h2, f, k, g);
}

//________________________________________________________
Bool_t AliTRDresolution::Process3D(ETRDresolutionPlot plot, Int_t idx, TF1 *f, Float_t k)
{
  //
  // Do the processing
  //

  if(!fContainer || !fGraphS || !fGraphM) return kFALSE;

  // retrive containers
  TH3S *h3(NULL);
  if(idx<0){
    if(!(h3= (TH3S*)(fContainer->At(plot)))) return kFALSE; 
  } else{
    TObjArray *a0(NULL);
    if(!(a0=(TObjArray*)(fContainer->At(plot)))) return kFALSE;
    if(!(h3=(TH3S*)a0->At(idx))) return kFALSE;
  }
  if(Int_t(h3->GetEntries())){ 
    AliDebug(2, Form("p[%d] idx[%d] h[%s] %s", plot, idx, h3->GetName(), h3->GetTitle()));
  } else {
    AliDebug(2, Form("p[%d] idx[%d] : Missing entries.", plot, idx));
    return kFALSE;
  }

  TObjArray *gm, *gs;
  if(!(gm = (TObjArray*)((TObjArray*)(fGraphM->At(plot)))->At(idx))) return kFALSE;
  if(!(gs = (TObjArray*)((TObjArray*)(fGraphS->At(plot)))->At(idx))) return kFALSE;
  TGraphErrors *g[2];

  TAxis *az = h3->GetZaxis();
  for(Int_t iz(0); iz<gm->GetEntriesFast(); iz++){
    if(!(g[0] = (TGraphErrors*)gm->At(iz))) return kFALSE;
    if(!(g[1] = (TGraphErrors*)gs->At(iz))) return kFALSE;
    az->SetRange(iz+1, iz+1);
    if(!Process((TH2*)h3->Project3D("yx"), f, k, g)) return kFALSE;
  }

  return kTRUE;
}


//________________________________________________________
Bool_t AliTRDresolution::Process3Dlinked(ETRDresolutionPlot plot, Int_t idx, TF1 *f, Float_t k)
{
  //
  // Do the processing
  //

  if(!fContainer || !fGraphS || !fGraphM) return kFALSE;

  // retrive containers
  TH3S *h3(NULL);
  if(idx<0){
    if(!(h3= (TH3S*)(fContainer->At(plot)))) return kFALSE; 
  } else{
    TObjArray *a0(NULL);
    if(!(a0=(TObjArray*)(fContainer->At(plot)))) return kFALSE;
    if(!(h3=(TH3S*)a0->At(idx))) return kFALSE;
  }
  if(Int_t(h3->GetEntries())){ 
    AliDebug(2, Form("p[%d] idx[%d] h[%s] %s", plot, idx, h3->GetName(), h3->GetTitle()));
  } else {
    AliDebug(2, Form("p[%d] idx[%d] : Missing entries.", plot, idx));
    return kFALSE;
  }

  TObjArray *gm, *gs;
  if(!(gm = (TObjArray*)((TObjArray*)(fGraphM->At(plot)))->At(idx))) return kFALSE;
  if(!(gs = (TObjArray*)((TObjArray*)(fGraphS->At(plot)))->At(idx))) return kFALSE;
  TGraphErrors *g[2];

  if(!(g[0] = (TGraphErrors*)gm->At(0))) return kFALSE;
  if(!(g[1] = (TGraphErrors*)gs->At(0))) return kFALSE;
  if(!Process((TH2*)h3->Project3D("yx"), f, k, g)) return kFALSE;

  if(!(g[0] = (TGraphErrors*)gm->At(1))) return kFALSE;
  if(!(g[1] = (TGraphErrors*)gs->At(1))) return kFALSE;
  if(!Process((TH2*)h3->Project3D("zx"), f, k, g)) return kFALSE;

  return kTRUE;
}


//________________________________________________________
Bool_t AliTRDresolution::Process3DL(ETRDresolutionPlot plot, Int_t idx, TF1 *f, Float_t k)
{
  //
  // Do the processing
  //

  if(!fContainer || !fGraphS || !fGraphM) return kFALSE;

  // retrive containers
  TH3S *h3 = (TH3S*)((TObjArray*)fContainer->At(plot))->At(idx);
  if(!h3) return kFALSE;
  AliDebug(2, Form("p[%d] idx[%d] h[%s] %s", plot, idx, h3->GetName(), h3->GetTitle()));

  TGraphAsymmErrors *gm; 
  TGraphErrors *gs;
  if(!(gm = (TGraphAsymmErrors*)((TObjArray*)fGraphM->At(plot))->At(0))) return kFALSE;
  if(!(gs = (TGraphErrors*)((TObjArray*)fGraphS->At(plot)))) return kFALSE;

  Float_t x, r, mpv, xM, xm;
  TAxis *ay = h3->GetYaxis();
  for(Int_t iy=1; iy<=h3->GetNbinsY(); iy++){
    ay->SetRange(iy, iy);
    x = ay->GetBinCenter(iy);
    TH2F *h2=(TH2F*)h3->Project3D("zx");
    TAxis *ax = h2->GetXaxis();
    for(Int_t ix=1; ix<=h2->GetNbinsX(); ix++){
      r = ax->GetBinCenter(ix);
      TH1D *h1 = h2->ProjectionY("py", ix, ix);
      if(h1->Integral()<50) continue;
      h1->Fit(f, "QN");

      GetLandauMpvFwhm(f, mpv, xm, xM);
      Int_t ip = gm->GetN();
      gm->SetPoint(ip, x, k*mpv);
      gm->SetPointError(ip, 0., 0., k*xm, k*xM);
      gs->SetPoint(ip, r, k*(xM-xm)/mpv);
      gs->SetPointError(ip, 0., 0.);
    }
  }

  return kTRUE;
}

//________________________________________________________
Bool_t AliTRDresolution::Process2Darray(ETRDresolutionPlot plot, Int_t idx, TF1 *f, Float_t k)
{
  //
  // Do the processing
  //

  if(!fContainer || !fGraphS || !fGraphM) return kFALSE;

  // retrive containers
  TObjArray *arr = (TObjArray*)(fContainer->At(plot));
  if(!arr) return kFALSE;
  AliDebug(2, Form("p[%d] idx[%d] arr[%s]", plot, idx, arr->GetName()));

  TObjArray *gm, *gs;
  if(!(gm = (TObjArray*)((TObjArray*)(fGraphM->At(plot)))->At(idx))) return kFALSE;
  if(!(gs = (TObjArray*)((TObjArray*)(fGraphS->At(plot)))->At(idx))) return kFALSE;

  TGraphErrors *g[2]; TH2I *h2(NULL); TObjArray *a0(NULL);
  for(Int_t ia(0); ia<arr->GetEntriesFast(); ia++){
    if(!(a0 = (TObjArray*)arr->At(ia))) continue;

    if(!(h2 = (TH2I*)a0->At(idx))) return kFALSE;
    if(Int_t(h2->GetEntries())){ 
      AliDebug(4, Form("   idx[%d] h[%s] %s", ia, h2->GetName(), h2->GetTitle()));
    } else {
      AliDebug(2, Form("   idx[%d] : Missing entries.", ia));
      continue;
    }

    if(!(g[0] = (TGraphErrors*)gm->At(ia))) return kFALSE;
    if(!(g[1] = (TGraphErrors*)gs->At(ia))) return kFALSE;
    if(!Process(h2, f, k, g)) return kFALSE;
  }

  return kTRUE;
}

//________________________________________________________
Bool_t AliTRDresolution::Process3Darray(ETRDresolutionPlot plot, Int_t idx, TF1 *f, Float_t k)
{
  //
  // Do the processing
  //

  if(!fContainer || !fGraphS || !fGraphM) return kFALSE;
  //printf("Process4D : processing plot[%d] idx[%d]\n", plot, idx);

  // retrive containers
  TObjArray *arr = (TObjArray*)(fContainer->At(plot));
  if(!arr) return kFALSE;
  AliDebug(2, Form("p[%d] idx[%d] arr[%s]", plot, idx, arr->GetName()));

  TObjArray *gm, *gs;
  if(!(gm = (TObjArray*)((TObjArray*)(fGraphM->At(plot)))->At(idx))) return kFALSE;
  if(!(gs = (TObjArray*)((TObjArray*)(fGraphS->At(plot)))->At(idx))) return kFALSE;

  TGraphErrors *g[2]; TH3S *h3(NULL); TObjArray *a0(NULL);
  Int_t in(0);
  for(Int_t ia(0); ia<arr->GetEntriesFast(); ia++){
    if(!(a0 = (TObjArray*)arr->At(ia))) continue;

    if(!(h3 = (TH3S*)a0->At(idx))) return kFALSE;
     if(Int_t(h3->GetEntries())){ 
       AliDebug(4, Form("   idx[%d] h[%s] %s", ia, h3->GetName(), h3->GetTitle()));
    } else {
      AliDebug(2, Form("   idx[%d] : Missing entries.", ia));
      continue;
    }
    TAxis *az = h3->GetZaxis();
    for(Int_t iz=1; iz<=az->GetNbins(); iz++, in++){
      if(in >= gm->GetEntriesFast()) break;
      if(!(g[0] = (TGraphErrors*)gm->At(in))) return kFALSE;
      if(!(g[1] = (TGraphErrors*)gs->At(in))) return kFALSE;
      az->SetRange(iz, iz);
      if(!Process((TH2*)h3->Project3D("yx"), f, k, g)) return kFALSE;
    }
  }
  AliDebug(2, Form("Projections [%d] from [%d]", in, gs->GetEntriesFast()));

  return kTRUE;
}

//________________________________________________________
Bool_t AliTRDresolution::Process3DlinkedArray(ETRDresolutionPlot plot, Int_t idx, TF1 *f, Float_t k)
{
  //
  // Do the processing
  //

  if(!fContainer || !fGraphS || !fGraphM) return kFALSE;
  //printf("Process4D : processing plot[%d] idx[%d]\n", plot, idx);

  // retrive containers
  TObjArray *arr = (TObjArray*)(fContainer->At(plot));
  if(!arr) return kFALSE;
  AliDebug(2, Form("p[%d] idx[%d] arr[%s]", plot, idx, arr->GetName()));

  TObjArray *gm, *gs;
  if(!(gm = (TObjArray*)((TObjArray*)(fGraphM->At(plot)))->At(idx))) return kFALSE;
  if(!(gs = (TObjArray*)((TObjArray*)(fGraphS->At(plot)))->At(idx))) return kFALSE;

  TGraphErrors *g[2]; TH3S *h3(NULL); TObjArray *a0(NULL);
  Int_t in(0);
  for(Int_t ia(0); ia<arr->GetEntriesFast(); ia++){
    if(!(a0 = (TObjArray*)arr->At(ia))) continue;
    if(!(h3 = (TH3S*)a0->At(idx))) return kFALSE;
    if(Int_t(h3->GetEntries())){     
      AliDebug(4, Form("   idx[%d] h[%s] %s", ia, h3->GetName(), h3->GetTitle()));
    } else {
      AliDebug(2, Form("   idx[%d] : Missing entries.", ia));
      continue;
    }
    if(!(g[0] = (TGraphErrors*)gm->At(in))) return kFALSE;
    if(!(g[1] = (TGraphErrors*)gs->At(in))) return kFALSE;
    if(!Process((TH2*)h3->Project3D("yx"), f, k, g)) return kFALSE;
    in++;

    if(!(g[0] = (TGraphErrors*)gm->At(in))) return kFALSE;
    if(!(g[1] = (TGraphErrors*)gs->At(in))) return kFALSE;
    if(!Process((TH2*)h3->Project3D("zx"), f, k, g)) return kFALSE;
    in++;
  }
  AliDebug(2, Form("Projections [%d] from [%d]", in, gs->GetEntriesFast()));

  return kTRUE;
}

//________________________________________________________
Bool_t AliTRDresolution::GetGraph(Float_t *bb, ETRDresolutionPlot ip, Int_t idx, Bool_t kLEG, const Char_t *explain)
{
  //
  // Get the graphs
  //

  if(!fGraphS || !fGraphM) return kFALSE;
  // axis titles look up
  Int_t nref = 0;
  for(Int_t jp=0; jp<(Int_t)ip; jp++) nref+=fgNproj[jp];
  UChar_t jdx = idx<0?0:idx;
  for(Int_t jc=0; jc<TMath::Min(jdx,fgNproj[ip]-1); jc++) nref++;
  Char_t **at = fAxTitle[nref];

  // build legends if requiered
  TLegend *leg(NULL);
  if(kLEG){
    leg=new TLegend(.65, .6, .95, .9);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
  }
  // build frame
  TH1S *h1(NULL);
  h1 = new TH1S(Form("h1TF_%02d", fIdxFrame++), Form("%s %s;%s;%s", at[0], explain?explain:"", at[1], at[2]), 2, bb[0], bb[2]);
  h1->SetMinimum(bb[1]);h1->SetMaximum(bb[3]);
  h1->SetLineColor(kBlack); h1->SetLineWidth(1);h1->SetLineStyle(2); 
  // axis range
  TAxis *ax = h1->GetXaxis();
  ax->CenterTitle();ax->SetMoreLogLabels();ax->SetTitleOffset(1.2);
  ax = h1->GetYaxis();
  ax->SetRangeUser(bb[1], bb[3]);
  ax->CenterTitle(); ax->SetTitleOffset(1.4);
  h1->Draw();
  // bounding box
  TBox *b = new TBox(-.15, bb[1], .15, bb[3]);
  b->SetFillStyle(3002);b->SetLineColor(0);
  b->SetFillColor(ip<=Int_t(kMCcluster)?kGreen:kBlue);
  b->Draw();

  TGraphErrors *gm = idx<0 ? (TGraphErrors*)fGraphM->At(ip) : (TGraphErrors*)((TObjArray*)(fGraphM->At(ip)))->At(idx);
  if(!gm) return kFALSE;
  TGraphErrors *gs = idx<0 ? (TGraphErrors*)fGraphS->At(ip) : (TGraphErrors*)((TObjArray*)(fGraphS->At(ip)))->At(idx);
  if(!gs) return kFALSE;

  Int_t n(0), nPlots(0);
  if((n=gm->GetN())) {
    nPlots++;
    gm->Draw("pl"); if(leg) leg->AddEntry(gm, gm->GetTitle(), "pl");
    PutTrendValue(Form("%s_%s", fgPerformanceName[ip], at[0]), gm->GetMean(2));
    PutTrendValue(Form("%s_%sRMS", fgPerformanceName[ip], at[0]), gm->GetRMS(2));
  }

  if((n=gs->GetN())){
    nPlots++;
    gs->Draw("pl"); if(leg) leg->AddEntry(gs, gs->GetTitle(), "pl");
    gs->Sort(&TGraph::CompareY);
    PutTrendValue(Form("%s_%sSigMin", fgPerformanceName[ip], at[0]), gs->GetY()[0]);
    PutTrendValue(Form("%s_%sSigMax", fgPerformanceName[ip], at[0]), gs->GetY()[n-1]);
    gs->Sort(&TGraph::CompareX); 
  }
  if(!nPlots) return kFALSE;
  if(leg) leg->Draw();
  return kTRUE;
}

//________________________________________________________
Bool_t AliTRDresolution::GetGraphArray(Float_t *bb, ETRDresolutionPlot ip, Int_t idx, Bool_t kLEG, Int_t n, Int_t *sel, const Char_t *explain)
{
  //
  // Get the graphs
  //

  if(!fGraphS || !fGraphM) return kFALSE;

  // axis titles look up
  Int_t nref(0);
  for(Int_t jp(0); jp<ip; jp++) nref+=fgNproj[jp];
  nref+=idx;
  Char_t **at = fAxTitle[nref];

  // build legends if requiered
  TLegend *legM(NULL), *legS(NULL);
  if(kLEG){
    legM=new TLegend(.35, .6, .65, .9);
    legM->SetHeader("Mean");
    legM->SetBorderSize(0);
    legM->SetFillStyle(0);
    legS=new TLegend(.65, .6, .95, .9);
    legS->SetHeader("Sigma");
    legS->SetBorderSize(0);
    legS->SetFillStyle(0);
  }
  // build frame
  TH1S *h1(NULL);
  h1 = new TH1S(Form("h1TF_%02d", fIdxFrame++), Form("%s %s;%s;%s", at[0], explain?explain:"", at[1], at[2]), 2, bb[0], bb[2]);
  h1->SetMinimum(bb[1]);h1->SetMaximum(bb[3]);
  h1->SetLineColor(kBlack); h1->SetLineWidth(1);h1->SetLineStyle(2); 
  // axis range
  TAxis *ax = h1->GetXaxis();
  ax->CenterTitle();ax->SetMoreLogLabels();ax->SetTitleOffset(1.2);
  ax = h1->GetYaxis();
  ax->SetRangeUser(bb[1], bb[3]);
  ax->CenterTitle(); ax->SetTitleOffset(1.4);
  h1->Draw();

  TGraphErrors *gm(NULL), *gs(NULL);
  TObjArray *a0(NULL), *a1(NULL);
  a0 = (TObjArray*)((TObjArray*)fGraphM->At(ip))->At(idx);
  a1 = (TObjArray*)((TObjArray*)fGraphS->At(ip))->At(idx);
  if(!n) n=a0->GetEntriesFast();
  AliDebug(4, Form("Graph : Ref[%d] Title[%s] Limits{x[%f %f] y[%f %f]} Comp[%d] Selection[%c]", nref, at[0], bb[0], bb[2], bb[1], bb[3], n, sel ? 'y' : 'n'));
  Int_t nn(0), nPlots(0);
  for(Int_t is(0), is0(0); is<n; is++){
    is0 = sel ? sel[is] : is;
    if(!(gs =  (TGraphErrors*)a1->At(is0))) return kFALSE;
    if(!(gm =  (TGraphErrors*)a0->At(is0))) return kFALSE;

    if((nn=gs->GetN())){
      nPlots++;
      gs->Draw("pc"); 
      if(legS){ 
        //printf("LegEntry %s [%s]%s\n", at[0], gs->GetName(), gs->GetTitle());
        legS->AddEntry(gs, gs->GetTitle(), "pl");
      }
      gs->Sort(&TGraph::CompareY);
      PutTrendValue(Form("%s_%sSigMin", fgPerformanceName[kMCtrack], at[0]), gs->GetY()[0]);
      PutTrendValue(Form("%s_%sSigMax", fgPerformanceName[kMCtrack], at[0]), gs->GetY()[nn-1]);
      gs->Sort(&TGraph::CompareX); 
    }
    if(gm->GetN()){
      nPlots++;
      gm->Draw("pc");
      if(legM){ 
        //printf("LegEntry %s [%s]%s\n", at[0], gm->GetName(), gm->GetTitle());
        legM->AddEntry(gm, gm->GetTitle(), "pl");
      }
      PutTrendValue(Form("%s_%s", fgPerformanceName[kMCtrack], at[0]), gm->GetMean(2));
      PutTrendValue(Form("%s_%sRMS", fgPerformanceName[kMCtrack], at[0]), gm->GetRMS(2));
    }
  }
  if(!nPlots) return kFALSE;
  if(kLEG){
    legM->Draw();
    legS->Draw();
  }
  return kTRUE;
}

//________________________________________________________
void AliTRDresolution::GetLandauMpvFwhm(TF1 * const f, Float_t &mpv, Float_t &xm, Float_t &xM)
{
  //
  // Get the most probable value and the full width half mean 
  // of a Landau distribution
  //

  const Float_t dx = 1.;
  mpv = f->GetParameter(1);
  Float_t fx, max = f->Eval(mpv);

  xm = mpv - dx;
  while((fx = f->Eval(xm))>.5*max){
    if(fx>max){ 
      max = fx;
      mpv = xm;
    }
    xm -= dx;
  }

  xM += 2*(mpv - xm);
  while((fx = f->Eval(xM))>.5*max) xM += dx;
}


//________________________________________________________
void AliTRDresolution::SetRecoParam(AliTRDrecoParam *r)
{

  fReconstructor->SetRecoParam(r);
}


//________________________________________________________
void AliTRDresolution::SetSegmentationLevel(Int_t l) 
{
// Setting the segmentation level to "l"
  fSegmentLevel = l;

  UShort_t const lNcomp[kNprojs] = {
    1,  1, //2, 
    fgkNresYsegm[fSegmentLevel], 2, //2, 
    2*fgkNresYsegm[fSegmentLevel], 2, 2, 2, 1, //5, 
    2*fgkNresYsegm[fSegmentLevel], 2, 2, 2, 1, //5,
    2*fgkNresYsegm[fSegmentLevel], 2, 2, 2, 1, //5,
  // MC
    fgkNresYsegm[fSegmentLevel], 2,          //2, 
    fgkNresYsegm[fSegmentLevel], 2, 2, 2, 1, //5, 
    fgkNresYsegm[fSegmentLevel], 2, 2, 2, 1, 1, 1, 1, 11, 11, 11, //11
    fgkNresYsegm[fSegmentLevel], 2, 2, 2, 1, 1, 1, 1, 11, 11, 11, //11
    6*fgkNresYsegm[fSegmentLevel], 6*2, 6*2, 6*2, 6, 6, 6, 6, 6*11, 6*11, 6*11  //11
  };
  memcpy(fNcomp, lNcomp, kNprojs*sizeof(UShort_t));

  Char_t const *lAxTitle[kNprojs][4] = {
    // Charge
    {"Impv", "x [cm]", "I_{mpv}", "x/x_{0}"}
  ,{"dI/Impv", "x/x_{0}", "#delta I/I_{mpv}", "x[cm]"}
    // Clusters to Kalman
  ,{"Cluster2Track residuals", "tg(#phi)", "y [#mum]", "#sigma_{y} [#mum]"}
  ,{"Cluster2Track  YZ pulls", fgkResYsegmName[fSegmentLevel], "y / z", "#sigma_{y}"}
    // TRD tracklet to Kalman fit
  ,{"Tracklet2Track Y residuals", "tg(#phi)", "y [#mum]", "#sigma_{y} [#mum]"}
  ,{"Tracklet2Track YZ pulls", fgkResYsegmName[fSegmentLevel], "y / z", "#sigma_{y}"}
  ,{"Tracklet2Track Z residuals", "tg(#theta)", "z [#mum]", "#sigma_{z} [#mum]"}
  ,{"Tracklet2Track Z pulls", "tg(#theta)", "z", "#sigma_{z}"}
  ,{"Tracklet2Track Phi residuals", "tg(#phi)", "#phi [mrad]", "#sigma_{#phi} [mrad]"}
    // TRDin 2 first TRD tracklet
  ,{"Tracklet2Track Y residuals @ TRDin", "tg(#phi)", "y [#mum]", "#sigma_{y} [#mum]"}
  ,{"Tracklet2Track YZ pulls @ TRDin", fgkResYsegmName[fSegmentLevel], "y / z", "#sigma_{y}"}
  ,{"Tracklet2Track Z residuals @ TRDin", "tg(#theta)", "z [#mum]", "#sigma_{z} [#mum]"}
  ,{"Tracklet2Track Z pulls @ TRDin", "tg(#theta)", "z", "#sigma_{z}"}
  ,{"Tracklet2Track Phi residuals @ TRDin", "tg(#phi)", "#phi [mrad]", "#sigma_{#phi} [mrad]"}
    // TRDout 2 first TRD tracklet
  ,{"Tracklet2Track Y residuals @ TRDout", "tg(#phi)", "y [#mum]", "#sigma_{y} [#mum]"}
  ,{"Tracklet2Track YZ pulls @ TRDout", fgkResYsegmName[fSegmentLevel], "y / z", "#sigma_{y}"}
  ,{"Tracklet2Track Z residuals @ TRDout", "tg(#theta)", "z [#mum]", "#sigma_{z} [#mum]"}
  ,{"Tracklet2Track Z pulls @ TRDout", "tg(#theta)", "z", "#sigma_{z}"}
  ,{"Tracklet2Track Phi residuals @ TRDout", "tg(#phi)", "#phi [mrad]", "#sigma_{#phi} [mrad]"}
    // MC cluster
  ,{"MC Cluster Y resolution", "tg(#phi)", "y [#mum]", "#sigma_{y} [#mum]"}
  ,{"MC Cluster YZ pulls", fgkResYsegmName[fSegmentLevel], "y / z", "#sigma_{y}"}
    // MC tracklet
  ,{"MC Tracklet Y resolution", "tg(#phi)", "y [#mum]",  "#sigma_{y}[#mum]"}
  ,{"MC Tracklet YZ pulls", fgkResYsegmName[fSegmentLevel], "y / z", "#sigma_{y}"}
  ,{"MC Tracklet Z resolution", "tg(#theta)", "z [#mum]", "#sigma_{z} [#mum]"}
  ,{"MC Tracklet Z pulls", "tg(#theta)", "z", "#sigma_{z}"}
  ,{"MC Tracklet Phi resolution", "tg(#phi)", "#phi [mrad]", "#sigma_{#phi} [mrad]"}
    // MC track TRDin
  ,{"MC Y resolution @ TRDin", "tg(#phi)", "y [#mum]", "#sigma_{y}[#mum]"}
  ,{"MC YZ pulls @ TRDin", fgkResYsegmName[fSegmentLevel], "y / z", "#sigma_{y}"}
  ,{"MC Z resolution @ TRDin", "tg(#theta)", "z [#mum]", "#sigma_{z} [#mum]"}
  ,{"MC Z pulls @ TRDin", "tg(#theta)", "z", "#sigma_{z}"}
  ,{"MC #Phi resolution @ TRDin", "tg(#phi)", "#phi [mrad]", "#sigma_{#phi} [mrad]"}
  ,{"MC SNP pulls @ TRDin", "tg(#phi)", "SNP", "#sigma_{snp}"}
  ,{"MC #Theta resolution @ TRDin", "tg(#theta)", "#theta [mrad]", "#sigma_{#theta} [mrad]"}
  ,{"MC TGL pulls @ TRDin", "tg(#theta)", "TGL", "#sigma_{tgl}"}
  ,{"MC P_{t} resolution @ TRDin", "p_{t}^{MC} [GeV/c]", "(p_{t}^{REC}-p_{t}^{MC})/p_{t}^{MC} [%]", "MC: #sigma^{TPC}(#Deltap_{t}/p_{t}^{MC}) [%]"}
  ,{"MC 1/P_{t} pulls @ TRDin", "1/p_{t}^{MC} [c/GeV]", "1/p_{t}^{REC}-1/p_{t}^{MC}", "MC PULL: #sigma_{1/p_{t}}^{TPC}"}
  ,{"MC P resolution @ TRDin", "p^{MC} [GeV/c]", "(p^{REC}-p^{MC})/p^{MC} [%]", "MC: #sigma^{TPC}(#Deltap/p^{MC}) [%]"}
    // MC track TRDout
  ,{"MC Y resolution @ TRDout", "tg(#phi)", "y [#mum]", "#sigma_{y}[#mum]"}
  ,{"MC YZ pulls @ TRDout", fgkResYsegmName[fSegmentLevel], "y / z", "#sigma_{y}"}
  ,{"MC Z resolution @ TRDout", "tg(#theta)", "z [#mum]", "#sigma_{z} [#mum]"}
  ,{"MC Z pulls @ TRDout", "tg(#theta)", "z", "#sigma_{z}"}
  ,{"MC #Phi resolution @ TRDout", "tg(#phi)", "#phi [mrad]", "#sigma_{#phi} [mrad]"}
  ,{"MC SNP pulls @ TRDout", "tg(#phi)", "SNP", "#sigma_{snp}"}
  ,{"MC #Theta resolution @ TRDout", "tg(#theta)", "#theta [mrad]", "#sigma_{#theta} [mrad]"}
  ,{"MC TGL pulls @ TRDout", "tg(#theta)", "TGL", "#sigma_{tgl}"}
  ,{"MC P_{t} resolution @ TRDout", "p_{t}^{MC} [GeV/c]", "(p_{t}^{REC}-p_{t}^{MC})/p_{t}^{MC} [%]", "MC: #sigma^{TPC}(#Deltap_{t}/p_{t}^{MC}) [%]"}
  ,{"MC 1/P_{t} pulls @ TRDout", "1/p_{t}^{MC} [c/GeV]", "1/p_{t}^{REC}-1/p_{t}^{MC}", "MC PULL: #sigma_{1/p_{t}}^{TPC}"}
  ,{"MC P resolution @ TRDout", "p^{MC} [GeV/c]", "(p^{REC}-p^{MC})/p^{MC} [%]", "MC: #sigma^{TPC}(#Deltap/p^{MC}) [%]"}
    // MC track in TRD
  ,{"MC Track Y resolution", "tg(#phi)", "y [#mum]", "#sigma_{y} [#mum]"}
  ,{"MC Track YZ pulls", fgkResYsegmName[fSegmentLevel], "y / z", "#sigma_{y}"}
  ,{"MC Track Z resolution", "tg(#theta)", "z [#mum]", "#sigma_{z} [#mum]"}
  ,{"MC Track Z pulls", "tg(#theta)", "z", "#sigma_{z}"}
  ,{"MC Track #Phi resolution", "tg(#phi)", "#phi [mrad]", "#sigma_{#phi} [mrad]"}
  ,{"MC Track SNP pulls", "tg(#phi)", "SNP", "#sigma_{snp}"}
  ,{"MC Track #Theta resolution", "tg(#theta)", "#theta [mrad]", "#sigma_{#theta} [mrad]"}
  ,{"MC Track TGL pulls", "tg(#theta)", "TGL", "#sigma_{tgl}"}
  ,{"MC P_{t} resolution", "p_{t} [GeV/c]", "(p_{t}^{REC}-p_{t}^{MC})/p_{t}^{MC} [%]", "#sigma(#Deltap_{t}/p_{t}^{MC}) [%]"}
  ,{"MC 1/P_{t} pulls", "1/p_{t}^{MC} [c/GeV]", "1/p_{t}^{REC} - 1/p_{t}^{MC}", "#sigma_{1/p_{t}}"}
  ,{"MC P resolution", "p [GeV/c]", "(p^{REC}-p^{MC})/p^{MC} [%]", "#sigma(#Deltap/p^{MC}) [%]"}
  };
  memcpy(fAxTitle, lAxTitle, 4*kNprojs*sizeof(Char_t*));
}
