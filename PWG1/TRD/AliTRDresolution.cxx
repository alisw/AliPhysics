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

UChar_t const AliTRDresolution::fgNhistos[kNviews] = {
  2, 2, 5, 5,
  2, 5, 12, 2, 14
};
UChar_t const AliTRDresolution::fgNproj[kNviews] = {
  2, 2, 5, 5,
  4, 7, 12, 2, 14
};
Char_t const * AliTRDresolution::fgPerformanceName[kNviews] = {
     "Charge"
    ,"Cluster2Track"
    ,"Tracklet2Track"
    ,"Tracklet2TPC" 
    ,"Cluster2MC"
    ,"Tracklet2MC"
    ,"TPC2MC"
    ,"TOF/HMPID2MC"
    ,"TRD2MC"
};
UChar_t const AliTRDresolution::fgNcomp[kNprojs] = {
  1, 1, //2, 
  1, 1, //2, 
  1, 1, 1, 1, 1, //5, 
  1, 1, 1, 1, 1, 1, 1, //5,
// MC
  1, 1, 1, 1,    //4, 
  1, 1, 1, 1, 1, //5, 
  1, 1, 1, 1, 1, 1, 1, 1, 11, 11, 11, 11, //12, 
  1, 1, //2, 
  1, 1, 1, 1, 1, 1, 1, 1, 66, 66, 66, 66, 66, 66 //14
};
Char_t const *AliTRDresolution::fgAxTitle[kNprojs][4] = {
  // Charge
  {"Impv", "x [cm]", "I_{mpv}", "x/x_{0}"}
 ,{"dI/Impv", "x/x_{0}", "#delta I/I_{mpv}", "x[cm]"}
  // Clusters to Kalman
 ,{"Cluster2Track residuals", "tg(#phi)", "#mu_{y} [#mum]", "#sigma_{y} [#mum]"}
 ,{"Cluster2Track  pulls", "tg(#phi)", "#mu_{y}", "#sigma_{y}"}
  // TRD tracklet to Kalman fit
 ,{"Tracklet2Track Y residuals", "tg(#phi)", "#mu_{y} [#mum]", "#sigma_{y} [#mum]"}
 ,{"Tracklet2Track Y pulls", "tg(#phi)", "#mu_{y}", "#sigma_{y}"}
 ,{"Tracklet2Track Z residuals", "tg(#theta)", "#mu_{z} [#mum]", "#sigma_{z} [#mum]"}
 ,{"Tracklet2Track Z pulls", "tg(#theta)", "#mu_{z}", "#sigma_{z}"}
 ,{"Tracklet2Track Phi residuals", "tg(#phi)", "#mu_{#phi} [mrad]", "#sigma_{#phi} [mrad]"}
  // TPC track 2 first TRD tracklet
 ,{"Tracklet2Track Y residuals @ TRDin", "tg(#phi)", "#mu_{y} [#mum]", "#sigma_{y} [#mum]"}
 ,{"Tracklet2Track Y pulls @ TRDin", "tg(#phi)", "#mu_{y}", "#sigma_{y}"}
 ,{"Tracklet2Track Z residuals @ TRDin", "tg(#theta)", "#mu_{z} [#mum]", "#sigma_{z} [#mum]"}
 ,{"Tracklet2Track Z pulls @ TRDin", "tg(#theta)", "#mu_{z}", "#sigma_{z}"}
 ,{"Tracklet2Track Phi residuals @ TRDin", "tg(#phi)", "#mu_{#phi} [mrad]", "#sigma_{#phi} [mrad]"}
  // MC cluster
 ,{"MC Cluster Y resolution (p_{t}<1 GeV/c)", "tg(#phi)", "#mu_{y} [#mum]", "#sigma_{y} [#mum]"}
 ,{"MC Cluster Y resolution (1<p_{t}<2 GeV/c)", "tg(#phi)", "#mu_{y} [#mum]", "#sigma_{y} [#mum]"}
 ,{"MC Cluster Y resolution (p_{t}>3 GeV/c)", "tg(#phi)", "#mu_{y} [#mum]", "#sigma_{y} [#mum]"}
 ,{"MC Cluster Y pulls", "tg(#phi)", "#mu_{y}", "#sigma_{y}"}
  // MC tracklet
 ,{"MC Tracklet Y resolution (p_{t}<1 GeV/c)", "tg(#phi)", "#mu_{y} [#mum]", "#sigma_{y} [#mum]"}
 ,{"MC Tracklet Y resolution (1<p_{t}<2 GeV/c)", "tg(#phi)", "#mu_{y} [#mum]", "#sigma_{y} [#mum]"}
 ,{"MC Tracklet Y resolution (p_{t}>3 GeV/c)", "tg(#phi)", "#mu_{y} [#mum]", "#sigma_{y}[#mum]"}
 ,{"MC Tracklet Y pulls", "tg(#phi)", "#mu_{y}", "#sigma_{y}"}
 ,{"MC Tracklet Cross Z resolution", "tg(#theta)", "#mu_{z} [#mum]", "#sigma_{z} [#mum]"}
 ,{"MC Tracklet Cross Z pulls", "tg(#theta)", "#mu_{z}", "#sigma_{z}"}
 ,{"MC Tracklet Phi resolution", "tg(#phi)", "#mu_{#phi} [mrad]", "#sigma_{#phi} [mrad]"}
  // MC track TPC
 ,{"Y resolution @ TRDin", "tg(#phi)", "#mu_{y} [#mum]", "#sigma_{y}[#mum]"}
 ,{"Y pulls @ TRDin", "tg(#phi)", "#mu_{y}", "#sigma_{y}"}
 ,{"Z resolution @ TRDin", "tg(#theta)", "#mu_{z} [#mum]", "#sigma_{z} [#mum]"}
 ,{"Z pulls @ TRDin", "tg(#theta)", "#mu_{z}", "#sigma_{z}"}
 ,{"Phi resolution @ TRDin", "tg(#phi)", "#mu_{#phi} [mrad]", "#sigma_{#phi} [mrad]"}
 ,{"SNP pulls @ TRDin", "tg(#phi)", "#mu_{snp}", "#sigma_{snp}"}
 ,{"Theta resolution @ TRDin", "tg(#theta)", "#mu_{#theta} [mrad]", "#sigma_{#theta} [mrad]"}
 ,{"TGL pulls @ TRDin", "tg(#theta)", "#mu_{tgl}", "#sigma_{tgl}"}
 ,{"P_{t} resolution @ TRDin", "p_{t}^{MC} [GeV/c]", "(p_{t}^{REC}-p_{t}^{MC})/p_{t}^{MC} [%]", "MC: #sigma^{TPC}(#Deltap_{t}/p_{t}^{MC}) [%]"}
 ,{"1/P_{t} pulls @ TRDin", "1/p_{t}^{MC} [c/GeV]", "1/p_{t}^{REC}-1/p_{t}^{MC}", "MC PULL: #sigma_{1/p_{t}}^{TPC}"}
 ,{"P resolution @ TRDin", "p^{MC} [GeV/c]", "(p^{REC}-p^{MC})/p^{MC} [%]", "MC: #sigma^{TPC}(#Deltap/p^{MC}) [%]"}
 ,{"P pulls @ TRDin", "p^{MC} [GeV/c]", "1/p^{REC}-1/p^{MC}", "MC PULL: #sigma^{TPC}(#Deltap/#sigma_{p})"}
  // MC track TOF
 ,{"PosZ", "tg(#theta)", "MC: #mu_{z}^{TOF} [#mum]", "MC: #sigma_{z}^{TOF} [#mum]"}
 ,{"PullsZ", "tg(#theta)", "MC PULL: #mu_{z}^{TOF}", "MC PULL: #sigma_{z}^{TOF}"}
  // MC track in TRD
 ,{"TRD track MC Y resolution", "tg(#phi)", "#mu_{y}^{Trk} [#mum]", "#sigma_{y}^{Trk} [#mum]"}
 ,{"TRD track MC Y pulls", "tg(#phi)", "#mu_{y}^{Trk}", "#sigma_{y}^{Trk}"}
 ,{"TRD track MC Z resolution", "tg(#theta)", "#mu_{z}^{Trk} [#mum]", "#sigma_{z}^{Trk} [#mum]"}
 ,{"TRD track MC Z pulls", "tg(#theta)", "#mu_{z}^{Trk}", "#sigma_{z}^{Trk}"}
 ,{"TRD track MC Phi resolution", "tg(#phi)", "#mu_{#phi}^{Trk} [mrad]", "#sigma_{#phi}^{Trk} [mrad]"}
 ,{"TRD track MC SNP pulls", "tg(#phi)", "#mu_{snp}^{Trk}", "#sigma_{snp}^{Trk}"}
 ,{"TRD track MC Theta resolution", "tg(#theta)", "#mu_{#theta}^{Trk} [mrad]", "#sigma_{#theta}^{Trk} [mrad]"}
 ,{"TRD track MC TGL pulls", "tg(#theta)", "#mu_{tgl}^{Trk}", "#sigma_{tgl}^{Trk}"}
 ,{"P_{t} resolution TRD Layer", "p_{t} [GeV/c]", "(p_{t}^{REC}-p_{t}^{MC})/p_{t}^{MC} [%]", "#sigma(#Deltap_{t}/p_{t}^{MC}) [%]"}
 ,{"1/P_{t} pulls TRD Layer", "1/p_{t}^{MC} [c/GeV]", "1/p_{t}^{REC} - 1/p_{t}^{MC}", "#sigma_{1/p_{t}}"}
 ,{"P resolution TRD Layer", "p [GeV/c]", "(p^{REC}-p^{MC})/p^{MC} [%]", "#sigma(#Deltap/p^{MC}) [%]"}
 ,{"[SA] P_{t} resolution TRD Layer", "p_{t}^{MC} [GeV/c]", "(p_{t}^{REC}-p_{t}^{MC})/p_{t}^{MC} [%]", "MC: #sigma^{Trk}(#Deltap_{t}/p_{t}^{MC}) [%]"}
 ,{"[SA] 1/P_{t} pulls TRD Layer", "1/p_{t}^{MC} [c/GeV]", "1/p_{t}^{REC}-1/p_{t}^{MC}", "MC PULL: #sigma_{1/p_{t}}^{Trk}"}
 ,{"[SA] P resolution TRD Layer", "p^{MC} [GeV/c]", "(p^{REC}-p^{MC})/p^{MC} [%]", "MC: #sigma^{Trk}(#Deltap/p^{MC}) [%]"}
};

//________________________________________________________
AliTRDresolution::AliTRDresolution()
  :AliTRDrecoTask()
  ,fStatus(0)
  ,fIdxPlot(0)
  ,fIdxFrame(0)
  ,fReconstructor(NULL)
  ,fGeo(NULL)
  ,fDBPDG(NULL)
  ,fGraphS(NULL)
  ,fGraphM(NULL)
  ,fCl(NULL)
  ,fTrklt(NULL)
  ,fMCcl(NULL)
  ,fMCtrklt(NULL)
{
  //
  // Default constructor
  //
  SetNameTitle("TRDresolution", "TRD spatial and momentum resolution");
}

//________________________________________________________
AliTRDresolution::AliTRDresolution(char* name)
  :AliTRDrecoTask(name, "TRD spatial and momentum resolution")
  ,fStatus(0)
  ,fIdxPlot(0)
  ,fIdxFrame(0)
  ,fReconstructor(NULL)
  ,fGeo(NULL)
  ,fDBPDG(NULL)
  ,fGraphS(NULL)
  ,fGraphM(NULL)
  ,fCl(NULL)
  ,fTrklt(NULL)
  ,fMCcl(NULL)
  ,fMCtrklt(NULL)
{
  //
  // Default constructor
  //

  fReconstructor = new AliTRDReconstructor();
  fReconstructor->SetRecoParam(AliTRDrecoParam::GetLowFluxParam());
  fGeo = new AliTRDgeometry();

  InitFunctorList();

  DefineOutput(kClToTrk, TObjArray::Class()); // cluster2track
  DefineOutput(kTrkltToTrk, TObjArray::Class()); // tracklet2track
  DefineOutput(kClToMC, TObjArray::Class()); // cluster2mc
  DefineOutput(kTrkltToMC, TObjArray::Class()); // tracklet2mc
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
  if(fTrklt){fTrklt->Delete(); delete fTrklt;}
  if(fMCcl){fMCcl->Delete(); delete fMCcl;}
  if(fMCtrklt){fMCtrklt->Delete(); delete fMCtrklt;}
}


//________________________________________________________
void AliTRDresolution::UserCreateOutputObjects()
{
  // spatial resolution
  OpenFile(1, "RECREATE");
  fContainer = Histos();

  fCl = new TObjArray();
  fCl->SetOwner(kTRUE);
  fTrklt = new TObjArray();
  fTrklt->SetOwner(kTRUE);
  fMCcl = new TObjArray();
  fMCcl->SetOwner(kTRUE);
  fMCtrklt = new TObjArray();
  fMCtrklt->SetOwner(kTRUE);
}

//________________________________________________________
void AliTRDresolution::UserExec(Option_t *opt)
{
  //
  // Execution part
  //

  fCl->Delete();
  fTrklt->Delete();
  fMCcl->Delete();
  fMCtrklt->Delete();
  AliTRDrecoTask::UserExec(opt);
  PostData(kClToTrk, fCl);
  PostData(kTrkltToTrk, fTrklt);
  PostData(kClToMC, fMCcl);
  PostData(kTrkltToMC, fMCtrklt);
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
  if(!(arr = ((TObjArray*)fContainer->At(kCharge)))){
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
  if(!(arr = ((TObjArray*)fContainer->At(kCluster)))){
    AliWarning("No output container defined.");
    return NULL;
  }
  ULong_t status = fkESD ? fkESD->GetStatus():0;

  Double_t covR[7], cov[3];
  Float_t pt, x0, y0, z0, dy, dz, dydx, dzdx;
  AliTRDseedV1 *fTracklet(NULL);  
  for(Int_t ily=0; ily<AliTRDgeometry::kNlayer; ily++){
    if(!(fTracklet = fkTrack->GetTracklet(ily))) continue;
    if(!fTracklet->IsOK()) continue;
    x0 = fTracklet->GetX0();
    pt = fTracklet->GetPt();

    // retrive the track angle with the chamber
    y0   = fTracklet->GetYref(0);
    z0   = fTracklet->GetZref(0);
    dydx = fTracklet->GetYref(1);
    dzdx = fTracklet->GetZref(1);
    fTracklet->GetCovRef(covR);
    Float_t tilt(fTracklet->GetTilt());
    AliTRDcluster *c = NULL;
    fTracklet->ResetClusterIter(kFALSE);
    while((c = fTracklet->PrevCluster())){
      Float_t xc = c->GetX();
      Float_t yc = c->GetY();
      Float_t zc = c->GetZ();
      Float_t dx = x0 - xc; 
      Float_t yt = y0 - dx*dydx;
      Float_t zt = z0 - dx*dzdx; 
      // calculate residuals using tilt correction
//       yc -= tilt*(zc-zt); 
//       dy = yt - yc;

      // calculate residuals using correct covariance matrix
      cov[0] = c->GetSigmaY2();
      cov[1] = c->GetSigmaYZ();
      cov[2] = c->GetSigmaZ2();
      // do rotation
      Double_t sy2(cov[0]), sz2(cov[2]);
      Double_t t2 = tilt*tilt;
      Double_t correction = 1./(1. + t2);
      cov[0] = (sy2+t2*sz2)*correction;
      cov[1] = tilt*(sz2 - sy2)*correction;
      cov[2] = (t2*sy2+sz2)*correction;
      // do inversion
      Double_t r00=cov[0]+covR[0], r01=cov[1]+covR[1], r11=cov[2]+covR[2];
      Double_t det=r00*r11 - r01*r01;
      Double_t tmp=r00; r00=r11/det; r11=tmp/det;
      dy = (yc - yt)*TMath::Sqrt(r00);
      dz = (zc - zt)*TMath::Sqrt(r11);

      ((TH3S*)arr->At(0))->Fill(dydx, dy/*, dz*/, pt);
      ((TH2I*)arr->At(1))->Fill(dydx, dy/TMath::Sqrt(cov[0] /*+ sx2*/ + sy2));
  
      if(DebugLevel()>=2){
        // Get z-position with respect to anode wire
        Int_t istk = fGeo->GetStack(c->GetDetector());
        AliTRDpadPlane *pp = fGeo->GetPadPlane(ily, istk);
        Float_t row0 = pp->GetRow0();
        Float_t d  =  row0 - zt + pp->GetAnodeWireOffset();
        d -= ((Int_t)(2 * d)) / 2.0;
        if (d > 0.25) d  = 0.5 - d;

        AliTRDclusterInfo *clInfo = new AliTRDclusterInfo;
        fCl->Add(clInfo);
        clInfo->SetCluster(c);
        Float_t covF[] = {cov[0], cov[1], cov[2]};
        clInfo->SetGlobalPosition(yt, zt, dydx, dzdx, covF);
        clInfo->SetResolution(dy);
        clInfo->SetAnisochronity(d);
        clInfo->SetDriftLength(dx);
        clInfo->SetTilt(tilt);
        (*DebugStream()) << "ClusterREC"
          <<"status="  << status
          <<"clInfo.=" << clInfo
          << "\n";
      }
    }
  }
  return (TH2I*)arr->At(0);
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
  if(!(arr = (TObjArray*)fContainer->At(kTrackTRD ))){
    AliWarning("No output container defined.");
    return NULL;
  }

  Double_t cov[3], covR[7]/*, sqr[3], inv[3]*/;
  Float_t x, dx, dy, dz;
  AliTRDseedV1 *fTracklet = NULL;  
  for(Int_t il=AliTRDgeometry::kNlayer; il--;){
    if(!(fTracklet = fkTrack->GetTracklet(il))) continue;
    if(!fTracklet->IsOK()) continue;
    x    = fTracklet->GetX();
    dx   = fTracklet->GetX0() - x;
    // compute dy^2 and dz^2
    dy   = fTracklet->GetYref(0)-dx*fTracklet->GetYref(1) - fTracklet->GetY();
    dz   = fTracklet->GetZref(0)-dx*fTracklet->GetZref(1) - fTracklet->GetZ();
    // compute covariance matrix
    fTracklet->GetCovAt(x, cov);
    fTracklet->GetCovRef(covR);
    cov[0] += covR[0]; cov[1] += covR[1]; cov[2] += covR[2]; 
/*  // Correct PULL calculation by considering off  
    // diagonal elements in the covariance matrix
    // compute square root matrix
    if(AliTRDseedV1::GetCovInv(cov, inv)==0.) continue;
    if(AliTRDseedV1::GetCovSqrt(inv, sqr)<0.) continue;
    Double_t y = sqr[0]*dy+sqr[1]*dz;
    Double_t z = sqr[1]*dy+sqr[2]*dz;
    ((TH3*)h)->Fill(y, z, fTracklet->GetYref(1));*/

    ((TH2I*)arr->At(0))->Fill(fTracklet->GetYref(1), dy);
    ((TH2I*)arr->At(1))->Fill(fTracklet->GetYref(1), dy/TMath::Sqrt(cov[0]));
    ((TH2I*)arr->At(4))->Fill(fTracklet->GetYref(1), TMath::ATan((fTracklet->GetYref(1)-fTracklet->GetYfit(1))/(1-fTracklet->GetYref(1)*fTracklet->GetYfit(1))));
    if(!fTracklet->IsRowCross()) continue;
    ((TH2I*)arr->At(2))->Fill(fTracklet->GetZref(1), dz);
    ((TH2I*)arr->At(3))->Fill(fTracklet->GetZref(1), dz/TMath::Sqrt(cov[2]));
  }


  return (TH2I*)arr->At(0);
}


//________________________________________________________
TH1* AliTRDresolution::PlotTrackTPC(const AliTRDtrackV1 *track)
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
  AliExternalTrackParam *tin = NULL;
  if(!(tin = fkTrack->GetTrackLow())){
    AliWarning("Track did not entered TRD fiducial volume.");
    return NULL;
  }
  TH1 *h = NULL;
  
  Double_t x = tin->GetX();
  AliTRDseedV1 *tracklet = NULL;  
  for(Int_t ily=0; ily<AliTRDgeometry::kNlayer; ily++){
    if(!(tracklet = fkTrack->GetTracklet(ily))) continue;
    break;
  }
  if(!tracklet || TMath::Abs(x-tracklet->GetX())>1.e-3){
    AliWarning("Tracklet did not match TRD entrance.");
    return NULL;
  }
  const Int_t kNPAR(5);
  Double_t parR[kNPAR]; memcpy(parR, tin->GetParameter(), kNPAR*sizeof(Double_t));
  Double_t covR[3*kNPAR]; memcpy(covR, tin->GetCovariance(), 3*kNPAR*sizeof(Double_t));
  Double_t cov[3]; tracklet->GetCovAt(x, cov);

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
  Double_t dy = parR[0] - tracklet->GetY(); 
  TObjArray *arr = (TObjArray*)fContainer->At(kTrackTPC);
  ((TH3S*)arr->At(0))->Fill(tracklet->GetYref(1), dy, 1./PAR[4]);
  ((TH2I*)arr->At(1))->Fill(tracklet->GetYref(1), dy/TMath::Sqrt(COV(0,0)+cov[0]));
  if(tracklet->IsRowCross()){
    Double_t dz = parR[1] - tracklet->GetZ(); 
    ((TH2I*)arr->At(2))->Fill(tracklet->GetZref(1), dz);
    ((TH2I*)arr->At(3))->Fill(tracklet->GetZref(1), dz/TMath::Sqrt(COV(1,1)+cov[2]));
  }
  Double_t dphi = TMath::ASin(PAR[2])-TMath::ATan(tracklet->GetYfit(1));  ((TH2I*)arr->At(4))->Fill(tracklet->GetYref(1), dphi);


  // register reference histo for mini-task
  h = (TH2I*)arr->At(0);

  if(DebugLevel()>=1){
    (*DebugStream()) << "trackIn"
      << "x="       << x
      << "P="       << &PAR
      << "C="       << &COV
      << "\n";

    Double_t y = tracklet->GetY(); 
    Double_t z = tracklet->GetZ(); 
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
  Float_t dx, pt0, x0=tracklet->GetX0(), y0, z0, dydx0, dzdx0;
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
  arr = (TObjArray*)fContainer->At(kMCtrackTPC);
  // y resolution/pulls
  ((TH3S*)arr->At(0))->Fill(dydx0, PARMC[0]-PAR[0], pt0);
  ((TH2I*)arr->At(1))->Fill(dydx0, (PARMC[0]-PAR[0])/TMath::Sqrt(COV(0,0)));
  // z resolution/pulls
  ((TH2I*)arr->At(2))->Fill(dzdx0, PARMC[1]-PAR[1]);
  ((TH2I*)arr->At(3))->Fill(dzdx0, (PARMC[1]-PAR[1])/TMath::Sqrt(COV(1,1)));
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
  if(DebugLevel()>=1){
    (*DebugStream()) << "trackInMC"
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
    AliWarning("No MC defined. Results will not be available.");
    return NULL;
  }
  if(track) fkTrack = track;
  if(!fkTrack){
    AliDebug(4, "No Track defined.");
    return NULL;
  }
  // retriev track characteristics
  Int_t pdg = fkMC->GetPDG(),
        sIdx(AliTRDpidUtil::Pdg2Pid(TMath::Abs(pdg))+1), // species index
        sign(0),
        det(-1),
        label(fkMC->GetLabel());
  if(!fDBPDG) fDBPDG=TDatabasePDG::Instance();
  TParticlePDG *ppdg(fDBPDG->GetParticle(pdg));
  if(ppdg) sign = ppdg->Charge() > 0. ? 1 : -1;

  TObjArray *arr(NULL);TH1 *h(NULL);
  UChar_t s;
  Double_t xAnode, x, y, z, pt, dydx, dzdx, dzdl;
  Float_t pt0, x0, y0, z0, dx, dy, dz, dydx0, dzdx0;
  Double_t covR[7]/*, cov[3]*/;

  if(DebugLevel()>=1){
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

    det = fTracklet->GetDetector();
    x0  = fTracklet->GetX0();
    //radial shift with respect to the MC reference (radial position of the pad plane)
    x= fTracklet->GetX();
    if(!fkMC->GetDirections(x0, y0, z0, dydx0, dzdx0, pt0, s)) continue;
    xAnode  = fTracklet->GetX0();

    // MC track position at reference radial position
    dx  = x0 - x;
    if(DebugLevel()>=1){
      (*DebugStream()) << "MC"
        << "det="     << det
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

    arr = (TObjArray*)((TObjArray*)fContainer->At(kMCtrackTRD))->At(ily);
    // y resolution/pulls
    ((TH3S*)arr->At(0))->Fill(dydx0, dy, pt0);
    ((TH2I*)arr->At(1))->Fill(dydx0, dy/TMath::Sqrt(covR[0]));
    // z resolution/pulls
    ((TH2I*)arr->At(2))->Fill(dzdx0, dz);
    ((TH2I*)arr->At(3))->Fill(dzdx0, dz/TMath::Sqrt(covR[2]));
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
    if(DebugLevel()>=1){
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
    tt.Fit(kTRUE, kTRUE);
    x= tt.GetX();y= tt.GetY();z= tt.GetZ();
    dydx = tt.GetYfit(1);
    dx = x0 - x;
    ymc = y0 - dx*dydx0;
    zmc = z0 - dx*dzdx0;
    Bool_t rc = tt.IsRowCross(); 
    
    // add tracklet residuals for y and dydx
    arr = (TObjArray*)fContainer->At(kMCtracklet);
    if(!rc){
      dy    = y-ymc;

      Float_t dphi  = (dydx - dydx0);
      dphi /= (1.- dydx*dydx0);

      ((TH3S*)arr->At(0))->Fill(dydx0, dy, pt0);
      if(tt.GetS2Y()>0.) ((TH2I*)arr->At(1))->Fill(dydx0, dy/TMath::Sqrt(tt.GetS2Y()));
      ((TH2I*)arr->At(4))->Fill(dydx0, TMath::ATan(dphi));
    } else {
      // add tracklet residuals for z
      dz = z-zmc;
      ((TH2I*)arr->At(2))->Fill(dzdl0, dz);
      if(tt.GetS2Z()>0.) ((TH2I*)arr->At(3))->Fill(dzdl0, dz/TMath::Sqrt(tt.GetS2Z()));
    }
  
    // Fill Debug stream for tracklet
    if(DebugLevel()>=1){
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

    Int_t istk = AliTRDgeometry::GetStack(det); 
    AliTRDpadPlane *pp = fGeo->GetPadPlane(ily, istk);
    Float_t zr0 = pp->GetRow0() + pp->GetAnodeWireOffset();
    Float_t tilt = fTracklet->GetTilt();
    //Double_t exb = AliTRDCommonParam::Instance()->GetOmegaTau(1.5);

    arr = (TObjArray*)fContainer->At(kMCcluster);
    AliTRDcluster *c = NULL;
    fTracklet->ResetClusterIter(kFALSE);
    while((c = fTracklet->PrevCluster())){
      Float_t  q = TMath::Abs(c->GetQ());
      x = c->GetX(); y = c->GetY();z = c->GetZ();
      dx = x0 - x; 
      ymc= y0 - dx*dydx0;
      zmc= z0 - dx*dzdx0;
      dy = (y - tilt*(z-zmc)) - ymc;
      // Fill Histograms
      if(q>20. && q<250.){ 
        ((TH3S*)arr->At(0))->Fill(dydx0, dy, pt0);
        ((TH2I*)arr->At(1))->Fill(dydx0, dy/TMath::Sqrt(c->GetSigmaY2()));
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
      clInfo->SetDriftLength(dx-.5*AliTRDgeometry::CamHght());
      clInfo->SetTilt(tilt);
      fMCcl->Add(clInfo);
      if(DebugLevel()>=2){ 
        if(!clInfoArr) clInfoArr=new TObjArray(AliTRDseedV1::kNclusters);
        clInfoArr->Add(clInfo);
      }
    }
    // Fill Debug Tree
    if(DebugLevel()>=2 && clInfoArr){
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
  Int_t first(2),     // first particle species to be drawn
        nspecies(7); // last particle species to be drawn
  TList *l = NULL; TVirtualPad *pad=NULL;
  switch(ifig){
  case kCharge:
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    ((TVirtualPad*)l->At(0))->cd();
    ((TGraphAsymmErrors*)((TObjArray*)fGraphM->At(kCharge))->At(0))->Draw("apl");
    ((TVirtualPad*)l->At(1))->cd();
    ((TGraphErrors*)((TObjArray*)fGraphS->At(kCharge))->At(0))->Draw("apl");
    break;
  case kCluster:
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    xy[0] = -.3; xy[1] = -200.; xy[2] = .3; xy[3] = 1000.;
    ((TVirtualPad*)l->At(0))->cd();
    if(!GetGraphPlot(&xy[0], kCluster, 0)) break;
    xy[0] = -.3; xy[1] = -0.5; xy[2] = .3; xy[3] = 2.5;
    ((TVirtualPad*)l->At(1))->cd();
    if(!GetGraphPlot(&xy[0], kCluster, 1)) break;
    return kTRUE;
  case kTrackTRD :
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    xy[0] = -.3; xy[1] = -500.; xy[2] = .3; xy[3] = 1500.;
    ((TVirtualPad*)l->At(0))->cd();
    if(!GetGraphPlot(&xy[0], kTrackTRD , 0)) break;
    xy[0] = -.3; xy[1] = -0.5; xy[2] = .3; xy[3] = 2.5;
    ((TVirtualPad*)l->At(1))->cd();
    if(!GetGraphPlot(&xy[0], kTrackTRD , 1)) break;
    return kTRUE;
  case 3: // kTrackTRD  z
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    xy[0] = -1.; xy[1] = -1000.; xy[2] = 1.; xy[3] = 4000.;
    ((TVirtualPad*)l->At(0))->cd();
    if(!GetGraphPlot(&xy[0], kTrackTRD , 2)) break;
    xy[0] = -1.; xy[1] = -0.5; xy[2] = 1.; xy[3] = 2.5;
    ((TVirtualPad*)l->At(1))->cd();
    if(!GetGraphPlot(&xy[0], kTrackTRD , 3)) break;
    return kTRUE;
  case 4: // kTrackTRD  phi
    xy[0] = -.3; xy[1] = -5.; xy[2] = .3; xy[3] = 50.;
    if(GetGraphPlot(&xy[0], kTrackTRD , 4)) return kTRUE;
    break;
  case 5: // kTrackTPC y
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    xy[0] = -.3; xy[1] = -500.; xy[2] = .3; xy[3] = 1500.;
    pad = ((TVirtualPad*)l->At(0)); pad->cd();
    pad->SetMargin(0.1, 0.1, 0.1, 0.01);
    if(!GetGraphPlot(&xy[0], kTrackTPC, 0)) break;
    xy[0] = -.3; xy[1] = -0.5; xy[2] = .3; xy[3] = 2.5;
    pad=((TVirtualPad*)l->At(1)); pad->cd();
    pad->SetMargin(0.1, 0.1, 0.1, 0.01);
    if(!GetGraphPlot(&xy[0], kTrackTPC, 1)) break;
    return kTRUE;
  case 6: // kTrackTPC z
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    xy[0] = -1.; xy[1] = -1000.; xy[2] = 1.; xy[3] = 4000.;
    pad = ((TVirtualPad*)l->At(0)); pad->cd();
    pad->SetMargin(0.1, 0.1, 0.1, 0.01);
    if(!GetGraphPlot(&xy[0], kTrackTPC, 2)) break;
    xy[0] = -1.; xy[1] = -0.5; xy[2] = 1.; xy[3] = 2.5;
    pad = ((TVirtualPad*)l->At(1)); pad->cd();
    pad->SetMargin(0.1, 0.1, 0.1, 0.01);
    if(!GetGraphPlot(&xy[0], kTrackTPC, 3)) break;
    return kTRUE;
  case 7: // kTrackTPC phi
    xy[0] = -.3; xy[1] = -5.; xy[2] = .3; xy[3] = 50.;
    if(GetGraphPlot(&xy[0], kTrackTPC, 4)) return kTRUE;
    break;
  case 8: // kMCcluster pt <2 GeV/c
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    xy[0]=-.3; xy[1]=-50.; xy[2]=.3; xy[3]=650.;
    ((TVirtualPad*)l->At(0))->cd();
    if(!GetGraphPlot(&xy[0], kMCcluster, 0)) break;
    ((TVirtualPad*)l->At(1))->cd();
    if(!GetGraphPlot(&xy[0], kMCcluster, 1)) break;
    return kTRUE;
  case 9: // kMCcluster pt > 3 GeV/c
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    xy[0]=-.3; xy[1]=-50.; xy[2]=.3; xy[3]=650.;
    ((TVirtualPad*)l->At(0))->cd();
    if(!GetGraphPlot(&xy[0], kMCcluster, 2)) break;
    xy[0] = -.3; xy[1] = -0.5; xy[2] = .3; xy[3] = 2.5;
    ((TVirtualPad*)l->At(1))->cd();
    if(!GetGraphPlot(&xy[0], kMCcluster, 3)) break;
    return kTRUE;
  case 10: //kMCtracklet [y] pt < 3 GeV/c
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    xy[0]=-.2; xy[1]=-50.; xy[2]=.2; xy[3] =250.;
    ((TVirtualPad*)l->At(0))->cd();
    if(!GetGraphPlot(&xy[0], kMCtracklet, 0)) break;
    ((TVirtualPad*)l->At(1))->cd();
    if(!GetGraphPlot(&xy[0], kMCtracklet, 1)) break;
    return kTRUE;
  case 11: //kMCtracklet [y] pt > 3 GeV/c
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    xy[0]=-.2; xy[1]=-50.; xy[2]=.2; xy[3] =250.;
    ((TVirtualPad*)l->At(0))->cd();
    if(!GetGraphPlot(&xy[0], kMCtracklet, 2)) break;
    xy[0] = -.2; xy[1] = -0.5; xy[2] = .2; xy[3] = 2.5;
    ((TVirtualPad*)l->At(1))->cd();
    if(!GetGraphPlot(&xy[0], kMCtracklet, 3)) break;
    return kTRUE;
  case 12: //kMCtracklet [z]
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    xy[0]=-1.; xy[1]=-100.; xy[2]=1.; xy[3] =2500.;
    ((TVirtualPad*)l->At(0))->cd();
    if(!GetGraphPlot(&xy[0], kMCtracklet, 4)) break;
    xy[0] = -1.; xy[1] = -0.5; xy[2] = 1.; xy[3] = 2.5;
    ((TVirtualPad*)l->At(1))->cd();
    if(!GetGraphPlot(&xy[0], kMCtracklet, 5)) break;
    return kTRUE;
  case 13: //kMCtracklet [phi]
    xy[0]=-.3; xy[1]=-3.; xy[2]=.3; xy[3] =25.;
    if(!GetGraphPlot(&xy[0], kMCtracklet, 6)) break;
    return kTRUE;
  case 14: //kMCtrackTRD [y]
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    xy[0]=-.2; xy[1]=-50.; xy[2]=.2; xy[3] =400.;
    ((TVirtualPad*)l->At(0))->cd();
    if(!GetGraphPlot(&xy[0], kMCtrackTRD, 0)) break;
    xy[0] = -.2; xy[1] = -0.5; xy[2] = .2; xy[3] = 3.5;
    ((TVirtualPad*)l->At(1))->cd();
    if(!GetGraphPlot(&xy[0], kMCtrackTRD, 1)) break;
    return kTRUE;
  case 15: //kMCtrackTRD [z]
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    xy[0]=-1.; xy[1]=-700.; xy[2]=1.; xy[3] =1500.;
    ((TVirtualPad*)l->At(0))->cd();
    if(!GetGraphPlot(&xy[0], kMCtrackTRD, 2)) break;
    xy[0] = -1.; xy[1] = -0.5; xy[2] = 1.; xy[3] = 2.5;
    ((TVirtualPad*)l->At(1))->cd();
    if(!GetGraphPlot(&xy[0], kMCtrackTRD, 3)) break;
    return kTRUE;
  case 16: //kMCtrackTRD [phi/snp]
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    xy[0]=-.2; xy[1]=-0.2; xy[2]=.2; xy[3] =2.;
    ((TVirtualPad*)l->At(0))->cd();
    if(!GetGraphPlot(&xy[0], kMCtrackTRD, 4)) break;
    xy[0] = -.2; xy[1] = -0.5; xy[2] = .2; xy[3] = 2.5;
    ((TVirtualPad*)l->At(1))->cd();
    if(!GetGraphPlot(&xy[0], kMCtrackTRD, 5)) break;
    return kTRUE;
  case 17: //kMCtrackTRD [theta/tgl]
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    xy[0]=-1.; xy[1]=-0.5; xy[2]=1.; xy[3] =5.;
    ((TVirtualPad*)l->At(0))->cd();
    if(!GetGraphPlot(&xy[0], kMCtrackTRD, 6)) break;
    xy[0] = -.2; xy[1] = -0.5; xy[2] = .2; xy[3] = 2.5;
    ((TVirtualPad*)l->At(1))->cd();
    if(!GetGraphPlot(&xy[0], kMCtrackTRD, 7)) break;
    return kTRUE;
  case 18: //kMCtrackTRD [pt]
    xy[0] = 0.2; xy[1] = -.7; xy[2] = 7.; xy[3] = 4.;
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    pad = (TVirtualPad*)l->At(0); pad->cd();
    pad->SetMargin(0.125, 0.015, 0.1, 0.015);
    if(!GetGraphTrack(&xy[0], 8, first, nspecies)) break;
    pad->Modified(); pad->Update(); pad->SetLogx();
    pad = (TVirtualPad*)l->At(1); pad->cd();
    pad->SetMargin(0.125, 0.015, 0.1, 0.015);
    if(!GetGraphTrack(&xy[0], 8, 55+first, nspecies, kTRUE)) break;
    pad->Modified(); pad->Update(); pad->SetLogx();
    return kTRUE;
  case 19: //kMCtrackTRD [1/pt] pulls
    xy[0] = 0.; xy[1] = -1.; xy[2] = 2.; xy[3] = 3.5;
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    pad = (TVirtualPad*)l->At(0); pad->cd();
    pad->SetMargin(0.125, 0.015, 0.1, 0.015);
    if(!GetGraphTrack(&xy[0], 9, first, nspecies)) break;
    pad = (TVirtualPad*)l->At(1); pad->cd();
    pad->SetMargin(0.125, 0.015, 0.1, 0.015);
    if(!GetGraphTrack(&xy[0], 9, 55+first, nspecies, kTRUE)) break;
    return kTRUE;
  case 20: //kMCtrackTRD [p]
    xy[0] = 0.2; xy[1] = -.7; xy[2] = 7.; xy[3] = 4.;
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    pad = (TVirtualPad*)l->At(0); pad->cd();
    pad->SetMargin(0.125, 0.015, 0.1, 0.015);
    if(!GetGraphTrack(&xy[0], 10, first, nspecies)) break;
    pad->Modified(); pad->Update(); pad->SetLogx();
    pad = (TVirtualPad*)l->At(1); pad->cd();
    pad->SetMargin(0.125, 0.015, 0.1, 0.015);
    if(!GetGraphTrack(&xy[0], 10, 55+first, nspecies, kTRUE)) break;
    pad->Modified(); pad->Update(); pad->SetLogx();
    return kTRUE;
  case 21: //kMCtrackTRD - SA [pt]
    xy[0] = 0.; xy[1] = -5.; xy[2] = 12.; xy[3] = 7.;
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    pad = (TVirtualPad*)l->At(0); pad->cd();
    pad->SetMargin(0.125, 0.015, 0.1, 0.015);
    if(!GetGraphTrack(&xy[0], 11, first, nspecies)) break;
    pad->SetLogx();
    pad = (TVirtualPad*)l->At(1); pad->cd();
    pad->SetMargin(0.125, 0.015, 0.1, 0.015);
    if(!GetGraphTrack(&xy[0], 11, 55+first, nspecies, kTRUE)) break;
    pad->SetLogx();
    return kTRUE;
  case 22: //kMCtrackTRD - SA [1/pt] pulls
    xy[0] = 0.; xy[1] = -1.5; xy[2] = 2.; xy[3] = 2.;
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    pad = (TVirtualPad*)l->At(0); pad->cd();
    pad->SetMargin(0.125, 0.015, 0.1, 0.015);
    if(!GetGraphTrack(&xy[0], 12, first, nspecies)) break;
    pad = (TVirtualPad*)l->At(1); pad->cd();
    pad->SetMargin(0.125, 0.015, 0.1, 0.015);
    if(!GetGraphTrack(&xy[0], 12, 55+first, nspecies, kTRUE)) break;
    return kTRUE;
  case 23: //kMCtrackTRD - SA [p]
    xy[0] = 0.; xy[1] = -7.5; xy[2] = 12.; xy[3] = 10.5;
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    pad = (TVirtualPad*)l->At(0); pad->cd();
    pad->SetMargin(0.125, 0.015, 0.1, 0.015);
    if(!GetGraphTrack(&xy[0], 13, first, nspecies)) break;
    pad->SetLogx();
    pad = (TVirtualPad*)l->At(1); pad->cd();
    pad->SetMargin(0.125, 0.015, 0.1, 0.015);
    if(!GetGraphTrack(&xy[0], 13, 55+first, nspecies, kTRUE)) break;
    pad->SetLogx();
    return kTRUE;
  case 24: // kMCtrackTPC [y]
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    xy[0]=-.25; xy[1]=-50.; xy[2]=.25; xy[3] =800.;
    ((TVirtualPad*)l->At(0))->cd();
    if(!GetGraphPlot(&xy[0], kMCtrackTPC, 0)) break;
    xy[0] = -.25; xy[1] = -0.5; xy[2] = .25; xy[3] = 2.5;
    ((TVirtualPad*)l->At(1))->cd();
    if(!GetGraphPlot(&xy[0], kMCtrackTPC, 1)) break;
    return kTRUE;
  case 25: // kMCtrackTPC [z]
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    xy[0]=-1.; xy[1]=-500.; xy[2]=1.; xy[3] =800.;
    ((TVirtualPad*)l->At(0))->cd();
    if(!GetGraphPlot(&xy[0], kMCtrackTPC, 2)) break;
    xy[0] = -1.; xy[1] = -0.5; xy[2] = 1.; xy[3] = 2.5;
    ((TVirtualPad*)l->At(1))->cd();
    if(!GetGraphPlot(&xy[0], kMCtrackTPC, 3)) break;
    return kTRUE;
  case 26: // kMCtrackTPC [phi|snp]
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    xy[0]=-.25; xy[1]=-0.5; xy[2]=.25; xy[3] =2.5;
    ((TVirtualPad*)l->At(0))->cd();
    if(!GetGraphPlot(&xy[0], kMCtrackTPC, 4)) break;
    xy[0] = -.25; xy[1] = -0.5; xy[2] = .25; xy[3] = 1.5;
    ((TVirtualPad*)l->At(1))->cd();
    if(!GetGraphPlot(&xy[0], kMCtrackTPC, 5)) break;
    return kTRUE;
  case 27: // kMCtrackTPC [theta|tgl]
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    xy[0]=-1.; xy[1]=-1.; xy[2]=1.; xy[3] =4.;
    ((TVirtualPad*)l->At(0))->cd();
    if(!GetGraphPlot(&xy[0], kMCtrackTPC, 6)) break;
    xy[0] = -1.; xy[1] = -0.5; xy[2] = 1.; xy[3] = 1.5;
    ((TVirtualPad*)l->At(1))->cd();
    if(!GetGraphPlot(&xy[0], kMCtrackTPC, 7)) break;
    return kTRUE;
  case 28: // kMCtrackTPC [pt]
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    xy[0] = 0.2; xy[1] = -.8; xy[2] = 7.; xy[3] = 6.;
    pad=(TVirtualPad*)l->At(0); pad->cd();
    pad->SetMargin(0.125, 0.015, 0.1, 0.015);
    if(!GetGraphTrackTPC(xy, 8, first, nspecies)) break;
    pad->SetLogx();
    xy[0]=0.; xy[1]=-0.5; xy[2]=2.; xy[3] =3.;
    pad = (TVirtualPad*)l->At(1); pad->cd();
    pad->SetMargin(0.125, 0.015, 0.1, 0.015);
    if(!GetGraphTrackTPC(xy, 9, first, nspecies, kTRUE)) break;
    return kTRUE;
  case 29: // kMCtrackTPC [p]
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    xy[0] = 0.2; xy[1] = -.8; xy[2] = 7.; xy[3] = 6.;
    pad = ((TVirtualPad*)l->At(0));pad->cd();
    pad->SetMargin(0.125, 0.015, 0.1, 0.015);
    if(!GetGraphTrackTPC(xy, 10, first, nspecies)) break;
    pad->SetLogx();
    xy[0]=0.; xy[1]=-0.5; xy[2]=8.; xy[3] =2.5;
    pad = ((TVirtualPad*)l->At(1)); pad->cd();
    pad->SetMargin(0.125, 0.015, 0.1, 0.015);
    if(!GetGraphTrackTPC(xy, 11, first, nspecies, kTRUE)) break;
    return kTRUE;
  case 30:  // kMCtrackTOF [z]
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
  Int_t nc(0);
  TGraph *gm= NULL, *gs= NULL;
  if(!fGraphS && !fGraphM){ 
    TObjArray *aM(NULL), *aS(NULL);
    Int_t n = fContainer->GetEntriesFast();
    fGraphS = new TObjArray(n); fGraphS->SetOwner();
    fGraphM = new TObjArray(n); fGraphM->SetOwner();
    for(Int_t ig=0; ig<n; ig++){
      fGraphM->AddAt(aM = new TObjArray(fgNproj[ig]), ig);
      fGraphS->AddAt(aS = new TObjArray(fgNproj[ig]), ig);

      for(Int_t ic=0; ic<fgNproj[ig]; ic++){
        //printf("g[%d] c[%d] n[%d]\n", ig, ic, fgNcomp[nc]);
        if(fgNcomp[nc]>1){
          TObjArray *agS(NULL), *agM(NULL);
          aS->AddAt(agS = new TObjArray(fgNcomp[nc]), ic); 
          aM->AddAt(agM = new TObjArray(fgNcomp[nc]), ic); 
          for(Int_t is=fgNcomp[nc]; is--;){
            agS->AddAt(gs = new TGraphErrors(), is);
            Int_t is0(is%11);
            gs->SetMarkerStyle(fgMarker[is0]);
            gs->SetMarkerColor(fgColorS[is0]);
            gs->SetLineColor(fgColorS[is0]);
            gs->SetNameTitle(Form("s_%d%02d%02d", ig, ic, is), fgParticle[is0]);

            agM->AddAt(gm = new TGraphErrors(), is);
            gm->SetMarkerStyle(fgMarker[is0]);
            gm->SetMarkerColor(fgColorM[is0]);
            gm->SetLineColor(fgColorM[is0]);
            gm->SetNameTitle(Form("m_%d%02d%02d", ig, ic, is), fgParticle[is0]);
          }
        } else {
          aS->AddAt(gs = new TGraphErrors(), ic);
          gs->SetMarkerStyle(23);
          gs->SetMarkerColor(kRed);
          gs->SetLineColor(kRed);
          gs->SetNameTitle(Form("s_%d%02d", ig, ic), "");
  
          aM->AddAt(gm = ig ? (TGraph*)new TGraphErrors() : (TGraph*)new TGraphAsymmErrors(), ic);
          gm->SetLineColor(kBlack);
          gm->SetMarkerStyle(7);
          gm->SetMarkerColor(kBlack);
          gm->SetNameTitle(Form("m_%d%02d", ig, ic), "");
        }
        nc+=fgNcomp[ic];
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
  Process2D(kCluster, 0, &fg, 1.e4); 
  Process2D(kCluster, 1, &fg); 
  fNRefFigures = 1;
  // Tracklet residual/pulls
  Process2D(kTrackTRD , 0, &fg, 1.e4); 
  Process2D(kTrackTRD , 1, &fg); 
  Process2D(kTrackTRD , 2, &fg, 1.e4); 
  Process2D(kTrackTRD , 3, &fg); 
  Process2D(kTrackTRD , 4, &fg, 1.e3); 
  fNRefFigures = 4;
  // TPC track residual/pulls
  Process2D(kTrackTPC, 0, &fg, 1.e4); 
  Process2D(kTrackTPC, 1, &fg); 
  Process2D(kTrackTPC, 2, &fg, 1.e4); 
  Process2D(kTrackTPC, 3, &fg); 
  Process2D(kTrackTPC, 4, &fg, 1.e3); 
  fNRefFigures = 7;

  if(!HasMCdata()) return kTRUE;


  //PROCESS MC RESIDUAL DISTRIBUTIONS

  // CLUSTER Y RESOLUTION/PULLS
  Process3Drange(kMCcluster, 0, 0, &fg, 1.e4, 1, 1);
  Process3Drange(kMCcluster, 0, 1, &fg, 1.e4, 2, 2);
  Process3Drange(kMCcluster, 0, 2, &fg, 1.e4, 3, 12);
  Process2D(kMCcluster, 1, &fg, 1., 3);
  fNRefFigures = 10;

  // TRACKLET RESOLUTION/PULLS
  Process3Drange(kMCtracklet, 0, 0, &fg, 1.e4, 1, 1); // y
  Process3Drange(kMCtracklet, 0, 1, &fg, 1.e4, 2, 2); // y
  Process3Drange(kMCtracklet, 0, 2, &fg, 1.e4, 3, 12); // y
  Process2D(kMCtracklet, 1, &fg, 1., 3);       // y pulls
  Process2D(kMCtracklet, 2, &fg, 1.e4, 4); // z
  Process2D(kMCtracklet, 3, &fg, 1., 5);       // z pulls
  Process2D(kMCtracklet, 4, &fg, 1.e3, 6); // phi
  fNRefFigures = 14;

  // TRACK RESOLUTION/PULLS
  Process2D(kMCtrackTRD, 0, &fg, 1.e4);   // y
  Process2D(kMCtrackTRD, 1, &fg);         // y PULL
  Process2D(kMCtrackTRD, 2, &fg, 1.e4);   // z
  Process2D(kMCtrackTRD, 3, &fg);         // z PULL
  Process2D(kMCtrackTRD, 4, &fg, 1.e3);   // phi
  Process2D(kMCtrackTRD, 5, &fg);         // snp PULL
  Process2D(kMCtrackTRD, 6, &fg, 1.e3);   // theta
  Process2D(kMCtrackTRD, 7, &fg);         // tgl PULL
  Process4D(kMCtrackTRD, 8, &fg, 1.e2);   // pt resolution
  Process4D(kMCtrackTRD, 8, &fpt, 1.e2, 4);// pt resolution e1- @ L0
  Process4D(kMCtrackTRD, 8, &fpt, 1.e2, 6);// pt resolution e1+ @ L0
  Process4D(kMCtrackTRD, 8, &fpt, 1.e2, 55+4);// pt resolution e1- @ L5
  Process4D(kMCtrackTRD, 8, &fpt, 1.e2, 55+6);// pt resolution e1+ @ L5
  Process4D(kMCtrackTRD, 9, &fg);         // 1/pt pulls
  Process4D(kMCtrackTRD, 10, &fg, 1.e2);  // p resolution
  Process4D(kMCtrackTRD, 10, &fpt, 1.e2, 4);// p resolution e1- @ L0
  Process4D(kMCtrackTRD, 10, &fpt, 1.e2, 6);// p resolution e1+ @ L0
  Process4D(kMCtrackTRD, 10, &fpt, 1.e2, 55+4);// p resolution e1- @ L5
  Process4D(kMCtrackTRD, 10, &fpt, 1.e2, 55+6);// p resolution e1+ @ L5
  Process4D(kMCtrackTRD, 11, &fg, 1.e2);   // pt resolution stand alone
  Process4D(kMCtrackTRD, 12, &fg);         // 1/pt pulls  stand alone
  Process4D(kMCtrackTRD, 13, &fg, 1.e2);  // p resolution  stand alone
  fNRefFigures = 24;

  // TRACK TPC RESOLUTION/PULLS
  Process2D(kMCtrackTPC, 0, &fg, 1.e4);// y resolution
  Process2D(kMCtrackTPC, 1, &fg);      // y pulls
  Process2D(kMCtrackTPC, 2, &fg, 1.e4);// z resolution
  Process2D(kMCtrackTPC, 3, &fg);      // z pulls
  Process2D(kMCtrackTPC, 4, &fg, 1.e3);// phi resolution
  Process2D(kMCtrackTPC, 5, &fg);      // snp pulls
  Process2D(kMCtrackTPC, 6, &fg, 1.e3);// theta resolution
  Process2D(kMCtrackTPC, 7, &fg);      // tgl pulls
  Process3D(kMCtrackTPC, 8, &fg, 1.e2);// pt resolution
  Process3D(kMCtrackTPC, 9, &fg);      // 1/pt pulls
  Process3D(kMCtrackTPC, 10, &fg, 1.e2);// p resolution
  Process3D(kMCtrackTPC, 11, &fg);      // p pulls
  fNRefFigures = 30;

  // TRACK HMPID RESOLUTION/PULLS
  Process2D(kMCtrackTOF, 0, &fg, 1.e4); // z towards TOF
  Process2D(kMCtrackTOF, 1, &fg);       // z towards TOF
  fNRefFigures = 31;

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
TObjArray* AliTRDresolution::BuildMonitorContainerCluster(const char* name)
{
// Build performance histograms for AliTRDcluster.vs TRD track or MC
//  - y reziduals/pulls

  TObjArray *arr = new TObjArray(2);
  arr->SetName(name); arr->SetOwner();
  TH1 *h(NULL); char hname[100], htitle[300];

  const Int_t kNpt(14), kNphi(48), kNdy(60);
  Float_t Phi=-.48, Dy=-.3, Pt=0.1;
  Float_t binsPhi[kNphi+1], binsDy[kNdy+1], binsPt[kNpt+1];
  for(Int_t i=0; i<kNphi+1; i++,Phi+=.02) binsPhi[i]=Phi;
  for(Int_t i=0; i<kNdy+1; i++,Dy+=.01) binsDy[i]=Dy;
  for(Int_t i=0;i<kNpt+1; i++,Pt=TMath::Exp(i*.15)-1.) binsPt[i]=Pt;

  // tracklet resolution/pull in y direction
  sprintf(hname, "%s_%s_Y", GetNameId(), name);
  sprintf(htitle, "Y res for \"%s\" @ %s;tg(#phi);#Delta y [cm];p_{t} [GeV/c]", GetNameId(), name);
  if(!(h = (TH3S*)gROOT->FindObject(hname))){
    h = new TH3S(hname, htitle, 
                 kNphi, binsPhi, kNdy, binsDy, kNpt, binsPt);
  } else h->Reset();
  arr->AddAt(h, 0);
  sprintf(hname, "%s_%s_Ypull", GetNameId(), name);
  sprintf(htitle, "Y pull for \"%s\" @ %s;tg(#phi);#Delta y  / #sigma_{y};entries", GetNameId(), name);
  if(!(h = (TH2I*)gROOT->FindObject(hname))){
    h = new TH2I(hname, htitle, 21, -.33, .33, 100, -4.5, 4.5);
  } else h->Reset();
  arr->AddAt(h, 1);

  return arr;
}

//________________________________________________________
TObjArray* AliTRDresolution::BuildMonitorContainerTracklet(const char* name)
{
// Build performance histograms for AliExternalTrackParam.vs TRD tracklet
//  - y reziduals/pulls
//  - z reziduals/pulls
//  - phi reziduals
  TObjArray *arr = BuildMonitorContainerCluster(name); 
  arr->Expand(5);
  TH1 *h(NULL); char hname[100], htitle[300];

  // tracklet resolution/pull in z direction
  sprintf(hname, "%s_%s_Z", GetNameId(), name);
  sprintf(htitle, "Z res for \"%s\" @ %s;tg(#theta);#Delta z [cm];entries", GetNameId(), name);
  if(!(h = (TH2I*)gROOT->FindObject(hname))){
    h = new TH2I(hname, htitle, 50, -1., 1., 100, -1.5, 1.5);
  } else h->Reset();
  arr->AddAt(h, 2);
  sprintf(hname, "%s_%s_Zpull", GetNameId(), name);
  sprintf(htitle, "Z pull for \"%s\" @ %s;tg(#theta);#Delta z  / #sigma_{z};entries", GetNameId(), name);
  if(!(h = (TH2I*)gROOT->FindObject(hname))){
    h = new TH2I(hname, htitle, 50, -1., 1., 100, -5.5, 5.5);
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
  } else h->Reset();
  arr->AddAt(h, 8);
  // 1/Pt pulls
  sprintf(hname, "%s_%s_1Pt", GetNameId(), name);
  sprintf(htitle, "1/P_{t} pull for \"%s\" @ %s;1/p_{t}^{MC} [c/GeV];#Delta(1/p_{t})/#sigma(1/p_{t});SPECIES", GetNameId(), name);
  if(!(h = (TH3S*)gROOT->FindObject(hname))){
    h = new TH3S(hname, htitle, 
                 kNpt, 0., 2., 100, -4., 4., kNspc, -5.5, 5.5);
  } else h->Reset();
  arr->AddAt(h, 9);
  // P resolution
  sprintf(hname, "%s_%s_P", GetNameId(), name);
  sprintf(htitle, "P res for \"%s\" @ %s;p [GeV/c];#Delta p/p^{MC};SPECIES", GetNameId(), name);
  if(!(h = (TH3S*)gROOT->FindObject(hname))){
    h = new TH3S(hname, htitle, 
                 kNpt, binsPt, kNdpt, binsDPt, kNspc, binsSpc);
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
  fContainer->AddAt(arr = new TObjArray(fgNhistos[kCharge]), kCharge);
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
  fContainer->AddAt(BuildMonitorContainerTracklet("Trk"), kTrackTRD);
  // tracklet to TRDin
  fContainer->AddAt(BuildMonitorContainerTracklet("TrkIN"), kTrackTPC);


  // Resolution histos
  if(!HasMCdata()) return fContainer;

  // cluster resolution 
  fContainer->AddAt(BuildMonitorContainerCluster("MCcl"),  kMCcluster);

  // tracklet resolution
  fContainer->AddAt(BuildMonitorContainerTracklet("MCtracklet"), kMCtracklet);

  // track resolution
  fContainer->AddAt(arr = new TObjArray(6/*fgNhistos[kMCtrackTRD]*/), kMCtrackTRD);
  arr->SetName("MCtrk");
  for(Int_t il(0); il<AliTRDgeometry::kNlayer; il++) arr->AddAt(BuildMonitorContainerTrack(Form("MCtrk_Ly%d", il)), il);

  // TRDin TRACK RESOLUTION
  fContainer->AddAt(BuildMonitorContainerTrack("MCtrkIN"), kMCtrackTPC);

  // TRDout TRACK RESOLUTION
  fContainer->AddAt(BuildMonitorContainerTrack("MCtrkOUT"), kMCtrackTOF);

  return fContainer;
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

  for(Int_t ibin = 1; ibin <= h2->GetNbinsX(); ibin++){
    Double_t x = h2->GetXaxis()->GetBinCenter(ibin);
    TH1D *h = h2->ProjectionY(pn, ibin, ibin);
    if(h->GetEntries()<100) continue;
    //AdjustF1(h, f);

    h->Fit(f, "QN");
    
    Int_t ip = g[0]->GetN();
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
  TH2I *h2 = idx<0 ? (TH2I*)(fContainer->At(plot)) : (TH2I*)((TObjArray*)(fContainer->At(plot)))->At(idx);
  if(!h2) return kFALSE;

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
  TH3S *h3 = idx<0 ? (TH3S*)(fContainer->At(plot)) : (TH3S*)((TObjArray*)(fContainer->At(plot)))->At(idx);
  if(!h3) return kFALSE;

  TObjArray *gm, *gs;
  if(!(gm = (TObjArray*)((TObjArray*)(fGraphM->At(plot)))->At(idx))) return kFALSE;
  if(!(gs = (TObjArray*)((TObjArray*)(fGraphS->At(plot)))->At(idx))) return kFALSE;
  TGraphErrors *g[2];

  TAxis *az = h3->GetZaxis();
  for(Int_t iz=1; iz<=az->GetNbins(); iz++){
    if(!(g[0] = (TGraphErrors*)gm->At(iz-1))) return kFALSE;
    if(!(g[1] = (TGraphErrors*)gs->At(iz-1))) return kFALSE;
    az->SetRange(iz, iz);
    if(!Process((TH2*)h3->Project3D("yx"), f, k, g)) return kFALSE;
  }

  return kTRUE;
}

//________________________________________________________
Bool_t AliTRDresolution::Process3Drange(ETRDresolutionPlot plot, Int_t hidx, Int_t gidx, TF1 *f, Float_t k, Int_t zbin0, Int_t zbin1)
{
  //
  // Do the processing
  //

  if(!fContainer || !fGraphS || !fGraphM) return kFALSE;

  // retrive containers
  TH3S *h3 = hidx<0 ? (TH3S*)(fContainer->At(plot)) : (TH3S*)((TObjArray*)(fContainer->At(plot)))->At(hidx);
  if(!h3) return kFALSE;

  TGraphErrors *g[2];
  if(!(g[0] = (TGraphErrors*)((TObjArray*)(fGraphM->At(plot)))->At(gidx))) return kFALSE;
  if(!(g[1] = (TGraphErrors*)((TObjArray*)(fGraphS->At(plot)))->At(gidx))) return kFALSE;

  TAxis *az = h3->GetZaxis();
  az->SetRange(zbin0, zbin1);
  if(!Process((TH2*)h3->Project3D("yx"), f, k, g)) return kFALSE;
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
Bool_t AliTRDresolution::Process4D(ETRDresolutionPlot plot, Int_t idx, TF1 *f, Float_t k, Int_t n)
{
  //
  // Do the processing
  //

  if(!fContainer || !fGraphS || !fGraphM) return kFALSE;
  //printf("Process4D : processing plot[%d] idx[%d]\n", plot, idx);

  // retrive containers
  TObjArray *arr = (TObjArray*)((TObjArray*)(fContainer->At(plot)))->At(idx);
  if(!arr) return kFALSE;

  TObjArray *gm, *gs;
  if(!(gm = (TObjArray*)((TObjArray*)(fGraphM->At(plot)))->At(idx))) return kFALSE;
  if(!(gs = (TObjArray*)((TObjArray*)(fGraphS->At(plot)))->At(idx))) return kFALSE;

  TGraphErrors *g[2];

  TH3S *h3(NULL);
  for(Int_t ix=0, in=0; ix<arr->GetEntriesFast(); ix++){
    if(!(h3 = (TH3S*)arr->At(ix))) return kFALSE;
    TAxis *az = h3->GetZaxis();
    //printf("  process ix[%d] bins[%d] in[%d]\n", ix, az->GetNbins(), in);
    for(Int_t iz=1; iz<=az->GetNbins(); iz++, in++){
      if(n>=0 && n!=in) continue;
      if(!(g[0] = (TGraphErrors*)gm->At(in))) return kFALSE;
      if(!(g[1] = (TGraphErrors*)gs->At(in))) return kFALSE;
      //printf("    g0[%s] g1[%s]\n", g[0]->GetName(), g[1]->GetName());
      az->SetRange(iz, iz);
      if(!Process((TH2*)h3->Project3D("yx"), f, k, g)) return kFALSE;
    }
  }

  return kTRUE;
}

//________________________________________________________
Bool_t AliTRDresolution::GetGraphPlot(Float_t *bb, ETRDresolutionPlot ip, Int_t idx)
{
  //
  // Get the graphs
  //

  if(!fGraphS || !fGraphM) return kFALSE;

  //printf("plotting task[%d] gidx[%d]\n", ip, idx);
  TGraphErrors *gm = idx<0 ? (TGraphErrors*)fGraphM->At(ip) : (TGraphErrors*)((TObjArray*)(fGraphM->At(ip)))->At(idx);
  if(!gm) return kFALSE;
  TGraphErrors *gs = idx<0 ? (TGraphErrors*)fGraphS->At(ip) : (TGraphErrors*)((TObjArray*)(fGraphS->At(ip)))->At(idx);
  if(!gs) return kFALSE;
  //printf("gs[%s] gm[%s]\n", gs->GetName(), gm->GetName());
  gs->Draw("apl"); gm->Draw("pl");
  //return kTRUE;
  // titles look up
  Int_t nref = 0;
  for(Int_t jp=0; jp<(Int_t)ip; jp++) nref+=fgNproj[jp];
  UChar_t jdx = idx<0?0:idx;
  for(Int_t jc=0; jc<TMath::Min(jdx,fgNproj[ip]-1); jc++) nref++;
  const Char_t **at = fgAxTitle[nref];

  Int_t n(0);
  if((n=gm->GetN())) {
    PutTrendValue(Form("%s_%s", fgPerformanceName[ip], at[0]), gm->GetMean(2));
    PutTrendValue(Form("%s_%sRMS", fgPerformanceName[ip], at[0]), gm->GetRMS(2));
  }

  if((n=gs->GetN())){
    gs->Sort(&TGraph::CompareY);
    PutTrendValue(Form("%s_%sSigMin", fgPerformanceName[ip], at[0]), gs->GetY()[0]);
    PutTrendValue(Form("%s_%sSigMax", fgPerformanceName[ip], at[0]), gs->GetY()[n-1]);
    gs->Sort(&TGraph::CompareX); 
  }

  // axis range
  TAxis *ax(NULL); TH1 *hf(NULL);
  hf = gs->GetHistogram();
  hf->SetTitle(at[0]);
  ax = hf->GetXaxis();
  ax->SetRangeUser(bb[0], bb[2]);
  ax->SetTitle(at[1]);ax->CenterTitle();

  ax = hf->GetYaxis();
  ax->SetRangeUser(bb[1], bb[3]);
  ax->SetTitleOffset(1.1);
  ax->SetTitle(at[2]);ax->CenterTitle();

  TGaxis *gax = NULL;
  gax = new TGaxis(bb[2], bb[1], bb[2], bb[3], bb[1], bb[3], 510, "+U");
  gax->SetLineColor(kRed);gax->SetLineWidth(2);gax->SetTextColor(kRed);
  //gax->SetVertical();
  gax->CenterTitle(); gax->SetTitleOffset(.7);
  gax->SetTitle(at[3]); gax->Draw();

  // bounding box
  TBox *b = new TBox(-.15, bb[1], .15, bb[3]);
  b->SetFillStyle(3002);b->SetLineColor(0);
  b->SetFillColor(ip<=Int_t(kMCcluster)?kGreen:kBlue);
  b->Draw();
  return kTRUE;
}

//________________________________________________________
Bool_t AliTRDresolution::GetGraphTrack(Float_t *bb, Int_t idx, Int_t first, Int_t n, Bool_t kLEG)
{
  //
  // Get the graphs
  //

  if(!fGraphS || !fGraphM) return kFALSE;

  // axis titles look up
  Int_t nref(0);
  for(Int_t jp=0; jp<Int_t(kMCtrackTRD); jp++) nref+=fgNproj[jp];
  nref+=idx;
  const Char_t **at = fgAxTitle[nref];
  //printf("nref[%d] ax[%s] x[%f %f] y[%f %f]\n", nref, at[0], bb[0], bb[2], bb[1], bb[3]);

  TLegend *legM(NULL), *legS(NULL);
  if(kLEG){
    legM=new TLegend(.35, .6, .65, .9);
    legM->SetHeader("Mean");
    legM->SetBorderSize(1);
    legM->SetFillColor(0);
    legS=new TLegend(.65, .6, .95, .9);
    legS->SetHeader("Sigma");
    legS->SetBorderSize(1);
    legS->SetFillColor(0);
  }
  Int_t layer(first/11);
  TH1S *h1(NULL);
  h1 = new TH1S(Form("h1TF_%02d", fIdxFrame++), Form("%s %d;%s;%s", at[0], layer, at[1], at[2]), 2, bb[0], bb[2]);
  h1->SetMinimum(bb[1]);h1->SetMaximum(bb[3]);
  h1->SetLineColor(kBlack); h1->SetLineWidth(1);h1->SetLineStyle(2); 
  // axis range
  TAxis *ax = h1->GetXaxis();
  ax->CenterTitle();ax->SetMoreLogLabels();ax->SetTitleOffset(1.2);
  ax = h1->GetYaxis();
  ax->SetRangeUser(bb[1], bb[3]);
  ax->CenterTitle(); ax->SetTitleOffset(1.4);
  h1->Draw();
//   TGaxis *gax = NULL;
//   gax = new TGaxis(bb[2], bb[1], bb[2], bb[3], bb[1], bb[3], 510, "+U");
//   gax->SetLineColor(kRed);gax->SetLineWidth(2);gax->SetTextColor(kRed);
//   //gax->SetVertical();
//   gax->CenterTitle(); //gax->SetTitleOffset(.5);gax->SetTitleSize(.06);
//   gax->SetTitle(at[3]); gax->Draw();

  TGraphErrors *gm = NULL, *gs = NULL;
  TObjArray *a0 = NULL, *a1 = NULL;
  a0 = (TObjArray*)((TObjArray*)fGraphM->At(kMCtrackTRD))->At(idx);
  a1 = (TObjArray*)((TObjArray*)fGraphS->At(kMCtrackTRD))->At(idx);
  Int_t nn(0), m(n/2);
  for(Int_t is=first, is0=0; is<first+n; is++, is0++){
    if(is0==m) continue; // no PID tracks
    if(is0==m-1||is0==m+1) continue; // electron tracks
    if(!(gs =  (TGraphErrors*)a1->At(is))) return kFALSE;
    if(!(gm =  (TGraphErrors*)a0->At(is))) return kFALSE;

    if((nn=gs->GetN())){
      gs->Draw("pc"); if(legS) legS->AddEntry(gs, gs->GetTitle(), "pl");
      gs->Sort(&TGraph::CompareY);
      PutTrendValue(Form("%s_%sSigMin%s", fgPerformanceName[kMCtrackTRD], at[0], AliPID::ParticleShortName(is0)), gs->GetY()[0]);
      PutTrendValue(Form("%s_%sSigMax%s", fgPerformanceName[kMCtrackTRD], at[0], AliPID::ParticleShortName(is0)), gs->GetY()[nn-1]);
      gs->Sort(&TGraph::CompareX); 
    }
    if(gm->GetN()){
      gm->Draw("pc");if(legM) legM->AddEntry(gm, gm->GetTitle(), "pl");
      PutTrendValue(Form("%s_%s_%s", fgPerformanceName[kMCtrackTRD], at[0], AliPID::ParticleShortName(is0)), gm->GetMean(2));
      PutTrendValue(Form("%s_%s_%sRMS", fgPerformanceName[kMCtrackTRD], at[0], AliPID::ParticleShortName(is0)), gm->GetRMS(2));
    }
  }
  if(kLEG){legM->Draw();legS->Draw();}
  return kTRUE;
}


//________________________________________________________
Bool_t AliTRDresolution::GetGraphTrackTPC(Float_t *bb, Int_t idx, Int_t ist, Int_t n, Bool_t kLEG)
{
  //
  // Get the graphs
  //

  if(!fGraphS || !fGraphM) return kFALSE;

  // axis titles look up
  Int_t nref = 0;
  for(Int_t jp=0; jp<Int_t(kMCtrackTPC); jp++) nref+=fgNproj[jp];
  nref+=idx;
  const Char_t **at = fgAxTitle[nref];

  TLegend *legM(NULL), *legS(NULL);
  if(kLEG){
    legM=new TLegend(.35, .6, .65, .9);
    legM->SetHeader("Mean");
    legM->SetBorderSize(1);
    legM->SetFillColor(0);
    legS=new TLegend(.65, .6, .95, .9);
    legS->SetHeader("Sigma");
    legS->SetBorderSize(1);
    legS->SetFillColor(0);
  }
  TH1S *h1(NULL);
  h1 = new TH1S(Form("h1TF_%02d", fIdxFrame++), at[0], 2, bb[0], bb[2]);
  h1->SetMinimum(bb[1]);h1->SetMaximum(bb[3]);
  h1->SetLineColor(kBlack); h1->SetLineWidth(1);h1->SetLineStyle(2); 
  // axis range
  TAxis *ax = h1->GetXaxis();
  ax->SetTitle(at[1]);ax->CenterTitle();ax->SetMoreLogLabels();ax->SetTitleOffset(1.2);
  ax = h1->GetYaxis();
  ax->SetRangeUser(bb[1], bb[3]);
  ax->SetTitleOffset(1.4);//ax->SetTitleSize(.06);
  ax->SetTitle(at[2]);ax->CenterTitle();
  h1->Draw();
//   TGaxis *gax = NULL;
//   gax = new TGaxis(bb[2], bb[1], bb[2], bb[3], bb[1], bb[3], 510, "+U");
//   gax->SetLineColor(kRed);gax->SetLineWidth(2);gax->SetTextColor(kRed);
//   //gax->SetVertical();
//   gax->CenterTitle(); //gax->SetTitleOffset(.5);gax->SetTitleSize(.06);
//   gax->SetTitle(at[3]); gax->Draw();

  Int_t nn(0), m(n/2);
  TGraphErrors *gm = NULL, *gs = NULL;
  TObjArray *a0 = NULL, *a1 = NULL;
  a0 = (TObjArray*)((TObjArray*)fGraphM->At(kMCtrackTPC))->At(idx);
  a1 = (TObjArray*)((TObjArray*)fGraphS->At(kMCtrackTPC))->At(idx);
  for(Int_t is=ist, is0=0; is<ist+n; is++, is0++){
    if(is0==m) continue; // no PID tracks
    //if(is0==m-1||is0==m+1) continue; // electron tracks
    if(!(gs =  (TGraphErrors*)a1->At(is))) return kFALSE;
    if(!(gm =  (TGraphErrors*)a0->At(is))) return kFALSE;
    if((nn=gs->GetN())){
      gs->Draw("pl");if(legS) legS->AddEntry(gs, gs->GetTitle(), "pl");
      gs->Sort(&TGraph::CompareY);
      PutTrendValue(Form("%s_%sSigMin%s", fgPerformanceName[kMCtrackTPC], at[0], AliPID::ParticleShortName(is0)), gs->GetY()[0]);
      PutTrendValue(Form("%s_%sSigMax%s", fgPerformanceName[kMCtrackTPC], at[0], AliPID::ParticleShortName(is0)), gs->GetY()[nn-1]);
      gs->Sort(&TGraph::CompareX); 
    }
    if(gm->GetN()){
      gm->Draw("pl"); if(legM) legM->AddEntry(gm, gm->GetTitle(), "pl");
      PutTrendValue(Form("%s_%s_%s", fgPerformanceName[kMCtrackTPC], at[0], AliPID::ParticleShortName(is0)), gm->GetMean(2));
      PutTrendValue(Form("%s_%s_%sRMS", fgPerformanceName[kMCtrackTPC], at[0], AliPID::ParticleShortName(is0)), gm->GetRMS(2));
    }
  }
  if(kLEG) {legM->Draw(); legS->Draw();}
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
