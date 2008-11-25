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

/* $Id: AliTRDtrackingResolution.cxx 27496 2008-07-22 08:35:45Z cblume $ */

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
//   AliTRDtrackingResolution *res = new AliTRDtrackingResolution();
//   //res->SetMCdata();
//   //res->SetVerbose();
//   //res->SetVisual();
//   res->Load("TRD.TaskResolution.root");
//   if(!res->PostProcess()) return;
//   res->GetRefFigure(0);
// }  
//
//  Authors:                                                              //
//    Alexandru Bercuci <A.Bercuci@gsi.de>                                //
//    Markus Fasel <M.Fasel@gsi.de>                                       //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include <cstring>

#include <TSystem.h>
#include <TObjArray.h>
#include <TH2.h>
#include <TH1.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TProfile.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include "TTreeStream.h"
#include "TGeoManager.h"

#include "AliAnalysisManager.h"
#include "AliTrackReference.h"
#include "AliTrackPointArray.h"
#include "AliCDBManager.h"

#include "AliTRDSimParam.h"
#include "AliTRDgeometry.h"
#include "AliTRDpadPlane.h"
#include "AliTRDcluster.h"
#include "AliTRDseedV1.h"
#include "AliTRDtrackV1.h"
#include "AliTRDtrackerV1.h"
#include "AliTRDReconstructor.h"
#include "AliTRDrecoParam.h"

#include "AliTRDtrackInfo/AliTRDtrackInfo.h"
#include "AliTRDtrackingResolution.h"

ClassImp(AliTRDtrackingResolution)

//________________________________________________________
AliTRDtrackingResolution::AliTRDtrackingResolution()
  :AliTRDrecoTask("Resolution", "Tracking Resolution")
  ,fStatus(0)
  ,fReconstructor(0x0)
  ,fGeo(0x0)
  ,fGraphS(0x0)
  ,fGraphM(0x0)
{
  fReconstructor = new AliTRDReconstructor();
  fReconstructor->SetRecoParam(AliTRDrecoParam::GetLowFluxParam());
  fGeo = new AliTRDgeometry();
  InitFunctorList();
}

//________________________________________________________
AliTRDtrackingResolution::~AliTRDtrackingResolution()
{
  if(fGraphS){fGraphS->Delete(); delete fGraphS;}
  if(fGraphM){fGraphM->Delete(); delete fGraphM;}
  delete fGeo;
  delete fReconstructor;
  if(gGeoManager) delete gGeoManager;
}


//________________________________________________________
void AliTRDtrackingResolution::CreateOutputObjects()
{
  // spatial resolution
  OpenFile(0, "RECREATE");

  fContainer = Histos();
}

// //________________________________________________________
// void AliTRDtrackingResolution::Exec(Option_t *)
// {
//   // spatial Resolution: res = pos_{Tracklet}(x = x_{Anode wire}) - pos_{TrackRef}(x = x_{Anode wire})
//   // angular Resolution: res = Tracklet angle - TrackRef Angle
// 
//   Int_t nTrackInfos = fTracks->GetEntriesFast();
//   if(fDebugLevel>=2 && nTrackInfos){ 
//     printf("Event[%d] TrackInfos[%d]\n", (Int_t)AliAnalysisManager::GetAnalysisManager()->GetCurrentEntry(), nTrackInfos);
//   }
//   const Int_t kNLayers = AliTRDgeometry::kNlayer;
//   Int_t pdg, nly, ncrs;
//   Double_t p, dy, theta/*, dphi, dymc, dzmc, dphimc*/;
//   Float_t fP[kNLayers], fX[kNLayers], fY[kNLayers], fZ[kNLayers], fPhi[kNLayers], fTheta[kNLayers]; // phi/theta angle per layer
//   Bool_t fMCMap[kNLayers], fLayerMap[kNLayers]; // layer map
// 
//   AliTRDpadPlane *pp = 0x0;
//   AliTrackPoint tr[kNLayers], tk[kNLayers];
//   AliExternalTrackParam *fOp = 0x0;
//   AliTRDtrackV1 *fTrack = 0x0;
//   AliTRDtrackInfo *fInfo = 0x0;
//   for(Int_t iTI = 0; iTI < nTrackInfos; iTI++){
//     // check if ESD and MC-Information are available
//     if(!(fInfo = dynamic_cast<AliTRDtrackInfo *>(fTracks->UncheckedAt(iTI)))) continue;
//     if(!(fTrack = fInfo->GetTrack())) continue;
//     if(!(fOp = fInfo->GetOuterParam())) continue;
//     pdg = fInfo->GetPDG();
//     nly = 0; ncrs = 0; theta = 0.;
// 
//     if(fDebugLevel>=3) printf("\tDoing track[%d] NTrackRefs[%d]\n", iTI, fInfo->GetNTrackRefs());
// 
//     p = fOp->P();
//     Int_t npts = 0;
//     memset(fP, 0, kNLayers*sizeof(Float_t));
//     memset(fX, 0, kNLayers*sizeof(Float_t));
//     memset(fY, 0, kNLayers*sizeof(Float_t));
//     memset(fZ, 0, kNLayers*sizeof(Float_t));
//     memset(fPhi, 0, kNLayers*sizeof(Float_t));
//     memset(fTheta, 0, kNLayers*sizeof(Float_t));
//     memset(fLayerMap, 0, kNLayers*sizeof(Bool_t));
//     memset(fMCMap, 0, kNLayers*sizeof(Bool_t));
//     AliTRDseedV1 *fTracklet = 0x0;
//     for(Int_t iplane = 0; iplane < kNLayers; iplane++){
//       if(!(fTracklet = fTrack->GetTracklet(iplane))) continue;
//       if(!fTracklet->IsOK()) continue;
// 
//       // Book point arrays
//       fLayerMap[iplane] = kTRUE;
//       tr[npts].SetXYZ(fTracklet->GetX0(), 0., 0.);
//       tk[npts].SetXYZ(fTracklet->GetX0(), fTracklet->GetYfit(0), fTracklet->GetZfit(0));
//       npts++;
// 
//       if(fDebugLevel>=4) printf("\t\tLy[%d] X0[%6.3f] Ncl[%d]\n", iplane, fTracklet->GetX0(), fTracklet->GetN());
// 
//       // define reference values
//       fP[iplane]   = p;
//       fX[iplane]   = fTracklet->GetX0();
//       fY[iplane]   = fTracklet->GetYref(0);
//       fZ[iplane]   = fTracklet->GetZref(0);
//       fPhi[iplane] = TMath::ATan(fTracklet->GetYref(1));
//       fTheta[iplane] = TMath::ATan(fTracklet->GetZref(1));
//       
// 
//       // RESOLUTION (compare to MC)
//       if(HasMCdata()){
//         if(fInfo->GetNTrackRefs() >= 2){ 
//           Double_t pmc, ymc, zmc, phiMC, thetaMC;
//           if(Resolution(fTracklet, fInfo, pmc, ymc, zmc, phiMC, thetaMC)){ 
//             fMCMap[iplane] = kTRUE;
//             fP[iplane]     = pmc;
//             fY[iplane]     = ymc;
//             fZ[iplane]     = zmc;
//             fPhi[iplane]   = phiMC;
//             fTheta[iplane] = thetaMC;
//           }
//         }
//       }
//       Float_t phi   = fPhi[iplane]*TMath::RadToDeg();
//       theta += fTheta[iplane]; nly++;
//       if(fTracklet->GetNChange()) ncrs++;
// 
//       // Do clusters residuals
//       if(!fTracklet->Fit(kFALSE)) continue;
//       AliTRDcluster *c = 0x0;
//       for(Int_t ic=AliTRDseed::knTimebins-1; ic>=0; ic--){
//         if(!(c = fTracklet->GetClusters(ic))) continue;
// 
//         dy = fTracklet->GetYat(c->GetX()) - c->GetY();
//         ((TH2I*)fContainer->At(kClusterYResidual))->Fill(phi, dy);
// 
//         if(fDebugLevel>=2){
//           Float_t q = c->GetQ();
//           // Get z-position with respect to anode wire
//           AliTRDSimParam    *simParam    = AliTRDSimParam::Instance();
//           Int_t det = c->GetDetector();
//           Float_t x = c->GetX();
//           Float_t z = fZ[iplane]-(fX[iplane]-x)*TMath::Tan(fTheta[iplane]);
//           Int_t stack = fGeo->GetStack(det);
//           pp = fGeo->GetPadPlane(iplane, stack);
//           Float_t row0 = pp->GetRow0();
//           Float_t d  =  row0 - z + simParam->GetAnodeWireOffset();
//           d -= ((Int_t)(2 * d)) / 2.0;
//           if (d > 0.25) d  = 0.5 - d;
//   
//           (*fDebugStream) << "ResidualClusters"
//             << "ly="   << iplane
//             << "stk="  << stack
//             << "pdg="  << pdg
//             << "phi="  << fPhi[iplane]
//             << "tht="  << fTheta[iplane]
//             << "q="    << q
//             << "x="    << x
//             << "z="    << z
//             << "d="    << d
//             << "dy="   << dy
//             << "\n";
//         }
//       }
//       pp = 0x0;
//     }
//     if(nly) theta /= nly; 
//     if(fDebugLevel>=1){
//       (*fDebugStream) << "TrackStatistics"
//         << "nly="   << nly
//         << "ncrs="  << ncrs
//         << "tht="   << theta
//         << "\n";
//     }
// 
// 
// //     // this protection we might drop TODO
// //     if(fTrack->GetNumberOfTracklets() < 6) continue;
// // 
// //     AliTRDtrackerV1::FitRiemanTilt(fTrack, 0x0, kTRUE, npts, tr);
// //     Int_t iref = 0;
// //     for(Int_t ip=0; ip<kNLayers; ip++){
// //       if(!fLayerMap[ip]) continue;
// //       fTracklet = fTrack->GetTracklet(ip);
// //       // recalculate fit based on the new tilt correction
// //       fTracklet->Fit();
// // 
// //       dy = fTracklet->GetYfit(0) - tr[iref].GetY();
// //       ((TH2I*)fContainer->At(kTrackletRiemanYResidual))->Fill(fPhi[ip]*TMath::RadToDeg(), dy);
// // 
// //       dphi = fTracklet->GetYfit(1)- fTracklet->GetYref(1);
// //       ((TH2I*)fContainer->At(kTrackletRiemanAngleResidual))->Fill(fPhi[ip]*TMath::RadToDeg(), dphi);
// // 
// //       if(HasMCdata()){
// //         dymc = fY[ip] - tr[iref].GetY();
// //         ((TH2I*)fContainer->At(kTrackRYResolution))->Fill(fPhi[ip]*TMath::RadToDeg(), dymc);
// // 
// //         dzmc = fZ[ip] - tr[iref].GetZ();
// //         ((TH2I*)fContainer->At(kTrackRZResolution))->Fill(fPhi[ip]*TMath::RadToDeg(), dzmc);
// // 
// //         dphimc = fPhi[ip] - fTracklet->GetYfit(1);
// //         ((TH2I*)fContainer->At(kTrackRAngleResolution))->Fill(fPhi[ip]*TMath::RadToDeg(), dphimc);
// //       }
// // 
// //       iref++;
// // 
// //       if(fDebugLevel>=1){
// //         (*fDebugStream) << "RiemannTrack"
// //           << "ly="    << ip
// //           << "mc="    << fMCMap[ip]
// //           << "p="     << fP[ip]
// //           << "phi="   << fPhi[ip]
// //           << "tht="   << fTheta[ip]
// //           << "dy="    << dy
// //           << "dphi="  << dphi
// //           << "dymc="  << dymc
// //           << "dzmc="  << dzmc
// //           << "dphimc="<< dphimc
// //           << "\n";
// //       }
// //     }
// 
// //  if(!gGeoManager) TGeoManager::Import("geometry.root");
// //     AliTRDtrackerV1::FitKalman(fTrack, 0x0, kFALSE, nc, tr);
// //     for(Int_t ip=0; ip<nc; ip++){
// //       dy = cl[ip].GetY() - tr[ip].GetY();
// //      ((TH2I*)fContainer->At(kTrackletKalmanYResidual))->Fill(phi*TMath::RadToDeg(), dy);
// //       dz = cl[ip].GetZ() - tr[ip].GetZ();
// //       if(fDebugLevel>=1){
// //         (*fDebugStream) << "KalmanTrack"
// //           << "dy="		  << dy
// //           << "dz="	 	  << dz
// // /*          << "phi="			<< phi
// //           << "theta="		<< theta
// //           << "dphi="		<< dphi*/
// //           << "\n";
// //       }
// //     }    
// 
// 
//   }
//   PostData(0, fContainer);
// }

//________________________________________________________
TH1* AliTRDtrackingResolution::PlotClusterResiduals(const AliTRDtrackV1 *track)
{
  if(track) fTrack = track;
  if(!fTrack){
    AliWarning("No Track defined.");
    return 0x0;
  }
  TH1 *h = 0x0;
  if(!(h = ((TH2I*)fContainer->At(kClusterResidual)))){
    AliWarning("No output histogram defined.");
    return 0x0;
  }

  Int_t pdg = (HasMCdata() && fMC) ? fMC->GetPDG() : 0;
  UChar_t s;
  Float_t x0, y0, z0, dy, dydx, dzdx;
  AliTRDseedV1 *fTracklet = 0x0;  
  for(Int_t ily=0; ily<AliTRDgeometry::kNlayer; ily++){
    if(!(fTracklet = fTrack->GetTracklet(ily))) continue;
    if(!fTracklet->IsOK()) continue;
    x0 = fTracklet->GetX0();

    // retrive the track angle with the chamber
    if(HasMCdata() && fMC){ 
      if(!fMC->GetDirections(x0, y0, z0, dydx, dzdx, s)) continue; 
    }else{ 
      y0   = fTracklet->GetYref(0);
      z0   = fTracklet->GetZref(0);
      dydx = fTracklet->GetYref(1);
      dzdx = fTracklet->GetZref(1);
    }

    AliTRDseedV1 trklt(*fTracklet);
    if(!trklt.Fit(kFALSE)) continue;

    AliTRDcluster *c = 0x0;
    fTracklet->ResetClusterIter(kFALSE);
    while((c = fTracklet->PrevCluster())){
      dy = trklt.GetYat(c->GetX()) - c->GetY();
      h->Fill(dydx, dy);
  
      if(fDebugLevel>=1){
        Float_t q = c->GetQ();
        // Get z-position with respect to anode wire
        AliTRDSimParam    *simParam    = AliTRDSimParam::Instance();
        Int_t det = c->GetDetector();
        Float_t x = c->GetX();
        Float_t z = z0-(x0-x)*dzdx;
        Int_t istk = fGeo->GetStack(det);
        AliTRDpadPlane *pp = fGeo->GetPadPlane(ily, istk);
        Float_t row0 = pp->GetRow0();
        Float_t d  =  row0 - z + simParam->GetAnodeWireOffset();
        d -= ((Int_t)(2 * d)) / 2.0;
        if (d > 0.25) d  = 0.5 - d;
  
        (*fDebugStream) << "ClusterResidual"
          << "ly="   << ily
          << "stk="  << istk
          << "pdg="  << pdg
          << "dydx=" << dydx
          << "dzdx=" << dzdx
          << "q="    << q
          << "x="    << x
          << "z="    << z
          << "d="    << d
          << "dy="   << dy
          << "\n";
      }
    }
  }
  return h;
}

//________________________________________________________
TH1* AliTRDtrackingResolution::PlotResolution(const AliTRDtrackV1 *track)
{
  if(!HasMCdata()){ 
    AliWarning("No MC defined. Results will not be available.");
    return 0x0;
  }
  if(track) fTrack = track;
  if(!fTrack){
    AliWarning("No Track defined.");
    return 0x0;
  }
  TH1 *h = 0x0;
  if(!(h = ((TH2I*)fContainer->At(kClusterResolution)))){
    AliWarning("No output histogram defined.");
    return 0x0;
  }
  if(!(h = ((TH2I*)fContainer->At(kTrackletYResolution)))){
    AliWarning("No output histogram defined.");
    return 0x0;
  }
  if(!(h = ((TH2I*)fContainer->At(kTrackletZResolution)))){
    AliWarning("No output histogram defined.");
    return 0x0;
  }
  if(!(h = ((TH2I*)fContainer->At(kTrackletAngleResolution)))){
    AliWarning("No output histogram defined.");
    return 0x0;
  }
  //printf("called for %d tracklets ... \n", fTrack->GetNumberOfTracklets());
  UChar_t s;
  Int_t pdg = fMC->GetPDG(), det=-1;
  Float_t x0, y0, z0, dy, dydx, dzdx;
  AliTRDseedV1 *fTracklet = 0x0;  
  for(Int_t ily=0; ily<AliTRDgeometry::kNlayer; ily++){
    if(!(fTracklet = fTrack->GetTracklet(ily))) continue;
    if(!fTracklet->IsOK()) continue;
    //printf("process layer[%d]\n", ily);
    // retrive the track position and direction within the chamber
    det = fTracklet->GetDetector();
    x0  = fTracklet->GetX0();
    if(!fMC->GetDirections(x0, y0, z0, dydx, dzdx, s)) continue; 

    // recalculate tracklet based on the MC info
    AliTRDseedV1 tt(*fTracklet);
    tt.SetZref(0, z0);
    tt.SetZref(1, dzdx); 
    if(!tt.Fit(kFALSE)) continue;
    //tt.Update();
    dy = tt.GetYfit(0) - y0;
    Float_t dphi   = TMath::ATan(tt.GetYfit(1)) - TMath::ATan(dydx);
    Float_t dz = 100.;
    Bool_t  cross = fTracklet->GetNChange();
    if(cross){
      Double_t *xyz = tt.GetCrossXYZ();
      dz = xyz[2] - (z0 - (x0 - xyz[0])*dzdx) ;
      ((TH2I*)fContainer->At(kTrackletZResolution))->Fill(dzdx, dz);
    }
  
    // Fill Histograms
    ((TH2I*)fContainer->At(kTrackletYResolution))->Fill(dydx, dy);
    ((TH2I*)fContainer->At(kTrackletAngleResolution))->Fill(dydx, dphi*TMath::RadToDeg());

    // Fill Debug stream
    if(fDebugLevel>=1){
      Float_t p = fMC->GetTrackRef() ? fMC->GetTrackRef()->P() : -1.;
      (*fDebugStream) << "TrkltResolution"
        << "det="	 	  << det
        << "pdg="     << pdg
        << "p="       << p
        << "ymc="     << y0
        << "zmc="     << z0
        << "dydx="    << dydx
        << "dzdx="    << dzdx
        << "cross="   << cross
        << "dy="		  << dy
        << "dz="	 	  << dz
        << "dphi="		<< dphi
        << "\n";
    }

    Int_t istk = AliTRDgeometry::GetStack(det); 
    AliTRDpadPlane *pp = fGeo->GetPadPlane(ily, istk);
    Float_t zr0 = pp->GetRow0() + AliTRDSimParam::Instance()->GetAnodeWireOffset();
    Float_t tilt = fTracklet->GetTilt();

    AliTRDcluster *c = 0x0;
    fTracklet->ResetClusterIter(kFALSE);
    while((c = fTracklet->PrevCluster())){
      Float_t  q = TMath::Abs(c->GetQ());
      Float_t xc = c->GetX();
      Float_t yc = c->GetY();
      Float_t zc = c->GetZ();
      Float_t dx = x0 - xc; 
      Float_t yt = y0 - dx*dydx;
      Float_t zt = z0 - dx*dzdx; 
      dy = yt - (yc - tilt*(zc-zt));

      // Fill Histograms
      if(q>100.) ((TH2I*)fContainer->At(kClusterResolution))->Fill(dydx, dy);
      
      // Fill Debug Tree
      if(fDebugLevel>=1){
        Float_t d = zr0 - zt;
        d -= ((Int_t)(2 * d)) / 2.0;
        if (d > 0.25) d  = 0.5 - d;
        Int_t label = fMC->GetLabel();
        (*fDebugStream) << "ClusterResolution"
          << "det="  << det
          << "pdg="  << pdg
          << "dydx=" << dydx
          << "dzdx=" << dzdx
          << "q="    << q
          << "d="    << d
          << "dy="   << dy
          << "xc="   << xc
          << "yc="   << yc
          << "zc="   << zc
          << "yt="   << yt
          << "zt="   << zt
          << "lbl="   << label
          << "\n";
      }
    }
  }
  return h;
}


//________________________________________________________
void AliTRDtrackingResolution::GetRefFigure(Int_t ifig)
{
  TAxis *ax = 0x0;
  TGraphErrors *g = 0x0;
  switch(ifig){
  case kClusterResidual:
    if(!(g = (TGraphErrors*)fGraphS->At(ifig))) break;
    g->Draw("apl");
    ax = g->GetHistogram()->GetYaxis();
    ax->SetRangeUser(-.5, 1.);
    ax->SetTitle("Clusters Y Residuals #sigma/#mu [mm]");
    ax = g->GetHistogram()->GetXaxis();
    ax->SetTitle("tg(#phi)");
    if(!(g = (TGraphErrors*)fGraphM->At(ifig))) break;
    g->Draw("pl");
    return;
  case kClusterResolution:
    if(!(g = (TGraphErrors*)fGraphS->At(ifig))) break;
    ax = g->GetHistogram()->GetYaxis();
    ax->SetRangeUser(-.5, 1.);
    ax->SetTitle("Cluster Y Resolution #sigma/#mu [mm]");
    ax = g->GetHistogram()->GetXaxis();
    ax->SetTitle("tg(#phi)");
    g->Draw("apl");
    if(!(g = (TGraphErrors*)fGraphM->At(ifig))) break;
    g->Draw("pl");
    return;
  case kTrackletYResolution:
    if(!(g = (TGraphErrors*)fGraphS->At(ifig))) break;
    ax = g->GetHistogram()->GetYaxis();
    ax->SetRangeUser(-.5, 1.);
    ax->SetTitle("Tracklet Y Resolution #sigma/#mu [mm]");
    ax = g->GetHistogram()->GetXaxis();
    ax->SetTitle("tg(#phi)");
    g->Draw("apl");
    if(!(g = (TGraphErrors*)fGraphM->At(ifig))) break;
    g->Draw("pl");
    return;
  case kTrackletZResolution:
    if(!(g = (TGraphErrors*)fGraphS->At(ifig))) break;
    ax = g->GetHistogram()->GetYaxis();
    ax->SetRangeUser(-.5, 1.);
    ax->SetTitle("Tracklet Z Resolution #sigma/#mu [mm]");
    ax = g->GetHistogram()->GetXaxis();
    ax->SetTitle("tg(#theta)");
    g->Draw("apl");
    if(!(g = (TGraphErrors*)fGraphM->At(ifig))) break;
    g->Draw("pl");
    return;
  case kTrackletAngleResolution:
    if(!(g = (TGraphErrors*)fGraphS->At(ifig))) break;
    ax = g->GetHistogram()->GetYaxis();
    ax->SetRangeUser(-.05, .2);
    ax->SetTitle("Tracklet Angular Resolution #sigma/#mu [deg]");
    ax = g->GetHistogram()->GetXaxis();
    ax->SetTitle("tg(#phi)");
    g->Draw("apl");
    if(!(g = (TGraphErrors*)fGraphM->At(ifig))) break;
    g->Draw("pl");
    return;
  default:
    AliInfo(Form("Reference plot [%d] not implemented yet", ifig));
    return;
  }
  AliInfo(Form("Reference plot [%d] missing result", ifig));
}


//________________________________________________________
Bool_t AliTRDtrackingResolution::PostProcess()
{
  //fContainer = dynamic_cast<TObjArray*>(GetOutputData(0));
  if (!fContainer) {
    Printf("ERROR: list not available");
    return kFALSE;
  }
  fNRefFigures = fContainer->GetEntriesFast();
  if(!fGraphS){ 
    fGraphS = new TObjArray(fNRefFigures);
    fGraphS->SetOwner();
  }
  if(!fGraphM){ 
    fGraphM = new TObjArray(fNRefFigures);
    fGraphM->SetOwner();
  }

  TH2I *h2 = 0x0;
  TH1D *h = 0x0;
  TGraphErrors *gm = 0x0, *gs = 0x0;

  // define models
  TF1 f("f1", "gaus", -.5, .5);  

  TF1 fb("fb", "[0]*exp(-0.5*((x-[1])/[2])**2)+[3]", -.5, .5);

  TF1 fc("fc", "[0]*exp(-0.5*((x-[1])/[2])**2)+[3]*exp(-0.5*((x-[4])/[5])**2)", -.5, .5);

  TCanvas *c = 0x0;
  if(IsVisual()) c = new TCanvas("c", Form("%s Visual", GetName()), 500, 500);
  char opt[5];
  sprintf(opt, "%c%c", IsVerbose() ? ' ' : 'Q', IsVisual() ? ' ': 'N');


  //PROCESS RESIDUAL DISTRIBUTIONS

  // Clusters residuals
  h2 = (TH2I *)(fContainer->At(kClusterResidual));
  gm = new TGraphErrors(h2->GetNbinsX());
  gm->SetLineColor(kBlue);
  gm->SetMarkerStyle(7);
  gm->SetMarkerColor(kBlue);
  gm->SetNameTitle("clm", "");
  fGraphM->AddAt(gm, kClusterResidual);
  gs = new TGraphErrors(h2->GetNbinsX());
  gs->SetLineColor(kRed);
  gs->SetMarkerStyle(23);
  gs->SetMarkerColor(kRed);
  gs->SetNameTitle("cls", "");
  fGraphS->AddAt(gs, kClusterResidual);
  for(Int_t ibin = 1; ibin <= h2->GetNbinsX(); ibin++){
    Double_t phi = h2->GetXaxis()->GetBinCenter(ibin);
    h = h2->ProjectionY("py", ibin, ibin);
    AdjustF1(h, &f);

    if(IsVisual()){c->cd(); c->SetLogy();}
    h->Fit(&f, opt, "", -0.5, 0.5);
    if(IsVisual()){c->Modified(); c->Update(); gSystem->Sleep(500);}
    
    gm->SetPoint(ibin - 1, phi, 10.*f.GetParameter(1));
    gm->SetPointError(ibin - 1, 0., 10.*f.GetParError(1));
    gs->SetPoint(ibin - 1, phi, 10.*f.GetParameter(2));
    gs->SetPointError(ibin - 1, 0., 10.*f.GetParError(2));
  }


  //PROCESS RESOLUTION DISTRIBUTIONS

  if(HasMCdata()){
    // cluster y resolution
    h2 = (TH2I*)fContainer->At(kClusterResolution);
    gm = new TGraphErrors(h2->GetNbinsX());
    gm->SetLineColor(kBlue);
    gm->SetMarkerStyle(7);
    gm->SetMarkerColor(kBlue);
    gm->SetNameTitle("clym", "");
    fGraphM->AddAt(gm, kClusterResolution);
    gs = new TGraphErrors(h2->GetNbinsX());
    gs->SetLineColor(kRed);
    gs->SetMarkerStyle(23);
    gs->SetMarkerColor(kRed);
    gs->SetNameTitle("clys", "");
    fGraphS->AddAt(gs, kClusterResolution);
    for(Int_t iphi=1; iphi<=h2->GetNbinsX(); iphi++){
      h = h2->ProjectionY("py", iphi, iphi);
      if(h->GetEntries()<100.) continue;
      AdjustF1(h, &f);

      if(IsVisual()){c->cd(); c->SetLogy();}
      h->Fit(&f, opt, "", -0.5, 0.5);
      if(IsVerbose()){
        printf("phi[%d] mean[%e] sigma[%e]\n\n", iphi, 10.*f.GetParameter(1), 10.*f.GetParameter(2));
      }
      if(IsVisual()){c->Modified(); c->Update(); gSystem->Sleep(500);}

      Double_t phi = h2->GetXaxis()->GetBinCenter(iphi);
      Int_t jphi = iphi -1;
      gm->SetPoint(jphi, phi, 10.*f.GetParameter(1));
      gm->SetPointError(jphi, 0., 10.*f.GetParError(1));
      gs->SetPoint(jphi, phi, 10.*f.GetParameter(2));
      gs->SetPointError(jphi, 0., 10.*f.GetParError(2));
    }
  
    // tracklet y resolution
    h2 = (TH2I*)fContainer->At(kTrackletYResolution);
    gm = new TGraphErrors(h2->GetNbinsX());
    gm->SetLineColor(kBlue);
    gm->SetMarkerStyle(7);
    gm->SetMarkerColor(kBlue);
    gm->SetNameTitle("trkltym", "");
    fGraphM->AddAt(gm, kTrackletYResolution);
    gs = new TGraphErrors(h2->GetNbinsX());
    gs->SetLineColor(kRed);
    gs->SetMarkerStyle(23);
    gs->SetMarkerColor(kRed);
    gs->SetNameTitle("trkltys", "");
    fGraphS->AddAt(gs, kTrackletYResolution);
    for(Int_t iphi=1; iphi<=h2->GetNbinsX(); iphi++){
      h = h2->ProjectionY("py", iphi, iphi);
      AdjustF1(h, &fb);

      if(IsVisual()){c->cd(); c->SetLogy();}
      h->Fit(&fb, opt, "", -0.5, 0.5);
      if(IsVisual()){c->Modified(); c->Update(); gSystem->Sleep(500);}

      Double_t phi = h2->GetXaxis()->GetBinCenter(iphi);
      Int_t jphi = iphi -1;
      gm->SetPoint(jphi, phi, 10.*fb.GetParameter(1));
      gm->SetPointError(jphi, 0., 10.*fb.GetParError(1));
      gs->SetPoint(jphi, phi, 10.*fb.GetParameter(2));
      gs->SetPointError(jphi, 0., 10.*fb.GetParError(2));
    }
  
    // tracklet z resolution
    h2 = (TH2I*)fContainer->At(kTrackletZResolution);
    gm = new TGraphErrors(h2->GetNbinsX());
    gm->SetLineColor(kBlue);
    gm->SetMarkerStyle(7);
    gm->SetMarkerColor(kBlue);
    gm->SetNameTitle("trkltzm", "");
    fGraphM->AddAt(gm, kTrackletZResolution);
    gs = new TGraphErrors(h2->GetNbinsX());
    gs->SetLineColor(kRed);
    gs->SetMarkerStyle(23);
    gs->SetMarkerColor(kRed);
    gs->SetNameTitle("trkltzs", "");
    fGraphS->AddAt(gs, kTrackletZResolution);
    for(Int_t iphi=1; iphi<=h2->GetNbinsX(); iphi++){
      h = h2->ProjectionY("py", iphi, iphi);
      AdjustF1(h, &fb);

      if(IsVisual()){c->cd(); c->SetLogy();}
      h->Fit(&fb, opt, "", -0.5, 0.5);
      if(IsVisual()){c->Modified(); c->Update(); gSystem->Sleep(500);}

      Double_t phi = h2->GetXaxis()->GetBinCenter(iphi);
      Int_t jphi = iphi -1;
      gm->SetPoint(jphi, phi, 10.*fb.GetParameter(1));
      gm->SetPointError(jphi, 0., 10.*fb.GetParError(1));
      gs->SetPoint(jphi, phi, 10.*fb.GetParameter(2));
      gs->SetPointError(jphi, 0., 10.*fb.GetParError(2));
    }
  
    // tracklet phi resolution
    h2 = (TH2I*)fContainer->At(kTrackletAngleResolution);
    gm = new TGraphErrors(h2->GetNbinsX());
    gm->SetLineColor(kBlue);
    gm->SetMarkerStyle(7);
    gm->SetMarkerColor(kBlue);
    gm->SetNameTitle("trkltym", "");
    fGraphM->AddAt(gm, kTrackletAngleResolution);
    gs = new TGraphErrors(h2->GetNbinsX());
    gs->SetLineColor(kRed);
    gs->SetMarkerStyle(23);
    gs->SetMarkerColor(kRed);
    gs->SetNameTitle("trkltys", "");
    fGraphS->AddAt(gs, kTrackletAngleResolution);
    for(Int_t iphi=1; iphi<=h2->GetNbinsX(); iphi++){
      h = h2->ProjectionY("py", iphi, iphi);

      if(IsVisual()){c->cd(); c->SetLogy();}
      h->Fit(&f, opt, "", -0.5, 0.5);
      if(IsVisual()){c->Modified(); c->Update(); gSystem->Sleep(500);}

      Double_t phi = h2->GetXaxis()->GetBinCenter(iphi);
      Int_t jphi = iphi -1;
      gm->SetPoint(jphi, phi, f.GetParameter(1));
      gm->SetPointError(jphi, 0., f.GetParError(1));
      gs->SetPoint(jphi, phi, f.GetParameter(2));
      gs->SetPointError(jphi, 0., f.GetParError(2));
    }
  }
  if(c) delete c;

  return kTRUE;
}


//________________________________________________________
void AliTRDtrackingResolution::Terminate(Option_t *)
{
  if(fDebugStream){ 
    delete fDebugStream;
    fDebugStream = 0x0;
    fDebugLevel = 0;
  }
  if(HasPostProcess()) PostProcess();
}

//________________________________________________________
void AliTRDtrackingResolution::AdjustF1(TH1 *h, TF1 *f)
{
// Helper function to avoid duplication of code
// Make first guesses on the fit parameters

  // find the intial parameters of the fit !! (thanks George)
  Int_t nbinsy = Int_t(.5*h->GetNbinsX());
  Double_t sum = 0.;
  for(Int_t jbin=nbinsy-4; jbin<=nbinsy+4; jbin++) sum+=h->GetBinContent(jbin); sum/=9.;
  f->SetParLimits(0, 0., 3.*sum);
  f->SetParameter(0, .9*sum);

  f->SetParLimits(1, -.2, .2);
  f->SetParameter(1, -0.1);

  f->SetParLimits(2, 0., 4.e-1);
  f->SetParameter(2, 2.e-2);
  if(f->GetNpar() <= 4) return;

  f->SetParLimits(3, 0., sum);
  f->SetParameter(3, .1*sum);

  f->SetParLimits(4, -.3, .3);
  f->SetParameter(4, 0.);

  f->SetParLimits(5, 0., 1.e2);
  f->SetParameter(5, 2.e-1);
}

//________________________________________________________
TObjArray* AliTRDtrackingResolution::Histos()
{
  if(fContainer) return fContainer;

  fContainer  = new TObjArray(5);

  TH1 *h = 0x0;
  // cluster to tracklet residuals [2]
  fContainer->AddAt(h = new TH2I("fYCl", "Clusters Residuals", 21, -.33, .33, 100, -.5, .5), kClusterResidual);
  h->GetXaxis()->SetTitle("tg(#phi)");
  h->GetYaxis()->SetTitle("#Delta y [cm]");
  h->GetZaxis()->SetTitle("entries");
//   // tracklet to Riemann fit residuals [2]
//   fContainer->AddAt(new TH2I("fYTrkltRRes", "Tracklet Riemann Residuals", 21, -21., 21., 100, -.5, .5), kTrackletRiemanYResidual);
//   fContainer->AddAt(new TH2I("fAngleTrkltRRes", "Tracklet Riemann Angular Residuals", 21, -21., 21., 100, -.5, .5), kTrackletRiemanAngleResidual);
//   fContainer->AddAt(new TH2I("fYTrkltKRes", "Tracklet Kalman Residuals", 21, -21., 21., 100, -.5, .5), kTrackletKalmanYResidual);
//   fContainer->AddAt(new TH2I("fAngleTrkltKRes", "Tracklet Kalman Angular Residuals", 21, -21., 21., 100, -.5, .5), kTrackletKalmanAngleResidual);

  // Resolution histos
  if(HasMCdata()){
    // cluster y resolution [0]
    fContainer->AddAt(h = new TH2I("fYClMC", "Cluster Resolution", 31, -.48, .48, 100, -.5, .5), kClusterResolution);
    h->GetXaxis()->SetTitle("tg(#phi)");
    h->GetYaxis()->SetTitle("#Delta y [cm]");
    h->GetZaxis()->SetTitle("entries");
    // tracklet y resolution [0]
    fContainer->AddAt(h = new TH2I("fYTrkltMC", "Tracklet Resolution (Y)", 31, -.48, .48, 100, -.5, .5), kTrackletYResolution);
    h->GetXaxis()->SetTitle("tg(#phi)");
    h->GetYaxis()->SetTitle("#Delta y [cm]");
    h->GetZaxis()->SetTitle("entries");
    // tracklet y resolution [0]
    fContainer->AddAt(h = new TH2I("fZTrkltMC", "Tracklet Resolution (Z)", 31, -.48, .48, 100, -.5, .5), kTrackletZResolution);
    h->GetXaxis()->SetTitle("tg(#theta)");
    h->GetYaxis()->SetTitle("#Delta z [cm]");
    h->GetZaxis()->SetTitle("entries");
    // tracklet angular resolution [1]
    fContainer->AddAt(h = new TH2I("fPhiTrkltMC", "Tracklet Resolution (Angular)", 31, -.48, .48, 100, -10., 10.), kTrackletAngleResolution);
    h->GetXaxis()->SetTitle("tg(#phi)");
    h->GetYaxis()->SetTitle("#Delta #phi [deg]");
    h->GetZaxis()->SetTitle("entries");

//     // Riemann track resolution [y, z, angular]
//     fContainer->AddAt(new TH2I("fYRT", "Track Riemann Y Resolution", 21, -21., 21., 100, -.5, .5), kTrackRYResolution);
//     fContainer->AddAt(new TH2I("fZRT", "Track Riemann Z Resolution", 21, -21., 21., 100, -.5, .5), kTrackRZResolution);
//     fContainer->AddAt(new TH2I("fPhiRT", "Track Riemann Angular Resolution", 21, -21., 21., 100, -10., 10.), kTrackRAngleResolution);
// 
//     Kalman track resolution [y, z, angular]
//     fContainer->AddAt(new TH2I("fYKT", "", 21, -21., 21., 100, -.5, .5), kTrackKYResolution);
//     fContainer->AddAt(new TH2I("fZKT", "", 21, -21., 21., 100, -.5, .5), kTrackKZResolution);
//     fContainer->AddAt(new TH2I("fPhiKT", "", 21, -21., 21., 100, -10., 10.), kTrackKAngleResolution);
  }
  return fContainer;
}


//________________________________________________________
void AliTRDtrackingResolution::SetRecoParam(AliTRDrecoParam *r)
{

  fReconstructor->SetRecoParam(r);
}
