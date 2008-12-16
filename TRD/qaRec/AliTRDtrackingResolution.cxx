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

#include <TROOT.h>
#include <TSystem.h>
#include <TObjArray.h>
#include <TH2.h>
#include <TH1.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TBox.h>
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

#include "AliTRDtrackInfo/AliTRDclusterInfo.h"
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
  ,fClResiduals(0x0)
  ,fTrkltResiduals(0x0)
  ,fTrkltPhiResiduals(0x0)
  ,fClResolution(0x0)
  ,fTrkltResolution(0x0)
{
  fReconstructor = new AliTRDReconstructor();
  fReconstructor->SetRecoParam(AliTRDrecoParam::GetLowFluxParam());
  fGeo = new AliTRDgeometry();

  InitFunctorList();

  DefineOutput(1+kClusterResidual, TObjArray::Class());
  DefineOutput(1+kTrackletYResidual, TObjArray::Class());
  DefineOutput(1+kTrackletPhiResidual, TObjArray::Class());
  DefineOutput(1+kClusterResolution, TObjArray::Class());
  DefineOutput(1+kTrackletYResolution, TObjArray::Class());
}

//________________________________________________________
AliTRDtrackingResolution::~AliTRDtrackingResolution()
{
  if(fGraphS){fGraphS->Delete(); delete fGraphS;}
  if(fGraphM){fGraphM->Delete(); delete fGraphM;}
  delete fGeo;
  delete fReconstructor;
  if(gGeoManager) delete gGeoManager;
  if(fClResiduals){fClResiduals->Delete(); delete fClResiduals;}
  if(fTrkltResiduals){fTrkltResiduals->Delete(); delete fTrkltResiduals;}
  if(fTrkltPhiResiduals){fTrkltPhiResiduals->Delete(); delete fTrkltPhiResiduals;}
  if(fClResolution){
    fClResolution->Delete(); 
    delete fClResolution;
  }
  if(fTrkltResolution){fTrkltResolution->Delete(); delete fTrkltResolution;}
}


//________________________________________________________
void AliTRDtrackingResolution::CreateOutputObjects()
{
  // spatial resolution
  OpenFile(0, "RECREATE");

  fContainer = Histos();

  fClResiduals = new TObjArray();
  fClResiduals->SetOwner(kTRUE);
  fTrkltResiduals = new TObjArray();
  fTrkltResiduals->SetOwner(kTRUE);
  fTrkltPhiResiduals = new TObjArray();
  fTrkltPhiResiduals->SetOwner(kTRUE);
  fClResolution = new TObjArray();
  fClResolution->SetOwner(kTRUE);
  fTrkltResolution = new TObjArray();
  fTrkltResolution->SetOwner(kTRUE);
}

//________________________________________________________
void AliTRDtrackingResolution::Exec(Option_t *opt)
{
  fClResiduals->Delete();
  fTrkltResiduals->Delete();
  fTrkltPhiResiduals->Delete();
  fClResolution->Delete();
  fTrkltResolution->Delete();

  AliTRDrecoTask::Exec(opt);

  PostData(1+kClusterResidual, fClResiduals);
  PostData(1+kTrackletYResidual, fTrkltResiduals);
  PostData(1+kTrackletPhiResidual, fTrkltPhiResiduals);
  PostData(1+kClusterResolution, fClResolution);
  PostData(1+kTrackletYResolution, fTrkltResolution);
}

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

  Float_t x0, y0, z0, dy, dydx, dzdx;
  AliTRDseedV1 *fTracklet = 0x0;  
  for(Int_t ily=0; ily<AliTRDgeometry::kNlayer; ily++){
    if(!(fTracklet = fTrack->GetTracklet(ily))) continue;
    if(!fTracklet->IsOK()) continue;
    x0 = fTracklet->GetX0();

    // retrive the track angle with the chamber
    y0   = fTracklet->GetYref(0);
    z0   = fTracklet->GetZref(0);
    dydx = fTracklet->GetYref(1);
    dzdx = fTracklet->GetZref(1);
    Float_t tilt = fTracklet->GetTilt();
    AliTRDcluster *c = 0x0;
    fTracklet->ResetClusterIter(kFALSE);
    while((c = fTracklet->PrevCluster())){
      Float_t xc = c->GetX();
      Float_t yc = c->GetY();
      Float_t zc = c->GetZ();
      Float_t dx = x0 - xc; 
      Float_t yt = y0 - dx*dydx;
      Float_t zt = z0 - dx*dzdx; 
      dy = yt - (yc - tilt*(zc-zt));

      //dy = trklt.GetYat(c->GetX()) - c->GetY();
      h->Fill(dydx, dy);
  
      if(fDebugLevel>=1){
        // Get z-position with respect to anode wire
        AliTRDSimParam    *simParam    = AliTRDSimParam::Instance();
        Int_t istk = fGeo->GetStack(c->GetDetector());
        AliTRDpadPlane *pp = fGeo->GetPadPlane(ily, istk);
        Float_t row0 = pp->GetRow0();
        Float_t d  =  row0 - zt + simParam->GetAnodeWireOffset();
        d -= ((Int_t)(2 * d)) / 2.0;
        if (d > 0.25) d  = 0.5 - d;

/*        AliTRDclusterInfo *clInfo = new AliTRDclusterInfo;
        fClResiduals->Add(clInfo);
        clInfo->SetCluster(c);
        clInfo->SetGlobalPosition(yt, zt, dydx, dzdx);
        clInfo->SetResolution(dy);
        clInfo->SetAnisochronity(d);
        clInfo->SetDriftLength(dx);
        (*fDebugStream) << "ClusterResiduals"
          <<"clInfo.=" << clInfo
          << "\n";*/
      }
    }
  }
  return h;
}


//________________________________________________________
TH1* AliTRDtrackingResolution::PlotTrackletResiduals(const AliTRDtrackV1 *track)
{
  if(track) fTrack = track;
  if(!fTrack){
    AliWarning("No Track defined.");
    return 0x0;
  }
  TH1 *h = 0x0;
  if(!(h = ((TH2I*)fContainer->At(kTrackletYResidual)))){
    AliWarning("No output histogram defined.");
    return 0x0;
  }

  AliTRDseedV1 *tracklet = 0x0;
  for(Int_t il=0; il<AliTRDgeometry::kNlayer; il++){
    if(!(tracklet = fTrack->GetTracklet(il))) continue;
    if(!tracklet->IsOK()) continue;
    h->Fill(tracklet->GetYref(1), tracklet->GetYref(0)-tracklet->GetYfit(0));
  }
  return h;

  // refit the track
//   AliRieman fRim(fTrack->GetNumberOfClusters());
//   Float_t x[AliTRDgeometry::kNlayer] = {-1., -1., -1., -1., -1., -1.}, y[AliTRDgeometry::kNlayer], dydx[AliTRDgeometry::kNlayer];

//     AliTRDcluster *c = 0x0;
//     tracklet->ResetClusterIter(kFALSE);
//     while((c = tracklet->PrevCluster())){
//       Float_t xc = c->GetX();
//       Float_t yc = c->GetY();
//       Float_t zc = c->GetZ();
//       Float_t zt = tracklet->GetZref(0) - (tracklet->GetX0()-xc)*tracklet->GetZref(1); 
//       yc -= tracklet->GetTilt()*(zc-zt);
//       fRim.AddPoint(xc, yc, zc, 1, 10);
//     }
//     tracklet->Fit(kTRUE);
// 
//     x[il] = tracklet->GetX0();
//     y[il] = tracklet->GetYfit(0)-tracklet->GetYfit(1)*(tracklet->GetX0()-x[il]);
//     dydx[il] = tracklet->GetYref(1);


//   fRim.Update();

//   for(Int_t il=0; il<AliTRDgeometry::kNlayer; il++){
//     if(x[il] < 0.) continue;
//     Float_t dy = y[il]-fRim.GetYat(x[il])/*/sigma_track*/;
//     h->Fill(dydx[il], dy);
// 
//     if(fDebugLevel>=1){
//       Double_t sigmay = fRim.GetErrY(x[il]);
//       Float_t yt = fRim.GetYat(x[il]);
//       (*fDebugStream) << "TrkltResiduals"
//         << "layer="  << il
//         << "x="      << x[il]
//         << "y="      << y[il]
//         << "yt="     << yt
//         << "dydx="   << dydx[il]
//         << "dy="     << dy
//         << "sigmay=" << sigmay
//         << "\n";
//     }
//   }
}

//________________________________________________________
TH1* AliTRDtrackingResolution::PlotTrackletPhiResiduals(const AliTRDtrackV1 *track)
{
  if(track) fTrack = track;
  if(!fTrack){
    AliWarning("No Track defined.");
    return 0x0;
  }
  TH1 *h = 0x0;
  if(!(h = ((TH2I*)fContainer->At(kTrackletPhiResidual)))){
    AliWarning("No output histogram defined.");
    return 0x0;
  }

  AliTRDseedV1 *tracklet = 0x0;
  for(Int_t il=0; il<AliTRDgeometry::kNlayer; il++){
    if(!(tracklet = fTrack->GetTracklet(il))) continue;
    if(!tracklet->IsOK()) continue;
    h->Fill(tracklet->GetYref(1), TMath::ATan(tracklet->GetYref(1))-TMath::ATan(tracklet->GetYfit(1)));
  }
  return h;

//   // refit the track
//   AliRieman fRim(fTrack->GetNumberOfClusters());
//   Float_t x[AliTRDgeometry::kNlayer] = {-1., -1., -1., -1., -1., -1.}, y[AliTRDgeometry::kNlayer], dydx[AliTRDgeometry::kNlayer];
//   Float_t dydx_ref=0, dydx_fit=0, phiref=0, phifit=0, phidiff=0;
//   AliTRDseedV1 *tracklet = 0x0;
//   for(Int_t il=0; il<AliTRDgeometry::kNlayer; il++){
//     if(!(tracklet = fTrack->GetTracklet(il))) continue;
//     if(!tracklet->IsOK()) continue;
//     AliTRDcluster *c = 0x0;
//     tracklet->ResetClusterIter(kFALSE);
//     while((c = tracklet->PrevCluster())){
//       Float_t xc = c->GetX();
//       Float_t yc = c->GetY();
//       Float_t zc = c->GetZ();
//       Float_t zt = tracklet->GetZref(0) - (tracklet->GetX0()-xc)*tracklet->GetZref(1); 
//       yc -= tracklet->GetTilt()*(zc-zt);
//       fRim.AddPoint(xc, yc, zc, 1, 10);
//     }
//     tracklet->Fit(kTRUE);
// 
//     x[il] = tracklet->GetX0();
//     y[il] = tracklet->GetYfit(0)-tracklet->GetYfit(1)*(tracklet->GetX0()-x[il]);
//     dydx[il] = tracklet->GetYref(1);
//   }
//   fRim.Update();
// 
//   for(Int_t il=0; il<AliTRDgeometry::kNlayer; il++){
//     if(x[il] < 0.) continue;
//     if(!(tracklet = fTrack->GetTracklet(il))) continue;
//     if(!tracklet->IsOK()) continue;
// 
//     dydx_ref = fRim.GetDYat(x[il]); 
//     dydx_fit = tracklet->GetYfit(1);
// 
//     phiref = TMath::ATan(dydx_ref);//*TMath::RadToDeg();
//     phifit = TMath::ATan(dydx_fit);//*TMath::RadToDeg();
//     
//     phidiff = phiref-phifit; /*/sigma_phi*/;
// 
//     h->Fill(dydx_ref, phidiff);
// 
// 
//     if(fDebugLevel>=1){
//       (*fDebugStream) << "TrkltPhiResiduals"
//         << "layer="  << il
//         << "dydx_fit="   << dydx_fit
//         << "dydx_ref="   << dydx_ref
//         << "phiref="     << phiref
//         << "phifit="     << phifit
//         << "phidiff="     << phidiff
//         << "\n";
//     }
//   }
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
  Int_t label = fMC->GetLabel();
  Float_t x0, y0, z0, dx, dy, dydx, dzdx;
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
    if(!tt.Fit(kTRUE)) continue;
    //tt.Update();

    dx = 0.;//x0 - tt.GetXref();
    Float_t yt = y0 - dx*dydx;
    Float_t yf = tt.GetYfit(0) - (fTracklet->GetX0() - x0)*tt.GetYfit(1);
    dy = yf-yt;
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
      dx = x0 - xc; 
      Float_t yt = y0 - dx*dydx;
      Float_t zt = z0 - dx*dzdx; 
      dy = yt - (yc - tilt*(zc-zt));

      // Fill Histograms
      if(q>20. && q<250.) ((TH2I*)fContainer->At(kClusterResolution))->Fill(dydx, dy);
      
      // Fill Debug Tree
      if(fDebugLevel>=1){
        Float_t d = zr0 - zt;
        d -= ((Int_t)(2 * d)) / 2.0;
        if (d > 0.25) d  = 0.5 - d;

        AliTRDclusterInfo *clInfo = new AliTRDclusterInfo;
        fClResolution->Add(clInfo);
        clInfo->SetCluster(c);
        clInfo->SetMC(pdg, label);
        clInfo->SetGlobalPosition(yt, zt, dydx, dzdx);
        clInfo->SetResolution(dy);
        clInfo->SetAnisochronity(d);
        clInfo->SetDriftLength(dx);
        clInfo->SetTilt(tilt);
        //clInfo->Print();
        (*fDebugStream) << "ClusterResolution"
          <<"clInfo.=" << clInfo
          << "\n";
      }
    }
  }
  return h;
}


//________________________________________________________
void AliTRDtrackingResolution::GetRefFigure(Int_t ifig)
{
  TBox *b = 0x0;
  TAxis *ax = 0x0;
  TGraphErrors *g = 0x0;
  switch(ifig){
  case kClusterResidual:
    if(!(g = (TGraphErrors*)fGraphS->At(ifig))) break;
    g->Draw("apl");
    ax = g->GetHistogram()->GetYaxis();
    ax->SetRangeUser(-.5, 2.5);
    ax->SetTitle("Clusters Y Residuals #sigma/#mu [mm]");
    ax = g->GetHistogram()->GetXaxis();
    ax->SetTitle("tg(#phi)");
    if(!(g = (TGraphErrors*)fGraphM->At(ifig))) break;
    g->Draw("pl");
    b = new TBox(-.15, -.5, .15, 1.);
    b->SetFillStyle(3002);b->SetFillColor(kGreen);
    b->SetLineColor(0); b->Draw();
    return;
  case kTrackletYResidual:
    if(!(g = (TGraphErrors*)fGraphS->At(ifig))) break;
    g->Draw("apl");
    ax = g->GetHistogram()->GetYaxis();
    ax->SetRangeUser(-.5, 3.);
    ax->SetTitle("Tracklet Y Residuals #sigma/#mu [mm]");
    ax = g->GetHistogram()->GetXaxis();
    ax->SetTitle("tg(#phi)");
    if(!(g = (TGraphErrors*)fGraphM->At(ifig))) break;
    g->Draw("pl");
    b = new TBox(-.15, -.5, .15, 1.);
    b->SetFillStyle(3002);b->SetFillColor(kGreen);
    b->SetLineColor(0); b->Draw();
    return;
  case kTrackletPhiResidual:
    if(!(g = (TGraphErrors*)fGraphS->At(ifig))) break;
    g->Draw("apl");
    ax = g->GetHistogram()->GetYaxis();
    ax->SetRangeUser(-.5, 2.);
    ax->SetTitle("Tracklet Phi Residuals #sigma/#mu [rad]");
    ax = g->GetHistogram()->GetXaxis();
    ax->SetTitle("tg(#phi)");
    if(!(g = (TGraphErrors*)fGraphM->At(ifig))) break;
    g->Draw("pl");
    b = new TBox(-.15, -.5, .15, 1.);
    b->SetFillStyle(3002);b->SetFillColor(kGreen);
    b->SetLineColor(0); b->Draw();
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
    b = new TBox(-.15, -.5, .15, 1.);
    b->SetFillStyle(3002);b->SetFillColor(kGreen);
    b->SetLineColor(0); b->Draw();
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
  gm = new TGraphErrors();
  gm->SetLineColor(kBlue);
  gm->SetMarkerStyle(7);
  gm->SetMarkerColor(kBlue);
  gm->SetNameTitle("clm", "");
  fGraphM->AddAt(gm, kClusterResidual);
  gs = new TGraphErrors();
  gs->SetLineColor(kRed);
  gs->SetMarkerStyle(23);
  gs->SetMarkerColor(kRed);
  gs->SetNameTitle("cls", "");
  fGraphS->AddAt(gs, kClusterResidual);
  for(Int_t ibin = 1; ibin <= h2->GetNbinsX(); ibin++){
    Double_t phi = h2->GetXaxis()->GetBinCenter(ibin);
    h = h2->ProjectionY("py", ibin, ibin);
    if(h->GetEntries()<100) continue;
    AdjustF1(h, &f);

    if(IsVisual()){c->cd(); c->SetLogy();}
    h->Fit(&f, opt, "", -0.5, 0.5);
    if(IsVisual()){c->Modified(); c->Update(); gSystem->Sleep(500);}
    
    Int_t ip = gm->GetN();
    gm->SetPoint(ip, phi, 10.*f.GetParameter(1));
    gm->SetPointError(ip, 0., 10.*f.GetParError(1));
    gs->SetPoint(ip, phi, 10.*f.GetParameter(2));
    gs->SetPointError(ip, 0., 10.*f.GetParError(2));
  }

  // Tracklet y residuals
  h2 = (TH2I *)(fContainer->At(kTrackletYResidual));
  gm = new TGraphErrors();
  gm->SetLineColor(kBlue);
  gm->SetMarkerStyle(7);
  gm->SetMarkerColor(kBlue);
  gm->SetNameTitle("tktm", "");
  fGraphM->AddAt(gm, kTrackletYResidual);
  gs = new TGraphErrors();
  gs->SetLineColor(kRed);
  gs->SetMarkerStyle(23);
  gs->SetMarkerColor(kRed);
  gs->SetNameTitle("tkts", "");
  fGraphS->AddAt(gs, kTrackletYResidual);
  for(Int_t ibin = 1; ibin <= h2->GetNbinsX(); ibin++){
    Double_t phi = h2->GetXaxis()->GetBinCenter(ibin);
    h = h2->ProjectionY("py", ibin, ibin);
    if(h->GetEntries()<100) continue;
    AdjustF1(h, &f);

    if(IsVisual()){c->cd(); c->SetLogy();}
    h->Fit(&f, opt, "", -0.5, 0.5);
    if(IsVisual()){c->Modified(); c->Update(); gSystem->Sleep(500);}
    
    Int_t ip = gm->GetN();
    gm->SetPoint(ip, phi, 10.*f.GetParameter(1));
    gm->SetPointError(ip, 0., 10.*f.GetParError(1));
    gs->SetPoint(ip, phi, 10.*f.GetParameter(2));
    gs->SetPointError(ip, 0., 10.*f.GetParError(2));
  }

  // Tracklet phi residuals
  h2 = (TH2I *)(fContainer->At(kTrackletPhiResidual));
  gm = new TGraphErrors();
  gm->SetLineColor(kBlue);
  gm->SetMarkerStyle(7);
  gm->SetMarkerColor(kBlue);
  gm->SetNameTitle("tktphim", "");
  fGraphM->AddAt(gm, kTrackletPhiResidual);
  gs = new TGraphErrors();
  gs->SetLineColor(kRed);
  gs->SetMarkerStyle(23);
  gs->SetMarkerColor(kRed);
  gs->SetNameTitle("tktphis", "");
  fGraphS->AddAt(gs, kTrackletPhiResidual);
  for(Int_t ibin = 1; ibin <= h2->GetNbinsX(); ibin++){
    Double_t phi = h2->GetXaxis()->GetBinCenter(ibin);
    h = h2->ProjectionY("py", ibin, ibin);
    if(h->GetEntries()<100) continue;
    AdjustF1(h, &f);

    if(IsVisual()){c->cd(); c->SetLogy();}
    h->Fit(&f, opt, "", -0.5, 0.5);
    if(IsVisual()){c->Modified(); c->Update(); gSystem->Sleep(500);}
    
    Int_t ip = gm->GetN();
    gm->SetPoint(ip, phi, 10.*f.GetParameter(1));
    gm->SetPointError(ip, 0., 10.*f.GetParError(1));
    gs->SetPoint(ip, phi, 10.*f.GetParameter(2));
    gs->SetPointError(ip, 0., 10.*f.GetParError(2));
  }


  //PROCESS RESOLUTION DISTRIBUTIONS

  if(HasMCdata()){
    // cluster y resolution
    h2 = (TH2I*)fContainer->At(kClusterResolution);
    gm = new TGraphErrors();
    gm->SetLineColor(kBlue);
    gm->SetMarkerStyle(7);
    gm->SetMarkerColor(kBlue);
    gm->SetNameTitle("clym", "");
    fGraphM->AddAt(gm, kClusterResolution);
    gs = new TGraphErrors();
    gs->SetLineColor(kRed);
    gs->SetMarkerStyle(23);
    gs->SetMarkerColor(kRed);
    gs->SetNameTitle("clys", "");
    fGraphS->AddAt(gs, kClusterResolution);
    for(Int_t iphi=1; iphi<=h2->GetNbinsX(); iphi++){
      h = h2->ProjectionY("py", iphi, iphi);
      if(h->GetEntries()<100) continue;
      AdjustF1(h, &f);

      if(IsVisual()){c->cd(); c->SetLogy();}
      h->Fit(&f, opt, "", -0.5, 0.5);
      if(IsVerbose()){
        printf("phi[%d] mean[%e] sigma[%e]\n\n", iphi, 10.*f.GetParameter(1), 10.*f.GetParameter(2));
      }
      if(IsVisual()){c->Modified(); c->Update(); gSystem->Sleep(500);}

      Double_t phi = h2->GetXaxis()->GetBinCenter(iphi);
      Int_t ip = gm->GetN();
      gm->SetPoint(ip, phi, 10.*f.GetParameter(1));
      gm->SetPointError(ip, 0., 10.*f.GetParError(1));
      gs->SetPoint(ip, phi, 10.*f.GetParameter(2));
      gs->SetPointError(ip, 0., 10.*f.GetParError(2));
    }
  
    // tracklet y resolution
    h2 = (TH2I*)fContainer->At(kTrackletYResolution);
    gm = new TGraphErrors();
    gm->SetLineColor(kBlue);
    gm->SetMarkerStyle(7);
    gm->SetMarkerColor(kBlue);
    gm->SetNameTitle("trkltym", "");
    fGraphM->AddAt(gm, kTrackletYResolution);
    gs = new TGraphErrors();
    gs->SetLineColor(kRed);
    gs->SetMarkerStyle(23);
    gs->SetMarkerColor(kRed);
    gs->SetNameTitle("trkltys", "");
    fGraphS->AddAt(gs, kTrackletYResolution);
    for(Int_t iphi=1; iphi<=h2->GetNbinsX(); iphi++){
      h = h2->ProjectionY("py", iphi, iphi);
      if(h->GetEntries()<100) continue;
      AdjustF1(h, &f);

      if(IsVisual()){c->cd(); c->SetLogy();}
      h->Fit(&f, opt, "", -0.5, 0.5);
      if(IsVisual()){c->Modified(); c->Update(); gSystem->Sleep(500);}

      Double_t phi = h2->GetXaxis()->GetBinCenter(iphi);
      Int_t ip = gm->GetN();
      gm->SetPoint(ip, phi, 10.*f.GetParameter(1));
      gm->SetPointError(ip, 0., 10.*f.GetParError(1));
      gs->SetPoint(ip, phi, 10.*f.GetParameter(2));
      gs->SetPointError(ip, 0., 10.*f.GetParError(2));
    }
  
    // tracklet z resolution
    h2 = (TH2I*)fContainer->At(kTrackletZResolution);
    gm = new TGraphErrors();
    gm->SetLineColor(kBlue);
    gm->SetMarkerStyle(7);
    gm->SetMarkerColor(kBlue);
    gm->SetNameTitle("trkltzm", "");
    fGraphM->AddAt(gm, kTrackletZResolution);
    gs = new TGraphErrors();
    gs->SetLineColor(kRed);
    gs->SetMarkerStyle(23);
    gs->SetMarkerColor(kRed);
    gs->SetNameTitle("trkltzs", "");
    fGraphS->AddAt(gs, kTrackletZResolution);
    for(Int_t iphi=1; iphi<=h2->GetNbinsX(); iphi++){
      h = h2->ProjectionY("py", iphi, iphi);
      if(h->GetEntries()<100) continue;
      AdjustF1(h, &fb);

      if(IsVisual()){c->cd(); c->SetLogy();}
      h->Fit(&fb, opt, "", -0.5, 0.5);
      if(IsVisual()){c->Modified(); c->Update(); gSystem->Sleep(500);}

      Double_t phi = h2->GetXaxis()->GetBinCenter(iphi);
      Int_t ip = gm->GetN();
      gm->SetPoint(ip, phi, 10.*fb.GetParameter(1));
      gm->SetPointError(ip, 0., 10.*fb.GetParError(1));
      gs->SetPoint(ip, phi, 10.*fb.GetParameter(2));
      gs->SetPointError(ip, 0., 10.*fb.GetParError(2));
    }
  
    //tracklet phi resolution
    h2 = (TH2I*)fContainer->At(kTrackletAngleResolution);
    gm = new TGraphErrors();
    gm->SetLineColor(kBlue);
    gm->SetMarkerStyle(7);
    gm->SetMarkerColor(kBlue);
    gm->SetNameTitle("trkltym", "");
    fGraphM->AddAt(gm, kTrackletAngleResolution);
    gs = new TGraphErrors();
    gs->SetLineColor(kRed);
    gs->SetMarkerStyle(23);
    gs->SetMarkerColor(kRed);
    gs->SetNameTitle("trkltys", "");
    fGraphS->AddAt(gs, kTrackletAngleResolution);
    for(Int_t iphi=1; iphi<=h2->GetNbinsX(); iphi++){
      h = h2->ProjectionY("py", iphi, iphi);
      if(h->GetEntries()<100) continue;

      if(IsVisual()){c->cd(); c->SetLogy();}
      h->Fit(&f, opt, "", -0.5, 0.5);
      if(IsVisual()){c->Modified(); c->Update(); gSystem->Sleep(500);}

      Double_t phi = h2->GetXaxis()->GetBinCenter(iphi);
      Int_t ip = gm->GetN();
      gm->SetPoint(ip, phi, f.GetParameter(1));
      gm->SetPointError(ip, 0., f.GetParError(1));
      gs->SetPoint(ip, phi, f.GetParameter(2));
      gs->SetPointError(ip, 0., f.GetParError(2));
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

  fContainer  = new TObjArray(7);
  //fContainer->SetOwner(kTRUE);

  TH1 *h = 0x0;
  // cluster to tracklet residuals [2]
  if(!(h = (TH2I*)gROOT->FindObject("fYCl"))){
    h = new TH2I("fYCl", "Clusters Residuals", 21, -.33, .33, 100, -.5, .5);
    h->GetXaxis()->SetTitle("tg(#phi)");
    h->GetYaxis()->SetTitle("#Delta y [cm]");
    h->GetZaxis()->SetTitle("entries");
  } else h->Reset();
  fContainer->AddAt(h, kClusterResidual);

  // tracklet to track residuals [2]
  if(!(h = (TH2I*)gROOT->FindObject("hTrkltYRez"))){ 
    h = new TH2I("hTrkltYRez", "Tracklets", 21, -.33, .33, 100, -.5, .5);
    h->GetXaxis()->SetTitle("tg(#phi)");
    h->GetYaxis()->SetTitle("#Delta y [cm]");
    h->GetZaxis()->SetTitle("entries");
  } else h->Reset();
  fContainer->AddAt(h, kTrackletYResidual);

  // tracklet to track residuals angular [2]
  if(!(h = (TH2I*)gROOT->FindObject("hTrkltPhiRez"))){ 
    h = new TH2I("hTrkltPhiRez", "Tracklets", 21, -.33, .33, 100, -.5, .5);
    h->GetXaxis()->SetTitle("tg(#phi)");
    h->GetYaxis()->SetTitle("#Delta phi [#circ]");
    h->GetZaxis()->SetTitle("entries");
  } else h->Reset();
  fContainer->AddAt(h, kTrackletPhiResidual);


  // Resolution histos
  if(HasMCdata()){
    // cluster y resolution [0]
    if(!(h = (TH2I*)gROOT->FindObject("fYClMC"))){
      h = new TH2I("fYClMC", "Cluster Resolution", 31, -.48, .48, 100, -.5, .5);
      h->GetXaxis()->SetTitle("tg(#phi)");
      h->GetYaxis()->SetTitle("#Delta y [cm]");
      h->GetZaxis()->SetTitle("entries");
    } else h->Reset();
    fContainer->AddAt(h, kClusterResolution);

    // tracklet y resolution [0]
    if(!(h = (TH2I*)gROOT->FindObject("fYTrkltMC"))){
      h = new TH2I("fYTrkltMC", "Tracklet Resolution (Y)", 31, -.48, .48, 100, -.5, .5);
      h->GetXaxis()->SetTitle("tg(#phi)");
      h->GetYaxis()->SetTitle("#Delta y [cm]");
      h->GetZaxis()->SetTitle("entries");
    } else h->Reset();
    fContainer->AddAt(h, kTrackletYResolution);

    // tracklet y resolution [0]
    if(!(h = (TH2I*)gROOT->FindObject("fZTrkltMC"))){
      h = new TH2I("fZTrkltMC", "Tracklet Resolution (Z)", 31, -.48, .48, 100, -.5, .5);
      h->GetXaxis()->SetTitle("tg(#theta)");
      h->GetYaxis()->SetTitle("#Delta z [cm]");
      h->GetZaxis()->SetTitle("entries");
    } else h->Reset();
    fContainer->AddAt(h, kTrackletZResolution);

    // tracklet angular resolution [1]
    if(!(h = (TH2I*)gROOT->FindObject("fPhiTrkltMC"))){
      h = new TH2I("fPhiTrkltMC", "Tracklet Resolution (Angular)", 31, -.48, .48, 100, -10., 10.);
      h->GetXaxis()->SetTitle("tg(#phi)");
      h->GetYaxis()->SetTitle("#Delta #phi [deg]");
      h->GetZaxis()->SetTitle("entries");
    } else h->Reset();
    fContainer->AddAt(h, kTrackletAngleResolution);
  }
  return fContainer;
}


//________________________________________________________
void AliTRDtrackingResolution::SetRecoParam(AliTRDrecoParam *r)
{

  fReconstructor->SetRecoParam(r);
}
