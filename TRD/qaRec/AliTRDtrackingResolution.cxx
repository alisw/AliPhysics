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
#include <TPDGCode.h>
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

  DefineOutput(1+kCluster, TObjArray::Class());
  DefineOutput(1+kTrackletY, TObjArray::Class());
  DefineOutput(1+kTrackletPhi, TObjArray::Class());
  DefineOutput(1+kMCcluster, TObjArray::Class());
  DefineOutput(1+kMCtrackletY, TObjArray::Class());
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

  PostData(1+kCluster, fClResiduals);
  PostData(1+kTrackletY, fTrkltResiduals);
  PostData(1+kTrackletPhi, fTrkltPhiResiduals);
  PostData(1+kMCcluster, fClResolution);
  PostData(1+kMCtrackletY, fTrkltResolution);
}

//________________________________________________________
TH1* AliTRDtrackingResolution::PlotCluster(const AliTRDtrackV1 *track)
{
  if(track) fTrack = track;
  if(!fTrack){
    AliWarning("No Track defined.");
    return 0x0;
  }
  TH1 *h = 0x0;
  if(!(h = ((TH2I*)fContainer->At(kCluster)))){
    AliWarning("No output histogram defined.");
    return 0x0;
  }

  Double_t cov[3];
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
    fTracklet->GetCovRef(cov);
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
      yc -= tilt*(zc-zt); // tilt correction
      dy = yt - yc;

      h->Fill(dydx, dy/TMath::Sqrt(cov[0] + c->GetSigmaY2()));
  
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
TH1* AliTRDtrackingResolution::PlotTracklet(const AliTRDtrackV1 *track)
{
  if(track) fTrack = track;
  if(!fTrack){
    AliWarning("No Track defined.");
    return 0x0;
  }
  TH1 *h = 0x0;
  if(!(h = ((TH2I*)fContainer->At(kTrackletY)))){
    AliWarning("No output histogram defined.");
    return 0x0;
  }

  Double_t cov[3], covR[3];
  Float_t xref, y0, dx, dy, dydx;
  AliTRDseedV1 *fTracklet = 0x0;  
  for(Int_t il=AliTRDgeometry::kNlayer; il--;){
    if(!(fTracklet = fTrack->GetTracklet(il))) continue;
    if(!fTracklet->IsOK()) continue;
    y0   = fTracklet->GetYref(0);
    dydx = fTracklet->GetYref(1);
    xref = fTracklet->GetXref();
    dx   = fTracklet->GetX0() - xref;
    dy   = y0-dx*dydx - fTracklet->GetYat(xref);
    fTracklet->GetCovAt(xref, cov);
    fTracklet->GetCovRef(covR);
    h->Fill(dydx, dy/TMath::Sqrt(cov[0] + covR[0]));
  }
  return h;
}

//________________________________________________________
TH1* AliTRDtrackingResolution::PlotTrackletPhi(const AliTRDtrackV1 *track)
{
  if(track) fTrack = track;
  if(!fTrack){
    AliWarning("No Track defined.");
    return 0x0;
  }
  TH1 *h = 0x0;
  if(!(h = ((TH2I*)fContainer->At(kTrackletPhi)))){
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
  UChar_t s;
  Int_t pdg = fMC->GetPDG(), det=-1;
  Int_t label = fMC->GetLabel();
  Float_t p, pt, x0, y0, z0, dx, dy, dz, dydx, dzdx;

  if(fDebugLevel>=1){
    Double_t DX[12], DY[12], DZ[12], DPt[12], COV[12][15];
    fMC->PropagateKalman(DX, DY, DZ, DPt, COV);
    (*fDebugStream) << "MCkalman"
      << "pdg="  << pdg
      << "dx0="  << DX[0]
      << "dx1="  << DX[1]
      << "dx2="  << DX[2]
      << "dy0="  << DY[0]
      << "dy1="  << DY[1]
      << "dy2="  << DY[2]
      << "dz0="  << DZ[0]
      << "dz1="  << DZ[1]
      << "dz2="  << DZ[2]
      << "dpt0=" << DPt[0]
      << "dpt1=" << DPt[1]
      << "dpt2=" << DPt[2]
      << "\n";
  }

  AliTRDseedV1 *fTracklet = 0x0;  
  for(Int_t ily=0; ily<AliTRDgeometry::kNlayer; ily++){
    if(!(fTracklet = fTrack->GetTracklet(ily)) ||
       !fTracklet->IsOK()) continue;

    det = fTracklet->GetDetector();
    x0  = fTracklet->GetX0();
    //radial shift with respect to the MC reference (radial position of the pad plane)
    dx  = x0 - fTracklet->GetXref();
    if(!fMC->GetDirections(x0, y0, z0, dydx, dzdx, pt, s)) continue; 
    // MC track position at reference radial position
    Float_t yt = y0 - dx*dydx;
    Float_t zt = z0 - dx*dzdx;
    p = pt*(1.+dzdx*dzdx); // pt -> p

    // add Kalman residuals for y, z and pt
    Float_t dxr= fTracklet->GetX0() - x0 + dx; 
    Float_t yr = fTracklet->GetYref(0) - dxr*fTracklet->GetYref(1);
    dy = yt - yr;
    Float_t zr = fTracklet->GetZref(0) - dxr*fTracklet->GetZref(1);
    dz = zt - zr;
    Float_t tgl = fTracklet->GetTgl();
    Float_t ptr = fTracklet->GetMomentum()/(1.+tgl*tgl);

    ((TH2I*)fContainer->At(kMCtrackY))->Fill(dydx, dy);
    ((TH2I*)fContainer->At(kMCtrackZ))->Fill(dzdx, dz);
    if(pdg!=kElectron && pdg!=kPositron) ((TH2I*)fContainer->At(kMCtrackPt))->Fill(1./pt, ptr-pt);
    // Fill Debug stream for Kalman track
    if(fDebugLevel>=1){
      Float_t dydxr = fTracklet->GetYref(1);
      (*fDebugStream) << "MCtrack"
        << "det="     << det
        << "pdg="     << pdg
        << "pt="      << pt
        << "yt="      << yt
        << "zt="      << zt
        << "dydx="    << dydx
        << "dzdx="    << dzdx
        << "ptr="     << ptr
        << "dy="      << dy
        << "dz="      << dz
        << "dydxr="   << dydxr
        << "dzdxr="   << tgl
        << "\n";
    }

    // recalculate tracklet based on the MC info
    AliTRDseedV1 tt(*fTracklet);
    tt.SetZref(0, z0);
    tt.SetZref(1, dzdx); 
    if(!tt.Fit(kTRUE)) continue;

    // add tracklet residuals for y and dydx
    Float_t yf = tt.GetYfit(0) - dxr*tt.GetYfit(1);
    dy = yt - yf;
    Float_t dphi   = (tt.GetYfit(1) - dydx);
    dphi /= 1.- tt.GetYfit(1)*dydx;
    ((TH2I*)fContainer->At(kMCtrackletY))->Fill(dydx, dy);
    ((TH2I*)fContainer->At(kMCtrackletPhi))->Fill(dydx, dphi*TMath::RadToDeg());

    Float_t dz = 100.;
    Bool_t rc = fTracklet->IsRowCross(); 
    if(rc){
      // add tracklet residuals for z
      Double_t *xyz = tt.GetCrossXYZ();
      dz = xyz[2] - (z0 - (x0 - xyz[0])*dzdx) ;
      ((TH2I*)fContainer->At(kMCtrackletZ))->Fill(dzdx, dz);
    }
  
    // Fill Debug stream for tracklet
    if(fDebugLevel>=1){
      (*fDebugStream) << "MCtracklet"
        << "det="     << det
        << "pdg="     << pdg
        << "p="       << p
        << "yt="      << yt
        << "zt="      << zt
        << "dydx="    << dydx
        << "dzdx="    << dzdx
        << "rowc="    << rc
        << "dy="      << dy
        << "dz="      << dz
        << "dphi="    << dphi
        << "\n";
    }

    Int_t istk = AliTRDgeometry::GetStack(det); 
    AliTRDpadPlane *pp = fGeo->GetPadPlane(ily, istk);
    Float_t zr0 = pp->GetRow0() + AliTRDSimParam::Instance()->GetAnodeWireOffset();
    Float_t tilt = fTracklet->GetTilt();

    Double_t x,y;
    AliTRDcluster *c = 0x0;
    fTracklet->ResetClusterIter(kFALSE);
    while((c = fTracklet->PrevCluster())){
      Float_t  q = TMath::Abs(c->GetQ());
      //AliTRDseedV1::GetClusterXY(c,x,y);
      x = c->GetX(); y = c->GetY();
      Float_t xc = x;
      Float_t yc = y;
      Float_t zc = c->GetZ();
      dx = x0 - xc; 
      Float_t yt = y0 - dx*dydx;
      Float_t zt = z0 - dx*dzdx; 
      dy = yt - (yc - tilt*(zc-zt));

      // Fill Histograms
      if(q>20. && q<250.) ((TH2I*)fContainer->At(kMCcluster))->Fill(dydx, dy);

      // Fill calibration container
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

      // Fill Debug Tree
      if(fDebugLevel>=2){
        //clInfo->Print();
        (*fDebugStream) << "MCcluster"
          <<"clInfo.=" << clInfo
          << "\n";
      }
    }
  }
  return h;
}


//________________________________________________________
Bool_t AliTRDtrackingResolution::GetRefFigure(Int_t ifig)
{
  Float_t y[2] = {0., 0.};
  TBox *b = 0x0;
  TAxis *ax = 0x0;
  TGraphErrors *g = 0x0;
  switch(ifig){
  case kCluster:
    if(!(g = (TGraphErrors*)fGraphS->At(ifig))) break;
    g->Draw("apl");
    ax = g->GetHistogram()->GetYaxis();
    y[0] = -0.5; y[1] = 2.5;
    ax->SetRangeUser(y[0], y[1]);
    ax->SetTitle("Cluster-Track Pools #sigma/#mu [mm]");
    ax = g->GetHistogram()->GetXaxis();
    ax->SetTitle("tg(#phi)");
    if(!(g = (TGraphErrors*)fGraphM->At(ifig))) break;
    g->Draw("pl");
    b = new TBox(-.15, y[0], .15, y[1]);
    b->SetFillStyle(3002);b->SetFillColor(kGreen);
    b->SetLineColor(0); b->Draw();
    return kTRUE;
  case kTrackletY:
    if(!(g = (TGraphErrors*)fGraphS->At(ifig))) break;
    g->Draw("apl");
    ax = g->GetHistogram()->GetYaxis();
    ax->SetRangeUser(-.5, 3.);
    ax->SetTitle("Tracklet-Track Y-Pools #sigma/#mu [mm]");
    ax = g->GetHistogram()->GetXaxis();
    ax->SetTitle("tg(#phi)");
    if(!(g = (TGraphErrors*)fGraphM->At(ifig))) break;
    g->Draw("pl");
    b = new TBox(-.15, y[0], .15, y[1]);
    b->SetFillStyle(3002);b->SetFillColor(kGreen);
    b->SetLineColor(0); b->Draw();
    return kTRUE;
  case kTrackletPhi:
    if(!(g = (TGraphErrors*)fGraphS->At(ifig))) break;
    g->Draw("apl");
    ax = g->GetHistogram()->GetYaxis();
    ax->SetRangeUser(-.5, 2.);
    ax->SetTitle("Tracklet-Track Phi-Residuals #sigma/#mu [rad]");
    ax = g->GetHistogram()->GetXaxis();
    ax->SetTitle("tg(#phi)");
    if(!(g = (TGraphErrors*)fGraphM->At(ifig))) break;
    g->Draw("pl");
    b = new TBox(-.15, y[0], .15, y[1]);
    b->SetFillStyle(3002);b->SetFillColor(kGreen);
    b->SetLineColor(0); b->Draw();
    return kTRUE;
  case kMCcluster:
    if(!(g = (TGraphErrors*)fGraphS->At(ifig))) break;
    ax = g->GetHistogram()->GetYaxis();
    y[0] = -.1; y[1] = 1.5;
    ax->SetRangeUser(y[0], y[1]);
    ax->SetTitle("Y_{cluster} #sigma/#mu [mm]");
    ax = g->GetHistogram()->GetXaxis();
    ax->SetTitle("tg(#phi)");
    g->Draw("apl");
    if(!(g = (TGraphErrors*)fGraphM->At(ifig))) break;
    g->Draw("pl");
    b = new TBox(-.15, y[0], .15, y[1]);
    b->SetFillStyle(3002);b->SetFillColor(kBlue);
    b->SetLineColor(0); b->Draw();
    return kTRUE;
  case kMCtrackletY:
    if(!(g = (TGraphErrors*)fGraphS->At(ifig))) break;
    ax = g->GetHistogram()->GetYaxis();
    y[0] = -.05; y[1] = 0.3;
    ax->SetRangeUser(y[0], y[1]);
    ax->SetTitle("Y_{tracklet} #sigma/#mu [mm]");
    ax = g->GetHistogram()->GetXaxis();
    ax->SetTitle("tg(#phi)");
    g->Draw("apl");
    if(!(g = (TGraphErrors*)fGraphM->At(ifig))) break;
    g->Draw("pl");
    b = new TBox(-.15, y[0], .15, y[1]);
    b->SetFillStyle(3002);b->SetFillColor(kBlue);
    b->SetLineColor(0); b->Draw();
    return kTRUE;
  case kMCtrackletZ:
    if(!(g = (TGraphErrors*)fGraphS->At(ifig))) break;
    ax = g->GetHistogram()->GetYaxis();
    ax->SetRangeUser(-.5, 1.);
    ax->SetTitle("Z_{tracklet} #sigma/#mu [mm]");
    ax = g->GetHistogram()->GetXaxis();
    ax->SetTitle("tg(#theta)");
    g->Draw("apl");
    if(!(g = (TGraphErrors*)fGraphM->At(ifig))) break;
    g->Draw("pl");
    return kTRUE;
  case kMCtrackletPhi:
    if(!(g = (TGraphErrors*)fGraphS->At(ifig))) break;
    ax = g->GetHistogram()->GetYaxis();
    y[0] = -.05; y[1] = .2;
    ax->SetRangeUser(y[0], y[1]);
    ax->SetTitle("#Phi_{tracklet} #sigma/#mu [deg]");
    ax = g->GetHistogram()->GetXaxis();
    ax->SetTitle("tg(#phi)");
    g->Draw("apl");
    if(!(g = (TGraphErrors*)fGraphM->At(ifig))) break;
    g->Draw("pl");
    return kTRUE;
  case kMCtrackY:
    if(!(g = (TGraphErrors*)fGraphS->At(ifig))) break;
    ax = g->GetHistogram()->GetYaxis();
    y[0] = -.05; y[1] = 0.2;
    ax->SetRangeUser(y[0], y[1]);
    ax->SetTitle("Y_{track} #sigma/#mu [mm]");
    ax = g->GetHistogram()->GetXaxis();
    ax->SetTitle("tg(#phi)");
    g->Draw("apl");
    if(!(g = (TGraphErrors*)fGraphM->At(ifig))) break;
    g->Draw("pl");
    b = new TBox(-.15, y[0], .15, y[1]);
    b->SetFillStyle(3002);b->SetFillColor(kBlue);
    b->SetLineColor(0); b->Draw();
    return kTRUE;
  case kMCtrackZ:
    if(!(g = (TGraphErrors*)fGraphS->At(ifig))) break;
    ax = g->GetHistogram()->GetYaxis();
    ax->SetRangeUser(-.5, 2.);
    ax->SetTitle("Z_{track} #sigma/#mu [mm]");
    ax = g->GetHistogram()->GetXaxis();
    ax->SetTitle("tg(#theta)");
    g->Draw("apl");
    if(!(g = (TGraphErrors*)fGraphM->At(ifig))) break;
    g->Draw("pl");
    return kTRUE;
  case kMCtrackPt:
    if(!(g = (TGraphErrors*)fGraphS->At(ifig))) break;
    ax = g->GetHistogram()->GetYaxis();
    ax->SetRangeUser(-.5, 2.);
    ax->SetTitle("#epsilon_{P_{t}}^{track} / #mu [%]");
    ax = g->GetHistogram()->GetXaxis();
    ax->SetTitle("1/p_{t}");
    g->Draw("apl");
    if(!(g = (TGraphErrors*)fGraphM->At(ifig))) break;
    g->Draw("pl");
    return kTRUE;
  }
  AliInfo(Form("Reference plot [%d] missing result", ifig));
  return kFALSE;
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
  TGraphErrors *gm = 0x0, *gs = 0x0;
  if(!fGraphS){ 
    fGraphS = new TObjArray(fNRefFigures);
    fGraphS->SetOwner();
    for(Int_t ig=0; ig<fNRefFigures; ig++){
      gs = new TGraphErrors();
      gs->SetLineColor(kRed);
      gs->SetMarkerStyle(23);
      gs->SetMarkerColor(kRed);
      gs->SetNameTitle(Form("s_%d", ig), "");
      fGraphS->AddAt(gs, ig);
    }
  }
  if(!fGraphM){ 
    fGraphM = new TObjArray(fNRefFigures);
    fGraphM->SetOwner();
    for(Int_t ig=0; ig<fNRefFigures; ig++){
      gm = new TGraphErrors();
      gm->SetLineColor(kBlue);
      gm->SetMarkerStyle(7);
      gm->SetMarkerColor(kBlue);
      gm->SetNameTitle(Form("m_%d", ig), "");
      fGraphM->AddAt(gm, ig);
    }
  }

  TH2I *h2 = 0x0;
  TH1D *h = 0x0;

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
  h2 = (TH2I *)(fContainer->At(kCluster));
  gm = (TGraphErrors*)fGraphM->At(kCluster);
  gs = (TGraphErrors*)fGraphS->At(kCluster);
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
  h2 = (TH2I *)(fContainer->At(kTrackletY));
  gm = (TGraphErrors*)fGraphM->At(kTrackletY);
  gs = (TGraphErrors*)fGraphS->At(kTrackletY);
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
  h2 = (TH2I *)(fContainer->At(kTrackletPhi));
  gm = (TGraphErrors*)fGraphM->At(kTrackletPhi);
  gs = (TGraphErrors*)fGraphS->At(kTrackletPhi);
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

  if(!HasMCdata()){
    if(c) delete c;
    return kTRUE;
  }

  //PROCESS MC RESIDUAL DISTRIBUTIONS

  // cluster y resolution
  h2 = (TH2I*)fContainer->At(kMCcluster);
  gm = (TGraphErrors*)fGraphM->At(kMCcluster);
  gs = (TGraphErrors*)fGraphS->At(kMCcluster);
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
  h2 = (TH2I*)fContainer->At(kMCtrackletY);
  gm = (TGraphErrors*)fGraphM->At(kMCtrackletY);
  gs = (TGraphErrors*)fGraphS->At(kMCtrackletY);
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
  h2 = (TH2I*)fContainer->At(kMCtrackletZ);
  gm = (TGraphErrors*)fGraphM->At(kMCtrackletZ);
  gs = (TGraphErrors*)fGraphS->At(kMCtrackletZ);
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
  h2 = (TH2I*)fContainer->At(kMCtrackletPhi);
  gm = (TGraphErrors*)fGraphM->At(kMCtrackletPhi);
  gs = (TGraphErrors*)fGraphS->At(kMCtrackletPhi);
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

  // track y resolution
  h2 = (TH2I*)fContainer->At(kMCtrackY);
  gm = (TGraphErrors*)fGraphM->At(kMCtrackY);
  gs = (TGraphErrors*)fGraphS->At(kMCtrackY);
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

  // track z resolution
  h2 = (TH2I*)fContainer->At(kMCtrackZ);
  gm = (TGraphErrors*)fGraphM->At(kMCtrackZ);
  gs = (TGraphErrors*)fGraphS->At(kMCtrackZ);
  for(Int_t iphi=1; iphi<=h2->GetNbinsX(); iphi++){
    h = h2->ProjectionY("pz", iphi, iphi);
    if(h->GetEntries()<70) continue;
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

  // track Pt resolution
  h2 = (TH2I*)fContainer->At(kMCtrackPt);
  TAxis *ax = h2->GetXaxis();
  gm = (TGraphErrors*)fGraphM->At(kMCtrackPt);
  gs = (TGraphErrors*)fGraphS->At(kMCtrackPt);
  TF1 fg("fg", "gaus", -1.5, 1.5);
  TF1 fl("fl", "landau", -4., 15.);
  TF1 fgl("fgl", "gaus(0)+landau(3)", -5., 20.);
  for(Int_t ip=1; ip<=ax->GetNbins(); ip++){
    h = h2->ProjectionY("ppt", ip, ip);
    if(h->GetEntries()<70) continue;

    h->Fit(&fg, "QN", "", -1.5, 1.5);
    fgl.SetParameter(0, fg.GetParameter(0));
    fgl.SetParameter(1, fg.GetParameter(1));
    fgl.SetParameter(2, fg.GetParameter(2));
    h->Fit(&fl, "QN", "", -4., 15.);
    fgl.SetParameter(3, fl.GetParameter(0));
    fgl.SetParameter(4, fl.GetParameter(1));
    fgl.SetParameter(5, fl.GetParameter(2));

    if(IsVisual()){c->cd(); c->SetLogy();}
    h->Fit(&fgl, opt, "", -5., 20.);
    if(IsVisual()){c->Modified(); c->Update(); gSystem->Sleep(500);}

    Float_t invpt = ax->GetBinCenter(ip);
    Int_t ip = gm->GetN();
    gm->SetPoint(ip, invpt, fgl.GetParameter(1));
    gm->SetPointError(ip, 0., fgl.GetParError(1));
    gs->SetPoint(ip, invpt, fgl.GetParameter(2)*invpt);
    gs->SetPointError(ip, 0., fgl.GetParError(2));
    // fgl.GetParameter(4) // Landau MPV
    // fgl.GetParameter(5) // Landau Sigma
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

  fContainer  = new TObjArray(10);
  //fContainer->SetOwner(kTRUE);

  TH1 *h = 0x0;
  // cluster to tracklet residuals [2]
  if(!(h = (TH2I*)gROOT->FindObject("hCl"))){
    h = new TH2I("hCl", "Clusters-Track Residuals", 21, -.33, .33, 100, -.5, .5);
    h->GetXaxis()->SetTitle("tg(#phi)");
    h->GetYaxis()->SetTitle("#Delta y [cm]");
    h->GetZaxis()->SetTitle("entries");
  } else h->Reset();
  fContainer->AddAt(h, kCluster);

  // tracklet to track residuals [2]
  if(!(h = (TH2I*)gROOT->FindObject("hTrkltY"))){
    h = new TH2I("hTrkltY", "Tracklets-Track Residuals (Y)", 21, -.33, .33, 100, -.5, .5);
    h->GetXaxis()->SetTitle("tg(#phi)");
    h->GetYaxis()->SetTitle("#Delta y [cm]");
    h->GetZaxis()->SetTitle("entries");
  } else h->Reset();
  fContainer->AddAt(h, kTrackletY);

  // tracklet to track residuals angular [2]
  if(!(h = (TH2I*)gROOT->FindObject("hTrkltPhi"))){
    h = new TH2I("hTrkltPhi", "Tracklets-Track Residuals (#Phi)", 21, -.33, .33, 100, -.5, .5);
    h->GetXaxis()->SetTitle("tg(#phi)");
    h->GetYaxis()->SetTitle("#Delta phi [#circ]");
    h->GetZaxis()->SetTitle("entries");
  } else h->Reset();
  fContainer->AddAt(h, kTrackletPhi);


  // Resolution histos
  if(!HasMCdata()) return fContainer;

  // cluster y resolution [0]
  if(!(h = (TH2I*)gROOT->FindObject("hMCcl"))){
    h = new TH2I("hMCcl", "Cluster Resolution", 31, -.48, .48, 100, -.5, .5);
    h->GetXaxis()->SetTitle("tg(#phi)");
    h->GetYaxis()->SetTitle("#Delta y [cm]");
    h->GetZaxis()->SetTitle("entries");
  } else h->Reset();
  fContainer->AddAt(h, kMCcluster);

  // tracklet y resolution [0]
  if(!(h = (TH2I*)gROOT->FindObject("hMCtrkltY"))){
    h = new TH2I("hMCtrkltY", "Tracklet Resolution (Y)", 31, -.48, .48, 100, -.5, .5);
    h->GetXaxis()->SetTitle("tg(#phi)");
    h->GetYaxis()->SetTitle("#Delta y [cm]");
    h->GetZaxis()->SetTitle("entries");
  } else h->Reset();
  fContainer->AddAt(h, kMCtrackletY);

  // tracklet y resolution [0]
  if(!(h = (TH2I*)gROOT->FindObject("hMCtrkltZ"))){
    h = new TH2I("hMCtrkltZ", "Tracklet Resolution (Z)", 31, -.48, .48, 100, -.5, .5);
    h->GetXaxis()->SetTitle("tg(#theta)");
    h->GetYaxis()->SetTitle("#Delta z [cm]");
    h->GetZaxis()->SetTitle("entries");
  } else h->Reset();
  fContainer->AddAt(h, kMCtrackletZ);

  // tracklet angular resolution [1]
  if(!(h = (TH2I*)gROOT->FindObject("hMCtrkltPhi"))){
    h = new TH2I("hMCtrkltPhi", "Tracklet Resolution (#Phi)", 31, -.48, .48, 100, -10., 10.);
    h->GetXaxis()->SetTitle("tg(#phi)");
    h->GetYaxis()->SetTitle("#Delta #phi [deg]");
    h->GetZaxis()->SetTitle("entries");
  } else h->Reset();
  fContainer->AddAt(h, kMCtrackletPhi);

  // Kalman track y resolution
  if(!(h = (TH2I*)gROOT->FindObject("hMCtrkY"))){
    h = new TH2I("hMCtrkY", "Kalman Track Resolution (Y)", 31, -.48, .48, 100, -.5, .5);
    h->GetXaxis()->SetTitle("tg(#phi)");
    h->GetYaxis()->SetTitle("#Delta y [cm]");
    h->GetZaxis()->SetTitle("entries");
  } else h->Reset();
  fContainer->AddAt(h, kMCtrackY);

  // Kalman track Z resolution
  if(!(h = (TH2I*)gROOT->FindObject("hMCtrkZ"))){
    h = new TH2I("hMCtrkZ", "Kalman Track Resolution (Z)", 20, -1., 1., 100, -1.5, 1.5);
    h->GetXaxis()->SetTitle("tg(#theta)");
    h->GetYaxis()->SetTitle("#Delta z [cm]");
    h->GetZaxis()->SetTitle("entries");
  } else h->Reset();
  fContainer->AddAt(h, kMCtrackZ);

  // Kalman track Pt resolution
  if(!(h = (TH2I*)gROOT->FindObject("hMCtrkPt"))){
    h = new TH2I("hMCtrkPt", "Kalman Track Resolution (Pt)", 100, 0., 2., 150, -5., 20.);
    h->GetXaxis()->SetTitle("1/p_{t}");
    h->GetYaxis()->SetTitle("#Delta p_{t} [GeV/c]");
    h->GetZaxis()->SetTitle("entries");
  } else h->Reset();
  fContainer->AddAt(h, kMCtrackPt);

  return fContainer;
}


//________________________________________________________
void AliTRDtrackingResolution::SetRecoParam(AliTRDrecoParam *r)
{

  fReconstructor->SetRecoParam(r);
}
