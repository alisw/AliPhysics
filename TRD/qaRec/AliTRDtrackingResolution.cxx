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
//  Reconstruction QA                                                     //
//                                                                        //
//  Authors:                                                              //
//    Markus Fasel <M.Fasel@gsi.de>                                       //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include <cstring>


#include <TObjArray.h>
#include <TH2.h>
#include <TH1.h>
#include <TF1.h>
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
  ,fReconstructor(0x0)
  ,fGeo(0x0)
{
  fReconstructor = new AliTRDReconstructor();
  fReconstructor->SetRecoParam(AliTRDrecoParam::GetLowFluxParam());
  fGeo = new AliTRDgeometry();
}

//________________________________________________________
AliTRDtrackingResolution::~AliTRDtrackingResolution()
{
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

  // cluster to tracklet residuals [2]
  fContainer->AddAt(new TH2I("fYClRes", "Clusters Residuals", 21, -21., 21., 100, -.5, .5), kClusterYResidual);
  // tracklet to Riemann fit residuals [2]
  fContainer->AddAt(new TH2I("fYTrkltRRes", "Tracklet Riemann Residuals", 21, -21., 21., 100, -.5, .5), kTrackletRiemanYResidual);
  fContainer->AddAt(new TH2I("fAngleTrkltRRes", "Tracklet Riemann Angular Residuals", 21, -21., 21., 100, -.5, .5), kTrackletRiemanAngleResidual);
  fContainer->AddAt(new TH2I("fYTrkltKRes", "Tracklet Kalman Residuals", 21, -21., 21., 100, -.5, .5), kTrackletKalmanYResidual);
  fContainer->AddAt(new TH2I("fAngleTrkltKRes", "Tracklet Kalman Angular Residuals", 21, -21., 21., 100, -.5, .5), kTrackletKalmanAngleResidual);

  // Resolution histos
  if(HasMCdata()){
    // tracklet resolution [0]
    fContainer->AddAt(new TH2I("fY", "Tracklet Resolution", 21, -21., 21., 100, -.5, .5), kTrackletYResolution);
    // tracklet angular resolution [1]
    fContainer->AddAt(new TH2I("fPhi", "Tracklet Angular Resolution", 21, -21., 21., 100, -10., 10.), kTrackletAngleResolution);

    // Riemann track resolution [y, z, angular]
    fContainer->AddAt(new TH2I("fYRT", "Track Riemann Y Resolution", 21, -21., 21., 100, -.5, .5), kTrackRYResolution);
    fContainer->AddAt(new TH2I("fZRT", "Track Riemann Z Resolution", 21, -21., 21., 100, -.5, .5), kTrackRZResolution);
    fContainer->AddAt(new TH2I("fPhiRT", "Track Riemann Angular Resolution", 21, -21., 21., 100, -10., 10.), kTrackRAngleResolution);

    // Kalman track resolution [y, z, angular]
    fContainer->AddAt(new TH2I("fYKT", "", 21, -21., 21., 100, -.5, .5), kTrackKYResolution);
    fContainer->AddAt(new TH2I("fZKT", "", 21, -21., 21., 100, -.5, .5), kTrackKZResolution);
    fContainer->AddAt(new TH2I("fPhiKT", "", 21, -21., 21., 100, -10., 10.), kTrackKAngleResolution);
  }

  // CREATE GRAPHS for DISPLAY
  
  // define iterator over graphs
  Int_t jgraph = (Int_t)kGraphStart;
  TH2I *h2 = (TH2I *)(fContainer->At(kClusterYResidual));
  // clusters tracklet residuals (mean-phi)
  TH1 *h = new TH1I("h", "", 100, -40., 40.);
  h->GetXaxis()->SetTitle("#Phi [deg]");
  h->GetYaxis()->SetTitle("Clusters Residuals : #sigma/#mu [mm]");
  h->GetYaxis()->SetRangeUser(-.05, 1.);
  fContainer->AddAt(h, jgraph++);

  TGraphErrors *g = new TGraphErrors(h2->GetNbinsX());
  g->SetLineColor(kGreen);
  g->SetMarkerStyle(22);
  g->SetMarkerColor(kGreen);
  g->SetNameTitle("clm", "Residuals Clusters-Tracklet Mean");
  fContainer->AddAt(g, jgraph++);

  // clusters tracklet residuals (sigma-phi)
  g = new TGraphErrors(h2->GetNbinsX());
  g->SetLineColor(kRed);
  g->SetMarkerStyle(23);
  g->SetMarkerColor(kRed);
  g->SetNameTitle("cls", "Residuals Clusters-Tracklet Sigma");
  fContainer->AddAt(g, jgraph++);

  if(HasMCdata()){
    // tracklet y resolution
    h2 = (TH2I*)fContainer->At(kTrackletYResolution);
    h = new TH1I("h", "", 100, -40., 40.);
    h->GetXaxis()->SetTitle("#Phi [deg]");
    h->GetYaxis()->SetTitle("Tracklet Resolution : #sigma/#mu [mm]");
    h->GetYaxis()->SetRangeUser(-.05, 1.);
    fContainer->AddAt(h, jgraph++);

    g = new TGraphErrors(h2->GetNbinsX());
    g->SetLineColor(kGreen);
    g->SetMarkerStyle(22);
    g->SetMarkerColor(kGreen);
    g->SetNameTitle("trkltym", "Resolution Tracklet Y Mean");
    fContainer->AddAt(g, jgraph++);
    g = new TGraphErrors(h2->GetNbinsX());
    g->SetLineColor(kRed);
    g->SetMarkerStyle(22);
    g->SetMarkerColor(kRed);
    g->SetNameTitle("trkltys", "Resolution Tracklet Y Sigma");
    fContainer->AddAt(g, jgraph++);

    // tracklet phi resolution
    h2 = (TH2I*)fContainer->At(kTrackletAngleResolution);
    h = new TH1I("h", "", 100, -40., 40.);
    h->GetXaxis()->SetTitle("#Phi [deg]");
    h->GetYaxis()->SetTitle("Tracklet Angular Resolution : #sigma/#mu [deg]");
    h->GetYaxis()->SetRangeUser(-.05, .2);
    fContainer->AddAt(h, jgraph++);

    g = new TGraphErrors(h2->GetNbinsX());
    g->SetLineColor(kGreen);
    g->SetMarkerStyle(22);
    g->SetMarkerColor(kGreen);
    g->SetNameTitle("trkltam", "Resolution Tracklet Y Mean");
    fContainer->AddAt(g, jgraph++);
    g = new TGraphErrors(h2->GetNbinsX());
    g->SetLineColor(kRed);
    g->SetMarkerStyle(22);
    g->SetMarkerColor(kRed);
    g->SetNameTitle("trkltas", "Angle Resolution Tracklet Sigma");
    fContainer->AddAt(g, jgraph++);
  }
}

//________________________________________________________
void AliTRDtrackingResolution::Exec(Option_t *)
{
  // spatial Resolution: res = pos_{Tracklet}(x = x_{Anode wire}) - pos_{TrackRef}(x = x_{Anode wire})
  // angular Resolution: res = Tracklet angle - TrackRef Angle

  Int_t nTrackInfos = fTracks->GetEntriesFast();
  if(fDebugLevel>=2) printf("Number of Histograms: %d\n", Histos()->GetEntries());

  Int_t pdg;
  Double_t p, dy, dphi, dymc, dzmc, dphimc;
  Float_t fP[kNLayers], fX[kNLayers], fY[kNLayers], fZ[kNLayers], fPhi[kNLayers], fTheta[kNLayers]; // phi/theta angle per layer
  Bool_t fMCMap[kNLayers], fLayerMap[kNLayers]; // layer map

  AliTRDpadPlane *pp = 0x0;
  AliTrackPoint tr[kNLayers], tk[kNLayers];
  AliExternalTrackParam *fOp = 0x0;
  AliTRDtrackV1 *fTrack = 0x0;
  AliTRDtrackInfo *fInfo = 0x0;
  if(fDebugLevel>=2) printf("Number of TrackInfos: %d\n", nTrackInfos);
  for(Int_t iTI = 0; iTI < nTrackInfos; iTI++){
    // check if ESD and MC-Information are available
    if(!(fInfo = dynamic_cast<AliTRDtrackInfo *>(fTracks->UncheckedAt(iTI)))) continue;
    if(!(fTrack = fInfo->GetTRDtrack())) continue;
    if(!(fOp = fInfo->GetOuterParam())) continue;
    pdg = fInfo->GetPDG();

    if(fDebugLevel>=3) printf("\tDoing track[%d] NTrackRefs[%d]\n", iTI, fInfo->GetNTrackRefs());

    p = fOp->P();
    Int_t npts = 0;
    memset(fP, 0, kNLayers*sizeof(Float_t));
    memset(fX, 0, kNLayers*sizeof(Float_t));
    memset(fY, 0, kNLayers*sizeof(Float_t));
    memset(fZ, 0, kNLayers*sizeof(Float_t));
    memset(fPhi, 0, kNLayers*sizeof(Float_t));
    memset(fTheta, 0, kNLayers*sizeof(Float_t));
    memset(fLayerMap, 0, kNLayers*sizeof(Bool_t));
    memset(fMCMap, 0, kNLayers*sizeof(Bool_t));
    AliTRDseedV1 *fTracklet = 0x0;
    for(Int_t iplane = 0; iplane < kNLayers; iplane++){
      if(!(fTracklet = fTrack->GetTracklet(iplane))) continue;
      if(!fTracklet->IsOK()) continue;

      // Book point arrays
      fLayerMap[iplane] = kTRUE;
      tr[npts].SetXYZ(fTracklet->GetX0(), 0., 0.);
      tk[npts].SetXYZ(fTracklet->GetX0(), fTracklet->GetYfit(0), fTracklet->GetZfit(0));
      npts++;

      if(fDebugLevel>=4) printf("\t\tLy[%d] X0[%6.3f] Ncl[%d]\n", iplane, fTracklet->GetX0(), fTracklet->GetN());

      // define reference values
      fP[iplane]   = p;
      fX[iplane]   = fTracklet->GetX0();
      fY[iplane]   = fTracklet->GetYref(0);
      fZ[iplane]   = fTracklet->GetZref(0);
      fPhi[iplane] = TMath::ATan(fTracklet->GetYref(1));
      fTheta[iplane] = TMath::ATan(fTracklet->GetZref(1));
      

      // RESOLUTION (compare to MC)
      if(HasMCdata()){
        if(fInfo->GetNTrackRefs() >= 2){ 
          Double_t pmc, ymc, zmc, phiMC, thetaMC;
          if(Resolution(fTracklet, fInfo, pmc, ymc, zmc, phiMC, thetaMC)){ 
            fMCMap[iplane] = kTRUE;
            fP[iplane]     = pmc;
            fY[iplane]     = ymc;
            fZ[iplane]     = zmc;
            fPhi[iplane]   = phiMC;
            fTheta[iplane] = thetaMC;
          }
        }
      }
      Float_t phi   = fPhi[iplane]*TMath::RadToDeg();
      //Float_t theta = fTheta[iplane]*TMath::RadToDeg();

      // Do clusters residuals
      if(!fTracklet->Fit(kFALSE)) continue;
      AliTRDcluster *c = 0x0;
      for(Int_t ic=AliTRDseed::knTimebins-1; ic>=0; ic--){
        if(!(c = fTracklet->GetClusters(ic))) continue;

        dy = fTracklet->GetYat(c->GetX()) - c->GetY();
        ((TH2I*)fContainer->At(kClusterYResidual))->Fill(phi, dy);

        if(fDebugLevel>=2){
          Float_t q = c->GetQ();
          // Get z-position with respect to anode wire
          AliTRDSimParam    *simParam    = AliTRDSimParam::Instance();
          Int_t det = c->GetDetector();
          Float_t x = c->GetX();
          Float_t z = fZ[iplane]-(fX[iplane]-x)*TMath::Tan(fTheta[iplane]);
          Int_t stack = fGeo->GetStack(det);
          pp = fGeo->GetPadPlane(iplane, stack);
          Float_t row0 = pp->GetRow0();
          Float_t d  =  row0 - z + simParam->GetAnodeWireOffset();
          d -= ((Int_t)(2 * d)) / 2.0;
          //if (d > 0.25) d  = 0.5 - d;
  
          (*fDebugStream) << "ResidualClusters"
            << "ly="   << iplane
            << "stk="  << stack
            << "pdg="  << pdg
            << "phi="  << fPhi[iplane]
            << "tht="  << fTheta[iplane]
            << "q="    << q
            << "x="    << x
            << "z="    << z
            << "d="    << d
            << "dy="   << dy
            << "\n";
        }
      }
      pp = 0x0;
    }


    // this protection we might drop TODO
    if(fTrack->GetNumberOfTracklets() < 6) continue;

    AliTRDtrackerV1::FitRiemanTilt(fTrack, 0x0, kTRUE, npts, tr);
    Int_t iref = 0;
    for(Int_t ip=0; ip<kNLayers; ip++){
      if(!fLayerMap[ip]) continue;
      fTracklet = fTrack->GetTracklet(ip);
      // recalculate fit based on the new tilt correction
      fTracklet->Fit();

      dy = fTracklet->GetYfit(0) - tr[iref].GetY();
      ((TH2I*)fContainer->At(kTrackletRiemanYResidual))->Fill(fPhi[ip]*TMath::RadToDeg(), dy);

      dphi = fTracklet->GetYfit(1)- fTracklet->GetYref(1);
      ((TH2I*)fContainer->At(kTrackletRiemanAngleResidual))->Fill(fPhi[ip]*TMath::RadToDeg(), dphi);

      if(HasMCdata()){
        dymc = fY[ip] - tr[iref].GetY();
        ((TH2I*)fContainer->At(kTrackRYResolution))->Fill(fPhi[ip]*TMath::RadToDeg(), dymc);

        dzmc = fZ[ip] - tr[iref].GetZ();
        ((TH2I*)fContainer->At(kTrackRZResolution))->Fill(fPhi[ip]*TMath::RadToDeg(), dzmc);

        dphimc = fPhi[ip] - fTracklet->GetYfit(1);
        ((TH2I*)fContainer->At(kTrackRAngleResolution))->Fill(fPhi[ip]*TMath::RadToDeg(), dphimc);
      }

      iref++;

      if(fDebugLevel>=2){
        (*fDebugStream) << "RiemannTrack"
          << "ly="    << ip
          << "mc="    << fMCMap[ip]
          << "p="     << fP[ip]
          << "phi="   << fPhi[ip]
          << "tht="   << fTheta[ip]
          << "dy="    << dy
          << "dphi="  << dphi
          << "dymc="  << dymc
          << "dzmc="  << dzmc
          << "dphimc="<< dphimc
          << "\n";
      }
    }

//  if(!gGeoManager) TGeoManager::Import("geometry.root");
//     AliTRDtrackerV1::FitKalman(fTrack, 0x0, kFALSE, nc, tr);
//     for(Int_t ip=0; ip<nc; ip++){
//       dy = cl[ip].GetY() - tr[ip].GetY();
//      ((TH2I*)fContainer->At(kTrackletKalmanYResidual))->Fill(phi*TMath::RadToDeg(), dy);
//       dz = cl[ip].GetZ() - tr[ip].GetZ();
//       if(fDebugLevel>=1){
//         (*fDebugStream) << "KalmanTrack"
//           << "dy="		  << dy
//           << "dz="	 	  << dz
// /*          << "phi="			<< phi
//           << "theta="		<< theta
//           << "dphi="		<< dphi*/
//           << "\n";
//       }
//     }    


  }
  PostData(0, fContainer);
}

//________________________________________________________
void AliTRDtrackingResolution::GetRefFigure(Int_t ifig, Int_t &first, Int_t &last, Option_t */*opt*/)
{
  switch(ifig){
  case 0:
    first = (Int_t)kGraphStart; last = first+3;
    break;
  case 1:
    first = (Int_t)kGraphStart+3; last = first+3;
    break;
  case 2:
    first = (Int_t)kGraphStart+6; last = first+3;
    break;
  default:
    first = (Int_t)kGraphStart; last = first;
    break;
  }
}


//________________________________________________________
Bool_t AliTRDtrackingResolution::Resolution(AliTRDseedV1 *tracklet, AliTRDtrackInfo *fInfo, Double_t &p, Double_t &ymc, Double_t &zmc, Double_t &phi, Double_t &theta)
{

  AliTrackReference *fTrackRefs[2] = {0x0, 0x0},   *tempTrackRef = 0x0;

  // check for 2 track ref where the radial position has a distance less than 3.7mm
  Int_t nFound = 0;
  for(Int_t itr = 0; itr < fInfo->GetNTrackRefs(); itr++){
    if(!(tempTrackRef = fInfo->GetTrackRef(itr))) continue;
    if(fDebugLevel>=5) printf("TrackRef %d: x = %f\n", itr, tempTrackRef->LocalX());
    if(TMath::Abs(tracklet->GetX0() - tempTrackRef->LocalX()) > 3.7) continue;
    fTrackRefs[nFound++] = tempTrackRef;
    if(nFound == 2) break;
  }
  if(nFound < 2){ 
    if(fDebugLevel>=4) printf("\t\tFound track crossing [%d] refX[%6.3f]\n", nFound, nFound>0 ? fTrackRefs[0]->LocalX() : 0.);
    return kFALSE;
  }
  // We found 2 track refs for the tracklet, get y and z at the anode wire by a linear approximation


  // RESOLUTION
  Double_t dx = fTrackRefs[1]->LocalX() - fTrackRefs[0]->LocalX();
  if(dx <= 0.){
    if(fDebugLevel>=4) printf("\t\ttrack ref in the wrong order refX0[%6.3f] refX1[%6.3f]\n", fTrackRefs[0]->LocalX(), fTrackRefs[1]->LocalX());
    return kFALSE;
  }
  Double_t dydx = (fTrackRefs[1]->LocalY() - fTrackRefs[0]->LocalY()) / dx;
  Double_t dzdx = (fTrackRefs[1]->Z() - fTrackRefs[0]->Z()) / dx;
  Double_t dx0 = fTrackRefs[1]->LocalX() - tracklet->GetX0();
  ymc =  fTrackRefs[1]->LocalY() - dydx*dx0;
  zmc =  fTrackRefs[1]->Z() - dzdx*dx0;
  
  // recalculate tracklet based on the MC info
  tracklet->SetZref(0, zmc);
  tracklet->SetZref(1, -dzdx); // TODO
  tracklet->Fit();
  Double_t dy = tracklet->GetYfit(0) - ymc;
  Double_t dz = tracklet->GetZfit(0) - zmc;
      
  //res_y *= 100; // in mm
  p = fTrackRefs[0]->P();

  phi   = TMath::ATan(dydx);
  theta = TMath::ATan(dzdx);
  Double_t dphi   = TMath::ATan(tracklet->GetYfit(1)) - phi;
  if(fDebugLevel>=4) printf("\t\tdx[%6.4f] dy[%6.4f] dz[%6.4f] dphi[%6.4f] \n", dx, dy, dz, dphi);
  
  // Fill Histograms
  if(TMath::Abs(dx-3.7)<1.E-3){
    ((TH2I*)fContainer->At(kTrackletYResolution))->Fill(phi*TMath::RadToDeg(), dy);
    ((TH2I*)fContainer->At(kTrackletAngleResolution))->Fill(phi*TMath::RadToDeg(), dphi*TMath::RadToDeg());
  }        
  // Fill Debug Tree
  if(fDebugLevel>=2){
    Int_t iplane = tracklet->GetPlane();
    Int_t pdg = fInfo->GetPDG();
    (*fDebugStream) << "ResolutionTrklt"
      << "ly="	 	  << iplane
      << "pdg="     << pdg
      << "p="       << p
      << "phi="			<< phi
      << "tht="		  << theta
      << "ymc="     << ymc
      << "zmc="     << zmc
      << "dx="      << dx
      << "dy="		  << dy
      << "dz="	 	  << dz
      << "dphi="		<< dphi
      << "\n";
  }

  return kTRUE;
}

//________________________________________________________
Bool_t AliTRDtrackingResolution::PostProcess()
{
  //fContainer = dynamic_cast<TObjArray*>(GetOutputData(0));
  fNRefFigures = 0;
  if (!fContainer) {
    Printf("ERROR: list not available");
    return kFALSE;
  }

  TH2I *h2 = 0x0;
  TH1D *h = 0x0;
  TGraphErrors *gm = 0x0, *gs = 0x0;
  TF1 f("f1", "gaus", -.5, .5);  
  // define iterator over graphs
  Int_t jgraph = (Int_t)kGraphStart;

  //PROCESS RESIDUAL DISTRIBUTIONS

  // Clusters residuals
  // define model
  TF1 fc("fc", "[0]*exp(-0.5*((x-[1])/[2])**2)+[3]*exp(-0.5*((x-[4])/[5])**2)", -.5, .5);
  h2 = (TH2I *)(fContainer->At(kClusterYResidual));
  jgraph++; //skip the frame histo 
  gm = (TGraphErrors*)fContainer->At(jgraph++);
  gs = (TGraphErrors*)fContainer->At(jgraph++);
  for(Int_t ibin = 1; ibin <= h2->GetNbinsX(); ibin++){
    Double_t phi = h2->GetXaxis()->GetBinCenter(ibin);
    Double_t dphi = h2->GetXaxis()->GetBinWidth(ibin)/2;
    h = h2->ProjectionY("py", ibin, ibin);
    Fit(h, &fc);
    gm->SetPoint(ibin - 1, phi, 10.*fc.GetParameter(1));
    gm->SetPointError(ibin - 1, dphi, 10.*fc.GetParError(1));
    gs->SetPoint(ibin - 1, phi, 10.*fc.GetParameter(2));
    gs->SetPointError(ibin - 1, dphi, 10.*fc.GetParError(2));
  }
  fNRefFigures++;


  //PROCESS RESOLUTION DISTRIBUTIONS
  if(HasMCdata()){
    // tracklet y resolution
    h2 = (TH2I*)fContainer->At(kTrackletYResolution);
    jgraph++; //skip the frame histo
    gm = (TGraphErrors*)fContainer->At(jgraph++);
    gs = (TGraphErrors*)fContainer->At(jgraph++);
    for(Int_t iphi=1; iphi<=h2->GetNbinsX(); iphi++){
      Double_t phi = h2->GetXaxis()->GetBinCenter(iphi);
      f.SetParameter(1, 0.);f.SetParameter(2, 2.e-2);
      h = h2->ProjectionY("py", iphi, iphi);
      Fit(h, &fc);
      Int_t jphi = iphi -1;
      gm->SetPoint(jphi, phi, 10.*f.GetParameter(1));
      gm->SetPointError(jphi, 0., 10.*f.GetParError(1));
      gs->SetPoint(jphi, phi, 10.*f.GetParameter(2));
      gs->SetPointError(jphi, 0., 10.*f.GetParError(2));
    }
    fNRefFigures++;
  
    // tracklet phi resolution
    h2 = (TH2I*)fContainer->At(kTrackletAngleResolution);
    jgraph++; //skip the frame histo
    gm = (TGraphErrors*)fContainer->At(jgraph++);
    gs = (TGraphErrors*)fContainer->At(jgraph++);
    for(Int_t iphi=1; iphi<=h2->GetNbinsX(); iphi++){
      Double_t phi = h2->GetXaxis()->GetBinCenter(iphi);
      f.SetParameter(1, 0.);f.SetParameter(2, 2.e-2);
      h = h2->ProjectionY("py", iphi, iphi);
      h->Fit(&f, "QN", "", -.5, .5);
      Int_t jphi = iphi -1;
      gm->SetPoint(jphi, phi, f.GetParameter(1));
      gm->SetPointError(jphi, 0., f.GetParError(1));
      gs->SetPoint(jphi, phi, f.GetParameter(2));
      gs->SetPointError(jphi, 0., f.GetParError(2));
    }
    fNRefFigures++;
  }

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
void AliTRDtrackingResolution::Fit(TH1 *h, TF1 *f)
{
// Helper function to avoid duplication of code
// Make first guesses on the fit parameters

  // find the intial parameters of the fit !! (thanks George)
  Int_t nbinsy = Int_t(.5*h->GetNbinsX());
  Double_t sum = 0.;
  for(Int_t jbin=nbinsy-4; jbin<=nbinsy+4; jbin++) sum+=h->GetBinContent(jbin); sum/=9.;
  f->SetParLimits(0, 0., 3.*sum);
  f->SetParameter(0, .9*sum);

  f->SetParLimits(1, -.1, .1);
  f->SetParameter(1, 0.);

  f->SetParLimits(2, 0., 1.e-1);
  f->SetParameter(2, 2.e-2);

  f->SetParLimits(3, 0., sum);
  f->SetParameter(3, .1*sum);

  f->SetParLimits(4, -.3, .3);
  f->SetParameter(4, 0.);

  f->SetParLimits(5, 0., 1.e2);
  f->SetParameter(5, 2.e-1);

  h->Fit(f, "QN", "", -0.5, 0.5);
}

//________________________________________________________
TObjArray* AliTRDtrackingResolution::Histos()
{
  if(!fContainer) fContainer  = new TObjArray(25);
  return fContainer;
}


//________________________________________________________
void AliTRDtrackingResolution::SetRecoParam(AliTRDrecoParam *r)
{

  fReconstructor->SetRecoParam(r);
}
