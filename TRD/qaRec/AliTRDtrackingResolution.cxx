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
#include <TList.h>
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
AliTRDtrackingResolution::AliTRDtrackingResolution(const char * name):
  AliAnalysisTask(name, "")
  ,fTrackObjects(0x0)
  ,fOutputHistograms(0x0)
  ,fYRes(0x0)
  ,fPhiRes(0x0)
  ,fReconstructor(0x0)
  ,fDebugLevel(0)
  ,fDebugStream(0x0)
{
  AliCDBManager *cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage("local://$ALICE_ROOT");
  cdb->SetRun(0);
  TGeoManager::Import("geometry.root");
  fReconstructor = new AliTRDReconstructor();
  fReconstructor->SetRecoParam(AliTRDrecoParam::GetLowFluxParam());
  AliTRDtrackerV1::SetNTimeBins(24);

  DefineInput(0, TObjArray::Class());
  DefineOutput(0, TList::Class());
}

//________________________________________________________
void AliTRDtrackingResolution::ConnectInputData(Option_t *){
  fTrackObjects = dynamic_cast<TObjArray *>(GetInputData(0));
}

//________________________________________________________
void AliTRDtrackingResolution::CreateOutputObjects()
{
  // spatial resolution
  OpenFile(0, "RECREATE");
  fOutputHistograms = new TList();
  fYRes = new TH2I("fY", "", 21, -21., 21., 100, -.5, .5);
  fOutputHistograms->Add(fYRes);
  fPhiRes = new TH2I("fPhi", "", 21, -21., 21., 100, -10., 10.);
  fOutputHistograms->Add(fPhiRes);
}

//________________________________________________________
void AliTRDtrackingResolution::Exec(Option_t *)
{
  // spatial Resolution: res = pos_{Tracklet}(x = x_{Anode wire}) - pos_{TrackRef}(x = x_{Anode wire})
  // angular Resolution: res = Tracklet angle - TrackRef Angle

  Int_t nTrackInfos = fTrackObjects->GetEntriesFast();
  if(fDebugLevel>=2) printf("Number of Histograms: %d\n", fOutputHistograms->GetEntries());

  Double_t dy, dz;
  AliTRDtrackV1 *fTrack = 0x0;
  AliTRDtrackInfo *fInfo = 0x0;
  if(fDebugLevel>=2) printf("Number of TrackInfos: %d\n", nTrackInfos);
  for(Int_t iTI = 0; iTI < nTrackInfos; iTI++){
    // check if ESD and MC-Information are available
    if(!(fInfo = dynamic_cast<AliTRDtrackInfo *>(fTrackObjects->UncheckedAt(iTI)))) continue;
    if(!(fTrack = fInfo->GetTRDtrack())) continue;
    if(fDebugLevel>=3) printf("\tDoing track[%d] NTrackRefs[%d]\n", iTI, fInfo->GetNTrackRefs());

    if(fInfo->GetNTrackRefs() < 2) continue; 
  
    AliTRDseedV1 *fTracklet = 0x0;
    AliTrackReference *fTrackRefs[2];
    for(Int_t iplane = 0; iplane < kNLayers; iplane++){
      if(!(fTracklet = fTrack->GetTracklet(iplane))) continue;
      if(!fTracklet->IsOK()) continue;

      if(fDebugLevel>=4) printf("\t\tLy[%d] X0[%6.3f] Ncl[%d]\n", iplane, fTracklet->GetX0(), fTracklet->GetN());

      // check for 2 track ref where the radial position has a distance less than 3.7mm
      Int_t nFound = 0;
      memset(fTrackRefs, 0, sizeof(AliTrackReference*) * 2);
      AliTrackReference *tempTrackRef = 0;
      for(Int_t itr = 0; itr < fInfo->GetNTrackRefs(); itr++){
        if(fDebugLevel>=5) printf("nFound = %d\n", nFound);
        if(nFound >= 2) break;
        tempTrackRef = fInfo->GetTrackRef(itr);
        if(!tempTrackRef) continue;
        if(fDebugLevel>=5) printf("TrackRef %d: x = %f\n", itr, tempTrackRef->LocalX());
        if(fTracklet->GetX0() - tempTrackRef->LocalX() > 3.7) continue;
        if(tempTrackRef->LocalX() - fTracklet->GetX0() > 3.7) break;
        if(fDebugLevel>=5) printf("accepted\n");
        if(nFound == 1)
          if(fTrackRefs[0]->LocalX() >= tempTrackRef->LocalX()) continue;
        fTrackRefs[nFound++] = tempTrackRef;
      }
      if(fDebugLevel>=4) printf("\t\tFound track crossing [%d]\n", nFound);
      if(nFound < 2) continue;
      // We found 2 track refs for the tracklet, get y and z at the anode wire by a linear approximation


      // RESOLUTION
      Double_t dx = fTrackRefs[1]->LocalX() - fTrackRefs[0]->LocalX();
      Double_t dydx = (fTrackRefs[1]->LocalY() - fTrackRefs[0]->LocalY()) / dx;
      Double_t dzdx = (fTrackRefs[1]->Z() - fTrackRefs[0]->Z()) / dx;
      Double_t dx0 = fTrackRefs[1]->LocalX() - fTracklet->GetX0();
      Double_t ymc =  fTrackRefs[1]->LocalY() - dydx*dx0;
      Double_t zmc =  fTrackRefs[1]->Z() - dzdx*dx0;
      
      // recalculate tracklet based on the MC info
      fTracklet->SetYref(0, zmc);
      fTracklet->SetYref(1, dzdx);
      fTracklet->Fit();
      dy = fTracklet->GetYfit(0) - ymc;
      dz = fTracklet->GetZfit(0) - zmc;
      //res_y *= 100; // in mm
      Double_t momentum = fTrackRefs[0]->P();

//       Double_t phi     = fTrackRefs[0]->Phi();
//       Double_t theta   = fTrackRefs[0]->Theta();
      Double_t phi   = TMath::ATan(dydx);
      Double_t theta = TMath::ATan(dzdx);
      Double_t dphi   = TMath::ATan(fTracklet->GetYfit(1)) - phi;
      
      // Fill Histograms
      if(TMath::Abs(dx-3.7)<1.E-3){
        fYRes->Fill(phi*TMath::RadToDeg(), dy);
        fPhiRes->Fill(phi*TMath::RadToDeg(), dphi*TMath::RadToDeg());
      }      

      if(fDebugLevel>=4) printf("\t\tdx[%6.4f] dy[%6.4f] dz[%6.4f] dphi[%6.4f] \n", dx, dy, dz, dphi);
      
      // Fill Debug Tree
      if(fDebugLevel>=1){
        (*fDebugStream) << "Resolution"
          << "plane="	 	<< iplane
          << "p="       << momentum
          << "dx="      << dx
          << "dy="		  << dy
          << "dz="	 	  << dz
          << "phi="			<< phi
          << "theta="		<< theta
          << "dphi="		<< dphi
          << "\n";
      }
    }
    if(fTrack->GetNumberOfTracklets() < 6) continue;

    // RESIDUALS
    Int_t nc = 0;
    AliTrackPoint tr[AliTRDtrackV1::kMAXCLUSTERSPERTRACK], cl[AliTRDtrackV1::kMAXCLUSTERSPERTRACK];
    // prepare
    for(Int_t iplane = 0; iplane < kNLayers; iplane++){
      if(!(fTracklet = fTrack->GetTracklet(iplane))) continue;
      if(!fTracklet->IsOK()) continue;
      AliTRDcluster *c = 0x0;
      for(Int_t ic=AliTRDseed::knTimebins-1; ic>=0; ic--){
        if(!(c = fTracklet->GetClusters(ic))) continue;
        tr[nc].SetXYZ(c->GetX(), 0., 0.);
        cl[nc].SetXYZ(c->GetX(), c->GetY(), c->GetZ());
        nc++;
      }
    }
    AliTRDtrackerV1::FitRiemanTilt(fTrack, 0x0, kTRUE, nc, tr);
    for(Int_t ip=0; ip<nc; ip++){
      dy = cl[ip].GetY() - tr[ip].GetY();
      dz = cl[ip].GetZ() - tr[ip].GetZ();
      if(fDebugLevel>=1){
        (*fDebugStream) << "ResidualsRT"
          << "dy="		  << dy
          << "dz="	 	  << dz
/*          << "phi="			<< phi
          << "theta="		<< theta
          << "dphi="		<< dphi*/
          << "\n";
      }
    }

//     AliTRDtrackerV1::FitKalman(fTrack, 0x0, kFALSE, nc, tr);
//     for(Int_t ip=0; ip<nc; ip++){
//       dy = cl[ip].GetY() - tr[ip].GetY();
//       dz = cl[ip].GetZ() - tr[ip].GetZ();
//       if(fDebugLevel>=1){
//         (*fDebugStream) << "ResidualsKF"
//           << "dy="		  << dy
//           << "dz="	 	  << dz
// /*          << "phi="			<< phi
//           << "theta="		<< theta
//           << "dphi="		<< dphi*/
//           << "\n";
//       }
//     }    
  }
  PostData(0, fOutputHistograms);
}

//________________________________________________________
void AliTRDtrackingResolution::Terminate(Option_t *)
{
  if(fDebugStream) delete fDebugStream;

  //process distributions
  fOutputHistograms = dynamic_cast<TList*>(GetOutputData(0));
  if (!fOutputHistograms) {
    Printf("ERROR: list not available");
    return;
  }
  TH2I *h2 = 0x0;
  TH1D *h = 0x0;

  // y resolution
  TF1 *f = new TF1("f1", "gaus", -.5, .5);  
  h2 = (TH2I*)fOutputHistograms->At(0);
  TGraphErrors *gm = new TGraphErrors(h2->GetNbinsX());
  gm->SetNameTitle("meany", "Mean dy");
  TGraphErrors *gs = new TGraphErrors(h2->GetNbinsX());
  gs->SetNameTitle("sigmy", "Sigma y");
  for(Int_t iphi=1; iphi<=h2->GetNbinsX(); iphi++){
    Double_t phi = h2->GetXaxis()->GetBinCenter(iphi);
    f->SetParameter(1, 0.);f->SetParameter(2, 2.e-2);
    h = h2->ProjectionY("py", iphi, iphi);
    h->Fit(f, "q", "goff", -.5, .5);
    Int_t jphi = iphi -1;
    gm->SetPoint(jphi, phi, f->GetParameter(1));
    gm->SetPointError(jphi, 0., f->GetParError(1));
    gs->SetPoint(jphi, phi, f->GetParameter(2));
    gs->SetPointError(jphi, 0., f->GetParError(2));
  }
  fOutputHistograms->Add(gm);
  fOutputHistograms->Add(gs);

  // phi resolution
  h2 = (TH2I*)fOutputHistograms->At(1);
  gm = new TGraphErrors(h2->GetNbinsX());
  gm->SetNameTitle("meanphi", "Mean Phi");
  gs = new TGraphErrors(h2->GetNbinsX());
  gs->SetNameTitle("sigmphi", "Sigma Phi");
  for(Int_t iphi=1; iphi<=h2->GetNbinsX(); iphi++){
    Double_t phi = h2->GetXaxis()->GetBinCenter(iphi);
    f->SetParameter(1, 0.);f->SetParameter(2, 2.e-2);
    h = h2->ProjectionY("py", iphi, iphi);
    h->Fit(f, "q", "goff", -.5, .5);
    Int_t jphi = iphi -1;
    gm->SetPoint(jphi, phi, f->GetParameter(1));
    gm->SetPointError(jphi, 0., f->GetParError(1));
    gs->SetPoint(jphi, phi, f->GetParameter(2));
    gs->SetPointError(jphi, 0., f->GetParError(2));
  }
  fOutputHistograms->Add(gm);
  fOutputHistograms->Add(gs);


  delete f;
}

//________________________________________________________
void AliTRDtrackingResolution::SetDebugLevel(Int_t level){
  fDebugLevel = level;
  if(!fDebugLevel) return;
  if(fDebugStream) return;
  fDebugStream = new TTreeSRedirector("TRD.Resolution.root");
}

//________________________________________________________
void AliTRDtrackingResolution::SetRecoParam(AliTRDrecoParam *r)
{
  fReconstructor->SetRecoParam(r);
}
