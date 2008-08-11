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

#include <TObjArray.h>
#include <TList.h>
#include <TH1F.h>
#include <TProfile.h>
#include <TMath.h>

#include "AliAnalysisManager.h"
#include "AliTRDseedV1.h"
#include "AliTrackReference.h"
#include "TTreeStream.h"

#include "AliTRDtrackInfo/AliTRDtrackInfo.h"
#include "AliTRDtrackingResolution.h"

#include <cstring>

ClassImp(AliTRDtrackingResolution)

//________________________________________________________
AliTRDtrackingResolution::AliTRDtrackingResolution(const char * name):
  AliAnalysisTask(name, ""),
  fTrackObjects(0x0),
  fOutputHistograms(0x0),
  fYres(0x0),
/*	fZres(0),
  fYresAngle(0),
  fPhiRes(0x0),
  fPhiResAngle(0x0),*/
  fDebugLevel(0),
  fDebugStream(0x0)
{
// 	memset(fYresLayer, 0, sizeof(TH1F *) * kNLayers);
// 	memset(fZresLayer, 0, sizeof(TH1F *) * kNLayers);
// 	memset(fYresLayerAngle, 0, sizeof(TProfile *) * kNLayers);
// 	memset(fPhiResLayer, 0, sizeof(TH1F *) * kNLayers);
// 	memset(fPhiResLayerAngle, 0, sizeof(TProfile *) * kNLayers);

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
  printf("Creating Histograms\n");
  OpenFile(0, "RECREATE");
  fOutputHistograms = new TList();
  fYres = new TH1F("fYres", "y-Resolution", 100, -1.5, 1.5);
  fOutputHistograms->Add(fYres);
// 	fZres = new TH1F("fZres", "z-Resolution", 100, -1.5, 1.5);
// 	fOutputHistograms->Add(fZres);
// 
// 	fYresAngle = new TProfile("fYresAngle", "y-Resolution - Angluar dependence", 80, -40, 40);
// 	fOutputHistograms->Add(fYresAngle);
// 
// 	// angular resolution
// 	fPhiRes = new TH1F("fPhiRes", "phi-resolution", 20, -10, 10);
// 	fOutputHistograms->Add(fPhiRes);
// 	
// 	fPhiResAngle = new TProfile("fPhiResAngle", "phi-resolution - Angular dependence", 80, -40, 40);
// 	fOutputHistograms->Add(fPhiResAngle);

/*	for(Int_t iplane = 0; iplane < kNLayers; iplane++){
    // spatial resolution
    fYresLayer[iplane] = new TH1F(Form("fYresLayer%d", iplane), Form("y-Resolution in Layer %d", iplane), 100, -1.5, 1.5);
    fOutputHistograms->Add(fYresLayer[iplane]);

    fZresLayer[iplane] = new TH1F(Form("fZresLayer%d", iplane), Form("z-Resolution in Layer %d", iplane), 100, -1.5, 1.5);
    fOutputHistograms->Add(fZresLayer[iplane]);

    fYresLayerAngle[iplane] = new TProfile(Form("fYresLayerAngle%d", iplane), Form("y-Resolution in Layer %d - Angluar dependence", iplane), 80, -40, 40);
    fOutputHistograms->Add(fYresLayerAngle[iplane]);

    // angular resolution
    fPhiResLayer[iplane] = new TH1F(Form("fPhiResLayer%d", iplane), Form("phi-resolution in Layer %d", iplane), 20, -10, 10);
    fOutputHistograms->Add(fPhiResLayer[iplane]);

    fPhiResLayerAngle[iplane] = new TProfile(Form("fPhiResAngle%d", iplane), Form("phi-resolution in Layer %d - Angular dependence", iplane), 80, -40, 40);
    fOutputHistograms->Add(fPhiResLayerAngle[iplane]);
  }*/
}

//________________________________________________________
void AliTRDtrackingResolution::Exec(Option_t *){
  // spatial Resolution: res = pos_{Tracklet}(x = x_{Anode wire}) - pos_{TrackRef}(x = x_{Anode wire})
  // angular Resolution: res = Tracklet angle - TrackRef Angle
  Int_t nTrackInfos = fTrackObjects->GetEntriesFast();
  if(fDebugLevel>=2) printf("Number of Histograms: %d\n", fOutputHistograms->GetEntries());
  AliTRDtrackInfo *fInfo = 0x0;
  if(fDebugLevel>=2) printf("Number of TrackInfos: %d\n", nTrackInfos);
  for(Int_t iTI = 0; iTI < nTrackInfos; iTI++){
    if(fDebugLevel>=2) printf("Doing Object %d\n", iTI);
    fInfo = dynamic_cast<AliTRDtrackInfo *>(fTrackObjects->UncheckedAt(iTI));
    // check if ESD and MC-Information are available
    if(!fInfo || !fInfo->GetTRDtrack() || fInfo->GetNTrackRefs() < 2) continue; 
    AliTRDseedV1*fTracklet = 0;
    AliTrackReference *fTrackRefs[2];
    for(Int_t iplane = 0; iplane < kNLayers; iplane++){
      if(fDebugLevel>=2) printf("plane %d\n", iplane);
      fTracklet = fInfo->GetTracklet(iplane);
      if(!fTracklet) continue;
      // check for 2 track ref where the radial position has a distance less than 3.7mm
      if(fDebugLevel>=2) printf("Find Track References for x = %f\n", fTracklet->GetX0());
      if(fDebugLevel>=2) printf("Number of Clusters: %d\n", fTracklet->GetN());
      Int_t nFound = 0;
      memset(fTrackRefs, 0, sizeof(AliTrackReference*) * 2);
      AliTrackReference *tempTrackRef = 0;
      for(Int_t itr = 0; itr < fInfo->GetNTrackRefs(); itr++){
        if(fDebugLevel>=2) printf("nFound = %d\n", nFound);
        if(nFound >= 2) break;
        tempTrackRef = fInfo->GetTrackRef(itr);
        if(!tempTrackRef) continue;
        if(fDebugLevel>=2) printf("TrackRef %d: x = %f\n", itr, tempTrackRef->LocalX());
        if(fTracklet->GetX0() - tempTrackRef->LocalX() > 3.7) continue;
        if(tempTrackRef->LocalX() - fTracklet->GetX0() > 3.7) break;
        if(fDebugLevel>=2) printf("accepted\n");
        if(nFound == 1)
          if(fTrackRefs[0]->LocalX() >= tempTrackRef->LocalX()) continue;
        fTrackRefs[nFound++] = tempTrackRef;
      }
      if(fDebugLevel>=2) printf("nFound = %d\n", nFound);
      if(nFound < 2) continue;
      // We found 2 track refs for the tracklet, get y and z at the anode wire by a linear approximation
      Double_t dx = fTrackRefs[1]->LocalX() - fTrackRefs[0]->LocalX();
      Double_t dydx = (fTrackRefs[1]->LocalY() - fTrackRefs[0]->LocalY()) / dx;
      Double_t dzdx = (fTrackRefs[1]->Z() - fTrackRefs[0]->Z()) / dx;
      Double_t dx0 = fTrackRefs[1]->LocalX() - fTracklet->GetX0();
      Double_t ymc =  fTrackRefs[1]->LocalY() - dydx*dx0;
      Double_t zmc =  fTrackRefs[1]->Z() - dzdx*dx0;
      
      Double_t dy = fTracklet->GetYfit(0) - ymc;
      Double_t dz = fTracklet->GetZfit(0) - zmc;
      //res_y *= 100; // in mm
      Double_t momentum = fTrackRefs[0]->P();
      
      // Fill Histograms
      if(fDebugLevel>=2) printf("dy = %f\n", dy);
      fYres->Fill(dy);
// 			fZres->Fill(res_z);
/*			fYresLayer[iplane]->Fill(res_y);
      fZresLayer[iplane]->Fill(res_z);*/

//       Double_t phi     = fTrackRefs[0]->Phi();
//       Double_t theta   = fTrackRefs[0]->Theta();
      Double_t phi   = TMath::ATan(dydx);
      Double_t theta = TMath::ATan(dzdx);
      
// 			fYresAngle->Fill(phi, res_y);
// 			fYresLayerAngle[iplane]->Fill(phi, res_y);
      
      Double_t dphi   = TMath::ATan(fTracklet->GetZfit(1)) - phi;
      
// 			fPhiRes->Fill(dphi);
// 			fPhiResLayer[iplane]->Fill(dphi);
// 			fPhiResAngle->Fill(phi, dphi);
// 			fPhiResLayerAngle[iplane]->Fill(phi, dphi);
      
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
  }
  PostData(0, fOutputHistograms);
}

//________________________________________________________
void AliTRDtrackingResolution::Terminate(Option_t *){
  //printf("Tracking Resolution: Terminate\n");
  if(fDebugStream) delete fDebugStream;
}

//________________________________________________________
void AliTRDtrackingResolution::SetDebugLevel(Int_t level){
  fDebugLevel = level;
  if(!fDebugLevel) return;
  if(fDebugStream) return;
  fDebugStream = new TTreeSRedirector("TRD.Resolution.root");
}
