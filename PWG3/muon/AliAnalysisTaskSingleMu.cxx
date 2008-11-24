/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
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

//-----------------------------------------------------------------------------
/// \class AliAnalysisTaskSingleMu
/// Analysis task for single muons in the spectrometer.
/// The output is a tree with:
///  - pt, y and phi of the muon
///  - z position of primary vertex
///  - transverse distance at vertex (DCA)
///  - matched trigger
///
/// \author Diego Stocco
//-----------------------------------------------------------------------------

//----------------------------------------------------------------------------
//    Implementation of the single muon analysis class
// An example of usage can be found in the macro RunSingleMuAnalysisFromAOD.C.
//----------------------------------------------------------------------------

#define AliAnalysisTaskSingleMu_cxx

// ROOT includes
#include "TChain.h"
#include "TROOT.h"
#include "TLorentzVector.h"
#include "TCanvas.h"

// STEER includes
#include "AliLog.h"

#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODVertex.h"
#include "AliAODInputHandler.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSingleMu.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskSingleMu) // Class implementation in ROOT context
/// \endcond

//________________________________________________________________________
AliAnalysisTaskSingleMu::AliAnalysisTaskSingleMu(const char *name) :
  AliAnalysisTask(name,""), 
  fAOD(0),
  fResults(0),
  fVarFloat(0),
  fVarInt(0),
  fFloatVarName(0),
  fIntVarName(0)
{
  //
  /// Constructor.
  //
  InitVariables();
  // Input slot #0 works with an Ntuple
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TTree container
  DefineOutput(0,  TTree::Class());
}

//___________________________________________________________________________
void AliAnalysisTaskSingleMu::ConnectInputData(Option_t *) 
{
  //
  /// Connect AOD here
  /// Called once
  //

  TTree* tree = dynamic_cast<TTree*> (GetInputData(0));
  if (!tree) {
    Printf("ERROR: Could not read chain from input slot 0");
  } else {
    /*
    // Disable all branches and enable only the needed ones
    // The next two lines are different when data produced as AliAODEvent is read
    tree->SetBranchStatus("*", kFALSE);
    tree->SetBranchStatus("fTracks.*", kTRUE);
    */

    AliAODInputHandler *aodH = dynamic_cast<AliAODInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

    if (!aodH) {
      Printf("ERROR: Could not get AODInputHandler");
    } else
      printf("   ConnectInputData of task %s\n", GetName());
    fAOD = aodH->GetEvent();
  }
}

//___________________________________________________________________________
void AliAnalysisTaskSingleMu::CreateOutputObjects() 
{
  //
  /// Create output histograms
  //
  printf("   CreateOutputObjects of task %s\n", GetName());

  // initialize tree
  if(!fResults) fResults = new TTree("Results", "Single mu selection results");

  TString baseName, suffixName;

  suffixName="/F";
  for(Int_t iVar=0; iVar<kNfloatVars; iVar++){
    baseName = fFloatVarName[iVar];
    if(iVar==0) baseName += suffixName;
    fResults->Branch(fFloatVarName[iVar].Data(), &fVarFloat[iVar], baseName.Data());
  }

  suffixName="/I";
  for(Int_t iVar=0; iVar<kNintVars; iVar++){
    baseName = fIntVarName[iVar];
    if(iVar==0) baseName += suffixName;
    fResults->Branch(fIntVarName[iVar].Data(), &fVarInt[iVar], baseName.Data());
  }
}

//________________________________________________________________________
void AliAnalysisTaskSingleMu::Exec(Option_t *) 
{
  //
  /// Main loop
  /// Called for each event
  //

  TTree *tinput = (TTree*)GetInputData(0);
  tinput->GetReadEntry();

  if (!fAOD) {
    Printf("ERROR: fAOD not available");
    return;
  }


  // Object declaration
  AliAODTrack *muonTrack = 0x0;

  Int_t nTracks = fAOD->GetNumberOfTracks();
  for (Int_t itrack = 0; itrack < nTracks; itrack++) {
    muonTrack = fAOD->GetTrack(itrack);

    // Apply cuts
    if(!FillTrackVariables(*muonTrack)) continue;
    fResults->Fill();
  }

  // Post final data. It will be written to a file with option "RECREATE"
  PostData(0, fResults);
}

//________________________________________________________________________
void AliAnalysisTaskSingleMu::Terminate(Option_t *) {
  //
  /// Draw some histogram at the end.
  //
  if (!gROOT->IsBatch()) {
    TCanvas *c1 = new TCanvas("c1","Vz vs Pt",10,10,310,310);
    c1->SetFillColor(10); c1->SetHighLightColor(10);
    c1->SetLeftMargin(0.15); c1->SetBottomMargin(0.15);  
    fResults->Draw("pt:vz","","COLZ");
  }
}

//________________________________________________________________________
void AliAnalysisTaskSingleMu::InitVariables() 
{
  //
  /// Reset histograms
  //

  fVarFloat = new Float_t[kNfloatVars];
  fVarInt = new Int_t[kNintVars];

  fFloatVarName = new TString[kNfloatVars];
  fFloatVarName[kVarPt]     = "pt";
  fFloatVarName[kVarY]      = "y";
  fFloatVarName[kVarPhi]    = "phi";
  fFloatVarName[kVarVz]     = "vz";
  fFloatVarName[kVarDCA]    = "dca";

  fIntVarName = new TString[kNintVars];
  fIntVarName[kVarTrig]     = "matchTrig";
}


//________________________________________________________________________
Bool_t AliAnalysisTaskSingleMu::FillTrackVariables(AliAODTrack &muonTrack)
{
  //
  /// Fill lorentz vector and check cuts
  //

  TLorentzVector lorVec;

  // Check if track is a muon
  if(muonTrack.GetMostProbablePID()!=AliAODTrack::kMuon) return kFALSE;

  // Check if track is triggered
  fVarInt[kVarTrig] = (muonTrack.GetMatchTrigger() && 0x3);
  
  // Fill track parameters
  Double_t px = muonTrack.Px();
  Double_t py = muonTrack.Py();
  Double_t pz = muonTrack.Pz();
  Double_t p  = muonTrack.P();

  const Double_t kMuonMass = 0.105658369;

  Double_t energy = TMath::Sqrt(p*p + kMuonMass*kMuonMass);
  lorVec.SetPxPyPzE(px,py,pz,energy);

  fVarFloat[kVarPt]  = lorVec.Pt();
  fVarFloat[kVarY]   = lorVec.Rapidity();
  fVarFloat[kVarPhi] = lorVec.Phi();

  fVarFloat[kVarVz] = fAOD->GetPrimaryVertex()->GetZ();

  Double_t xDca = muonTrack.XAtDCA();
  Double_t yDca = muonTrack.YAtDCA();

  fVarFloat[kVarDCA] = TMath::Sqrt(xDca*xDca + yDca*yDca);

  return kTRUE;
}
