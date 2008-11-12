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

//----------------------------------------------------------------------------
//    Implementation of the single muon analysis class
// An example of usage can be found in the macro runSingleMuAnalysis.C.
//----------------------------------------------------------------------------

#define AliAnalysisTaskSingleMu_cxx

// ROOT includes
#include "Riostream.h"
#include "TChain.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TParticle.h"
#include "TLorentzVector.h"

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

ClassImp(AliAnalysisTaskSingleMu)

//________________________________________________________________________
AliAnalysisTaskSingleMu::AliAnalysisTaskSingleMu(const char *name) :
  AliAnalysisTask(name,""), 
  fAOD(0),
  fOutputContainer(0)
{
  //
  /// Constructor.
  //
  ResetHistos();
  // Input slot #0 works with an Ntuple
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TObjArray container
  DefineOutput(0,  TObjArray::Class());
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

  Int_t ptBins = 60;
  Float_t ptLow = 0., ptHigh = 30.;
  const Char_t *ptName = "P_{t} (GeV/c)";

  Int_t vzBins = 40;
  Float_t vzLow = -20., vzHigh = 20.;
  const Char_t *vzName = "Vz (cm)";

  TString baseName, histoName;
  fOutputContainer = new TObjArray(fgkNhistos*fgkNTrigCuts);
  fOutputContainer->SetName("SingleMuAnalysisContainer");
  Int_t iHisto = 0;

  for(Int_t iTrig=0; iTrig<fgkNTrigCuts; iTrig++){
    
    // 2D histos
    if(!fVzVsPt[iTrig]){
      baseName = "fVzVsPt";
      histoName = baseName + trigName[iTrig];
      fVzVsPt[iTrig] = new TH2F(histoName, histoName,
				ptBins, ptLow, ptHigh,
				vzBins, vzLow, vzHigh);
      fVzVsPt[iTrig]->GetXaxis()->SetTitle(ptName);
      fVzVsPt[iTrig]->GetYaxis()->SetTitle(vzName);

      fOutputContainer->AddAt(fVzVsPt[iTrig], iHisto);
      iHisto++;
    }
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
  TLorentzVector lorVec;
  Int_t trigMatch = -1;

  Int_t nTracks = fAOD->GetNumberOfTracks();
  for (Int_t itrack = 0; itrack < nTracks; itrack++) {
    muonTrack = fAOD->GetTrack(itrack);

    // Apply cuts
    if(!MuonPassesCuts(*muonTrack, lorVec, trigMatch)) continue;

    for(Int_t iTrig=0; iTrig<=trigMatch; iTrig++){
      fVzVsPt[iTrig]->Fill(lorVec.Pt(), fAOD->GetPrimaryVertex()->GetZ());
    }
  }

  // Post final data. It will be written to a file with option "RECREATE"
  PostData(0, fOutputContainer);
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
    c1->Divide(2,2);
    for(Int_t iTrig=0; iTrig<fgkNTrigCuts; iTrig++){
      c1->cd(iTrig+1);
      fVzVsPt[iTrig]->DrawCopy("COLZ");
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskSingleMu::ResetHistos() 
{
  //
  /// Reset histograms
  //
  for(Int_t iTrig=0; iTrig<fgkNTrigCuts; iTrig++){
    fVzVsPt[iTrig] = 0x0;
  }
  trigName[kNoMatchTrig] = "NoMatchTrig";
  trigName[kAllPtTrig]   = "AllPtTrig";
  trigName[kLowPtTrig]   = "LowPtTrig";
  trigName[kHighPtTrig]  = "HighPtTrig";
}


//________________________________________________________________________
Bool_t AliAnalysisTaskSingleMu::MuonPassesCuts(AliAODTrack &muonTrack,
					       TLorentzVector &lorVec,
					       Int_t &trigMatch)
{
  //
  /// Fill lorentz vector and check cuts
  //

  // Check if track is a muon
  if(muonTrack.GetMostProbablePID()!=AliAODTrack::kMuon) return kFALSE;

  // Check if track is triggered
  trigMatch = kNoMatchTrig;
  if (muonTrack.MatchTriggerHighPt()) {
    trigMatch = kHighPtTrig;
  } else if (muonTrack.MatchTriggerLowPt()) {
    trigMatch = kLowPtTrig;
  } else if (muonTrack.MatchTriggerAnyPt()){
    trigMatch = kAllPtTrig;
  }
  
  // Fill track parameters
  Double_t px = muonTrack.Px();
  Double_t py = muonTrack.Py();
  Double_t pz = muonTrack.Pz();
  Double_t p  = muonTrack.P();

  const Double_t kMuonMass = 0.105658369;

  Double_t energy = TMath::Sqrt(p*p + kMuonMass*kMuonMass);
  lorVec.SetPxPyPzE(px,py,pz,energy);

  return kTRUE;
}
