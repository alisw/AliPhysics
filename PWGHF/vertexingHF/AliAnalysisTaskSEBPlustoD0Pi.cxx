

/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appeuear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

//
//
//   Base class for (BPlus -> D0 pi-> K pi pi) Analysis
//
//
//             Cuts are centralized in AliRDHFCutsBPlustoD0Pi
//             Like sign + track rotation backgrounds are implemented in the macro
//
//-----------------------------------------------------------------------
//
//                 Author Lennart van Doremalen
//           Utrecht University - l.v.r.vandoremalen@uu.nl
//
//-----------------------------------------------------------------------

#include <TSystem.h>
#include <TChain.h>
#include <TParticle.h>
#include "TROOT.h"
#include <TDatabasePDG.h>
#include <AliAnalysisDataSlot.h>
#include <AliAnalysisDataContainer.h>
#include "AliRDHFCutsBPlustoD0Pi.h"
#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliAnalysisManager.h"
#include "AliAODMCHeader.h"
#include "AliAODHandler.h"
#include "AliLog.h"
#include "AliVertex.h"
#include "AliVVertex.h"
#include "AliESDVertex.h"
#include "AliAODVertex.h"
#include "AliVertexerTracks.h"
#include "AliExternalTrackParam.h"
#include "AliNeutralTrackParam.h"
#include "AliAODRecoDecay.h"
#include "AliAODRecoDecayHF.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliAnalysisVertexingHF.h"
#include "AliVertexingHFUtils.h"
#include "AliESDtrack.h"
#include "AliAODMCParticle.h"
#include "AliAODEvent.h"
#include "AliAnalysisTaskSEBPlustoD0Pi.h"
#include "AliAODInputHandler.h"
#include <vector>
#include <TMatrix.h>
#include <TVector3.h>
#include <TArrayI.h>
#include <bitset>
#include <TH3F.h>

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskSEBPlustoD0Pi);
/// \endcond

//__________________________________________________________________________
AliAnalysisTaskSEBPlustoD0Pi::AliAnalysisTaskSEBPlustoD0Pi():
  AliAnalysisTaskSE(),
  fEvents(0),
  fUseMCInfo(kFALSE),
  fShowMask(0),
  fShowRejection(0),
  fQuickSignalAnalysis(0),
  fHistMassWindow(0.125),
  fDegreePerRotation(0),
  fNumberOfRotations(0),
  fPerformCutOptimization(0),
  fRemoveInjected(0),
  fOutput(0),
  fListCuts(0),
  fOutputBPlusResults(0),
  fOutputD0FirstDaughter(0),
  fOutputD0SecondDaughter(0),
  fOutputBPlusPion(0),
  fOutputD0(0),
  fOutputBPlus(0),
  fOutputD0_D0Pt(0),
  fCuts(0),
  fCEvents(0),
  fCTrigger(0),
  fCRejected(0),
  fBPlusPionTracks(0x0),
  fD0Tracks(0x0),
  fnPtBins(0),
  fnPtBinLimits(0),
  fPtBinLimits(0x0),
  fnPtBinsD0forD0ptbin(0),
  fnPtBinsD0forD0ptbinLimits(0),
  fPtBinLimitsD0forD0ptbin(0x0),
  fCutVariableValueArray(),
  fDaughterHistogramArray(),
  fDaughterHistogramArray2D(),
  fDaughterHistogramArrayExtra(),
  fMotherHistogramArray(),
  fMotherHistogramArray2D(),
  fMotherHistogramArrayExtra(),
  fResultsHistogramArray(),
  fResultsHistogramArray2D()
{
  //
  /// Default constructor
  //

}
//___________________________________________________________________________
AliAnalysisTaskSEBPlustoD0Pi::AliAnalysisTaskSEBPlustoD0Pi(const Char_t* name, AliRDHFCutsBPlustoD0Pi* cuts) :
  AliAnalysisTaskSE(name),
  fEvents(0),
  fUseMCInfo(kFALSE),
  fShowMask(0),
  fShowRejection(0),
  fQuickSignalAnalysis(0),
  fHistMassWindow(0.125),
  fDegreePerRotation(0),
  fNumberOfRotations(0),
  fPerformCutOptimization(0),
  fRemoveInjected(0),
  fOutput(0),
  fListCuts(0),
  fOutputBPlusResults(0),
  fOutputD0FirstDaughter(0),
  fOutputD0SecondDaughter(0),
  fOutputBPlusPion(0),
  fOutputD0(0),
  fOutputBPlus(0),
  fOutputD0_D0Pt(0),
  fCuts(0),
  fCEvents(0),
  fCTrigger(0),
  fCRejected(0),
  fBPlusPionTracks(0x0),
  fD0Tracks(0x0),
  fnPtBins(0),
  fnPtBinLimits(0),
  fPtBinLimits(0x0),
  fnPtBinsD0forD0ptbin(0),
  fnPtBinsD0forD0ptbinLimits(0),
  fPtBinLimitsD0forD0ptbin(0x0),
  fCutVariableValueArray(),
  fDaughterHistogramArray(),
  fDaughterHistogramArray2D(),
  fDaughterHistogramArrayExtra(),
  fMotherHistogramArray(),
  fMotherHistogramArray2D(),
  fMotherHistogramArrayExtra(),
  fResultsHistogramArray(),
  fResultsHistogramArray2D()
{
  //
  /// Constructor. Initialization of Inputs and Outputs
  //

  Info("AliAnalysisTaskSEBPlustoD0Pi", "Calling Constructor");

  fCuts = cuts;

  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class()); //conters
  DefineOutput(2, TList::Class());  //My private output
  DefineOutput(3, TList::Class()); // D0 pion output
  DefineOutput(4, TList::Class()); // D0 kaon output
  DefineOutput(5, TList::Class()); // BPlus pion output
  DefineOutput(6, TList::Class()); // D0 output
  DefineOutput(7, TList::Class());  // BPlus output
  DefineOutput(8, TList::Class());  // BPlus output
  DefineOutput(9, TList::Class());  // BPlusMC output
}

//___________________________________________________________________________
AliAnalysisTaskSEBPlustoD0Pi::~AliAnalysisTaskSEBPlustoD0Pi() {
  //
  /// destructor
  //
  Info("~AliAnalysisTaskSEBPlustoD0Pi", "Calling Destructor");

  delete fOutput;
  delete fOutputD0FirstDaughter;
  delete fOutputD0SecondDaughter;
  delete fOutputBPlusPion;
  delete fOutputD0;
  delete fOutputBPlus;
  delete fOutputD0_D0Pt;
  delete fOutputBPlusResults;
  delete fCuts;
  delete fCEvents;
  delete fCTrigger;
  delete fCRejected;
  delete fD0Tracks;
  delete fBPlusPionTracks;
  delete fListCuts;
  delete fPtBinLimits;
  delete fPtBinLimitsD0forD0ptbin;
  for (Int_t i = 0; i < 4; i++) {
    for (Int_t j = 0; j < 6; j++) {
      for (Int_t k = 0; k < 15; k++) {
        delete fDaughterHistogramArray[i][j][k]; fDaughterHistogramArray[i][j][k] = nullptr;
      }
    }
  }
  for (Int_t i = 0; i < 4; i++) {
    for (Int_t j = 0; j < 6; j++) {
      delete fDaughterHistogramArray2D[i][j]; fDaughterHistogramArray2D[i][j] = nullptr;
    }
  }
  for (Int_t i = 0; i < 4; i++) {
    for (Int_t j = 0; j < 6; j++) {
      delete fDaughterHistogramArrayExtra[i][j]; fDaughterHistogramArrayExtra[i][j] = nullptr;
    }
  }
  for (Int_t i = 0; i < 6; i++) {
    for (Int_t j = 0; j < 99; j++) {
      for (Int_t k = 0; k < 60; k++) {
        delete fMotherHistogramArray[i][j][k]; fMotherHistogramArray[i][j][k] = nullptr;
      }
    }
  }
  for (Int_t i = 0; i < 6; i++) {
    for (Int_t j = 0; j < 99; j++) {
      for (Int_t k = 0; k < 60; k++) {
        delete fMotherHistogramArray2D[i][j][k]; fMotherHistogramArray2D[i][j][k] = nullptr;
      }
    }
  }
  for (Int_t i = 0; i < 7; i++) {
    for (Int_t j = 0; j < 10; j++) {
      delete fMotherHistogramArrayExtra[i][j]; fMotherHistogramArrayExtra[i][j] = nullptr;
    }
  }
  for (Int_t i = 0; i < 20; i++) {
    for (Int_t j = 0; j < 99; j++) {
      delete fResultsHistogramArray[i][j]; fResultsHistogramArray[i][j] = nullptr;
    }
  }
  for (Int_t i = 0; i < 20; i++) {
    for (Int_t j = 0; j < 99; j++) {
      delete fResultsHistogramArray2D[i][j]; fResultsHistogramArray2D[i][j] = nullptr;
    }
  }
}
//_________________________________________________
void AliAnalysisTaskSEBPlustoD0Pi::Init() {
  //
  /// Initialization
  //

  if (fDebug > 1) printf("AliAnalysisTaskSEBPlustoD0Pi::Init() \n");

  return;
}
//_________________________________________________
void AliAnalysisTaskSEBPlustoD0Pi::UserExec(Option_t *) {

  //==================================================================================
  //  USER EXECUTION FUNCTION - start
  //==================================================================================
  //
  // This is the main function that has to be altered to do heavy flavour analysis.
  //
  //==================================================================================

  if (!fInputEvent) 
  {
    Error("UserExec", "NO EVENT FOUND!");
    return;
  }

  if (fEvents % 500 == 0) 
  {
    std::cout << "\r" << "Analysing event number: " << fEvents << std::endl;
  }

  fEvents++;

  // Show trigger mask
  if (fShowMask)
  {
    std::bitset<32> maskEV(((AliAODInputHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected());
    std::cout << "Event mask:   " << maskEV << std::endl;
    std::cout << "Trigger mask: " << std::bitset<32>(fCuts->GetTriggerMask()) << std::endl;
  }

  //save trigger information
  for (int i = 0; i < 36; ++i)
  {
    ULong64_t compareValue = BIT(i);
    // std::cout << "Compare mask: " << std::bitset<32>(compareValue) << std::endl;
    Bool_t trigger1 = kFALSE;
    Bool_t trigger2 = kFALSE;

    if (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & compareValue)
    {
      // std::cout << "Match!" << std::endl;
      fCTrigger->Fill(i+1);
      if(i==1) trigger1 = kTRUE;
      if(i==10) trigger2 = kTRUE;
    }
    if(trigger1 && trigger2) fCTrigger->Fill(35);
  }

  

  //==================================================================================
  //  EVENT INITIALIZATION - start
  //==================================================================================


  AliAODEvent* aodEvent = dynamic_cast<AliAODEvent*>(fInputEvent);
  TClonesArray * D0TracksFromFriendFile = 0;
  fCEvents->Fill(1);

  if (!aodEvent && AODEvent() && IsStandardAOD()) 
  {
    // In case there is an AOD handler writing a standard AOD, use the AOD
    // event in memory rather than the input (ESD) event.
    aodEvent = dynamic_cast<AliAODEvent*> (AODEvent());
    // in this case the braches in the deltaAOD (AliAOD.VertexingHF.root)
    // have to be taken from the AOD event hold by the AliAODExtension
    AliAODHandler* aodHandler = (AliAODHandler*)((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
    if (aodHandler->GetExtensions()) 
    {
      AliAODExtension *ext = (AliAODExtension*)aodHandler->GetExtensions()->FindObject("AliAOD.VertexingHF.root");
      AliAODEvent *aodFromExt = ext->GetAOD();
      D0TracksFromFriendFile = (TClonesArray*)aodFromExt->GetList()->FindObject("D0toKpi");
    }
  } else {
    D0TracksFromFriendFile = (TClonesArray*)aodEvent->GetList()->FindObject("D0toKpi");
  }

  // fix for temporary bug in ESDfilter
  // the AODs with null vertex pointer didn't pass the PhysSel
  if (!aodEvent->GetPrimaryVertex() || TMath::Abs(aodEvent->GetMagneticField()) < 0.001) return;
  fCEvents->Fill(2);

  // trigger class for PbPb C0SMH-B-NOPF-ALLNOTRD
  TString trigclass = aodEvent->GetFiredTriggerClasses();
  if (trigclass.Contains("C0SMH-B-NOPF-ALLNOTRD") || trigclass.Contains("C0SMH-B-NOPF-ALL")) fCEvents->Fill(5);

  fCRejected->Fill(0);
  if (!fCuts->IsEventSelected(aodEvent))
  {
    if (fShowRejection) std::cout << "Event rejected by code: " << fCuts->GetWhyRejection() << std::endl;
    fCRejected->Fill(fCuts->GetWhyRejection() + 1);
    return;
  }

  fCEvents->Fill(3);

  //get the magnetic field
  Double_t bz = (Double_t)aodEvent->GetMagneticField();

  // AOD primary vertex
  AliAODVertex *primaryVertex = (AliAODVertex*)aodEvent->GetPrimaryVertex();
  if (!primaryVertex) return;
  if (primaryVertex->GetNContributors() < 1) return;
  fCEvents->Fill(4);

  if (!D0TracksFromFriendFile)
  {
    AliInfo("Could not find array of HF vertices, skipping the event");
    return;
  }
  else AliDebug(2, Form("Found %d vertices", D0TracksFromFriendFile->GetEntriesFast()));

  AliAODMCHeader *mcHeader = 0;
  if (fUseMCInfo)
  {
    mcHeader = (AliAODMCHeader*)aodEvent->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if (!mcHeader) {
      printf(" MC header branch not found!\n");
      return;
    }
  }


  //==================================================================================
  //  EVENT INITIALIZATION - end
  //==================================================================================
  //  GET TRUE MC INFO - start
  //==================================================================================

  TClonesArray *mcTrackArray = nullptr;
  if (fUseMCInfo) mcTrackArray = dynamic_cast<TClonesArray*>(aodEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  if (fUseMCInfo && !mcTrackArray) return;

  //Here we create an array and vectors that we use to look how many BPluss can be reconstructed after each cut
  TMatrix * BPlustoD0PiLabelMatrix = new TMatrix(0, 5);

  //we fill the array with all BPlus->D0Pi tracks
  if (fUseMCInfo) {
    BPlustoD0PiSignalTracksInMC(mcTrackArray, aodEvent, BPlustoD0PiLabelMatrix, fOutputBPlusResults);
  }

  //==================================================================================
  //  GET TRUE MC INFO - end
  //==================================================================================
  //  PARTICLE SELECTION LOOP - start
  //==================================================================================
  //
  // Here we select and reconstruct the particles for the BPlus -> D0 + Pion decay.
  //
  //==================================================================================

  D0Selection(aodEvent, primaryVertex, bz, mcTrackArray, BPlustoD0PiLabelMatrix, D0TracksFromFriendFile, mcHeader);

  if((Int_t)fD0Tracks->size() > 0)
  {
    BPlusPionSelection(aodEvent, primaryVertex, bz, mcTrackArray, BPlustoD0PiLabelMatrix, mcHeader);

    if((Int_t)fBPlusPionTracks->size() > 0)
    {
      BPlusSelection(aodEvent, primaryVertex, bz, mcTrackArray, BPlustoD0PiLabelMatrix, D0TracksFromFriendFile, mcHeader);
    }
  }

  // Clear arrays and memory management:
  fD0Tracks->erase(fD0Tracks->begin(), fD0Tracks->end());
  fBPlusPionTracks->erase(fBPlusPionTracks->begin(), fBPlusPionTracks->end());

  delete BPlustoD0PiLabelMatrix; BPlustoD0PiLabelMatrix = NULL;

  //==================================================================================
  //  PARTICLE SELECTION LOOP - end
  //==================================================================================

  PostData(1, fOutput);
  PostData(3, fOutputD0FirstDaughter);
  PostData(4, fOutputD0SecondDaughter);
  PostData(5, fOutputBPlusPion);
  PostData(6, fOutputD0);
  PostData(7, fOutputBPlus);
  PostData(8, fOutputD0_D0Pt);
  PostData(9, fOutputBPlusResults);


  //==================================================================================
  //  USER EXECUTION FUNCTION - end
  //==================================================================================

}
//___________________________________ terminate ___________________________
void AliAnalysisTaskSEBPlustoD0Pi::Terminate(Option_t*) {
  /// The Terminate() function is the last function to be called during
  /// a query. It always runs on the client, it can be used to present
  /// the results graphically or save the results to file.

  //Info("Terminate","");
  AliAnalysisTaskSE::Terminate();

  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput) {
    printf("ERROR: fOutput not available\n");
    return;
  }

  fCEvents        = dynamic_cast<TH1F*>(fOutput->FindObject("fCEvents"));
  fCTrigger       = dynamic_cast<TH1F*>(fOutput->FindObject("fCTrigger"));
  fCRejected      = dynamic_cast<TH1F*>(fOutput->FindObject("fCRejected"));


  fListCuts = dynamic_cast<TList*> (GetOutputData(2));
  if (!fListCuts) {
    printf("ERROR: fListCuts not available\n");
    return;
  }
  fOutputD0FirstDaughter = dynamic_cast<TList*> (GetOutputData(3));
  if (!fOutputD0FirstDaughter) {
    printf("ERROR: fOutputD0FirstDaughter not available\n");
    return;
  }
  fOutputD0SecondDaughter = dynamic_cast<TList*> (GetOutputData(4));
  if (!fOutputD0SecondDaughter) {
    printf("ERROR: fOutputD0SecondDaughter not available\n");
    return;
  }
  fOutputBPlusPion = dynamic_cast<TList*> (GetOutputData(5));
  if (!fOutputBPlusPion) {
    printf("ERROR: fOutputBPlusPion not available\n");
    return;
  }
  fOutputD0 = dynamic_cast<TList*> (GetOutputData(6));
  if (!fOutputD0) {
    printf("ERROR: fOutputD0 not available\n");
    return;
  }
  fOutputBPlus = dynamic_cast<TList*> (GetOutputData(7));
  if (!fOutputBPlus) {
    printf("ERROR: fOutputBPlus not available\n");
    return;
  }
  fOutputD0_D0Pt = dynamic_cast<TList*> (GetOutputData(8));
  if (!fOutputD0_D0Pt) {
    printf("ERROR: fOutputD0_D0Pt not available\n");
    return;
  }
  fOutputBPlusResults = dynamic_cast<TList*> (GetOutputData(9));
  if (!fOutputBPlusResults) {
    printf("ERROR: fOutputBPlusResults not available\n");
    return;
  }
  return;
}

//___________________________________________________________________________
void AliAnalysisTaskSEBPlustoD0Pi::UserCreateOutputObjects() {
  /// output
  Info("UserCreateOutputObjects", "CreateOutputObjects of task %s\n", GetName());

  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("General_information");

  fOutputD0FirstDaughter = new TList();
  fOutputD0FirstDaughter->SetOwner();
  fOutputD0FirstDaughter->SetName("listD0Pion");

  fOutputD0SecondDaughter = new TList();
  fOutputD0SecondDaughter->SetOwner();
  fOutputD0SecondDaughter->SetName("listD0Kaon");

  fOutputBPlusPion = new TList();
  fOutputBPlusPion->SetOwner();
  fOutputBPlusPion->SetName("listBPlusPion");

  fOutputD0 = new TList();
  fOutputD0->SetOwner();
  fOutputD0->SetName("listD0");

  fOutputBPlus = new TList();
  fOutputBPlus->SetOwner();
  fOutputBPlus->SetName("listBPlus");

  fOutputD0_D0Pt = new TList();
  fOutputD0_D0Pt->SetOwner();
  fOutputD0_D0Pt->SetName("listD0_D0Pt");

  fOutputBPlusResults = new TList();
  fOutputBPlusResults->SetOwner();
  fOutputBPlusResults->SetName("listBPlusResults");

  // we prepare vectors that will save the positions of the daughter tracks in the track list during the reconstruction
  fBPlusPionTracks = new std::vector<Int_t>;
  fD0Tracks = new std::vector<Int_t>;

  // we get information on the pt bins
  fnPtBins = fCuts->GetNPtBins();
  fnPtBinLimits = fnPtBins + 1;
  fPtBinLimits = fCuts->GetPtBinLimits();

  fnPtBinsD0forD0ptbin = fCuts->GetNPtBinsD0forD0ptbin();
  fnPtBinsD0forD0ptbinLimits = fnPtBinsD0forD0ptbin + 1;
  fPtBinLimitsD0forD0ptbin = fCuts->GetPtBinLimitsD0forD0ptbin();

  std::cout << "Nr. of BPlus meson bins: " <<  fCuts->GetNPtBins() << " limits: " << std::endl;
  for (int i = 0; i < fnPtBinLimits; ++i)
  {
    std::cout << fPtBinLimits[i] << " " << std::endl;
  }
  std::cout << std::endl;
  std::cout << "Nr. of D0 meson bins: " <<  fCuts->GetNPtBinsD0forD0ptbin() << " limits: " << std::endl;
  for (int i = 0; i < fnPtBinsD0forD0ptbinLimits; ++i)
  {
    std::cout << fPtBinLimitsD0forD0ptbin[i] << " " << std::endl;
  }
  std::cout << std::endl;

  fListCuts = new TList();
  fListCuts->SetOwner();
  fListCuts->SetName("Cuts");
  AliRDHFCutsBPlustoD0Pi* copyfCuts = new AliRDHFCutsBPlustoD0Pi(*fCuts);
  fListCuts->Add(copyfCuts);

  // define histograms
  DefineHistograms();

  // Post the data
  PostData(1, fOutput);
  PostData(2, fListCuts);
  PostData(3, fOutputD0FirstDaughter);
  PostData(4, fOutputD0SecondDaughter);
  PostData(5, fOutputBPlusPion);
  PostData(6, fOutputD0);
  PostData(7, fOutputBPlus);
  PostData(8, fOutputD0_D0Pt);
  PostData(9, fOutputBPlusResults);

  return;
}
//___________________________________ histograms _______________________________________
void  AliAnalysisTaskSEBPlustoD0Pi::DefineHistograms() {

  /// Create histograms

  fCEvents = new TH1F("fCEvents", "Event counter", 13, 0, 13);
  fCEvents->SetStats(kTRUE);
  fCEvents->GetXaxis()->SetTitle("");
  fCEvents->GetYaxis()->SetTitle("Number");
  fCEvents->GetXaxis()->SetBinLabel(2, "total nr. of events");
  fCEvents->GetXaxis()->SetBinLabel(3, "good prim vtx and B field");
  fCEvents->GetXaxis()->SetBinLabel(4, "total events selected");
  fCEvents->GetXaxis()->SetBinLabel(5, "nr. vtx contributors");
  fCEvents->GetXaxis()->SetBinLabel(6, "trigger for PbPb");
  fCEvents->GetXaxis()->SetBinLabel(7, "no z vtx");
  fCEvents->GetXaxis()->SetBinLabel(8, "");
  fCEvents->GetXaxis()->SetBinLabel(9, "");
  fCEvents->GetXaxis()->SetBinLabel(10, "");
  fCEvents->GetXaxis()->SetBinLabel(11, "");
  fCEvents->GetXaxis()->SetBinLabel(12, "nr. of D0 fail to be rec");
  fOutput->Add(fCEvents);

  fCTrigger= new TH1F("fCTrigger", "Trigger counter", 36, 0, 36);
  fCTrigger->SetStats(kTRUE);
  fCTrigger->GetXaxis()->SetTitle("");
  fCTrigger->GetYaxis()->SetTitle("Number");
  fCTrigger->GetXaxis()->SetBinLabel(1, "total nr. of events");
  fCTrigger->GetXaxis()->SetBinLabel(2, "Bit 0 - kMB/kINT1");
  fCTrigger->GetXaxis()->SetBinLabel(3, "Bit 1 - kINT7");
  fCTrigger->GetXaxis()->SetBinLabel(4, "Bit 2 - kMUON");
  fCTrigger->GetXaxis()->SetBinLabel(5, "Bit 3 - kHighMult/kHighMultSPD");
  fCTrigger->GetXaxis()->SetBinLabel(6, "Bit 4 - kEMC1");
  fCTrigger->GetXaxis()->SetBinLabel(7, "Bit 5 - kCINT5/kINT5");
  fCTrigger->GetXaxis()->SetBinLabel(8, "Bit 6 - kCMUS5/kMUSPB/kINT7inMUON");
  fCTrigger->GetXaxis()->SetBinLabel(9, "Bit 7 - kMuonSingleHighPt7/kMUSH7/kMUSHPB");
  fCTrigger->GetXaxis()->SetBinLabel(10, "Bit 8 - kMuonLikeLowPt7/kMUL7/kMuonLikePB");
  fCTrigger->GetXaxis()->SetBinLabel(11, "Bit 9 - kMuonUnlikeLowPt7/kMUU7/kMuonUnlikePB");
  fCTrigger->GetXaxis()->SetBinLabel(12, "Bit 10 - kEMC7/kEMC8");
  fCTrigger->GetXaxis()->SetBinLabel(13, "Bit 11 - kMUS7/kMuonSingleLowPt7");
  fCTrigger->GetXaxis()->SetBinLabel(14, "Bit 12 - kPHI1");
  fCTrigger->GetXaxis()->SetBinLabel(15, "Bit 13 - kPHI7/kPHI8/kPHOSPb");
  fCTrigger->GetXaxis()->SetBinLabel(16, "Bit 14 - kEMCEJE");
  fCTrigger->GetXaxis()->SetBinLabel(17, "Bit 15 - kEMCEGA");
  fCTrigger->GetXaxis()->SetBinLabel(18, "Bit 16 - kHighMultV0/kCentral");
  fCTrigger->GetXaxis()->SetBinLabel(19, "Bit 17 - kSemiCentral");
  fCTrigger->GetXaxis()->SetBinLabel(20, "Bit 18 - kDG/kDG5");
  fCTrigger->GetXaxis()->SetBinLabel(21, "Bit 19 - kZED");
  fCTrigger->GetXaxis()->SetBinLabel(22, "Bit 20 - kSPI7/kSPI");
  fCTrigger->GetXaxis()->SetBinLabel(23, "Bit 21 - kINT8");
  fCTrigger->GetXaxis()->SetBinLabel(24, "Bit 22 - kMuonSingleLowPt8");
  fCTrigger->GetXaxis()->SetBinLabel(25, "Bit 23 - kMuonSingleHighPt8");
  fCTrigger->GetXaxis()->SetBinLabel(26, "Bit 24 - kMuonLikeLowPt8");
  fCTrigger->GetXaxis()->SetBinLabel(27, "Bit 25 - kMuonUnlikeLowPt8");
  fCTrigger->GetXaxis()->SetBinLabel(28, "Bit 26 - kMuonUnlikeLowPt0");
  fCTrigger->GetXaxis()->SetBinLabel(29, "Bit 27 - kUserDefined");
  fCTrigger->GetXaxis()->SetBinLabel(30, "Bit 28 - kTRD");
  fCTrigger->GetXaxis()->SetBinLabel(31, "Bit 29 - ...");
  fCTrigger->GetXaxis()->SetBinLabel(32, "Bit 30 - kFastOnly");
  fCTrigger->GetXaxis()->SetBinLabel(33, "Bit 31 - ...");
  fCTrigger->GetXaxis()->SetBinLabel(34, "Bit 32 - ...");
  fCTrigger->GetXaxis()->SetBinLabel(35, "kINT7 && kEMC7");
  fOutput->Add(fCTrigger);

  fCRejected= new TH1F("fCRejected", "Rejected counter", 35, 0, 35);
  fCRejected->SetStats(kTRUE);
  fCRejected->GetXaxis()->SetTitle("");
  fCRejected->GetYaxis()->SetTitle("Number");
  fCRejected->GetXaxis()->SetBinLabel(1, "total nr. of events");
  fCRejected->GetXaxis()->SetBinLabel(2, "0");
  fCRejected->GetXaxis()->SetBinLabel(3, "1");
  fCRejected->GetXaxis()->SetBinLabel(4, "2");
  fCRejected->GetXaxis()->SetBinLabel(5, "3");
  fCRejected->GetXaxis()->SetBinLabel(6, "4");
  fCRejected->GetXaxis()->SetBinLabel(7, "5");
  fCRejected->GetXaxis()->SetBinLabel(8, "6");
  fCRejected->GetXaxis()->SetBinLabel(9, "7");
  fCRejected->GetXaxis()->SetBinLabel(10, "8");
  fCRejected->GetXaxis()->SetBinLabel(11, "9");
  fCRejected->GetXaxis()->SetBinLabel(12, "10");
  fCRejected->GetXaxis()->SetBinLabel(13, "11");
  fCRejected->GetXaxis()->SetBinLabel(14, "12");
  fCRejected->GetXaxis()->SetBinLabel(15, "13");
  fCRejected->GetXaxis()->SetBinLabel(16, "14");
  fCRejected->GetXaxis()->SetBinLabel(17, "15");
  fCRejected->GetXaxis()->SetBinLabel(18, "16");
  fCRejected->GetXaxis()->SetBinLabel(19, "17");
  fCRejected->GetXaxis()->SetBinLabel(20, "18");
  fCRejected->GetXaxis()->SetBinLabel(21, "19");
  fCRejected->GetXaxis()->SetBinLabel(22, "20");
  fCRejected->GetXaxis()->SetBinLabel(23, "21");
  fCRejected->GetXaxis()->SetBinLabel(24, "22");
  fCRejected->GetXaxis()->SetBinLabel(25, "23");
  fCRejected->GetXaxis()->SetBinLabel(26, "24");
  fCRejected->GetXaxis()->SetBinLabel(27, "25");
  fCRejected->GetXaxis()->SetBinLabel(28, "26");
  fCRejected->GetXaxis()->SetBinLabel(29, "27");
  fCRejected->GetXaxis()->SetBinLabel(30, "28");
  fCRejected->GetXaxis()->SetBinLabel(31, "29");
  fCRejected->GetXaxis()->SetBinLabel(32, "30");
  fCRejected->GetXaxis()->SetBinLabel(33, "31");
  fCRejected->GetXaxis()->SetBinLabel(34, "32");
  fOutput->Add(fCRejected);

  //====================================================

  TString name_mc_BPlus_pt = "mc_BPlus_pt";
  TH1F* hist_mc_BPlus_pt = new TH1F(name_mc_BPlus_pt.Data(), "Pt monte carlo BPlus in BPlus->D*#pi; p_{T} [GeV/c]; Entries", 400, 0, 20);
  hist_mc_BPlus_pt->Sumw2();
  hist_mc_BPlus_pt->SetLineColor(6);
  hist_mc_BPlus_pt->SetMarkerStyle(20);
  hist_mc_BPlus_pt->SetMarkerSize(0.6);
  hist_mc_BPlus_pt->SetMarkerColor(6);
  TH1F* histogram_mc_BPlus_pt = (TH1F*)hist_mc_BPlus_pt->Clone();
  fOutputBPlusResults->Add(histogram_mc_BPlus_pt);
  fResultsHistogramArray[0][0] = histogram_mc_BPlus_pt;

  TString name_mc_BPlus_pion_pt = "mc_BPlus_pion_pt";
  TH1F* hist_mc_BPlus_pion_pt = new TH1F(name_mc_BPlus_pion_pt.Data(), "Pt monte carlo pion of BPlus in BPlus->D*#pi; p_{T} [GeV/c]; Entries", 400, 0, 20);
  hist_mc_BPlus_pion_pt->Sumw2();
  hist_mc_BPlus_pion_pt->SetLineColor(6);
  hist_mc_BPlus_pion_pt->SetMarkerStyle(20);
  hist_mc_BPlus_pion_pt->SetMarkerSize(0.6);
  hist_mc_BPlus_pion_pt->SetMarkerColor(6);
  TH1F* histogram_mc_BPlus_pion_pt = (TH1F*)hist_mc_BPlus_pion_pt->Clone();
  fOutputBPlusResults->Add(histogram_mc_BPlus_pion_pt);
  fResultsHistogramArray[0][1] = histogram_mc_BPlus_pion_pt;

  TString name_mc_D0_pt = "mc_D0_pt";
  TH1F* hist_mc_D0_pt = new TH1F(name_mc_D0_pt.Data(), "Pt monte carlo D0 in BPlus->D*#pi; p_{T} [GeV/c]; Entries", 400, 0, 20);
  hist_mc_D0_pt->Sumw2();
  hist_mc_D0_pt->SetLineColor(6);
  hist_mc_D0_pt->SetMarkerStyle(20);
  hist_mc_D0_pt->SetMarkerSize(0.6);
  hist_mc_D0_pt->SetMarkerColor(6);
  TH1F* histogram_mc_D0_pt = (TH1F*)hist_mc_D0_pt->Clone();
  fOutputBPlusResults->Add(histogram_mc_D0_pt);
  fResultsHistogramArray[0][2] = histogram_mc_D0_pt;

  TString name_mc_D0_pion_pt = "mc_D0_pion_pt";
  TH1F* hist_mc_D0_pion_pt = new TH1F(name_mc_D0_pion_pt.Data(), "Pt monte carlo pion of D0 in BPlus->D*#pi; p_{T} [GeV/c]; Entries", 400, 0, 20);
  hist_mc_D0_pion_pt->Sumw2();
  hist_mc_D0_pion_pt->SetLineColor(6);
  hist_mc_D0_pion_pt->SetMarkerStyle(20);
  hist_mc_D0_pion_pt->SetMarkerSize(0.6);
  hist_mc_D0_pion_pt->SetMarkerColor(6);
  TH1F* histogram_mc_D0_pion_pt = (TH1F*)hist_mc_D0_pion_pt->Clone();
  fOutputBPlusResults->Add(histogram_mc_D0_pion_pt);
  fResultsHistogramArray[0][3] = histogram_mc_D0_pion_pt;

  TString name_mc_D0_kaon_pt = "mc_D0_kaon_pt";
  TH1F* hist_mc_D0_kaon_pt = new TH1F(name_mc_D0_kaon_pt.Data(), "Pt monte carlo kaon of D0 in BPlus->D*#pi; p_{T} [GeV/c]; Entries", 400, 0, 20);
  hist_mc_D0_kaon_pt->Sumw2();
  hist_mc_D0_kaon_pt->SetLineColor(6);
  hist_mc_D0_kaon_pt->SetMarkerStyle(20);
  hist_mc_D0_kaon_pt->SetMarkerSize(0.6);
  hist_mc_D0_kaon_pt->SetMarkerColor(6);
  TH1F* histogram_mc_D0_kaon_pt = (TH1F*)hist_mc_D0_kaon_pt->Clone();
  fOutputBPlusResults->Add(histogram_mc_D0_kaon_pt);
  fResultsHistogramArray[0][4] = histogram_mc_D0_kaon_pt;

  TString name_mc_BPlus_rapidity_true = "mc_BPlus_rapidity_true";
  TH1F* hist_mc_BPlus_rapidity_true = new TH1F(name_mc_BPlus_rapidity_true.Data(), "rapidity_true monte carlo BPlus in BPlus->D*#pi; Y; Entries", 5000, -20, 20);
  hist_mc_BPlus_rapidity_true->Sumw2();
  hist_mc_BPlus_rapidity_true->SetLineColor(6);
  hist_mc_BPlus_rapidity_true->SetMarkerStyle(20);
  hist_mc_BPlus_rapidity_true->SetMarkerSize(0.6);
  hist_mc_BPlus_rapidity_true->SetMarkerColor(6);
  TH1F* histogram_mc_BPlus_rapidity_true = (TH1F*)hist_mc_BPlus_rapidity_true->Clone();
  fOutputBPlusResults->Add(histogram_mc_BPlus_rapidity_true);
  fResultsHistogramArray[0][5] = histogram_mc_BPlus_rapidity_true;

  TString name_mc_BPlus_pion_rapidity_true = "mc_BPlus_pion_rapidity_true";
  TH1F* hist_mc_BPlus_pion_rapidity_true = new TH1F(name_mc_BPlus_pion_rapidity_true.Data(), "rapidity_true monte carlo pion of BPlus in BPlus->D*#pi; Y; Entries", 5000, -20, 20);
  hist_mc_BPlus_pion_rapidity_true->Sumw2();
  hist_mc_BPlus_pion_rapidity_true->SetLineColor(6);
  hist_mc_BPlus_pion_rapidity_true->SetMarkerStyle(20);
  hist_mc_BPlus_pion_rapidity_true->SetMarkerSize(0.6);
  hist_mc_BPlus_pion_rapidity_true->SetMarkerColor(6);
  TH1F* histogram_mc_BPlus_pion_rapidity_true = (TH1F*)hist_mc_BPlus_pion_rapidity_true->Clone();
  fOutputBPlusResults->Add(histogram_mc_BPlus_pion_rapidity_true);
  fResultsHistogramArray[0][6] = histogram_mc_BPlus_pion_rapidity_true;

  TString name_mc_D0_rapidity_true = "mc_D0_rapidity_true";
  TH1F* hist_mc_D0_rapidity_true = new TH1F(name_mc_D0_rapidity_true.Data(), "rapidity_true monte carlo D0 in BPlus->D*#pi; Y; Entries", 5000, -20, 20);
  hist_mc_D0_rapidity_true->Sumw2();
  hist_mc_D0_rapidity_true->SetLineColor(6);
  hist_mc_D0_rapidity_true->SetMarkerStyle(20);
  hist_mc_D0_rapidity_true->SetMarkerSize(0.6);
  hist_mc_D0_rapidity_true->SetMarkerColor(6);
  TH1F* histogram_mc_D0_rapidity_true = (TH1F*)hist_mc_D0_rapidity_true->Clone();
  fOutputBPlusResults->Add(histogram_mc_D0_rapidity_true);
  fResultsHistogramArray[0][7] = histogram_mc_D0_rapidity_true;

  TString name_mc_D0_pion_rapidity_true = "mc_D0_pion_rapidity_true";
  TH1F* hist_mc_D0_pion_rapidity_true = new TH1F(name_mc_D0_pion_rapidity_true.Data(), "rapidity_true monte carlo pion of D0 in BPlus->D*#pi; Y; Entries", 5000, -20, 20);
  hist_mc_D0_pion_rapidity_true->Sumw2();
  hist_mc_D0_pion_rapidity_true->SetLineColor(6);
  hist_mc_D0_pion_rapidity_true->SetMarkerStyle(20);
  hist_mc_D0_pion_rapidity_true->SetMarkerSize(0.6);
  hist_mc_D0_pion_rapidity_true->SetMarkerColor(6);
  TH1F* histogram_mc_D0_pion_rapidity_true = (TH1F*)hist_mc_D0_pion_rapidity_true->Clone();
  fOutputBPlusResults->Add(histogram_mc_D0_pion_rapidity_true);
  fResultsHistogramArray[0][8] = histogram_mc_D0_pion_rapidity_true;

  TString name_mc_D0_kaon_rapidity_true = "mc_D0_kaon_rapidity_true";
  TH1F* hist_mc_D0_kaon_rapidity_true = new TH1F(name_mc_D0_kaon_rapidity_true.Data(), "rapidity_true monte carlo kaon of D0 in BPlus->D*#pi; Y; Entries", 5000, -20, 20);
  hist_mc_D0_kaon_rapidity_true->Sumw2();
  hist_mc_D0_kaon_rapidity_true->SetLineColor(6);
  hist_mc_D0_kaon_rapidity_true->SetMarkerStyle(20);
  hist_mc_D0_kaon_rapidity_true->SetMarkerSize(0.6);
  hist_mc_D0_kaon_rapidity_true->SetMarkerColor(6);
  TH1F* histogram_mc_D0_kaon_rapidity_true = (TH1F*)hist_mc_D0_kaon_rapidity_true->Clone();
  fOutputBPlusResults->Add(histogram_mc_D0_kaon_rapidity_true);
  fResultsHistogramArray[0][9] = histogram_mc_D0_kaon_rapidity_true;

  TString name_mc_BPlus_pseudorapidity_true = "mc_BPlus_pseudorapidity_true";
  TH1F* hist_mc_BPlus_pseudorapidity_true = new TH1F(name_mc_BPlus_pseudorapidity_true.Data(), "pseudorapidity_true monte carlo BPlus in BPlus->D*#pi; #eta; Entries", 5000, -20, 20);
  hist_mc_BPlus_pseudorapidity_true->Sumw2();
  hist_mc_BPlus_pseudorapidity_true->SetLineColor(6);
  hist_mc_BPlus_pseudorapidity_true->SetMarkerStyle(20);
  hist_mc_BPlus_pseudorapidity_true->SetMarkerSize(0.6);
  hist_mc_BPlus_pseudorapidity_true->SetMarkerColor(6);
  TH1F* histogram_mc_BPlus_pseudorapidity_true = (TH1F*)hist_mc_BPlus_pseudorapidity_true->Clone();
  fOutputBPlusResults->Add(histogram_mc_BPlus_pseudorapidity_true);
  fResultsHistogramArray[0][10] = histogram_mc_BPlus_pseudorapidity_true;

  TString name_mc_BPlus_pion_pseudorapidity_true = "mc_BPlus_pion_pseudorapidity_true";
  TH1F* hist_mc_BPlus_pion_pseudorapidity_true = new TH1F(name_mc_BPlus_pion_pseudorapidity_true.Data(), "pseudorapidity_true monte carlo pion of BPlus in BPlus->D*#pi; #eta; Entries", 5000, -20, 20);
  hist_mc_BPlus_pion_pseudorapidity_true->Sumw2();
  hist_mc_BPlus_pion_pseudorapidity_true->SetLineColor(6);
  hist_mc_BPlus_pion_pseudorapidity_true->SetMarkerStyle(20);
  hist_mc_BPlus_pion_pseudorapidity_true->SetMarkerSize(0.6);
  hist_mc_BPlus_pion_pseudorapidity_true->SetMarkerColor(6);
  TH1F* histogram_mc_BPlus_pion_pseudorapidity_true = (TH1F*)hist_mc_BPlus_pion_pseudorapidity_true->Clone();
  fOutputBPlusResults->Add(histogram_mc_BPlus_pion_pseudorapidity_true);
  fResultsHistogramArray[0][11] = histogram_mc_BPlus_pion_pseudorapidity_true;

  TString name_mc_D0_pseudorapidity_true = "mc_D0_pseudorapidity_true";
  TH1F* hist_mc_D0_pseudorapidity_true = new TH1F(name_mc_D0_pseudorapidity_true.Data(), "pseudorapidity_true monte carlo D0 in BPlus->D*#pi; #eta; Entries", 5000, -20, 20);
  hist_mc_D0_pseudorapidity_true->Sumw2();
  hist_mc_D0_pseudorapidity_true->SetLineColor(6);
  hist_mc_D0_pseudorapidity_true->SetMarkerStyle(20);
  hist_mc_D0_pseudorapidity_true->SetMarkerSize(0.6);
  hist_mc_D0_pseudorapidity_true->SetMarkerColor(6);
  TH1F* histogram_mc_D0_pseudorapidity_true = (TH1F*)hist_mc_D0_pseudorapidity_true->Clone();
  fOutputBPlusResults->Add(histogram_mc_D0_pseudorapidity_true);
  fResultsHistogramArray[0][12] = histogram_mc_D0_pseudorapidity_true;

  TString name_mc_D0_pion_pseudorapidity_true = "mc_D0_pion_pseudorapidity_true";
  TH1F* hist_mc_D0_pion_pseudorapidity_true = new TH1F(name_mc_D0_pion_pseudorapidity_true.Data(), "pseudorapidity_true monte carlo pion of D0 in BPlus->D*#pi; #eta; Entries", 5000, -20, 20);
  hist_mc_D0_pion_pseudorapidity_true->Sumw2();
  hist_mc_D0_pion_pseudorapidity_true->SetLineColor(6);
  hist_mc_D0_pion_pseudorapidity_true->SetMarkerStyle(20);
  hist_mc_D0_pion_pseudorapidity_true->SetMarkerSize(0.6);
  hist_mc_D0_pion_pseudorapidity_true->SetMarkerColor(6);
  TH1F* histogram_mc_D0_pion_pseudorapidity_true = (TH1F*)hist_mc_D0_pion_pseudorapidity_true->Clone();
  fOutputBPlusResults->Add(histogram_mc_D0_pion_pseudorapidity_true);
  fResultsHistogramArray[0][13] = histogram_mc_D0_pion_pseudorapidity_true;

  TString name_mc_D0_kaon_pseudorapidity_true = "mc_D0_kaon_pseudorapidity_true";
  TH1F* hist_mc_D0_kaon_pseudorapidity_true = new TH1F(name_mc_D0_kaon_pseudorapidity_true.Data(), "pseudorapidity_true monte carlo kaon of D0 in BPlus->D*#pi; #eta; Entries", 5000, -20, 20);
  hist_mc_D0_kaon_pseudorapidity_true->Sumw2();
  hist_mc_D0_kaon_pseudorapidity_true->SetLineColor(6);
  hist_mc_D0_kaon_pseudorapidity_true->SetMarkerStyle(20);
  hist_mc_D0_kaon_pseudorapidity_true->SetMarkerSize(0.6);
  hist_mc_D0_kaon_pseudorapidity_true->SetMarkerColor(6);
  TH1F* histogram_mc_D0_kaon_pseudorapidity_true = (TH1F*)hist_mc_D0_kaon_pseudorapidity_true->Clone();
  fOutputBPlusResults->Add(histogram_mc_D0_kaon_pseudorapidity_true);
  fResultsHistogramArray[0][14] = histogram_mc_D0_kaon_pseudorapidity_true;

  if (fPerformCutOptimization)
  {
    Int_t nCuts = fCuts->GetnCutsForOptimization();
    Int_t nVariables = fCuts->GetnVariablesForCutOptimization();
    Int_t nCutOptimizationBins = TMath::Power(nCuts, nVariables);
    Int_t nSigmaBins = fCuts->GetNumberOfSigmaBinsForCutOptimization();

    for (Int_t k = 0; k < fnPtBins; ++k)
    {
      TString ptBinMother = "";
      ptBinMother += "_ptbin_";
      ptBinMother += fPtBinLimits[k];
      ptBinMother += "_to_";
      ptBinMother += fPtBinLimits[k + 1];

      TString name_cut_optimization_signal = "cut_optimization_signal";
      name_cut_optimization_signal += ptBinMother;
      TH2F* hist_cut_optimization_signal = new TH2F(name_cut_optimization_signal.Data(), "Total signal for different cuts; Cut number; Entries", nCutOptimizationBins, 0, nCutOptimizationBins, 2 * nSigmaBins, -nSigmaBins, nSigmaBins);
      hist_cut_optimization_signal->Sumw2();
      hist_cut_optimization_signal->SetLineColor(6);
      hist_cut_optimization_signal->SetMarkerStyle(20);
      hist_cut_optimization_signal->SetMarkerSize(0.6);
      hist_cut_optimization_signal->SetMarkerColor(6);
      TH2F* histogram_cut_optimization_signal = (TH2F*)hist_cut_optimization_signal->Clone();
      fOutputBPlusResults->Add(histogram_cut_optimization_signal);
      fResultsHistogramArray2D[0][k] = histogram_cut_optimization_signal;

      TString name_cut_optimization_background = "cut_optimization_background";
      name_cut_optimization_background += ptBinMother;
      TH2F* hist_cut_optimization_background = new TH2F(name_cut_optimization_background.Data(), "Total background for different cuts; Cut number; Entries", nCutOptimizationBins, 0, nCutOptimizationBins, 2 * nSigmaBins, -nSigmaBins, nSigmaBins);
      hist_cut_optimization_background->Sumw2();
      hist_cut_optimization_background->SetLineColor(6);
      hist_cut_optimization_background->SetMarkerStyle(20);
      hist_cut_optimization_background->SetMarkerSize(0.6);
      hist_cut_optimization_background->SetMarkerColor(6);
      TH2F* histogram_cut_optimization_background = (TH2F*)hist_cut_optimization_background->Clone();
      fOutputBPlusResults->Add(histogram_cut_optimization_background);
      fResultsHistogramArray2D[1][k] = histogram_cut_optimization_background;
    }
  }


  //==================================================

  TString name_BPluss_in_analysis = "BPlus_in_analysis";
  TH1F* hist_BPluss_in_analysis = new TH1F(name_BPluss_in_analysis.Data(), "Number of BPluss to kpipipi in the Analysis; Entries", 10, 0, 10);
  hist_BPluss_in_analysis->Sumw2();
  hist_BPluss_in_analysis->SetLineColor(6);
  hist_BPluss_in_analysis->SetMarkerStyle(20);
  hist_BPluss_in_analysis->SetMarkerSize(0.6);
  hist_BPluss_in_analysis->SetMarkerColor(6);
  hist_BPluss_in_analysis->SetStats(kTRUE);
  hist_BPluss_in_analysis->GetXaxis()->SetBinLabel(1, "no. of BPlus");
  hist_BPluss_in_analysis->GetXaxis()->SetBinLabel(2, "no. of BPlus to kpipi");
  hist_BPluss_in_analysis->GetXaxis()->SetBinLabel(3, "no. with all tracks in event");
  hist_BPluss_in_analysis->GetXaxis()->SetBinLabel(4, "no. ...");
  hist_BPluss_in_analysis->GetXaxis()->SetBinLabel(5, "no. ...");
  hist_BPluss_in_analysis->GetXaxis()->SetBinLabel(6, "no. ...");
  hist_BPluss_in_analysis->GetXaxis()->SetBinLabel(7, "no. ...");
  hist_BPluss_in_analysis->GetXaxis()->SetBinLabel(8, "no. ...");
  TH1F* hist_BPluss_in_analysis_mc = (TH1F*)hist_BPluss_in_analysis->Clone();
  fOutputBPlusResults->Add(hist_BPluss_in_analysis_mc);
  fResultsHistogramArray[3][0] = hist_BPluss_in_analysis_mc;

  TString name_BPlus_per_bin = "BPlus_per_bin";
  TH1F* hist_BPlus_per_bin = new TH1F(name_BPlus_per_bin.Data(), "Number of BPlus to kpipi in the Analysis per bin; Entries", fnPtBins, 0, fnPtBins);
  for (Int_t i = 0; i < fnPtBins; ++i)
  {
    TString bin_name = "";
    bin_name += fPtBinLimits[i];
    bin_name += "-";
    bin_name += fPtBinLimits[i + 1];
    hist_BPlus_per_bin->GetXaxis()->SetBinLabel(i + 1, bin_name);
  }
  TH1F* hist_BPlus_per_bin_mc = (TH1F*)hist_BPlus_per_bin->Clone();
  fOutputBPlusResults->Add(hist_BPlus_per_bin_mc);
  fResultsHistogramArray[3][1] = hist_BPlus_per_bin_mc;

  TString name_BPlus_per_bin_in_Acc = "BPlus_per_bin_in_Acc";
  TH1F* hist_BPlus_per_bin_in_Acc = new TH1F(name_BPlus_per_bin_in_Acc.Data(), "Number of BPlus to kpipi in the Analysis per bin with all daughters in acceptance; Entries", fnPtBins, 0, fnPtBins);
  for (Int_t i = 0; i < fnPtBins; ++i)
  {
    TString bin_name = "";
    bin_name += fPtBinLimits[i];
    bin_name += "-";
    bin_name += fPtBinLimits[i + 1];
    hist_BPlus_per_bin_in_Acc->GetXaxis()->SetBinLabel(i + 1, bin_name);
  }
  TH1F* hist_BPlus_per_bin_in_Acc_mc = (TH1F*)hist_BPlus_per_bin_in_Acc->Clone();
  fOutputBPlusResults->Add(hist_BPlus_per_bin_in_Acc_mc);
  fResultsHistogramArray[3][2] = hist_BPlus_per_bin_in_Acc_mc;

  //======================================================================================================================================================

  //we make the histograms for the Pions and Kaon
  for (Int_t i = 0; i < 3; i++) {

    TString add_name = "";
    TList * listout = 0x0;
    if (i == 0) listout = fOutputD0FirstDaughter;
    if (i == 1) listout = fOutputD0SecondDaughter;
    if (i == 2) listout = fOutputBPlusPion;

    for (Int_t j = 0; j < 6; j++) {
      if (j == 0) add_name = "";
      if (j == 1) add_name = "Signal";
      if (j == 2) add_name = "Cut";
      if (j == 3) add_name = "SignalCut";
      if (j == 4) add_name = "Result";
      if (j == 5) add_name = "SignalResult";

      TString name_Histogram = "";
      TString discription_Histogram = "";
      Int_t numberOfBins = 0;
      Double_t lowerBound = 0.0;
      Double_t upperBound = 0.0;

      for (Int_t k = 0; k < 10; ++k)
      {
        if (k == 0) {name_Histogram = "ptTrack"; discription_Histogram = "pt track; p_{T} [GeV/c]; Entries"; numberOfBins = 600; lowerBound = 0; upperBound = 30;}
        if (k == 1) {name_Histogram = "momentumTrack"; discription_Histogram = "momentum track; p [GeV/c]; Entries"; numberOfBins = 600; lowerBound = 0; upperBound = 30;}
        if (k == 2) {name_Histogram = "numberOfITS"; discription_Histogram = "Number of ITS clusters track; [#]; Entries"; numberOfBins = 10; lowerBound = -0.5; upperBound = 9.5;}
        if (k == 3) {name_Histogram = "numberOfTPC"; discription_Histogram = "Number of TPC clusters track; [#]; Entries"; numberOfBins = 601; lowerBound = -0.5; upperBound = 600.5;}
        if (k == 4) {name_Histogram = "pointsOnITS"; discription_Histogram = "Number of ITS clusters track per layer; [#]; Entries"; numberOfBins = 10; lowerBound = -0.5; upperBound = 9.5;}
        if (k == 5) {name_Histogram = "nSigmaTPC"; discription_Histogram = "n sigma TPC for track PID; sigma; Entries"; numberOfBins = 500; lowerBound = -5; upperBound = 5;}
        if (k == 6) {name_Histogram = "nSigmaTOF"; discription_Histogram = "n sigma TOF for track PID; sigma; Entries"; numberOfBins = 500; lowerBound = -5; upperBound = 5;}
        if (k == 7) {name_Histogram = "nSigmaTPCandTOF"; discription_Histogram = "n sigma TPC and TOF for track PID; a.u.; Entries"; numberOfBins = 1000; lowerBound = 0; upperBound = 10;}
        if (k == 8) {name_Histogram = "impactParameter"; discription_Histogram = "Impact Parameter track;  [cm]; Entries"; numberOfBins = 2000; lowerBound = 0; upperBound = 0.5;}
        if (k == 9) {name_Histogram = "EtaTrack"; discription_Histogram = "Eta Track; Entries"; numberOfBins = 1000; lowerBound = -3; upperBound = 3;}

        name_Histogram += add_name;
        TH1F* histogram = new TH1F(name_Histogram.Data(), discription_Histogram.Data(), numberOfBins, lowerBound, upperBound);
        histogram->Sumw2();
        if (j % 2 == 0) histogram->SetLineColor(6);
        if (j % 2 == 1) histogram->SetLineColor(4);
        histogram->SetMarkerStyle(20);
        histogram->SetMarkerSize(0.6);
        if (j % 2 == 0) histogram->SetMarkerColor(6);
        if (j % 2 == 1) histogram->SetMarkerColor(4);
        TH1F* histogram_Clone = (TH1F*)histogram->Clone();
        listout->Add(histogram_Clone);
        fDaughterHistogramArray[i][j][k] = histogram_Clone;
      }

      TString numberofparticlesperevent = "numberofparticlesperevent";
      numberofparticlesperevent += add_name;
      TH1F* hist_numberofparticlesperevent = new TH1F(numberofparticlesperevent.Data(), "Number of particles per event; number of particles in one event; Entries", 100, 0, 100);
      hist_numberofparticlesperevent->Sumw2();
      hist_numberofparticlesperevent->SetLineColor(6);
      hist_numberofparticlesperevent->SetMarkerStyle(20);
      hist_numberofparticlesperevent->SetMarkerSize(0.6);
      hist_numberofparticlesperevent->SetMarkerColor(6);
      TH1F* histogram_numberofparticlesperevent = (TH1F*)hist_numberofparticlesperevent->Clone();
      listout->Add(histogram_numberofparticlesperevent);
      fDaughterHistogramArray[i][j][10] = histogram_numberofparticlesperevent;
    }

    TH1F * effectOfCuts = new TH1F("effectOfCuts", "Removal counter", 18, 0, 18);
    effectOfCuts->SetStats(kTRUE);
    effectOfCuts->GetXaxis()->SetTitle("Cut number");
    effectOfCuts->GetYaxis()->SetTitle("Particles cut");
    effectOfCuts->GetXaxis()->SetBinLabel(1, "total");
    effectOfCuts->GetXaxis()->SetBinLabel(2, "1");
    effectOfCuts->GetXaxis()->SetBinLabel(3, "2");
    effectOfCuts->GetXaxis()->SetBinLabel(4, "3");
    effectOfCuts->GetXaxis()->SetBinLabel(5, "4");
    effectOfCuts->GetXaxis()->SetBinLabel(6, "5");
    effectOfCuts->GetXaxis()->SetBinLabel(7, "6");
    effectOfCuts->GetXaxis()->SetBinLabel(8, "7");
    effectOfCuts->GetXaxis()->SetBinLabel(9, "8");
    effectOfCuts->GetXaxis()->SetBinLabel(10, "9");
    effectOfCuts->GetXaxis()->SetBinLabel(11, "10");
    effectOfCuts->GetXaxis()->SetBinLabel(12, "11");
    effectOfCuts->GetXaxis()->SetBinLabel(13, "12");
    effectOfCuts->GetXaxis()->SetBinLabel(14, "13");
    effectOfCuts->GetXaxis()->SetBinLabel(15, "14");
    effectOfCuts->GetXaxis()->SetBinLabel(16, "15");
    effectOfCuts->GetXaxis()->SetBinLabel(17, "16");
    effectOfCuts->GetXaxis()->SetBinLabel(18, "17");
    listout->Add(effectOfCuts);
    fDaughterHistogramArrayExtra[i][0] = effectOfCuts;

    TH1F * effectOfCutsSignal = new TH1F("effectOfCutsSignal", "Removal counter", 18, 0, 18);
    effectOfCutsSignal->SetStats(kTRUE);
    effectOfCutsSignal->GetXaxis()->SetTitle("Cut number");
    effectOfCutsSignal->GetYaxis()->SetTitle("Particles cut");
    effectOfCutsSignal->GetXaxis()->SetBinLabel(1, "total");
    effectOfCutsSignal->GetXaxis()->SetBinLabel(2, "1");
    effectOfCutsSignal->GetXaxis()->SetBinLabel(3, "2");
    effectOfCutsSignal->GetXaxis()->SetBinLabel(4, "3");
    effectOfCutsSignal->GetXaxis()->SetBinLabel(5, "4");
    effectOfCutsSignal->GetXaxis()->SetBinLabel(6, "5");
    effectOfCutsSignal->GetXaxis()->SetBinLabel(7, "6");
    effectOfCutsSignal->GetXaxis()->SetBinLabel(8, "7");
    effectOfCutsSignal->GetXaxis()->SetBinLabel(9, "8");
    effectOfCutsSignal->GetXaxis()->SetBinLabel(10, "9");
    effectOfCutsSignal->GetXaxis()->SetBinLabel(11, "10");
    effectOfCutsSignal->GetXaxis()->SetBinLabel(12, "11");
    effectOfCutsSignal->GetXaxis()->SetBinLabel(13, "12");
    effectOfCutsSignal->GetXaxis()->SetBinLabel(14, "13");
    effectOfCutsSignal->GetXaxis()->SetBinLabel(15, "14");
    effectOfCutsSignal->GetXaxis()->SetBinLabel(16, "15");
    effectOfCutsSignal->GetXaxis()->SetBinLabel(17, "16");
    effectOfCutsSignal->GetXaxis()->SetBinLabel(18, "17");
    listout->Add(effectOfCutsSignal);
    fDaughterHistogramArrayExtra[i][1] = effectOfCutsSignal;

    TString name_particle_pdg = "particle_pdg";
    TH1F* hist_particle_pdg = new TH1F(name_particle_pdg.Data(), "Pdg code particle; pdg code; Entries", 2000, -0.5, 1999.5);
    hist_particle_pdg->Sumw2();
    hist_particle_pdg->SetLineColor(6);
    hist_particle_pdg->SetMarkerStyle(20);
    hist_particle_pdg->SetMarkerSize(0.6);
    hist_particle_pdg->SetMarkerColor(6);
    TH1F* histogram_particle_pdg = (TH1F*)hist_particle_pdg->Clone();
    listout->Add(histogram_particle_pdg);
    fDaughterHistogramArrayExtra[i][2] = histogram_particle_pdg;

    TString name_particle_mother_pdg = "particle_mother_pdg";
    TH1F* hist_particle_mother_pdg = new TH1F(name_particle_mother_pdg.Data(), "Pdg code particle mother; pdg code; Entries", 2000, -0.5, 1999.5);
    hist_particle_mother_pdg->Sumw2();
    hist_particle_mother_pdg->SetLineColor(6);
    hist_particle_mother_pdg->SetMarkerStyle(20);
    hist_particle_mother_pdg->SetMarkerSize(0.6);
    hist_particle_mother_pdg->SetMarkerColor(6);
    TH1F* histogram_particle_mother_pdg = (TH1F*)hist_particle_mother_pdg->Clone();
    listout->Add(histogram_particle_mother_pdg);
    fDaughterHistogramArrayExtra[i][3] = histogram_particle_mother_pdg;

    TString name_ptBPlus_vs_ptTrack = "ptBPlus_vs_ptTrackBackground";
    TH2F* hist_ptBPlus_vs_ptTrack = new TH2F(name_ptBPlus_vs_ptTrack.Data(), "Pt BPlus vs Pt ; p_{T} BPlus [GeV/c]; p_{T} track [GeV/c]", 100, 0, 30, 100, 0, 30);
    hist_ptBPlus_vs_ptTrack->Sumw2();
    hist_ptBPlus_vs_ptTrack->SetLineColor(6);
    hist_ptBPlus_vs_ptTrack->SetMarkerStyle(20);
    hist_ptBPlus_vs_ptTrack->SetMarkerSize(0.6);
    hist_ptBPlus_vs_ptTrack->SetMarkerColor(6);
    TH2F* histogram_ptBPlus_vs_ptTrack = (TH2F*)hist_ptBPlus_vs_ptTrack->Clone();
    listout->Add(histogram_ptBPlus_vs_ptTrack);
    fDaughterHistogramArray2D[i][4] = histogram_ptBPlus_vs_ptTrack;

    TString name_ptBPlus_vs_ptTrackSignal = "ptBPlus_vs_ptTrackSignal";
    TH2F* hist_ptBPlus_vs_ptTrackSignal = new TH2F(name_ptBPlus_vs_ptTrackSignal.Data(), "Pt BPlus vs Pt ; p_{T} BPlus [GeV/c]; p_{T} track [GeV/c]", 100, 0, 30, 100, 0, 30);
    hist_ptBPlus_vs_ptTrackSignal->Sumw2();
    hist_ptBPlus_vs_ptTrackSignal->SetLineColor(4);
    hist_ptBPlus_vs_ptTrackSignal->SetMarkerStyle(20);
    hist_ptBPlus_vs_ptTrackSignal->SetMarkerSize(0.6);
    hist_ptBPlus_vs_ptTrackSignal->SetMarkerColor(6);
    TH2F* histogram_ptBPlus_vs_ptTrackSignal = (TH2F*)hist_ptBPlus_vs_ptTrackSignal->Clone();
    listout->Add(histogram_ptBPlus_vs_ptTrackSignal);
    fDaughterHistogramArray2D[i][5] = histogram_ptBPlus_vs_ptTrackSignal;
  }

  //we make the histograms for the reconstructed particles
  for (Int_t i = 0; i < 3; i++) {

    TString add_name = "";
    TList * listout = 0x0;
    Int_t nHistogramSets = 0;
    if (i == 0) {listout = fOutputD0; nHistogramSets = 6 + 2 * fnPtBins;}
    if (i == 1) {listout = fOutputBPlus; nHistogramSets = 6 + 2 * fnPtBins;}
    if (i == 2) {listout = fOutputD0_D0Pt; nHistogramSets = 2 * fnPtBinsD0forD0ptbin;}


    for (Int_t j = 0; j < nHistogramSets; j++) {
      if (i < 2)
      {
        if (j == 0) add_name = "Uncut";
        if (j == 1) add_name = "SignalUncut";
        if (j == 2) add_name = "Cut";
        if (j == 3) add_name = "SignalCut";
        if (j == 4) add_name = "Result";
        if (j == 5) add_name = "SignalResult";
        if (j % 2 == 0 && j > 5) {add_name = "_ptbin_"; add_name += fPtBinLimits[(j - 6) / 2]; add_name += "_to_"; add_name += fPtBinLimits[(j - 6) / 2 + 1];}
        if (j % 2 == 1 && j > 5) {add_name = "Signal_ptbin_"; add_name += fPtBinLimits[(j - 7) / 2]; add_name += "_to_"; add_name += fPtBinLimits[(j - 7) / 2 + 1];}
      }
      if (i == 2)
      {
        if (j % 2 == 0) {add_name = "_ptbin_"; add_name += fPtBinLimitsD0forD0ptbin[j / 2]; add_name += "_to_"; add_name += fPtBinLimitsD0forD0ptbin[1 + j / 2];}
        if (j % 2 == 1) {add_name = "Signal_ptbin_"; add_name += fPtBinLimitsD0forD0ptbin[(j - 1) / 2]; add_name += "_to_"; add_name += fPtBinLimitsD0forD0ptbin[1 + (j - 1) / 2];}
      }

      TList * listCandidates = new TList();
      listCandidates->SetOwner();
      listCandidates->SetName(add_name);
      listout->Add(listCandidates);

      TString name_Histogram = "";
      TString discription_Histogram  = "";
      Int_t numberOfBins = 0;
      Double_t lowerBound = 0.0;
      Double_t upperBound = 0.0;
      Int_t numberOfBinsTwo = 0;
      Double_t lowerBoundTwo = 0.0;
      Double_t upperBoundTwo = 0.0;

      for (Int_t k = 0; k < 49; ++k)
      {
        if (k == 0) {name_Histogram = "ptMother"; discription_Histogram = "pt mother; p_{T} [GeV/c]; Entries"; numberOfBins = 300; lowerBound = 0; upperBound = 30;}
        if (k == 1) {name_Histogram = "ptFirstDaughter"; discription_Histogram = "pt first daughter; p_{T} [GeV/c]; Entries"; numberOfBins = 300; lowerBound = 0; upperBound = 30;}
        if (k == 2) {name_Histogram = "ptSecondDaughter"; discription_Histogram = "pt second daughter; p_{T} [GeV/c]; Entries"; numberOfBins = 300; lowerBound = 0; upperBound = 30;}
        if (k == 3) {name_Histogram = "etaMother"; discription_Histogram = "eta mother; #eta; Entries"; numberOfBins = 100; lowerBound = -2; upperBound = 2;}
        if (k == 4) {name_Histogram = "phiMother"; discription_Histogram = "phi mother; #phi; Entries"; numberOfBins = 25; lowerBound = 0; upperBound = 2 * TMath::Pi();}
        if (k == 5) {name_Histogram = "d0Mother"; discription_Histogram = "d0 mother;  [cm]; Entries"; numberOfBins = 500; lowerBound = 0; upperBound = 1.0;}
        if (k == 6) {name_Histogram = "d0FirstDaughter"; discription_Histogram = "d0 first daughter;  [cm]; Entries"; numberOfBins = 500; lowerBound = 0; upperBound = 1;}

        if (k == 7) {name_Histogram = "d0SecondDaughter"; discription_Histogram = "d0 second daughter;  [cm]; Entries"; numberOfBins = 500; lowerBound = 0; upperBound = 1;}

        if (k == 8) {name_Histogram = "pointingAngleMother"; discription_Histogram = "pointing angle;  [Cos(#theta)]; Entries"; numberOfBins = 100; lowerBound = -1; upperBound = 1;}
        if (k == 9) {name_Histogram = "impactProduct"; discription_Histogram = "impact product; [cm^{2}]; Entries"; numberOfBins = 500; lowerBound = -0.01; upperBound = 0.01;}
        if (k == 10) {name_Histogram = "impactProductXY"; discription_Histogram = "impact product XY; [cm^{2}]; Entries"; numberOfBins = 200; lowerBound = 0; upperBound = 0.5;}
        if (k == 11) {name_Histogram = "invariantMassMother"; discription_Histogram = "mass mother candidate; m [GeV/c^{2}]; Entries"; numberOfBins = 10000; lowerBound = 0; upperBound = 10;}
        if (k == 12) {name_Histogram = "deltaMassMother"; discription_Histogram = "mass mother candidate; m [GeV/c^{2}]; Entries"; numberOfBins = 10000; lowerBound = 0; upperBound = 10;}
        if (k == 13) {name_Histogram = "dcaMother"; discription_Histogram = "dca mother; distance [cm]; Entries"; numberOfBins = 500; lowerBound = 0; upperBound = 0.25;}
        if (k == 14) {name_Histogram = "vertexDistance"; discription_Histogram = "vertex distance between mother and primary vertex; distance [cm]; Entries"; numberOfBins = 100; lowerBound = 0; upperBound = 1;}
        if (k == 15) {name_Histogram = "normDecayLength"; discription_Histogram = "Normalized decay length w.r.t primary vertex; [cm]; Entries"; numberOfBins = 100; lowerBound = 0; upperBound = 50;}
        if (k == 16) {name_Histogram = "pseudoProperDecayTime"; discription_Histogram = "Pseudo Proper Decay Time w.r.t primary vertex; [a.u.]; Entries"; numberOfBins = 1000; lowerBound = -10; upperBound = 10;}
        if (k == 17) {name_Histogram = "DecayTime"; discription_Histogram = "Decay Time w.r.t primary vertex; [a.u.]; Entries"; numberOfBins = 100; lowerBound = 0; upperBound = 0.00000001;}
        if (k == 18) {name_Histogram = "normDecayTime"; discription_Histogram = "Normalized Decay Time w.r.t primary vertex; [a.u.]; Entries"; numberOfBins = 100; lowerBound = 0; upperBound = 0.0000001;}
        if (k == 19) {name_Histogram = "angleMotherFirstDaughter"; discription_Histogram = "flight angle mother and first daughter; [Cos(#phi)]; Entries"; numberOfBins = 100; lowerBound = -1; upperBound = 1;}
        if (k == 20) {name_Histogram = "angleMotherSecondDaughter"; discription_Histogram = "flight angle mother and second daughter; [Cos(#phi)]; Entries"; numberOfBins = 100; lowerBound = 0.5; upperBound = 1;}
        if (k == 21) {name_Histogram = "angleBetweenBothDaughters"; discription_Histogram = "angle between both daughters; [Cos(#phi)]; Entries"; numberOfBins = 100; lowerBound = -1; upperBound = 1;}
        if (k == 22) {name_Histogram = "cosThetaStar"; discription_Histogram = "cosThetaStar; [Cos(#theta*)]; Entries"; numberOfBins = 200; lowerBound = -2; upperBound = 2;}
        if (k == 23) {name_Histogram = "vertexX"; discription_Histogram = "Vertex position; [cm]; Entries"; numberOfBins = 500; lowerBound = -5; upperBound = 5;}
        if (k == 24) {name_Histogram = "vertexY"; discription_Histogram = "Vertex position; [cm]; Entries"; numberOfBins = 500; lowerBound = -5; upperBound = 5;}
        if (k == 25) {name_Histogram = "vertexZ"; discription_Histogram = "Vertex position; [cm]; Entries"; numberOfBins = 500; lowerBound = -20; upperBound = 20;}


        if (k == 26) {
          if (i == 0) {name_Histogram = "pointingAngleToBPlus"; discription_Histogram = "Pointing angle w.r.t. BPlus decay vertex; [Cos(#theta)]; Entries"; numberOfBins = 200; lowerBound = -1; upperBound = 1;}
          else continue;
        }
        if (k == 27) {
          if (i == 0) {name_Histogram = "d0MotherToBPlus"; discription_Histogram = "d0 Mother w.r.t. BPlus decay vertex; [cm]; Entries"; numberOfBins = 500; lowerBound = 0; upperBound = 1;}
          else continue;
        }
        if (k == 28) {
          if (i == 0) {name_Histogram = "d0FirstDaughterToBPlus"; discription_Histogram = "d0 first daughter w.r.t. BPlus decay vertex; [cm]; Entries"; numberOfBins = 500; lowerBound = 0; upperBound = 1;}
          else continue;
        }
        if (k == 29) {
          if (i == 0) {name_Histogram = "d0SecondDaughterToBPlus"; discription_Histogram = "d0 second daughter w.r.t. BPlus decay vertex; [cm]; Entries"; numberOfBins = 500; lowerBound = 0; upperBound = 1;}
          else continue;
        }
        if (k == 30) {
          if (i == 0) {name_Histogram = "impactProductToBPlus"; discription_Histogram = "impact product w.r.t. BPlus decay vertex; [cm]; Entries"; numberOfBins = 400; lowerBound = -0.02; upperBound = 0.02;}
          else continue;
        }
        if (k == 31) {
          if (i == 0) {name_Histogram = "impactProductXYToBPlus"; discription_Histogram = "impact product XY w.r.t. BPlus decay vertex; [cm^{2}]; Entries"; numberOfBins = 100; lowerBound = 0; upperBound = 0.5;}
          else continue;
        }
        if (k == 32) {
          if (i == 0) {name_Histogram = "normDecayLengthToBPlus"; discription_Histogram = "Normalized decay length w.r.t. BPlus decay vertex; [cm]; Entries"; numberOfBins = 100; lowerBound = 0; upperBound = 50;}
          else continue;
        }
        if (k == 33) {
          if (i == 0) {name_Histogram = "pseudoProperDecayTimeToBPlus"; discription_Histogram = "Pseudo Proper Decay Time w.r.t BPlus vertex; [a.u.]; Entries"; numberOfBins = 1000; lowerBound = -1; upperBound = 1;}
          else continue;
        }
        if (k == 34) {
          if (i == 0) {name_Histogram = "DecayTimeToBPlus"; discription_Histogram = "Decay Time w.r.t BPlus vertex; [a.u.]; Entries"; numberOfBins = 100; lowerBound = 0; upperBound = 0.00000001;}
          else continue;
        }
        if (k == 35) {
          if (i == 0) {name_Histogram = "normDecayTimeToBPlus"; discription_Histogram = "Normalized Decay Time w.r.t BPlus vertex; [a.u.]; Entries"; numberOfBins = 100; lowerBound = 0; upperBound = 0.00000001;}
          else continue;
        }

        if (k == 36) {name_Histogram = "topomaticFirstDaughter"; discription_Histogram = "topomatic d0 first daughter; [cm]; Entries"; numberOfBins = 500; lowerBound = 0; upperBound = 20;}
        if (k == 37) {name_Histogram = "topomaticSecondDaughter"; discription_Histogram = "topomatic d0 second daughter; [cm]; Entries"; numberOfBins = 500; lowerBound = 0; upperBound = 20;}
        if (k == 38) {name_Histogram = "topomaticMax"; discription_Histogram = "Max topomatic; [cm]; Entries"; numberOfBins = 500; lowerBound = 0; upperBound = 20;}
        if (k == 39) {name_Histogram = "topomaticMin"; discription_Histogram = "Min topomatic; [cm]; Entries"; numberOfBins = 500; lowerBound = 0; upperBound = 20;}
        if (k == 40) {name_Histogram = "pointingAngleMotherXY"; discription_Histogram = "pointing angle XY;  [Cos(#theta)]; Entries"; numberOfBins = 1000; lowerBound = -1; upperBound = 1;}
        if (k == 41) {name_Histogram = "vertexDistanceXY"; discription_Histogram = "vertex distance between mother and primary vertex XY; distance [cm]; Entries"; numberOfBins = 1000; lowerBound = 0; upperBound = 10;}
        if (k == 42) {name_Histogram = "normDecayLengthXY"; discription_Histogram = "Normalized decay length w.r.t primary vertex XY; [cm]; Entries"; numberOfBins = 100; lowerBound = 0; upperBound = 50;}
        if (k == 43) {name_Histogram = "vertexChi2"; discription_Histogram = "#Chi^{2} Vertex; [#Chi^{2}]; Entries"; numberOfBins = 1000; lowerBound = 0; upperBound = 50;}
        if (k == 44) {name_Histogram = "vertexChi2NDF"; discription_Histogram = "#Chi^{2} per NDF Vertex; [#Chi^{2}/NDF]; Entries"; numberOfBins = 1000; lowerBound = 0; upperBound = 50;}
        
        if (k == 45) {name_Histogram = "Normd0FirstDaughter"; discription_Histogram = "norm d0 first daughter;  [cm]; Entries"; numberOfBins = 500; lowerBound = 0; upperBound = 100;}
        if (k == 46) {name_Histogram = "Normd0SecondDaughter"; discription_Histogram = "norm d0 second daughter;  [cm]; Entries"; numberOfBins = 500; lowerBound = 0; upperBound = 100;}
        if (k == 47) {name_Histogram = "Normd0Mother"; discription_Histogram = "norm d0 mother; [cm]; Entries"; numberOfBins = 500; lowerBound = 0; upperBound = 100;}
        if (k == 48) {name_Histogram = "NormimpactProduct"; discription_Histogram = "norm impact product; [cm^{2}]; Entries"; numberOfBins = 500; lowerBound = -200; upperBound = 200;}

        name_Histogram += add_name;
        TH1F* histogram = new TH1F(name_Histogram.Data(), discription_Histogram.Data(), numberOfBins, lowerBound, upperBound);
        histogram->Sumw2();
        if (j % 2 == 0) histogram->SetLineColor(6);
        if (j % 2 == 1) histogram->SetLineColor(4);
        histogram->SetMarkerStyle(20);
        histogram->SetMarkerSize(0.6);
        if (j % 2 == 0) histogram->SetMarkerColor(6);
        if (j % 2 == 1) histogram->SetMarkerColor(4);
        TH1F* histogram_Clone = (TH1F*)histogram->Clone();
        listCandidates->Add(histogram_Clone);
        fMotherHistogramArray[i][j][k] = histogram_Clone;
      }


      name_Histogram = "";
      discription_Histogram  = "";
      numberOfBins = 0;
      lowerBound = 0.0;
      upperBound = 0.0;
      numberOfBinsTwo = 0;
      lowerBoundTwo = 0.0;
      upperBoundTwo = 0.0;

      //we make the 2D histograms for the reconstructed particles
      Int_t nFirst = 0;
      Int_t nSecond = 1;
      Int_t nVariables = 10;
      Int_t nHistograms = nVariables * (nVariables - 1) / 2;

      TList * list2D = new TList();
      list2D->SetOwner();
      TString name2D = "2D_Histograms";
      name2D += add_name;
      list2D->SetName(name2D);
      listCandidates->Add(list2D);

      for (Int_t k = 0; k < nHistograms; ++k)
      {
        numberOfBins = 50; numberOfBinsTwo = 50;
        if (nFirst == 0) {name_Histogram = "d0FirstDaughter"; discription_Histogram = "d0 first daughter [cm];"; lowerBound = 0; upperBound = 1;}
        if (nFirst == 1) {name_Histogram = "d0SecondDaughter"; discription_Histogram = "d0 second daughter [cm];"; lowerBound = 0; upperBound = 1;}
        if (nFirst == 2) {name_Histogram = "d0Mother"; discription_Histogram = "d0 mother [cm];"; lowerBound = 0; upperBound = 1;}
        if (nFirst == 3) {name_Histogram = "pointingAngleMother"; discription_Histogram = "pointing angle [Cos(#theta)];"; lowerBound = -1; upperBound = 1;}
        if (nFirst == 4) {name_Histogram = "impactProduct"; discription_Histogram = "impact product [cm^{2}];"; lowerBound = -0.01; upperBound = 0.01;}
        if (nFirst == 5) {name_Histogram = "impactProductXY"; discription_Histogram = "impact product XY [cm^{2}];"; lowerBound = 0; upperBound = 0.5;}
        if (nFirst == 6) {name_Histogram = "vertexDistance"; discription_Histogram = "vertex distance between mother and primary vertex  [cm];"; lowerBound = 0; upperBound = 1;}
        if (nFirst == 7) {name_Histogram = "normDecayLength"; discription_Histogram = "Normalized decay length w.r.t primary vertex [cm];"; lowerBound = 0; upperBound = 50;}
        if (nFirst == 8) {name_Histogram = "pointingAngleMotherXY"; discription_Histogram = "pointing angle XY [Cos(#theta)];"; lowerBound = -1; upperBound = 1;}
        if (nFirst == 9) {name_Histogram = "vertexDistanceXY"; discription_Histogram = "vertex distance between mother and primary vertex XY [cm];"; lowerBound = 0; upperBound = 1;}
        if (nFirst == 10) {name_Histogram = "normDecayLengthXY"; discription_Histogram = "Normalized decay length w.r.t primary vertex XY [cm];"; lowerBound = 0; upperBound = 50;}

        if (nSecond == 0) {name_Histogram += "d0FirstDaughter"; discription_Histogram += "d0 first daughter [cm];"; lowerBoundTwo = 0; upperBoundTwo = 1;}
        if (nSecond == 1) {name_Histogram += "d0SecondDaughter"; discription_Histogram += "d0 second daughter [cm];"; lowerBoundTwo = 0; upperBoundTwo = 1;}
        if (nSecond == 2) {name_Histogram += "d0Mother"; discription_Histogram += "d0 mother [cm];"; lowerBoundTwo = 0; upperBoundTwo = 1;}
        if (nSecond == 3) {name_Histogram += "pointingAngleMother"; discription_Histogram += "pointing angle [Cos(#theta)];"; lowerBoundTwo = -1; upperBoundTwo = 1;}
        if (nSecond == 4) {name_Histogram += "impactProduct"; discription_Histogram += "impact product [cm^{2}];"; lowerBoundTwo = -0.01; upperBoundTwo = 0.01;}
        if (nSecond == 5) {name_Histogram += "impactProductXY"; discription_Histogram += "impact product XY [cm^{2}];"; lowerBoundTwo = 0; upperBoundTwo = 0.5;}
        if (nSecond == 6) {name_Histogram += "vertexDistance"; discription_Histogram += "vertex distance between mother and primary vertex  [cm];"; lowerBoundTwo = 0; upperBoundTwo = 1;}
        if (nSecond == 7) {name_Histogram += "normDecayLength"; discription_Histogram += "Normalized decay length w.r.t primary vertex [cm];"; lowerBoundTwo = 0; upperBoundTwo = 50;}
        if (nSecond == 8) {name_Histogram += "_pointingAngleMotherXY"; discription_Histogram += "pointing angle XY [Cos(#theta)];"; lowerBoundTwo = -1; upperBoundTwo = 1;}
        if (nSecond == 9) {name_Histogram += "_vertexDistanceXY"; discription_Histogram += "vertex distance between mother and primary vertex XY [cm];"; lowerBoundTwo = 0; upperBoundTwo = 1;}
        if (nSecond == 10) {name_Histogram += "_normDecayLengthXY"; discription_Histogram += "Normalized decay length w.r.t primary vertex XY [cm];"; lowerBoundTwo = 0; upperBoundTwo = 50;}

        name_Histogram += add_name;
        TH2F* histogram = new TH2F(name_Histogram.Data(), discription_Histogram.Data(), numberOfBins, lowerBound, upperBound, numberOfBinsTwo, lowerBoundTwo, upperBoundTwo);
        histogram->Sumw2();
        if (j % 2 == 0) histogram->SetLineColor(6);
        if (j % 2 == 1) histogram->SetLineColor(4);
        histogram->SetMarkerStyle(20);
        histogram->SetMarkerSize(0.6);
        histogram->SetMarkerColor(6);
        TH2F* histogram_Clone = (TH2F*)histogram->Clone();
        list2D->Add(histogram_Clone);
        fMotherHistogramArray2D[i][j][k] = histogram_Clone;

        nSecond++;
        if (nSecond > nVariables)
        {
          nFirst++;
          nSecond = nFirst + 1;
        }
      }
    }

    TH1F * effectOfCuts = new TH1F("effectOfCuts", "Removal counter", 100, 0, 100);
    effectOfCuts->SetStats(kTRUE);
    effectOfCuts->GetXaxis()->SetTitle("Cut number");
    effectOfCuts->GetYaxis()->SetTitle("Particles cut");
    effectOfCuts->GetXaxis()->SetBinLabel(1, "total");
    for (Int_t i = 1; i < 100; ++i)
    {
      TString integerText = "";
      integerText += i - 1;
      effectOfCuts->GetXaxis()->SetBinLabel(i + 1, integerText);
    }
    listout->Add(effectOfCuts);
    fMotherHistogramArrayExtra[i][0] = effectOfCuts;

    TH1F * effectOfCutsSignal = new TH1F("effectOfCutsSignal", "Removal counter", 100, 0, 100);
    effectOfCutsSignal->SetStats(kTRUE);
    effectOfCutsSignal->GetXaxis()->SetTitle("Cut number");
    effectOfCutsSignal->GetYaxis()->SetTitle("Particles cut");
    effectOfCutsSignal->GetXaxis()->SetBinLabel(1, "total");
    for (Int_t i = 1; i < 100; ++i)
    {
      TString integerText = "";
      integerText += i - 1;
      effectOfCutsSignal->GetXaxis()->SetBinLabel(i + 1, integerText);
    }
    listout->Add(effectOfCutsSignal);
    fMotherHistogramArrayExtra[i][1] = effectOfCutsSignal;

    TString name_particle_pdg = "particle_pdg";
    TH1F* hist_particle_pdg = new TH1F(name_particle_pdg.Data(), "Pdg code particle; pdg code; Entries", 2000, -0.5, 1999.5);
    hist_particle_pdg->Sumw2();
    hist_particle_pdg->SetLineColor(6);
    hist_particle_pdg->SetMarkerStyle(20);
    hist_particle_pdg->SetMarkerSize(0.6);
    hist_particle_pdg->SetMarkerColor(6);
    TH1F* histogram_particle_pdg = (TH1F*)hist_particle_pdg->Clone();
    listout->Add(histogram_particle_pdg);
    fMotherHistogramArrayExtra[i][2] = histogram_particle_pdg;

    TString name_particle_mother_pdg = "particle_mother_pdg";
    TH1F* hist_particle_mother_pdg = new TH1F(name_particle_mother_pdg.Data(), "Pdg code particle mother; pdg code; Entries", 2000, -0.5, 1999.5);
    hist_particle_mother_pdg->Sumw2();
    hist_particle_mother_pdg->SetLineColor(6);
    hist_particle_mother_pdg->SetMarkerStyle(20);
    hist_particle_mother_pdg->SetMarkerSize(0.6);
    hist_particle_mother_pdg->SetMarkerColor(6);
    TH1F* histogram_particle_mother_pdg = (TH1F*)hist_particle_mother_pdg->Clone();
    listout->Add(histogram_particle_mother_pdg);
    fMotherHistogramArrayExtra[i][3] = histogram_particle_mother_pdg;

    TString name_distance_vertex_from_real = "distance_vertex_from_real";
    TH1F* hist_distance_vertex_from_real = new TH1F(name_distance_vertex_from_real.Data(), "Distance reconstructed vertex from real vertex; distance [cm]; Entries", 500, 0, 0.5);
    hist_distance_vertex_from_real->Sumw2();
    hist_distance_vertex_from_real->SetLineColor(6);
    hist_distance_vertex_from_real->SetMarkerStyle(20);
    hist_distance_vertex_from_real->SetMarkerSize(0.6);
    hist_distance_vertex_from_real->SetMarkerColor(6);
    TH1F* histogram_distance_vertex_from_real = (TH1F*)hist_distance_vertex_from_real->Clone();
    listout->Add(histogram_distance_vertex_from_real);
    fMotherHistogramArrayExtra[i][4] = histogram_distance_vertex_from_real;

    TString name_distance_vertex_from_real_new = "distance_vertex_from_real_new";
    TH1F* hist_distance_vertex_from_real_new = new TH1F(name_distance_vertex_from_real_new.Data(), "Distance reconstructed vertex from real vertex; distance [cm]; Entries", 500, 0, 0.5);
    hist_distance_vertex_from_real_new->Sumw2();
    hist_distance_vertex_from_real_new->SetLineColor(6);
    hist_distance_vertex_from_real_new->SetMarkerStyle(20);
    hist_distance_vertex_from_real_new->SetMarkerSize(0.6);
    hist_distance_vertex_from_real_new->SetMarkerColor(6);
    TH1F* histogram_distance_vertex_from_real_new = (TH1F*)hist_distance_vertex_from_real_new->Clone();
    listout->Add(histogram_distance_vertex_from_real_new);
    fMotherHistogramArrayExtra[i][5] = histogram_distance_vertex_from_real_new;

    TString name_momentum_resolution = "momentum_resolution";
    TH1F* hist_momentum_resolution = new TH1F(name_momentum_resolution.Data(), "Momentum resolution; difference between real and reconstructed momentum [GeV/c]; Entries", 1000, 0, 1);
    hist_momentum_resolution->Sumw2();
    hist_momentum_resolution->SetLineColor(6);
    hist_momentum_resolution->SetMarkerStyle(20);
    hist_momentum_resolution->SetMarkerSize(0.6);
    hist_momentum_resolution->SetMarkerColor(6);
    TH1F* histogram_momentum_resolution = (TH1F*)hist_momentum_resolution->Clone();
    listout->Add(histogram_momentum_resolution);
    fMotherHistogramArrayExtra[i][6] = histogram_momentum_resolution;
  }



  //we make the histograms for the same sign method histograms and the pt bins
  for (Int_t i = 0; i < 7; ++i) {
    TString typeListName = "";
    if (i == 0) typeListName = "MassPlots";
    if (i == 1) typeListName = "MassPlots_SameSign";
    if (i == 2) typeListName = "MassPlots_SignSum";
    if (i == 3) typeListName = "MassPlots_Background_rotation";
    if (i == 4) typeListName = "MassPlots_HIJING_Background";
    if (i == 5) typeListName = "MassPlots_HIJING_Signal";
    if (i == 6) typeListName = "MassPlots_HIJING_Background_rotation";

    TList * listMassPlots = new TList();
    listMassPlots->SetOwner();
    listMassPlots->SetName(typeListName);
    fOutputBPlusResults->Add(listMassPlots);

    for (Int_t k = 0; k < fnPtBins + 3; ++k) {
      TString ptBinMother = "";
      if (k == 0) ptBinMother = "";
      if (k == 1) ptBinMother = "_ptbin_3_to_inf";
      if (k == 2) ptBinMother = "_ptbin_6_to_inf";
      if (k > 2) {ptBinMother += "_ptbin_"; ptBinMother += fPtBinLimits[k - 3]; ptBinMother += "_to_"; ptBinMother += fPtBinLimits[k - 2];}

      TString name_invariantMassMother = "invariantMassBPlus";
      name_invariantMassMother += ptBinMother;
      TH1F* hist_invariantMassMother = new TH1F(name_invariantMassMother.Data(), "mass mother candidate; m [GeV/c^2]; Entries", 2000, 0, 20);
      hist_invariantMassMother->Sumw2();
      hist_invariantMassMother->SetLineColor(6);
      hist_invariantMassMother->SetMarkerStyle(20);
      hist_invariantMassMother->SetMarkerSize(0.6);
      hist_invariantMassMother->SetMarkerColor(6);
      TH1F* histogram_invariantMassMother = (TH1F*)hist_invariantMassMother->Clone();
      listMassPlots->Add(histogram_invariantMassMother); 
      fResultsHistogramArray[4 + i][k] = histogram_invariantMassMother;
    }
    for (Int_t k = 0; k < fnPtBins + 3; ++k) {
      TString ptBinMother = "";
      if (k == 0) ptBinMother = "";
      if (k == 1) ptBinMother = "_ptbin_3_to_inf";
      if (k == 2) ptBinMother = "_ptbin_6_to_inf";
      if (k > 2) {ptBinMother += "_ptbin_"; ptBinMother += fPtBinLimits[k - 3]; ptBinMother += "_to_"; ptBinMother += fPtBinLimits[k - 2];}

      TString name_deltainvariantMassMother = "deltainvariantMassBPlus";
      name_deltainvariantMassMother += ptBinMother;
      TH1F* hist_deltainvariantMassMother = new TH1F(name_deltainvariantMassMother.Data(), "delta mass mother candidate; m [GeV/c^2]; Entries", 2000, 0, 20);
      hist_deltainvariantMassMother->Sumw2();
      hist_deltainvariantMassMother->SetLineColor(6);
      hist_deltainvariantMassMother->SetMarkerStyle(20);
      hist_deltainvariantMassMother->SetMarkerSize(0.6);
      hist_deltainvariantMassMother->SetMarkerColor(6);
      TH1F* histogram_deltainvariantMassMother = (TH1F*)hist_deltainvariantMassMother->Clone();
      listMassPlots->Add(histogram_deltainvariantMassMother);
      fResultsHistogramArray[4 + i][k + fnPtBins + 3] = histogram_deltainvariantMassMother;
    }
  }


  TString name_cutEffectBackground = "cutEffectBackground";
  TH2F* hist_cutEffectBackground = new TH2F(name_cutEffectBackground.Data(), "Effect of Cuts on background; cut number; cut number", 99, 0, 99, 99, 0, 99);
  for (int i = 0; i < 99; ++i)
  {
    TString integerText = "";
    integerText += i;
    hist_cutEffectBackground->GetXaxis()->SetBinLabel(i + 1, integerText);
    hist_cutEffectBackground->GetYaxis()->SetBinLabel(i + 1, integerText);
  }
  TH2F* histogram_cutEffectBackground = (TH2F*)hist_cutEffectBackground->Clone();
  fOutputBPlusResults->Add(histogram_cutEffectBackground);
  fResultsHistogramArray2D[2][0] = histogram_cutEffectBackground;


  TString name_cutEffectSignal = "cutEffectSignal";
  TH2F* hist_cutEffectSignal = new TH2F(name_cutEffectSignal.Data(), "Effect of Cuts on Signal; cut number; cut number", 99, 0, 99, 99, 0, 99);
  for (Int_t i = 0; i < 99; ++i)
  {
    TString integerText = "";
    integerText += i;
    hist_cutEffectSignal->GetXaxis()->SetBinLabel(i + 1, integerText);
    hist_cutEffectSignal->GetYaxis()->SetBinLabel(i + 1, integerText);
  }
  TH2F* histogram_cutEffectSignal = (TH2F*)hist_cutEffectSignal->Clone();
  fOutputBPlusResults->Add(histogram_cutEffectSignal);
  fResultsHistogramArray2D[2][1] = histogram_cutEffectSignal;

  TString name_cutEffectUniqueBackground = "cutEffectUniqueBackground";
  TH1F* hist_cutEffectUniqueBackground = new TH1F(name_cutEffectUniqueBackground.Data(), "Effect of Cuts on Signal; cut number; cut number", 99, 0, 99);
  for (Int_t i = 0; i < 99; ++i)
  {
    TString integerText = "";
    integerText += i;
    hist_cutEffectUniqueBackground->GetXaxis()->SetBinLabel(i + 1, integerText);
  }
  TH1F* histogram_cutEffectUniqueBackground = (TH1F*)hist_cutEffectUniqueBackground->Clone();
  fOutputBPlusResults->Add(histogram_cutEffectUniqueBackground);
  fResultsHistogramArray[13][0] = histogram_cutEffectUniqueBackground;

  TString name_cutEffectUniqueSignal = "cutEffectUniqueSignal";
  TH1F* hist_cutEffectUniqueSignal = new TH1F(name_cutEffectUniqueSignal.Data(), "Effect of Cuts on Signal; cut number; cut number", 99, 0, 99);
  for (Int_t i = 0; i < 99; ++i)
  {
    TString integerText = "";
    integerText += i;
    hist_cutEffectUniqueSignal->GetXaxis()->SetBinLabel(i + 1, integerText);
  }
  TH1F* histogram_cutEffectUniqueSignal = (TH1F*)hist_cutEffectUniqueSignal->Clone();
  fOutputBPlusResults->Add(histogram_cutEffectUniqueSignal);
  fResultsHistogramArray[13][1] = histogram_cutEffectUniqueSignal;

  TString name_totalITSBackground = "totalITSBackground";
  TH1F* hist_totalITSBackground = new TH1F(name_totalITSBackground.Data(), "Total nr. of ITS hits for the daughters; number [#]; Entries", 30, 0, 30);
  hist_totalITSBackground->Sumw2();
  hist_totalITSBackground->SetLineColor(6);
  hist_totalITSBackground->SetMarkerStyle(20);
  hist_totalITSBackground->SetMarkerSize(0.6);
  hist_totalITSBackground->SetMarkerColor(6);
  TH1F* histogram_totalITSBackground = (TH1F*)hist_totalITSBackground->Clone();
  fOutputBPlusResults->Add(histogram_totalITSBackground);
  fResultsHistogramArray[11][0] = histogram_totalITSBackground;

  TString name_totalITSSignal = "totalITSSignal";
  TH1F* hist_totalITSSignal = new TH1F(name_totalITSSignal.Data(), "Total nr. of ITS hits for the daughters; number [#]; Entries", 30, 0, 30);
  hist_totalITSSignal->Sumw2();
  hist_totalITSSignal->SetLineColor(6);
  hist_totalITSSignal->SetMarkerStyle(20);
  hist_totalITSSignal->SetMarkerSize(0.6);
  hist_totalITSSignal->SetMarkerColor(6);
  TH1F* histogram_totalITSSignal = (TH1F*)hist_totalITSSignal->Clone();
  fOutputBPlusResults->Add(histogram_totalITSSignal);
  fResultsHistogramArray[11][1] = histogram_totalITSSignal;

  TString name_totalTPCBackground = "totalTPCBackground";
  TH1F* hist_totalTPCBackground = new TH1F(name_totalTPCBackground.Data(), "Total nr. of TPC hits for the daughters; number [#]; Entries", 1000, 0, 1000);
  hist_totalTPCBackground->Sumw2();
  hist_totalTPCBackground->SetLineColor(6);
  hist_totalTPCBackground->SetMarkerStyle(20);
  hist_totalTPCBackground->SetMarkerSize(0.6);
  hist_totalTPCBackground->SetMarkerColor(6);
  TH1F* histogram_totalTPCBackground = (TH1F*)hist_totalTPCBackground->Clone();
  fOutputBPlusResults->Add(histogram_totalTPCBackground);
  fResultsHistogramArray[11][2] = histogram_totalTPCBackground;

  TString name_totalTPCSignal = "totalTPCSignal";
  TH1F* hist_totalTPCSignal = new TH1F(name_totalTPCSignal.Data(), "Total nr. of TPC hits for the daughters; number [#]; Entries", 1000, 0, 1000);
  hist_totalTPCSignal->Sumw2();
  hist_totalTPCSignal->SetLineColor(6);
  hist_totalTPCSignal->SetMarkerStyle(20);
  hist_totalTPCSignal->SetMarkerSize(0.6);
  hist_totalTPCSignal->SetMarkerColor(6);
  TH1F* histogram_totalTPCSignal = (TH1F*)hist_totalTPCSignal->Clone();
  fOutputBPlusResults->Add(histogram_totalTPCSignal);
  fResultsHistogramArray[11][3] = histogram_totalTPCSignal;

  TString name_totalSigmaPIDBackground = "totalSigmaPIDBackground";
  TH1F* hist_totalSigmaPIDBackground = new TH1F(name_totalSigmaPIDBackground.Data(), "Total sigma of TPC and TOF PID for the daughters; number [#]; Entries", 1000, 0, 100);
  hist_totalSigmaPIDBackground->Sumw2();
  hist_totalSigmaPIDBackground->SetLineColor(6);
  hist_totalSigmaPIDBackground->SetMarkerStyle(20);
  hist_totalSigmaPIDBackground->SetMarkerSize(0.6);
  hist_totalSigmaPIDBackground->SetMarkerColor(6);
  TH1F* histogram_totalSigmaPIDBackground = (TH1F*)hist_totalSigmaPIDBackground->Clone();
  fOutputBPlusResults->Add(histogram_totalSigmaPIDBackground);
  fResultsHistogramArray[11][4] = histogram_totalSigmaPIDBackground;

  TString name_totalSigmaPIDSignal = "totalSigmaPIDSignal";
  TH1F* hist_totalSigmaPIDSignal = new TH1F(name_totalSigmaPIDSignal.Data(), "Total sigma of TPC and TOF PID for the daughters; number [#]; Entries", 1000, 0, 100);
  hist_totalSigmaPIDSignal->Sumw2();
  hist_totalSigmaPIDSignal->SetLineColor(6);
  hist_totalSigmaPIDSignal->SetMarkerStyle(20);
  hist_totalSigmaPIDSignal->SetMarkerSize(0.6);
  hist_totalSigmaPIDSignal->SetMarkerColor(6);
  TH1F* histogram_totalSigmaPIDSignal = (TH1F*)hist_totalSigmaPIDSignal->Clone();
  fOutputBPlusResults->Add(histogram_totalSigmaPIDSignal);
  fResultsHistogramArray[11][5] = histogram_totalSigmaPIDSignal;

  for (int i = 0; i < 3; ++i)
  {
    TString name_Histogram;
    TString discription_Histogram;
    if (i == 0) {name_Histogram = "invmassD0PionBPlusPion"; discription_Histogram = "Invariant mass D0 Pion and BPlus Pion; inv. mass [GeV/c^{2}]; Entries";}
    if (i == 1) {name_Histogram = "invmassD0KaonBPlusPion"; discription_Histogram = "Invariant mass D0 Kaon and BPlus Pion; inv. mass [GeV/c^{2}]; Entries";}
    if (i == 2) {name_Histogram = "invmassD0PionD0KaonBPlusPion"; discription_Histogram = "Invariant mass D0 Pion, D0 Kaon, and BPlus Pion; inv. mass [GeV/c^{2}]; Entries";}

    for (int j = 0; j < 2; ++j)
    {
      TString add_name = "";
      if (j == 1) add_name = "_Signal";
      name_Histogram += add_name;
      TH1F* histogram = new TH1F(name_Histogram.Data(), discription_Histogram.Data(), 1000, 0, 30);
      histogram->Sumw2();
      if (j % 2 == 0) histogram->SetLineColor(6);
      if (j % 2 == 1) histogram->SetLineColor(4);
      histogram->SetMarkerStyle(20);
      histogram->SetMarkerSize(0.6);
      if (j % 2 == 0) histogram->SetMarkerColor(6);
      if (j % 2 == 1) histogram->SetMarkerColor(4);
      TH1F* histogram_Clone = (TH1F*)histogram->Clone();
      fOutputBPlusResults->Add(histogram_Clone);
      fResultsHistogramArray[12][2*i + j] = histogram_Clone;
    }
  }



  return;
}
//-------------------------------------------------------------------------------------
AliAODVertex* AliAnalysisTaskSEBPlustoD0Pi::RecalculateVertex(const AliVVertex *primary, TObjArray *tracks, Double_t bField, Double_t dispersion, Bool_t optUseFitter, Bool_t optPropagate, Bool_t optUseDiamondConstraint) {
  //
  // Helper function to recalculate a vertex.
  //

  AliESDVertex *vertexESD = 0;
  AliAODVertex *vertexAOD = 0;

  AliVertexerTracks vertexer;
  vertexer.SetFieldkG(bField);

  vertexer.SetVtxStart((AliESDVertex*)primary); //primary vertex
  vertexESD = (AliESDVertex*)vertexer.VertexForSelectedESDTracks(tracks, optUseFitter, optPropagate, optUseDiamondConstraint);

  // delete vertexer; vertexer=NULL;

  if (!vertexESD) return vertexAOD;


  if (vertexESD->GetNContributors() != tracks->GetEntriesFast())
  {
    delete vertexESD; vertexESD = nullptr;
    return vertexAOD;
  }

  // convert to AliAODVertex
  Double_t pos[3], cov[6], chi2perNDF;
  for (Int_t a = 0; a < 3; a++)pos[a] = 0.;
  for (Int_t b = 0; b < 6; b++)cov[b] = 0.;
  chi2perNDF = 0;

  vertexESD->GetXYZ(pos); // position
  vertexESD->GetCovMatrix(cov); //covariance matrix


  Double_t vertRadius2 = pos[0] * pos[0] + pos[1] * pos[1];
  if (vertRadius2 > 8.) //(2.82)^2 radius beam pipe
  {
    delete vertexESD; vertexESD = nullptr;
    return vertexAOD;
  }

  chi2perNDF = vertexESD->GetChi2toNDF();
  dispersion = vertexESD->GetDispersion();
  delete vertexESD; vertexESD = nullptr;
  Int_t nprongs = 2; //tracks->GetEntriesFast();
  vertexAOD = new AliAODVertex(pos, cov, chi2perNDF, 0x0, -1, AliAODVertex::kUndef, nprongs);

  return vertexAOD;
}
//-------------------------------------------------------------------------------------
void AliAnalysisTaskSEBPlustoD0Pi::BPlustoD0PiSignalTracksInMC(TClonesArray * mcTrackArray, AliAODEvent*  /*aodevent*/, TMatrix * BPlustoD0PiLabelMatrix, TList *listout) {

  TMatrix &particleMatrix = *BPlustoD0PiLabelMatrix;
  for (Int_t i = 0; i < mcTrackArray->GetEntriesFast(); i++) {

    Int_t mcLabelPionBPlus = 0;
    Int_t mcLabelPionD0 = 0;
    Int_t mcLabelKaon = 0;
    Int_t mcLabelD0 = 0;
    Int_t mcLabelBPlus = 0;

    Double_t ptMC[5] = {0.0};
    Double_t yMC[5] = {0.0};
    Double_t pseudoYMC[5] = {0.0};

    Bool_t mcPionBPlusPresent = kFALSE;
    Bool_t mcPionD0Present = kFALSE;
    Bool_t mcKaonPresent = kFALSE;


    AliAODMCParticle *mcTrackParticle = dynamic_cast< AliAODMCParticle*>(mcTrackArray->At(i));
    if (!mcTrackParticle) {std::cout << "no particle" << std::endl; continue;}
    Int_t pdgCodeMC = TMath::Abs(mcTrackParticle->GetPdgCode());

    if (pdgCodeMC == 521)
    { //if the track is a BPlus we look at its daughters

      mcLabelBPlus = i;
      Int_t nDaughterBPlus = mcTrackParticle->GetNDaughters();
      ptMC[0] = mcTrackParticle->Pt();
      yMC[0] = mcTrackParticle->Y();
      pseudoYMC[0] = mcTrackParticle->Eta();

      // TString fillthis = "BPlus_in_analysis";
      ((TH1F*)fResultsHistogramArray[3][0])->Fill(0);

      if (nDaughterBPlus == 2)
      {
        for (Int_t iDaughterBPlus = 0; iDaughterBPlus < 2; iDaughterBPlus++)
        {
          AliAODMCParticle* daughterBPlus = (AliAODMCParticle*)mcTrackArray->At(mcTrackParticle->GetDaughterLabel(iDaughterBPlus));
          if (!daughterBPlus) break;
          Int_t pdgCodeDaughterBPlus = TMath::Abs(daughterBPlus->GetPdgCode());

          if (pdgCodeDaughterBPlus == 211) //if the track is a pion we save its monte carlo label
          {
            mcLabelPionBPlus = mcTrackParticle->GetDaughterLabel(iDaughterBPlus);
            mcPionBPlusPresent = kTRUE;
            ptMC[1] = daughterBPlus->Pt();
            yMC[1] = daughterBPlus->Y();
            pseudoYMC[1] = daughterBPlus->Eta();
          } else if (pdgCodeDaughterBPlus == 421) //if the track is a D0 we look at its daughters
          {
            mcLabelD0 = mcTrackParticle->GetDaughterLabel(iDaughterBPlus);
            Int_t nDaughterD0 = daughterBPlus->GetNDaughters();
            ptMC[2] = daughterBPlus->Pt();
            yMC[2] = daughterBPlus->Y();
            pseudoYMC[2] = daughterBPlus->Eta();

            if (nDaughterD0 == 2)
            {
              for (Int_t iDaughterD0 = 0; iDaughterD0 < 2; iDaughterD0++)
              {
                AliAODMCParticle* daughterD0 = (AliAODMCParticle*)mcTrackArray->At(daughterBPlus->GetDaughterLabel(iDaughterD0));
                if (!daughterD0) break;
                Int_t pdgCodeDaughterD0 = TMath::Abs(daughterD0->GetPdgCode());
                if (pdgCodeDaughterD0 == 211) //if the track is a pion we save its monte carlo label
                {
                  mcLabelPionD0 = daughterBPlus->GetDaughterLabel(iDaughterD0);
                  ptMC[3] = daughterD0->Pt();
                  yMC[3] = daughterD0->Y();
                  pseudoYMC[3] = daughterD0->Eta();
                  mcPionD0Present = kTRUE;
                } else if (pdgCodeDaughterD0 == 321) //if the track is a kaon we save its monte carlo label
                {
                  mcLabelKaon = daughterBPlus->GetDaughterLabel(iDaughterD0);
                  mcKaonPresent = kTRUE;
                  ptMC[4] = daughterD0->Pt();
                  yMC[4] = daughterD0->Y();
                  pseudoYMC[4] = daughterD0->Eta();
                } else break;
              }
            }
          } else break;
        }
      }
    }
    if (mcPionBPlusPresent && mcPionD0Present && mcKaonPresent) {

      // TString fillthis = "BPlus_in_analysis";
      ((TH1F*)fResultsHistogramArray[3][0])->Fill(1);

      
      for (Int_t j = 0; j < fnPtBins; ++j)
      { 
        // fillthis = "BPlus_per_bin";
        if (fPtBinLimits[j] < ptMC[0] && ptMC[0] < fPtBinLimits[j + 1]) {((TH1F*)fResultsHistogramArray[3][1])->Fill(j); break;}
      }

      ((TH1F*)fResultsHistogramArray[0][0])->Fill(ptMC[0]);
      ((TH1F*)fResultsHistogramArray[0][1])->Fill(ptMC[1]);
      ((TH1F*)fResultsHistogramArray[0][2])->Fill(ptMC[2]);
      ((TH1F*)fResultsHistogramArray[0][3])->Fill(ptMC[3]);
      ((TH1F*)fResultsHistogramArray[0][4])->Fill(ptMC[4]);

      ((TH1F*)fResultsHistogramArray[0][5])->Fill(yMC[0]);
      ((TH1F*)fResultsHistogramArray[0][6])->Fill(yMC[1]);
      ((TH1F*)fResultsHistogramArray[0][7])->Fill(yMC[2]);
      ((TH1F*)fResultsHistogramArray[0][8])->Fill(yMC[3]);
      ((TH1F*)fResultsHistogramArray[0][9])->Fill(yMC[4]);

      ((TH1F*)fResultsHistogramArray[0][10])->Fill(pseudoYMC[0]);
      ((TH1F*)fResultsHistogramArray[0][11])->Fill(pseudoYMC[1]);
      ((TH1F*)fResultsHistogramArray[0][12])->Fill(pseudoYMC[2]);
      ((TH1F*)fResultsHistogramArray[0][13])->Fill(pseudoYMC[3]);
      ((TH1F*)fResultsHistogramArray[0][14])->Fill(pseudoYMC[4]);

      // We check if the tracks are in acceptance
      if (ptMC[1] < 0.1 || TMath::Abs(pseudoYMC[1]) > 0.8 ) continue;
      if (ptMC[3] < 0.1 || TMath::Abs(pseudoYMC[3]) > 0.8 ) continue;
      if (ptMC[4] < 0.1 || TMath::Abs(pseudoYMC[4]) > 0.8 ) continue;

      // We check if the BPlus is in the fiducial region
      if (TMath::Abs(yMC[0]) > 0.8) continue;

      Int_t rows = BPlustoD0PiLabelMatrix->GetNrows();

      BPlustoD0PiLabelMatrix->ResizeTo(rows + 1, 5);
      particleMatrix(rows, 0) = mcLabelPionBPlus;
      particleMatrix(rows, 1) = mcLabelPionD0;
      particleMatrix(rows, 2) = mcLabelKaon;
      particleMatrix(rows, 3) = mcLabelD0;
      particleMatrix(rows, 4) = mcLabelBPlus;

      // fillthis = "BPlus_in_analysis";
      ((TH1F*)fResultsHistogramArray[3][0])->Fill(2);


      
      for (Int_t j = 0; j < fnPtBins; ++j)
      {
        // fillthis = "BPlus_per_bin_in_Acc";
        if (fPtBinLimits[j] < ptMC[0] && ptMC[0] < fPtBinLimits[j + 1]) {((TH1F*)fResultsHistogramArray[3][2])->Fill(j); break;}
      }
    }
  }


  return;
}
//-------------------------------------------------------------------------------------
Bool_t AliAnalysisTaskSEBPlustoD0Pi::D0FirstDaughterSelection(AliAODTrack* aodTrack, AliAODVertex *primaryVertex, Double_t bz, TClonesArray * mcTrackArray, TMatrix * BPlustoD0PiLabelMatrix, AliAODMCHeader * header, AliAODEvent* aodEvent) {

  // we select the D0 pion and save its information
  if (!aodTrack) AliFatal("Not a standard AOD");

  //quick quality cut
  if (aodTrack->GetITSNcls() < 1) return kFALSE;
  // if (aodTrack->GetTPCNcls() < 1) return kFALSE;
  if (aodTrack->GetStatus()&AliESDtrack::kITSpureSA) return kFALSE;
  if (!(aodTrack->GetStatus()&AliESDtrack::kITSin)) return kFALSE;
  if (aodTrack->GetID() < 0) return kFALSE;
  Double_t covtest[21];
  if (!aodTrack->GetCovarianceXYZPxPyPz(covtest)) return kFALSE;

  Int_t mcLabelParticle = -1;
  mcLabelParticle = aodTrack->GetLabel();

  // we fill histograms with information of the track
  Double_t pt_track = aodTrack->Pt();
  Double_t momentum_track = aodTrack->P();
  Int_t numberOfITS = aodTrack->GetITSNcls();
  Int_t numberOfTPC = aodTrack->GetTPCNcls();

  AliExternalTrackParam particleTrack;
  particleTrack.CopyFromVTrack(aodTrack);
  Double_t d0[2], covd0[3];
  particleTrack.PropagateToDCA(primaryVertex, bz, 100., d0, covd0);

  //we check if the particle is a signal track, we look at both daughter options therefore the signal will be too high in the D0 daughter signal histograms
  Bool_t isDesiredCandidate = kFALSE;
  if (fUseMCInfo) {
    TMatrix &particleMatrix = *BPlustoD0PiLabelMatrix;
    for (Int_t k = 0; k < BPlustoD0PiLabelMatrix->GetNrows(); ++k) {
      if (mcLabelParticle == (Int_t)particleMatrix(k, 1) || mcLabelParticle == (Int_t)particleMatrix(k, 2)) {
        isDesiredCandidate = kTRUE;
        break;
      }
    }
  }

  if (fUseMCInfo) {
    if (IsTrackInjected(aodTrack, header, mcTrackArray) && !isDesiredCandidate && fQuickSignalAnalysis == 2) return kFALSE;
  }

  Int_t daughterType = 0;


  Int_t histType = 0;
  ((TH1F*)fDaughterHistogramArray[daughterType][histType][0])->Fill(pt_track);
  ((TH1F*)fDaughterHistogramArray[daughterType][histType][1])->Fill(momentum_track);
  ((TH1F*)fDaughterHistogramArray[daughterType][histType][2])->Fill(numberOfITS);
  ((TH1F*)fDaughterHistogramArray[daughterType][histType][3])->Fill(numberOfTPC);

  for (Int_t j = 0; j < 10; ++j)
  {
    if (aodTrack->HasPointOnITSLayer(j)) ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(j);

  }
  ((TH1F*)fDaughterHistogramArray[daughterType][histType][8])->Fill(d0[0]);

  if (isDesiredCandidate)
  {
    histType = 1;
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][0])->Fill(pt_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][1])->Fill(momentum_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][2])->Fill(numberOfITS);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][3])->Fill(numberOfTPC);

    for (Int_t j = 0; j < 10; ++j)
    {
      if (aodTrack->HasPointOnITSLayer(j)) ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(j);

    }
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][8])->Fill(d0[0]);
  }

  //we apply a number of cuts on the particle
  Bool_t bCut = kFALSE;

  //We do not apply a PID cut at this stage since we don't know if we are dealing with a kaon or a pion

  if (aodTrack->GetITSNcls() < fCuts->GetMinITSNclsD0FirstDaughter()) {
    if (isDesiredCandidate) {
      ((TH1F*)fDaughterHistogramArrayExtra[0][1])->Fill(3);
    } else ((TH1F*)fDaughterHistogramArrayExtra[0][0])->Fill(3);
    bCut = kTRUE;
  }

  // TPC cluster cut turned off, cut is done with crossed row method instead using the IsThisDaughterSelected function
  // if (aodTrack->GetTPCNcls() < fCuts->GetMinTPCNclsD0FirstDaughter()) {
  if (!fCuts->IsThisDaughterSelected(aodTrack,primaryVertex,aodEvent)) {
    if (isDesiredCandidate) {
      ((TH1F*)fDaughterHistogramArrayExtra[0][1])->Fill(4);
    } else ((TH1F*)fDaughterHistogramArrayExtra[0][0])->Fill(4);
    bCut = kTRUE;
  }

  if (fCuts->UseITSRefitD0FirstDaughter() == kTRUE) {
    if (!(aodTrack->GetStatus()&AliESDtrack::kITSrefit)) {
      if (isDesiredCandidate) {
        ((TH1F*)fDaughterHistogramArrayExtra[0][1])->Fill(5);
      } else ((TH1F*)fDaughterHistogramArrayExtra[0][0])->Fill(5);
      bCut = kTRUE;
    }
  }

  if (fCuts->UseTPCRefitD0FirstDaughter() == kTRUE) {
    if ((!(aodTrack->GetStatus()&AliESDtrack::kTPCrefit))) {
      if (isDesiredCandidate) {
        ((TH1F*)fDaughterHistogramArrayExtra[0][1])->Fill(6);
      } else ((TH1F*)fDaughterHistogramArrayExtra[0][0])->Fill(6);
      bCut = kTRUE;
    }
  }

  if (fCuts->UseFilterBitD0FirstDaughter() == kTRUE) {
    if (!(aodTrack->TestFilterMask(BIT(fCuts->GetFilterBitD0FirstDaughter())))) {
      if (isDesiredCandidate) {
        ((TH1F*)fDaughterHistogramArrayExtra[0][1])->Fill(7);
      } else ((TH1F*)fDaughterHistogramArrayExtra[0][0])->Fill(7);
      bCut = kTRUE;
    }
  }

  if (aodTrack->Pt() < fCuts->GetMinPtD0FirstDaughter()) {
    if (isDesiredCandidate) {
      ((TH1F*)fDaughterHistogramArrayExtra[0][1])->Fill(8);
    } else ((TH1F*)fDaughterHistogramArrayExtra[0][0])->Fill(8);
    bCut = kTRUE;
  }

  if (TMath::Abs(d0[0]) < fCuts->GetMind0D0FirstDaughter()) {
    if (isDesiredCandidate) {
      ((TH1F*)fDaughterHistogramArrayExtra[0][1])->Fill(12);
    } else ((TH1F*)fDaughterHistogramArrayExtra[0][0])->Fill(12);
    bCut = kTRUE;
  }

  if (TMath::Abs(aodTrack->Eta()) > fCuts->GetMaxAbsEtaD0FirstDaughter()) {
    if (isDesiredCandidate) {
      ((TH1F*)fDaughterHistogramArrayExtra[0][1])->Fill(9);
    } else ((TH1F*)fDaughterHistogramArrayExtra[0][0])->Fill(9);
    bCut = kTRUE;
  }

  Bool_t bHardSelectionArrayITS[7] = {kFALSE};
  fCuts->GetHardSelectionArrayITSD0FirstDaughter(bHardSelectionArrayITS);
  Bool_t bSoftSelectionArrayITS[7] = {kFALSE};
  fCuts->GetSoftSelectionArrayITSD0FirstDaughter(bSoftSelectionArrayITS);

  Bool_t bHardITSPass = kTRUE;
  for (Int_t j = 0; j < 7; ++j)
  {
    if (bHardSelectionArrayITS[j])
    {
      if (!aodTrack->HasPointOnITSLayer(j)) bHardITSPass = kFALSE;
    }
  }

  Int_t nCounterSoftSelection = 0;
  Bool_t bSoftITSPass = kTRUE;
  for (Int_t j = 0; j < 7; ++j)
  {
    if (bSoftSelectionArrayITS[j])
    {
      if (aodTrack->HasPointOnITSLayer(j)) nCounterSoftSelection++;
    }
  }
  if (nCounterSoftSelection < fCuts->GetNSoftITSCutD0FirstDaughter()) bSoftITSPass = kFALSE;

  if (!bHardITSPass) {
    if (isDesiredCandidate) {
      ((TH1F*)fDaughterHistogramArrayExtra[0][1])->Fill(10);
    } else ((TH1F*)fDaughterHistogramArrayExtra[0][0])->Fill(10);
    bCut = kTRUE;
  }

  if (!bSoftITSPass) {
    if (isDesiredCandidate) {
      ((TH1F*)fDaughterHistogramArrayExtra[0][1])->Fill(11);
    } else ((TH1F*)fDaughterHistogramArrayExtra[0][0])->Fill(11);
    bCut = kTRUE;
  }

  if (!isDesiredCandidate && fQuickSignalAnalysis == 1) {
    if (isDesiredCandidate) {
      ((TH1F*)fDaughterHistogramArrayExtra[0][1])->Fill(13);
    } else ((TH1F*)fDaughterHistogramArrayExtra[0][0])->Fill(13);
    bCut = kTRUE;
  }

  if (bCut) {
    if (isDesiredCandidate) {
      ((TH1F*)fDaughterHistogramArrayExtra[0][1])->Fill(0);
    } else ((TH1F*)fDaughterHistogramArrayExtra[0][0])->Fill(0);
    return kFALSE;
  }

  //we fill histograms with track information of the tracks that pass the cuts
  histType = 2;
  ((TH1F*)fDaughterHistogramArray[daughterType][histType][0])->Fill(pt_track);
  ((TH1F*)fDaughterHistogramArray[daughterType][histType][1])->Fill(momentum_track);
  ((TH1F*)fDaughterHistogramArray[daughterType][histType][2])->Fill(numberOfITS);
  ((TH1F*)fDaughterHistogramArray[daughterType][histType][3])->Fill(numberOfTPC);

  for (Int_t j = 0; j < 10; ++j)
  {
    if (aodTrack->HasPointOnITSLayer(j)) ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(j);

  }
  ((TH1F*)fDaughterHistogramArray[daughterType][histType][8])->Fill(d0[0]);

  if (isDesiredCandidate)
  {
    histType = 3;
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][0])->Fill(pt_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][1])->Fill(momentum_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][2])->Fill(numberOfITS);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][3])->Fill(numberOfTPC);

    for (Int_t j = 0; j < 10; ++j)
    {
      if (aodTrack->HasPointOnITSLayer(j)) ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(j);

    }
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][8])->Fill(d0[0]);
  }

  return kTRUE;
}
//-------------------------------------------------------------------------------------
Bool_t AliAnalysisTaskSEBPlustoD0Pi::D0SecondDaughterSelection(AliAODTrack* aodTrack, AliAODVertex *primaryVertex, Double_t bz, TClonesArray * mcTrackArray, TMatrix * BPlustoD0PiLabelMatrix, AliAODMCHeader * header, AliAODEvent* aodEvent) {

  // we select the D0 pion and save its information
  if (!aodTrack) AliFatal("Not a standard AOD");

  //quick quality cut
  if (aodTrack->GetITSNcls() < 1) return kFALSE;
  // if (aodTrack->GetTPCNcls() < 1) return kFALSE;
  if (aodTrack->GetStatus()&AliESDtrack::kITSpureSA) return kFALSE;
  if (!(aodTrack->GetStatus()&AliESDtrack::kITSin)) return kFALSE;
  if (aodTrack->GetID() < 0) return kFALSE;
  Double_t covtest[21];
  if (!aodTrack->GetCovarianceXYZPxPyPz(covtest)) return kFALSE;

  Int_t mcLabelParticle = -1;
  mcLabelParticle = aodTrack->GetLabel();

  // we fill histograms with information of the track
  Double_t pt_track = aodTrack->Pt();
  Double_t momentum_track = aodTrack->P();
  Int_t numberOfITS = aodTrack->GetITSNcls();
  Int_t numberOfTPC = aodTrack->GetTPCNcls();

  AliExternalTrackParam particleTrack;
  particleTrack.CopyFromVTrack(aodTrack);
  Double_t d0[2], covd0[3];
  particleTrack.PropagateToDCA(primaryVertex, bz, 100., d0, covd0);

  //we check if the particle is a signal track, we look at both daughter options therefore the signal will be too high in the D0 daughter signal histograms
  Bool_t isDesiredCandidate = kFALSE;
  if (fUseMCInfo) {
    TMatrix &particleMatrix = *BPlustoD0PiLabelMatrix;
    for (Int_t k = 0; k < BPlustoD0PiLabelMatrix->GetNrows(); ++k) {
      if (mcLabelParticle == (Int_t)particleMatrix(k, 1) || mcLabelParticle == (Int_t)particleMatrix(k, 2)) {
        isDesiredCandidate = kTRUE;
        break;
      }
    }
  }

  if (fUseMCInfo) {
    if (IsTrackInjected(aodTrack, header, mcTrackArray) && !isDesiredCandidate && fQuickSignalAnalysis == 2) return kFALSE;
  }

  Int_t daughterType = 1;


  Int_t histType = 0;
  ((TH1F*)fDaughterHistogramArray[daughterType][histType][0])->Fill(pt_track);
  ((TH1F*)fDaughterHistogramArray[daughterType][histType][1])->Fill(momentum_track);
  ((TH1F*)fDaughterHistogramArray[daughterType][histType][2])->Fill(numberOfITS);
  ((TH1F*)fDaughterHistogramArray[daughterType][histType][3])->Fill(numberOfTPC);

  for (Int_t j = 0; j < 10; ++j)
  {
    if (aodTrack->HasPointOnITSLayer(j)) ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(j);

  }
  ((TH1F*)fDaughterHistogramArray[daughterType][histType][8])->Fill(d0[0]);

  if (isDesiredCandidate)
  {
    histType = 1;
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][0])->Fill(pt_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][1])->Fill(momentum_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][2])->Fill(numberOfITS);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][3])->Fill(numberOfTPC);

    for (Int_t j = 0; j < 10; ++j)
    {
      if (aodTrack->HasPointOnITSLayer(j)) ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(j);
    }
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][8])->Fill(d0[0]);
  }

  //we apply a number of cuts on the particle
  Bool_t bCut = kFALSE;

  //We do not apply a PID cut at this stage since we don't know if we are dealing with a kaon or a pion

  if (aodTrack->GetITSNcls() < fCuts->GetMinITSNclsD0SecondDaughter()) {
    if (isDesiredCandidate) {
      ((TH1F*)fDaughterHistogramArrayExtra[1][1])->Fill(3);
    } else ((TH1F*)fDaughterHistogramArrayExtra[1][0])->Fill(3);
    bCut = kTRUE;
  }

  // TPC cluster cut turned off, cut is done with crossed row method instead using the IsThisDaughterSelected function
  // if (aodTrack->GetTPCNcls() < fCuts->GetMinTPCNclsD0SecondDaughter()) {
  if (!fCuts->IsThisDaughterSelected(aodTrack,primaryVertex,aodEvent)) {
    if (isDesiredCandidate) {
      ((TH1F*)fDaughterHistogramArrayExtra[1][1])->Fill(4);
    } else ((TH1F*)fDaughterHistogramArrayExtra[1][0])->Fill(4);
    bCut = kTRUE;
  }

  if (fCuts->UseITSRefitD0SecondDaughter() == kTRUE) {
    if (!(aodTrack->GetStatus()&AliESDtrack::kITSrefit)) {
      if (isDesiredCandidate) {
        ((TH1F*)fDaughterHistogramArrayExtra[1][1])->Fill(5);
      } else ((TH1F*)fDaughterHistogramArrayExtra[1][0])->Fill(5);
      bCut = kTRUE;
    }
  }

  if (fCuts->UseTPCRefitD0SecondDaughter() == kTRUE) {
    if ((!(aodTrack->GetStatus()&AliESDtrack::kTPCrefit))) {
      if (isDesiredCandidate) {
        ((TH1F*)fDaughterHistogramArrayExtra[1][1])->Fill(6);
      } else ((TH1F*)fDaughterHistogramArrayExtra[1][0])->Fill(6);
      bCut = kTRUE;
    }
  }

  if (fCuts->UseFilterBitD0SecondDaughter() == kTRUE) {
    if (!(aodTrack->TestFilterMask(BIT(fCuts->GetFilterBitD0SecondDaughter())))) {
      if (isDesiredCandidate) {
        ((TH1F*)fDaughterHistogramArrayExtra[1][1])->Fill(7);
      } else ((TH1F*)fDaughterHistogramArrayExtra[1][0])->Fill(7);
      bCut = kTRUE;
    }
  }

  if (aodTrack->Pt() < fCuts->GetMinPtD0SecondDaughter()) {
    if (isDesiredCandidate) {
      ((TH1F*)fDaughterHistogramArrayExtra[1][1])->Fill(8);
    } else ((TH1F*)fDaughterHistogramArrayExtra[1][0])->Fill(8);
    bCut = kTRUE;
  }

  if (TMath::Abs(d0[0]) < fCuts->GetMind0D0SecondDaughter()) {
    if (isDesiredCandidate) {
      ((TH1F*)fDaughterHistogramArrayExtra[1][1])->Fill(12);
    } else ((TH1F*)fDaughterHistogramArrayExtra[1][0])->Fill(12);
    bCut = kTRUE;
  }

  if (TMath::Abs(aodTrack->Eta()) > fCuts->GetMaxAbsEtaD0SecondDaughter()) {
    if (isDesiredCandidate) {
      ((TH1F*)fDaughterHistogramArrayExtra[1][1])->Fill(9);
    } else ((TH1F*)fDaughterHistogramArrayExtra[1][0])->Fill(9);
    bCut = kTRUE;
  }

  Bool_t bHardSelectionArrayITS[7] = {kFALSE};
  fCuts->GetHardSelectionArrayITSD0SecondDaughter(bHardSelectionArrayITS);
  Bool_t bSoftSelectionArrayITS[7] = {kFALSE};
  fCuts->GetSoftSelectionArrayITSD0SecondDaughter(bSoftSelectionArrayITS);

  Bool_t bHardITSPass = kTRUE;
  for (Int_t j = 0; j < 7; ++j)
  {
    if (bHardSelectionArrayITS[j])
    {
      if (!aodTrack->HasPointOnITSLayer(j)) bHardITSPass = kFALSE;
    }
  }

  Int_t nCounterSoftSelection = 0;
  Bool_t bSoftITSPass = kTRUE;
  for (Int_t j = 0; j < 7; ++j)
  {
    if (bSoftSelectionArrayITS[j])
    {
      if (aodTrack->HasPointOnITSLayer(j)) nCounterSoftSelection++;
    }
  }
  if (nCounterSoftSelection < fCuts->GetNSoftITSCutD0SecondDaughter()) bSoftITSPass = kFALSE;

  if (!bHardITSPass) {
    if (isDesiredCandidate) {
      ((TH1F*)fDaughterHistogramArrayExtra[1][1])->Fill(10);
    } else ((TH1F*)fDaughterHistogramArrayExtra[1][0])->Fill(10);
    bCut = kTRUE;
  }

  if (!bSoftITSPass) {
    if (isDesiredCandidate) {
      ((TH1F*)fDaughterHistogramArrayExtra[1][1])->Fill(11);
    } else ((TH1F*)fDaughterHistogramArrayExtra[1][0])->Fill(11);
    bCut = kTRUE;
  }

  if (!isDesiredCandidate && fQuickSignalAnalysis == 1) bCut = kTRUE;

  if (bCut) {
    if (isDesiredCandidate) {
      ((TH1F*)fDaughterHistogramArrayExtra[1][1])->Fill(0);
    } else ((TH1F*)fDaughterHistogramArrayExtra[1][0])->Fill(0);
    return kFALSE;
  }

  //we fill histograms with track information of the tracks that pass the cuts
  histType = 2;
  ((TH1F*)fDaughterHistogramArray[daughterType][histType][0])->Fill(pt_track);
  ((TH1F*)fDaughterHistogramArray[daughterType][histType][1])->Fill(momentum_track);
  ((TH1F*)fDaughterHistogramArray[daughterType][histType][2])->Fill(numberOfITS);
  ((TH1F*)fDaughterHistogramArray[daughterType][histType][3])->Fill(numberOfTPC);

  for (Int_t j = 0; j < 10; ++j)
  {
    if (aodTrack->HasPointOnITSLayer(j)) ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(j);

  }
  ((TH1F*)fDaughterHistogramArray[daughterType][histType][8])->Fill(d0[0]);

  if (isDesiredCandidate)
  {
    histType = 3;
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][0])->Fill(pt_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][1])->Fill(momentum_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][2])->Fill(numberOfITS);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][3])->Fill(numberOfTPC);

    for (Int_t j = 0; j < 10; ++j)
    {
      if (aodTrack->HasPointOnITSLayer(j)) ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(j);

    }
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][8])->Fill(d0[0]);
  }

  return kTRUE;
}
//-------------------------------------------------------------------------------------
void AliAnalysisTaskSEBPlustoD0Pi::BPlusPionSelection(AliAODEvent* aodEvent, AliAODVertex *primaryVertex, Double_t bz, TClonesArray * mcTrackArray, TMatrix * BPlustoD0PiLabelMatrix, AliAODMCHeader * header) {

  //we keep track of the number of particles we could use and how many we actually use after cuts
  Int_t numberofparticles = 0;
  Int_t numberofparticlesused = 0;

  TString fillthis = "";

  //we loop over all tracks in the event
  for (Int_t i = 0; i < aodEvent->GetNumberOfTracks(); i++) {
    AliAODTrack* aodTrack = dynamic_cast<AliAODTrack*>(aodEvent->GetTrack(i));
    if (!aodTrack) AliFatal("Not a standard AOD");

    //quick quality cut
    if (aodTrack->GetITSNcls() < 1) continue;
    // if (aodTrack->GetTPCNcls() < 1) continue;
    if (aodTrack->GetStatus()&AliESDtrack::kITSpureSA) continue;
    if (!(aodTrack->GetStatus()&AliESDtrack::kITSin)) continue;
    if (aodTrack->GetID() < 0) continue;
    Double_t covtest[21];
    if (!aodTrack->GetCovarianceXYZPxPyPz(covtest)) continue;

    // Double_t pos[3],cov[6];
    // primaryVertex->GetXYZ(pos);
    // primaryVertex->GetCovarianceMatrix(cov);
    // const AliESDVertex vESD(pos,cov,100.,100);
    // if(!fCuts->IsDaughterSelected(aodTrack,&vESD,fCuts->GetTrackCuts(),aodEvent)) continue;

    Int_t mcLabelParticle = -1;
    mcLabelParticle = aodTrack->GetLabel();

    numberofparticles++;

    // we fill histograms with information of the track
    Double_t pt_track = aodTrack->Pt();
    Double_t momentum_track = aodTrack->P();
    Int_t numberOfITS = aodTrack->GetITSNcls();
    Int_t numberOfTPC = aodTrack->GetTPCNcls();
    Double_t nSigmaTPC = 0;
    Double_t nSigmaTOF = 0;
    Int_t pionPIDnumber = 2;
    Int_t TPCok = 0;
    Int_t TOFok = 0;

    AliAODPidHF* trackPIDHF = (AliAODPidHF*)fCuts->GetPidHF();
    if(trackPIDHF) TPCok = trackPIDHF->GetnSigmaTPC(aodTrack, pionPIDnumber, nSigmaTPC);
    if(trackPIDHF) TOFok = trackPIDHF->GetnSigmaTOF(aodTrack, pionPIDnumber, nSigmaTOF);

    AliExternalTrackParam particleTrack;
    particleTrack.CopyFromVTrack(aodTrack);
    Double_t d0[2], covd0[3];
    particleTrack.PropagateToDCA(primaryVertex, bz, 100., d0, covd0);


    //we check if the particle is a signal track
    Bool_t isDesiredCandidate = kFALSE;
    if (fUseMCInfo) {
      TMatrix &particleMatrix = *BPlustoD0PiLabelMatrix;
      for (Int_t k = 0; k < BPlustoD0PiLabelMatrix->GetNrows(); ++k) {
        if (mcLabelParticle == (Int_t)particleMatrix(k, 0)) {
          isDesiredCandidate = kTRUE;
          break;
        }
      }
    }

    if (fUseMCInfo) {
      if (IsTrackInjected(aodTrack, header, mcTrackArray) && !isDesiredCandidate && fQuickSignalAnalysis == 2) continue;
    }

    Int_t daughterType = 2;

    Int_t histType = 0;
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][0])->Fill(pt_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][1])->Fill(momentum_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][2])->Fill(numberOfITS);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][3])->Fill(numberOfTPC);

    for (Int_t j = 0; j < 10; ++j)
    {
      if (aodTrack->HasPointOnITSLayer(j)) ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(j);

    }

    if (TPCok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(nSigmaTPC);
    if (TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][5])->Fill(nSigmaTOF);
    if (TPCok != -1 && TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][6])->Fill(TMath::Sqrt(nSigmaTPC * nSigmaTPC + nSigmaTOF * nSigmaTOF));
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][8])->Fill(d0[0]);

    if (isDesiredCandidate)
    {
      histType = 1;
      ((TH1F*)fDaughterHistogramArray[daughterType][histType][0])->Fill(pt_track);
      ((TH1F*)fDaughterHistogramArray[daughterType][histType][1])->Fill(momentum_track);
      ((TH1F*)fDaughterHistogramArray[daughterType][histType][2])->Fill(numberOfITS);
      ((TH1F*)fDaughterHistogramArray[daughterType][histType][3])->Fill(numberOfTPC);

      for (Int_t j = 0; j < 10; ++j)
      {
        if (aodTrack->HasPointOnITSLayer(j)) ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(j);

      }

      if (TPCok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(nSigmaTPC);
      if (TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][5])->Fill(nSigmaTOF);
      if (TPCok != -1 && TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][6])->Fill(TMath::Sqrt(nSigmaTPC * nSigmaTPC + nSigmaTOF * nSigmaTOF));
      ((TH1F*)fDaughterHistogramArray[daughterType][histType][8])->Fill(d0[0]);
    }


    //we apply a number of cuts on the particle
    Bool_t bCut = kFALSE;


    //we apply a PID cut, 2 is used to indicate we look for a pion
    if (!(fCuts->SelectPID(aodTrack, 2))) {
      if (isDesiredCandidate) {
        ((TH1F*)fDaughterHistogramArrayExtra[2][1])->Fill(2);
      } else ((TH1F*)fDaughterHistogramArrayExtra[2][0])->Fill(2);
      bCut = kTRUE;
    }

    if (aodTrack->GetITSNcls() < fCuts->GetMinITSNclsBPlusPion()) {
      if (isDesiredCandidate) {
        ((TH1F*)fDaughterHistogramArrayExtra[2][1])->Fill(3);
      } else ((TH1F*)fDaughterHistogramArrayExtra[2][0])->Fill(3);
      bCut = kTRUE;
    }

    // TPC cluster cut turned off, cut is done with crossed row method instead using the IsThisDaughterSelected function
    // if (aodTrack->GetTPCNcls() < fCuts->GetMinTPCNclsBPlusPion()) {
    if (!fCuts->IsThisDaughterSelected(aodTrack,primaryVertex,aodEvent)) {
      if (isDesiredCandidate) {
        ((TH1F*)fDaughterHistogramArrayExtra[2][1])->Fill(4);
      } else ((TH1F*)fDaughterHistogramArrayExtra[2][0])->Fill(4);
      bCut = kTRUE;
    }

    if (fCuts->UseITSRefitBPlusPion() == kTRUE) {
      if (!(aodTrack->GetStatus()&AliESDtrack::kITSrefit)) {
        if (isDesiredCandidate) {
          ((TH1F*)fDaughterHistogramArrayExtra[2][1])->Fill(5);
        } else ((TH1F*)fDaughterHistogramArrayExtra[2][0])->Fill(5);
        bCut = kTRUE;
      }
    }

    if (fCuts->UseTPCRefitBPlusPion() == kTRUE) {
      if ((!(aodTrack->GetStatus()&AliESDtrack::kTPCrefit))) {
        if (isDesiredCandidate) {
          ((TH1F*)fDaughterHistogramArrayExtra[2][1])->Fill(6);
        } else ((TH1F*)fDaughterHistogramArrayExtra[2][0])->Fill(6);
        bCut = kTRUE;
      }
    }

    if (fCuts->UseFilterBitBPlusPion() == kTRUE) {
      if (!(aodTrack->TestFilterMask(BIT(fCuts->GetFilterBitBPlusPion())))) {
        if (isDesiredCandidate) {
          ((TH1F*)fDaughterHistogramArrayExtra[2][1])->Fill(7);
        } else ((TH1F*)fDaughterHistogramArrayExtra[2][0])->Fill(7);
        bCut = kTRUE;
      }
    }


    if (aodTrack->Pt() < fCuts->GetMinPtBPlusPion()) {
      if (isDesiredCandidate) {
        ((TH1F*)fDaughterHistogramArrayExtra[2][1])->Fill(8);
      } else ((TH1F*)fDaughterHistogramArrayExtra[2][0])->Fill(8);
      bCut = kTRUE;
    }


    if (TMath::Abs(d0[0]) < fCuts->GetMind0BPlusPion()) {
      if (isDesiredCandidate) {
        ((TH1F*)fDaughterHistogramArrayExtra[2][1])->Fill(12);
      } else ((TH1F*)fDaughterHistogramArrayExtra[2][0])->Fill(12);
      bCut = kTRUE;
    }

    if (TMath::Abs(aodTrack->Eta()) > fCuts->GetMaxAbsEtaBPlusPion()) {
      if (isDesiredCandidate) {
        ((TH1F*)fDaughterHistogramArrayExtra[2][1])->Fill(9);
      } else ((TH1F*)fDaughterHistogramArrayExtra[2][0])->Fill(9);
      bCut = kTRUE;
    }

    Bool_t bHardSelectionArrayITS[7] = {kFALSE};
    fCuts->GetHardSelectionArrayITSBPlusPion(bHardSelectionArrayITS);
    Bool_t bSoftSelectionArrayITS[7] = {kFALSE};
    fCuts->GetSoftSelectionArrayITSBPlusPion(bSoftSelectionArrayITS);

    Bool_t bHardITSPass = kTRUE;
    for (Int_t j = 0; j < 7; ++j)
    {
      if (bHardSelectionArrayITS[j])
      {
        if (!aodTrack->HasPointOnITSLayer(j)) bHardITSPass = kFALSE;
      }
    }

    Int_t nCounterSoftSelection = 0;
    Bool_t bSoftITSPass = kTRUE;
    for (Int_t j = 0; j < 7; ++j)
    {
      if (bSoftSelectionArrayITS[j])
      {
        if (aodTrack->HasPointOnITSLayer(j)) nCounterSoftSelection++;
      }
    }
    if (nCounterSoftSelection < fCuts->GetNSoftITSCutBPlusPion()) bSoftITSPass = kFALSE;

    if (!bHardITSPass) {
      if (isDesiredCandidate) {
        ((TH1F*)fDaughterHistogramArrayExtra[2][1])->Fill(10);
      } else ((TH1F*)fDaughterHistogramArrayExtra[2][0])->Fill(10);
      bCut = kTRUE;
    }

    if (!bSoftITSPass) {
      if (isDesiredCandidate) {
        ((TH1F*)fDaughterHistogramArrayExtra[2][1])->Fill(11);
      } else ((TH1F*)fDaughterHistogramArrayExtra[2][0])->Fill(11);
      bCut = kTRUE;
    }


    if (!isDesiredCandidate && fQuickSignalAnalysis == 1) bCut = kTRUE;

    if (bCut) {
      if (isDesiredCandidate) {
        ((TH1F*)fDaughterHistogramArrayExtra[2][1])->Fill(0);
      } else ((TH1F*)fDaughterHistogramArrayExtra[2][0])->Fill(0);
      continue;
    }

    //we fill histograms with track information of the tracks that pass the cuts
    histType = 2;
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][0])->Fill(pt_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][1])->Fill(momentum_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][2])->Fill(numberOfITS);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][3])->Fill(numberOfTPC);

    for (Int_t j = 0; j < 10; ++j)
    {
      if (aodTrack->HasPointOnITSLayer(j)) ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(j);

    }

    if (TPCok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(nSigmaTPC);
    if (TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][5])->Fill(nSigmaTOF);
    if (TPCok != -1 && TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][6])->Fill(TMath::Sqrt(nSigmaTPC * nSigmaTPC + nSigmaTOF * nSigmaTOF));
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][8])->Fill(d0[0]);

    if (isDesiredCandidate)
    {
      histType = 3;
      ((TH1F*)fDaughterHistogramArray[daughterType][histType][0])->Fill(pt_track);
      ((TH1F*)fDaughterHistogramArray[daughterType][histType][1])->Fill(momentum_track);
      ((TH1F*)fDaughterHistogramArray[daughterType][histType][2])->Fill(numberOfITS);
      ((TH1F*)fDaughterHistogramArray[daughterType][histType][3])->Fill(numberOfTPC);

      for (Int_t j = 0; j < 10; ++j)
      {
        if (aodTrack->HasPointOnITSLayer(j)) ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(j);

      }

      if (TPCok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(nSigmaTPC);
      if (TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][5])->Fill(nSigmaTOF);
      if (TPCok != -1 && TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][6])->Fill(TMath::Sqrt(nSigmaTPC * nSigmaTPC + nSigmaTOF * nSigmaTOF));
      ((TH1F*)fDaughterHistogramArray[daughterType][histType][8])->Fill(d0[0]);
    }

    fBPlusPionTracks->push_back(i);
    numberofparticlesused++;
  }

  // cout << "BPlus pions used: " << numberofparticlesused << endl;

  ((TH1F*)fDaughterHistogramArray[2][0][10])->Fill(numberofparticles);
  ((TH1F*)fDaughterHistogramArray[2][1][10])->Fill(numberofparticlesused);
  return;
}
//-------------------------------------------------------------------------------------
void AliAnalysisTaskSEBPlustoD0Pi::D0Selection(AliAODEvent* aodEvent, AliAODVertex *primaryVertex, Double_t bz, TClonesArray * mcTrackArray, TMatrix * BPlustoD0PiLabelMatrix, TClonesArray * D0TracksFromFriendFile, AliAODMCHeader * header) {

  TString fillthis = "";

  AliAnalysisVertexingHF *vHF = new AliAnalysisVertexingHF();

  //next we loop over all the D0 candidates
  for (Int_t j = 0; j < D0TracksFromFriendFile->GetEntriesFast(); j++)
  {

    //we get the track of the D0
    AliAODRecoDecayHF2Prong * trackD0 = (AliAODRecoDecayHF2Prong*)(D0TracksFromFriendFile->At(j));
    if (!trackD0) {std::cout << "found none" << std::endl; continue;}
    if (trackD0 == nullptr) {std::cout << "found nullptr" << std::endl; continue;}

    if (!(vHF->FillRecoCand(aodEvent, trackD0))) //Fill the data members of the candidate only if they are empty.
    {
      fCEvents->Fill(12); //monitor how often this fails
      continue;
    }

    AliAODTrack * trackFirstDaughter = (AliAODTrack*)(trackD0->GetDaughter(0));
    AliAODTrack * trackSecondDaughter = (AliAODTrack*)(trackD0->GetDaughter(1));
    if (!D0FirstDaughterSelection(trackFirstDaughter, primaryVertex, bz, mcTrackArray, BPlustoD0PiLabelMatrix, header, aodEvent)) continue;
    if (!D0SecondDaughterSelection(trackSecondDaughter, primaryVertex, bz, mcTrackArray, BPlustoD0PiLabelMatrix, header, aodEvent)) continue;


    AliAODVertex *vertexMother = (AliAODVertex*)trackD0->GetSecondaryVtx();

    //we check if the track is a desired candidate
    Int_t pdgCodeMother = -1;
    Bool_t isDesiredCandidate = kFALSE;
    Int_t motherType, histType;
    motherType = 0;
    Int_t mcLabelD0 = -1;

    if (fUseMCInfo)
    {
      mcLabelD0 = MatchCandidateToMonteCarlo(421, trackD0, mcTrackArray, BPlustoD0PiLabelMatrix, kTRUE);

      if (mcLabelD0 >= 0)
      {
        isDesiredCandidate = kTRUE;

        Int_t mcLabelFirstTrack = -1;
        mcLabelFirstTrack = trackFirstDaughter->GetLabel();

        if (mcLabelFirstTrack >= 0)
        {
          AliAODMCParticle *mcParticleFirstTrack = (AliAODMCParticle*)mcTrackArray->At(mcLabelFirstTrack);
          AliAODMCParticle *mcMotherParticle = (AliAODMCParticle*)mcTrackArray->At(mcLabelD0);

          if (mcParticleFirstTrack && mcMotherParticle)
          {
            pdgCodeMother = mcMotherParticle->GetPdgCode();

            Double_t vertex_distance = TMath::Sqrt((vertexMother->GetX() - mcParticleFirstTrack->Xv()) * (vertexMother->GetX() - mcParticleFirstTrack->Xv()) + (vertexMother->GetY() - mcParticleFirstTrack->Yv()) * (vertexMother->GetY() - mcParticleFirstTrack->Yv()) + (vertexMother->GetZ() - mcParticleFirstTrack->Zv()) * (vertexMother->GetZ() - mcParticleFirstTrack->Zv()));
            ((TH1F*)fMotherHistogramArrayExtra[motherType][4])->Fill(vertex_distance);

            // Double_t vertex_distance_new = TMath::Sqrt((vertexMotherNew->GetX() - mcParticleFirstTrack->Xv()) * (vertexMotherNew->GetX() - mcParticleFirstTrack->Xv()) + (vertexMotherNew->GetY() - mcParticleFirstTrack->Yv()) * (vertexMotherNew->GetY() - mcParticleFirstTrack->Yv()) + (vertexMotherNew->GetZ() - mcParticleFirstTrack->Zv()) * (vertexMotherNew->GetZ() - mcParticleFirstTrack->Zv()));
            // ((TH1F*)fMotherHistogramArrayExtra[motherType][5])->Fill(vertex_distance_new);

            Double_t momentum_resolution = TMath::Sqrt((trackD0->Px() - mcMotherParticle->Px()) * (trackD0->Px() - mcMotherParticle->Px()) + (trackD0->Py() - mcMotherParticle->Py()) * (trackD0->Py() - mcMotherParticle->Py()) + (trackD0->Pz() - mcMotherParticle->Pz()) * (trackD0->Pz() - mcMotherParticle->Pz()));
            ((TH1F*)fMotherHistogramArrayExtra[motherType][6])->Fill(momentum_resolution);
          }
        }
      }
    }

    // We fill the histograms
    histType = 0;
    FillD0Histograms(trackD0, primaryVertex, bz, motherType, histType);
    if (isDesiredCandidate && fUseMCInfo) {
      histType = 1;
      FillD0Histograms(trackD0, primaryVertex, bz, motherType, histType, pdgCodeMother);
    }

    // Here we apply cuts on the particle
    Bool_t cutMother = kFALSE;

    Bool_t bCutArray[29] = {0};
    Int_t cutReturnValue = fCuts->IsD0forD0ptbinSelected(trackD0, 0, aodEvent, bCutArray);
    if (cutReturnValue == -1) cutMother = kTRUE;
    if (cutReturnValue == 0) cutMother = kTRUE;
    for (Int_t k = 0; k < 29; ++k)
    {
      if (bCutArray[k] == kTRUE) {
        if (isDesiredCandidate) {
          ((TH1F*)fMotherHistogramArrayExtra[motherType][1])->Fill(k + 1);
        } else ((TH1F*)fMotherHistogramArrayExtra[motherType][0])->Fill(k + 1);
        cutMother = kTRUE;
      }
    }

    if (!fCuts->AreDaughtersSelected(trackD0, aodEvent)) cutMother = kTRUE;

    if (!isDesiredCandidate && fQuickSignalAnalysis == 1) cutMother = kTRUE;

    if (cutMother) {
      if (isDesiredCandidate) {
        ((TH1F*)fMotherHistogramArrayExtra[motherType][1])->Fill(0);
      } else ((TH1F*)fMotherHistogramArrayExtra[motherType][0])->Fill(0);
      continue;
    }

    // We fill the cut histograms
    histType = 2;
    FillD0Histograms(trackD0, primaryVertex, bz, motherType, histType);
    if (isDesiredCandidate && fUseMCInfo) {
      histType = 3;
      FillD0Histograms(trackD0, primaryVertex, bz, motherType, histType, pdgCodeMother);
    }

    //we save the location of the D0 candidate
    fD0Tracks->push_back(j);
  }

  delete vHF; vHF = nullptr;
  return;
}
//-------------------------------------------------------------------------------------
void AliAnalysisTaskSEBPlustoD0Pi::BPlusSelection(AliAODEvent* aodEvent, AliAODVertex *primaryVertex, Double_t bz, TClonesArray * mcTrackArray, TMatrix * BPlustoD0PiLabelMatrix, TClonesArray * D0TracksFromFriendFile, AliAODMCHeader * header) {

  TString fillthis = "";

  //we loop over all the D0 candidates
  for (Int_t j = 0; j < (Int_t)fD0Tracks->size(); j++)
  {

    //Save current Object count
    Int_t ObjectNumber = TProcessID::GetObjectCount();

    //we get the track of the D0
    AliAODRecoDecayHF2Prong * trackD0 = (AliAODRecoDecayHF2Prong*)(D0TracksFromFriendFile->At(fD0Tracks->at(j)));
    if (!trackD0) {std::cout << "found none" << std::endl; continue;}
    if (trackD0 == nullptr) {std::cout << "found nullptr" << std::endl; continue;}

    //we loop over all the BPlus pion candidates
    for (Int_t i = 0; i < (Int_t)fBPlusPionTracks->size(); i++)
    {

      //we get the track of the BPlus daughter
      AliAODTrack * trackBPlusPion = dynamic_cast<AliAODTrack*>(aodEvent->GetTrack(fBPlusPionTracks->at(i)));
      if (!trackBPlusPion) continue;

      //we check if the IDs of the tracks are different
      AliAODTrack* twoProngdaughter0 = (AliAODTrack*)trackD0->GetDaughter(0);
      AliAODTrack* twoProngdaughter1 = (AliAODTrack*)trackD0->GetDaughter(1);
      UShort_t idProng0 = twoProngdaughter0->GetID();
      UShort_t idProng1 = twoProngdaughter1->GetID();

      if (trackBPlusPion->GetID() == idProng0 || trackBPlusPion->GetID() == idProng1) continue;

      Int_t pdgD0 = 421;
      if (trackBPlusPion->Charge() == 1) pdgD0 = -421;

      //we check if the pions have the opposite charge
      //this only works if pid is turned on
      Bool_t bWrongSign = kFALSE;
      if (trackBPlusPion->Charge() == -1)
      {
        if ((fCuts->SelectPID(((AliAODTrack*)trackD0->GetDaughter(0)), 2)) && (fCuts->SelectPID(((AliAODTrack*)trackD0->GetDaughter(1)), 3))) bWrongSign = kFALSE;
        else if ((fCuts->SelectPID(((AliAODTrack*)trackD0->GetDaughter(0)), 3)) && (fCuts->SelectPID(((AliAODTrack*)trackD0->GetDaughter(1)), 2))) bWrongSign = kTRUE;
        else continue;
      } else if (trackBPlusPion->Charge() == 1) {
        pdgD0 = -421;
        if ((fCuts->SelectPID(((AliAODTrack*)trackD0->GetDaughter(0)), 3)) && (fCuts->SelectPID(((AliAODTrack*)trackD0->GetDaughter(1)), 2))) bWrongSign = kFALSE;
        else if ((fCuts->SelectPID(((AliAODTrack*)trackD0->GetDaughter(0)), 2)) && (fCuts->SelectPID(((AliAODTrack*)trackD0->GetDaughter(1)), 3))) bWrongSign = kTRUE;
        else continue;
      }

      //location BPlus pion rotation around PV
      for (Int_t iRot = 0; iRot < fNumberOfRotations + 1; ++iRot)
      {
        //we create a copy of the track that we will rotate
        AliAODTrack * trackBPlusPionRotated = new AliAODTrack(*trackBPlusPion);

        //for iRot == 0, we use the original unrotated track. For iRot > 0 we rotate the track and set the label to -1
        if (iRot != 0)
        {
          //should still check if track is already at PV
          Double_t dPhiRotated = trackBPlusPionRotated->Phi() + TMath::Pi() - (TMath::Pi() * fDegreePerRotation * fNumberOfRotations / (180.0 * 2.0)) + (TMath::Pi() * fDegreePerRotation * iRot / 180.0);
          trackBPlusPionRotated->SetPhi(dPhiRotated);
        }

        //we use the BPlus pion and D0 track to reconstruct the vertex for the BPlus
        AliExternalTrackParam firstTrack;
        firstTrack.CopyFromVTrack(trackBPlusPionRotated);
        AliExternalTrackParam secondTrack;
        secondTrack.CopyFromVTrack(trackD0);
        
        UInt_t prongs[2];
        prongs[1] = 421;
        prongs[0] = 211;

        UShort_t id[2];
        id[0] = firstTrack.GetID();
        id[1] = secondTrack.GetID();

        // we calculate the vertex of the mother candidate using two tracks 
        TObjArray daughterTracks;

        daughterTracks.Add(&firstTrack);
        daughterTracks.Add(&secondTrack);


        Double_t dispersion = 0;
        AliAODVertex *vertexMother = RecalculateVertex(primaryVertex, &daughterTracks, bz, dispersion, kTRUE, kTRUE, kFALSE);
        if (!vertexMother)
        {
          delete vertexMother; vertexMother = nullptr;
          delete trackBPlusPionRotated; trackBPlusPionRotated = nullptr;
          continue;
        }

        //use the new vertex to create the BPlus candidate
        Double_t xdummy = 0., ydummy = 0., dca;
        Double_t d0z0[2], covd0z0[3], d0[2], d0err[2];


        firstTrack.PropagateToDCA(vertexMother, bz, 100., d0z0, covd0z0);
        secondTrack.PropagateToDCA(vertexMother, bz, 100., d0z0, covd0z0);

        //we reconstruct the mother decay prong
        Double_t px[2], py[2], pz[2];
        px[0] = firstTrack.Px();
        py[0] = firstTrack.Py();
        pz[0] = firstTrack.Pz();
        px[1] = secondTrack.Px();
        py[1] = secondTrack.Py();
        pz[1] = secondTrack.Pz();

        firstTrack.PropagateToDCA(primaryVertex, bz, 100., d0z0, covd0z0);
        d0[0] = d0z0[0];
        d0err[0] = TMath::Sqrt(covd0z0[0]);
        secondTrack.PropagateToDCA(primaryVertex, bz, 100., d0z0, covd0z0);
        d0[1] = d0z0[0];
        d0err[1] = TMath::Sqrt(covd0z0[0]);


        dca = secondTrack.GetDCA(&firstTrack, bz, xdummy, ydummy);


        Short_t chargeMother = trackD0->Charge() + trackBPlusPionRotated->Charge();

        AliAODRecoDecayHF2Prong trackBPlus(vertexMother, px, py, pz, d0, d0err, dca);

        trackBPlus.SetCharge(chargeMother);
        trackBPlus.GetSecondaryVtx()->AddDaughter(trackBPlusPionRotated);
        trackBPlus.GetSecondaryVtx()->AddDaughter(trackD0);
        trackBPlus.SetPrimaryVtxRef((AliAODVertex*)aodEvent->GetPrimaryVertex());
        trackBPlus.SetProngIDs(2, id);

        // Fiducial cut
        if (TMath::Abs(trackBPlus.Y(521)) > 0.8) {
          delete vertexMother; vertexMother = nullptr;
          delete trackBPlusPionRotated; trackBPlusPionRotated = nullptr;
          continue;
        }

        // We check if the signal is injected, optionally we can reject injected signals
        Bool_t bIsInjected = kFALSE;
        if (fUseMCInfo) {
          bIsInjected = IsCandidateInjected(&trackBPlus, header, mcTrackArray);
          if (fRemoveInjected && bIsInjected) {
            delete vertexMother; vertexMother = nullptr;
            delete trackBPlusPionRotated; trackBPlusPionRotated = nullptr;
            continue;
          }
        }

        // We check if the BPlus candidate is a true signal in Monte Carlo
        Bool_t isDesiredCandidate = kFALSE;
        Int_t mcLabelBPlus = -1;
        fillthis = "";
        Int_t motherType, histType;
        motherType = 1;

        if (fUseMCInfo)
        {
          mcLabelBPlus = MatchCandidateToMonteCarlo(521, &trackBPlus, mcTrackArray, BPlustoD0PiLabelMatrix);
          if (mcLabelBPlus >= 0 && trackBPlusPionRotated->GetLabel() >= 0 && iRot == 0)
          {
            isDesiredCandidate = kTRUE;
          }
        }

        histType = 0;

        if (isDesiredCandidate)
        {
          AliAODMCParticle *mcTrackBPlusPion = (AliAODMCParticle*)mcTrackArray->At(trackBPlusPion->GetLabel());
          AliAODMCParticle *mcTrackBPlus = (AliAODMCParticle*)mcTrackArray->At(mcLabelBPlus);

          // Double_t vertex_distance = TMath::Sqrt((vertexMotherTemp->GetX() - mcTrackBPlusPion->Xv()) * (vertexMotherTemp->GetX() - mcTrackBPlusPion->Xv()) + (vertexMotherTemp->GetY() - mcTrackBPlusPion->Yv()) * (vertexMotherTemp->GetY() - mcTrackBPlusPion->Yv()) + (vertexMotherTemp->GetZ() - mcTrackBPlusPion->Zv()) * (vertexMotherTemp->GetZ() - mcTrackBPlusPion->Zv()));
          // ((TH1F*)fMotherHistogramArrayExtra[motherType][4])->Fill(vertex_distance);

          Double_t vertex_distance_new = TMath::Sqrt((vertexMother->GetX() - mcTrackBPlusPion->Xv()) * (vertexMother->GetX() - mcTrackBPlusPion->Xv()) + (vertexMother->GetY() - mcTrackBPlusPion->Yv()) * (vertexMother->GetY() - mcTrackBPlusPion->Yv()) + (vertexMother->GetZ() - mcTrackBPlusPion->Zv()) * (vertexMother->GetZ() - mcTrackBPlusPion->Zv()));
          ((TH1F*)fMotherHistogramArrayExtra[motherType][5])->Fill(vertex_distance_new);

          Double_t momentum_resolution = TMath::Sqrt((trackBPlus.Px() - mcTrackBPlus->Px()) * (trackBPlus.Px() - mcTrackBPlus->Px()) + (trackBPlus.Py() - mcTrackBPlus->Py()) * (trackBPlus.Py() - mcTrackBPlus->Py()) + (trackBPlus.Pz() - mcTrackBPlus->Pz()) * (trackBPlus.Pz() - mcTrackBPlus->Pz()));
          ((TH1F*)fMotherHistogramArrayExtra[motherType][6])->Fill(momentum_resolution);
        }

        if (!bWrongSign)
        {
          // We fill the histograms
          histType = 0;
          FillBPlusHistograms(&trackBPlus, primaryVertex, bz, motherType, histType);

          if (isDesiredCandidate)
          {
            histType = 1;
            FillBPlusHistograms(&trackBPlus, primaryVertex, bz, motherType, histType);
          }
        }


        // We apply cuts
        Bool_t cutMother = kFALSE;

        Bool_t bCutArray[75] = {0};
        Int_t numberOfCuts = 75;
        Int_t cutReturnValue = fCuts->IsSelected(&trackBPlus, 0, aodEvent, bCutArray);
        if (cutReturnValue == -1) cutMother = kTRUE;
        if (cutReturnValue == 0) cutMother = kTRUE;


        // We save information about the cuts
        TString histName = "";
        Double_t invariantMassMother = trackBPlus.InvMass(2, prongs);
        Double_t pdgMassMother = TDatabasePDG::Instance()->GetParticle(521)->Mass();
        Double_t massWindow = fHistMassWindow; //GeV/c^2

        for (Int_t n = 0; n < 75; ++n)
        {
          if (bCutArray[n] == kTRUE) {
            if (isDesiredCandidate) {
              ((TH1F*)fMotherHistogramArrayExtra[motherType][1])->Fill(n + 1);
            } else ((TH1F*)fMotherHistogramArrayExtra[motherType][0])->Fill(n + 1);
            cutMother = kTRUE;
          }
        }

        if (TMath::Abs(invariantMassMother - pdgMassMother) < massWindow) {
          for (Int_t l = 0; l < numberOfCuts; ++l) //total
          {
            if (bCutArray[l] == kFALSE) continue;
            for (Int_t j = 0; j < numberOfCuts; ++j)
            {
              if (bCutArray[j] == kFALSE) continue;
              if (isDesiredCandidate == kFALSE) ((TH2F*)(fResultsHistogramArray2D[2][0]))->Fill(l, j);
              if (isDesiredCandidate == kTRUE) ((TH2F*)(fResultsHistogramArray2D[2][1]))->Fill(l, j);
            }
          }

          for (Int_t l = 0; l < numberOfCuts; ++l) //unique
          {
            if (bCutArray[l] == kFALSE) continue;
            Bool_t bFill = kTRUE;
            for (Int_t j = 0; j < numberOfCuts; ++j)
            {
              if (l == j) continue;
              if (bCutArray[j] == kTRUE)
              {
                bFill = kFALSE;
                break;
              }

            }
            if (bFill == kTRUE)
            {
              if (isDesiredCandidate == kFALSE) ((TH1F*)(fResultsHistogramArray[13][0]))->Fill(l);
              if (isDesiredCandidate == kTRUE) ((TH1F*)(fResultsHistogramArray[13][1]))->Fill(l);
            }
          }
        }


        if (!isDesiredCandidate && fQuickSignalAnalysis == 1) cutMother = kTRUE;

        if (cutMother)
        {
          if (isDesiredCandidate)
          {
            ((TH1F*)fMotherHistogramArrayExtra[motherType][1])->Fill(0);
          } else ((TH1F*)fMotherHistogramArrayExtra[motherType][0])->Fill(0);
          delete vertexMother; vertexMother = nullptr;
          delete trackBPlusPionRotated; trackBPlusPionRotated = nullptr;
          continue;
        }

        if (!bWrongSign)
        {
          histType = 2;
          FillBPlusHistograms(&trackBPlus, primaryVertex, bz, motherType, histType);
          if (fUseMCInfo && isDesiredCandidate)
          {
            //fill mc histograms
            histType = 3;
            FillBPlusHistograms(&trackBPlus, primaryVertex, bz, motherType, histType);
          }
        }

        if (TMath::Abs(invariantMassMother - pdgMassMother) < massWindow)
        {
          if (!bWrongSign)
          {
            FillFinalTrackHistograms(&trackBPlus, primaryVertex, bz, isDesiredCandidate, mcTrackArray);
            if (!isDesiredCandidate)
            {
              motherType = 0; histType = 4; FillD0Histograms(trackD0, primaryVertex, bz, motherType, histType, pdgD0);
              motherType = 1; histType = 4; FillBPlusHistograms(&trackBPlus, primaryVertex, bz, motherType, histType);
            }
            if (isDesiredCandidate)
            {
              motherType = 0; histType = 5; FillD0Histograms(trackD0, primaryVertex, bz, motherType, histType, pdgD0);
              motherType = 1; histType = 5; FillBPlusHistograms(&trackBPlus, primaryVertex, bz, motherType, histType);
            }
          }
        }


        // Here we fill the histograms per pt bin and apply the same sign method
        TString ptBinMother = "";
        Int_t ptBin = fCuts->PtBin(trackBPlus.Pt());
        ptBinMother += "_ptbin_"; ptBinMother += fPtBinLimits[ptBin]; ptBinMother += "_to_"; ptBinMother += fPtBinLimits[ptBin + 1];
        histType = 6 + 2 * ptBin;

        Int_t d0PtBin = fCuts->PtBinD0forD0ptbin(trackD0->Pt());
        Int_t histTypeD0 = 2 * d0PtBin;


        if (TMath::Abs(invariantMassMother - pdgMassMother) < massWindow)
        {
          if (!bWrongSign && histType > 5)
          {
            if (!isDesiredCandidate)
            {
              motherType = 0; FillD0Histograms(trackD0, primaryVertex, bz, motherType, histType, pdgD0);
              motherType = 1; FillBPlusHistograms(&trackBPlus, primaryVertex, bz, motherType, histType);
              motherType = 2; FillD0Histograms(trackD0, primaryVertex, bz, motherType, histTypeD0, pdgD0);
            }

            if (isDesiredCandidate)
            {
              histType += 1;
              motherType = 0; FillD0Histograms(trackD0, primaryVertex, bz, motherType, histType, pdgD0);
              motherType = 1; FillBPlusHistograms(&trackBPlus, primaryVertex, bz, motherType, histType);
              motherType = 2; FillD0Histograms(trackD0, primaryVertex, bz, motherType, histTypeD0 + 1, pdgD0);
            }
          }
        }

        // after the (loose) cut on candidates we can get data for cut optimization
        if (!bWrongSign && fPerformCutOptimization)
        {
          Double_t sigmaWindowForCutOptimization = fCuts->GetSigmaForCutOptimization(ptBin);
          Int_t nSigmaBins = fCuts->GetNumberOfSigmaBinsForCutOptimization();
          Double_t massDifference = TMath::Abs(invariantMassMother - pdgMassMother);

          if (massDifference < nSigmaBins * sigmaWindowForCutOptimization)
          {
            Int_t nSigmaBin = 0;
            if (invariantMassMother < pdgMassMother)
            {
              for (Int_t iSigma = 0; iSigma < nSigmaBins; ++iSigma)
              {
                if ((massDifference > iSigma * sigmaWindowForCutOptimization) && (massDifference < (iSigma + 1) * sigmaWindowForCutOptimization)) {nSigmaBin = -(iSigma + 1); break;}
              }
            }
            if (invariantMassMother > pdgMassMother)
            {
              for (Int_t iSigma = 0; iSigma < nSigmaBins; ++iSigma)
              {
                if ((massDifference > iSigma * sigmaWindowForCutOptimization) && (massDifference < (iSigma + 1) * sigmaWindowForCutOptimization)) {nSigmaBin = iSigma; break;}
              }
            }

            Int_t nStartVariable = 0;
            Int_t nStartFillNumber = 0;
            Int_t nVariables = fCuts->GetnVariablesForCutOptimization();
            Int_t nCuts = fCuts->GetnCutsForOptimization();
            CutOptimizationVariableValues(&trackBPlus, aodEvent);
            CutOptimizationLoop(nStartVariable, nVariables, nCuts, ptBin, nStartFillNumber, isDesiredCandidate, nSigmaBin);
          }
        }


        Double_t invmassDelta = DeltaInvMassBPlusKpipi(&trackBPlus);

        // Fill invariant mass histograms
        for (Int_t m = 0; m < 7; ++m) 
        {
          Int_t fillfactor = 1;
          if (bWrongSign && !(m == 1 || m == 2)) continue;
          if (!bWrongSign && m == 1) continue;
          if (bWrongSign && m == 2) fillfactor = -1;
          if (iRot > 0 && !(m == 3 || m == 6)) continue;
          if (iRot == 0 && (m == 3 || m == 6)) continue;
          if (m > 3 && bIsInjected) continue;
          if (m == 4 && isDesiredCandidate) continue;
          if (m == 5 && !isDesiredCandidate) continue;
          if (m == 6 && isDesiredCandidate) continue;

          ((TH1F*)fResultsHistogramArray[4 + m][0])->Fill(invariantMassMother,fillfactor);
          ((TH1F*)fResultsHistogramArray[4 + m][0 + fnPtBins + 3])->Fill(invmassDelta,fillfactor);

          if (ptBin > 0)
          {
            ((TH1F*)fResultsHistogramArray[4 + m][1])->Fill(invariantMassMother,fillfactor);
            ((TH1F*)fResultsHistogramArray[4 + m][1 + fnPtBins + 3])->Fill(invmassDelta,fillfactor);
          }
          if (ptBin > 1)
          {
            ((TH1F*)fResultsHistogramArray[4 + m][2])->Fill(invariantMassMother,fillfactor);
            ((TH1F*)fResultsHistogramArray[4 + m][2 + fnPtBins + 3])->Fill(invmassDelta,fillfactor);
          }

          ((TH1F*)fResultsHistogramArray[4 + m][ptBin + 3])->Fill(invariantMassMother,fillfactor);
          ((TH1F*)fResultsHistogramArray[4 + m][ptBin + 3 + fnPtBins + 3])->Fill(invmassDelta,fillfactor);
        }

        delete vertexMother; vertexMother = NULL;
        delete trackBPlusPionRotated; trackBPlusPionRotated = nullptr;
      }
    }

    //Restore Object count
    //To save space in the table keeping track of all referenced objects,
    //we reset the object count to what it was at the beginning of the loop.
    TProcessID::SetObjectCount(ObjectNumber);
  }
  return;
}
//-------------------------------------------------------------------------------------
void AliAnalysisTaskSEBPlustoD0Pi::FillFinalTrackHistograms(AliAODRecoDecayHF2Prong * selectedBPlus, AliAODVertex *primaryVertex, Double_t bz, Bool_t isDesiredCandidate, TClonesArray * mcTrackArray) {

  //In this function we fill histograms with the properties of all the daughters of our selected signal candidate

  AliAODTrack* selectedBPlusPion = (AliAODTrack*)selectedBPlus->GetDaughter(0);
  AliAODRecoDecayHF2Prong* selectedD0 = (AliAODRecoDecayHF2Prong*)selectedBPlus->GetDaughter(1);

  AliAODTrack* selectedD0Pion = 0x0;
  AliAODTrack* selectedD0Kaon = 0x0;

  if (selectedBPlus->Charge() == 1) selectedD0Pion = (AliAODTrack*)selectedD0->GetDaughter(0);
  if (selectedBPlus->Charge() == -1) selectedD0Pion = (AliAODTrack*)selectedD0->GetDaughter(1);

  if (selectedBPlus->Charge() == 1) selectedD0Kaon = (AliAODTrack*)selectedD0->GetDaughter(1);
  if (selectedBPlus->Charge() == -1) selectedD0Kaon = (AliAODTrack*)selectedD0->GetDaughter(0);

  Double_t d0BPluspion = TMath::Abs(selectedBPlus->Getd0Prong(0));
  Double_t d0D0pion = 0;
  Double_t d0D0kaon = 0;

  if (selectedBPlus->Charge() == 1) d0D0pion = selectedD0->Getd0Prong(0);
  if (selectedBPlus->Charge() == -1) d0D0pion = selectedD0->Getd0Prong(1);

  if (selectedBPlus->Charge() == 1) d0D0kaon = selectedD0->Getd0Prong(1);
  if (selectedBPlus->Charge() == -1) d0D0kaon = selectedD0->Getd0Prong(0);

  Double_t pt_track = 0;
  Double_t momentum_track = 0;
  Int_t numberOfITS = 0;
  Int_t numberOfTPC = 0;
  Int_t daughterType, histType;
  Int_t totalNumberOfITS = 0;
  Int_t totalNumberOfTPC = 0;
  Double_t nSigmaTPC = 0;
  Double_t nSigmaTOF = 0;
  Double_t nSigmaTPCtotal = 0;
  Double_t nSigmaTOFtotal = 0;
  Int_t pionPIDnumber = 2;
  Int_t kaonPIDnumber = 3;
  Int_t TPCok = 0;
  Int_t TOFok = 0;

  AliAODPidHF* trackPIDHF = (AliAODPidHF*)fCuts->GetPidHF();

  //fill the D0 pion info
  pt_track = selectedD0Pion->Pt();
  momentum_track = selectedD0Pion->P();
  numberOfITS = selectedD0Pion->GetITSNcls();
  numberOfTPC = selectedD0Pion->GetTPCNcls();
  totalNumberOfITS += numberOfITS;
  totalNumberOfTPC += numberOfTPC;
  if(trackPIDHF) TPCok = trackPIDHF->GetnSigmaTPC(selectedD0Pion, pionPIDnumber, nSigmaTPC);
  if(trackPIDHF) TOFok = trackPIDHF->GetnSigmaTOF(selectedD0Pion, pionPIDnumber, nSigmaTOF);
  if (TPCok != -1) nSigmaTPCtotal += nSigmaTPC * nSigmaTPC;
  if (TOFok != -1) nSigmaTOFtotal += nSigmaTOF * nSigmaTOF;

  Double_t ptBPlus = selectedBPlus->Pt();

  daughterType = 0;
  histType = 4;
  if (!isDesiredCandidate)
  {
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][0])->Fill(pt_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][1])->Fill(momentum_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][2])->Fill(numberOfITS);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][3])->Fill(numberOfTPC);

    for (Int_t j = 0; j < 10; ++j)
    {
      if (selectedD0Pion->HasPointOnITSLayer(j)) ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(j);

    }

    if (TPCok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(nSigmaTPC);
    if (TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][5])->Fill(nSigmaTOF);
    if (TPCok != -1 && TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][6])->Fill(TMath::Sqrt(nSigmaTPC * nSigmaTPC + nSigmaTOF * nSigmaTOF));
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][8])->Fill(d0D0pion);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][9])->Fill(selectedD0Pion->Eta());
    ((TH2F*)fDaughterHistogramArray2D[daughterType][4])->Fill(ptBPlus, pt_track);
  }

  if (isDesiredCandidate)
  {
    histType = 5;
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][0])->Fill(pt_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][1])->Fill(momentum_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][2])->Fill(numberOfITS);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][3])->Fill(numberOfTPC);

    for (Int_t j = 0; j < 10; ++j)
    {
      if (selectedD0Pion->HasPointOnITSLayer(j)) ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(j);

    }

    if (TPCok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(nSigmaTPC);
    if (TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][5])->Fill(nSigmaTOF);
    if (TPCok != -1 && TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][6])->Fill(TMath::Sqrt(nSigmaTPC * nSigmaTPC + nSigmaTOF * nSigmaTOF));
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][8])->Fill(d0D0pion);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][9])->Fill(selectedD0Pion->Eta());
    ((TH2F*)fDaughterHistogramArray2D[daughterType][5])->Fill(ptBPlus, pt_track);
  }

  //we save the pdgcode of the used particle and its mother to check PID efficiency
  if (fUseMCInfo)
  {
    Float_t pdgCodeParticle = -1;
    Float_t pdgCodeParticleMother = -1;
    Int_t mcLabelParticle = -1;
    Int_t mcLabelParticleMother = -1;
    mcLabelParticle = selectedD0Pion->GetLabel();

    if (mcLabelParticle >= 0) {

      AliAODMCParticle *mcTrackParticle = (AliAODMCParticle*)mcTrackArray->At(mcLabelParticle);
      pdgCodeParticle = TMath::Abs(mcTrackParticle->GetPdgCode());
      ((TH1F*)fDaughterHistogramArrayExtra[0][2])->Fill(pdgCodeParticle);
      mcLabelParticleMother = mcTrackParticle->GetMother();

      if (mcLabelParticleMother >= 0) {
        AliAODMCParticle *mcTrackParticleMother = (AliAODMCParticle*)mcTrackArray->At(mcLabelParticleMother);
        pdgCodeParticleMother = TMath::Abs(mcTrackParticleMother->GetPdgCode());
        ((TH1F*)fDaughterHistogramArrayExtra[0][3])->Fill(pdgCodeParticleMother);
      }
    }
  }

  //fill the D0 kaon info
  pt_track = selectedD0Kaon->Pt();
  momentum_track = selectedD0Kaon->P();
  numberOfITS = selectedD0Kaon->GetITSNcls();
  numberOfTPC = selectedD0Kaon->GetTPCNcls();
  totalNumberOfITS += numberOfITS;
  totalNumberOfTPC += numberOfTPC;
  if(trackPIDHF) TPCok = trackPIDHF->GetnSigmaTPC(selectedD0Kaon, kaonPIDnumber, nSigmaTPC);
  if(trackPIDHF) TOFok = trackPIDHF->GetnSigmaTOF(selectedD0Kaon, kaonPIDnumber, nSigmaTOF);
  if (TPCok != -1) nSigmaTPCtotal += nSigmaTPC * nSigmaTPC;
  if (TOFok != -1) nSigmaTOFtotal += nSigmaTOF * nSigmaTOF;

  daughterType = 1;
  histType = 4;
  if (!isDesiredCandidate)
  {
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][0])->Fill(pt_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][1])->Fill(momentum_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][2])->Fill(numberOfITS);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][3])->Fill(numberOfTPC);

    for (Int_t j = 0; j < 10; ++j)
    {
      if (selectedD0Kaon->HasPointOnITSLayer(j)) ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(j);

    }

    if (TPCok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(nSigmaTPC);
    if (TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][5])->Fill(nSigmaTOF);
    if (TPCok != -1 && TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][6])->Fill(TMath::Sqrt(nSigmaTPC * nSigmaTPC + nSigmaTOF * nSigmaTOF));
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][8])->Fill(d0D0kaon);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][9])->Fill(selectedD0Kaon->Eta());
    ((TH2F*)fDaughterHistogramArray2D[daughterType][4])->Fill(ptBPlus, pt_track);
  }

  if (isDesiredCandidate)
  {
    histType = 5;
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][0])->Fill(pt_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][1])->Fill(momentum_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][2])->Fill(numberOfITS);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][3])->Fill(numberOfTPC);

    for (Int_t j = 0; j < 10; ++j)
    {
      if (selectedD0Kaon->HasPointOnITSLayer(j)) ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(j);

    }

    if (TPCok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(nSigmaTPC);
    if (TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][5])->Fill(nSigmaTOF);
    if (TPCok != -1 && TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][6])->Fill(TMath::Sqrt(nSigmaTPC * nSigmaTPC + nSigmaTOF * nSigmaTOF));
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][8])->Fill(d0D0kaon);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][9])->Fill(selectedD0Kaon->Eta());
    ((TH2F*)fDaughterHistogramArray2D[daughterType][5])->Fill(ptBPlus, pt_track);
  }

  //we save the pdgcode of the used particle and its mother to check PID efficiency
  if (fUseMCInfo)
  {
    Float_t pdgCodeParticle = -1;
    Float_t pdgCodeParticleMother = -1;
    Int_t mcLabelParticle = -1;
    Int_t mcLabelParticleMother = -1;
    mcLabelParticle = selectedD0Kaon->GetLabel();

    if (mcLabelParticle >= 0) {

      AliAODMCParticle *mcTrackParticle = (AliAODMCParticle*)mcTrackArray->At(mcLabelParticle);
      pdgCodeParticle = TMath::Abs(mcTrackParticle->GetPdgCode());
      ((TH1F*)fDaughterHistogramArrayExtra[1][2])->Fill(pdgCodeParticle);
      mcLabelParticleMother = mcTrackParticle->GetMother();

      if (mcLabelParticleMother >= 0) {
        AliAODMCParticle *mcTrackParticleMother = (AliAODMCParticle*)mcTrackArray->At(mcLabelParticleMother);
        pdgCodeParticleMother = TMath::Abs(mcTrackParticleMother->GetPdgCode());
        ((TH1F*)fDaughterHistogramArrayExtra[1][3])->Fill(pdgCodeParticleMother);
      }
    }
  }

  //fill the BPlus pion info
  pt_track = selectedBPlusPion->Pt();
  momentum_track = selectedBPlusPion->P();
  numberOfITS = selectedBPlusPion->GetITSNcls();
  numberOfTPC = selectedBPlusPion->GetTPCNcls();
  totalNumberOfITS += numberOfITS;
  totalNumberOfTPC += numberOfTPC;
  if(trackPIDHF) TPCok = trackPIDHF->GetnSigmaTPC(selectedBPlusPion, pionPIDnumber, nSigmaTPC);
  if(trackPIDHF) TOFok = trackPIDHF->GetnSigmaTOF(selectedBPlusPion, pionPIDnumber, nSigmaTOF);
  if (TPCok != -1) nSigmaTPCtotal += nSigmaTPC * nSigmaTPC;
  if (TOFok != -1) nSigmaTOFtotal += nSigmaTOF * nSigmaTOF;

  daughterType = 2;
  histType = 4;
  if (!isDesiredCandidate)
  {
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][0])->Fill(pt_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][1])->Fill(momentum_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][2])->Fill(numberOfITS);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][3])->Fill(numberOfTPC);

    for (Int_t j = 0; j < 10; ++j)
    {
      if (selectedBPlusPion->HasPointOnITSLayer(j)) ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(j);

    }

    if (TPCok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(nSigmaTPC);
    if (TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][5])->Fill(nSigmaTOF);
    if (TPCok != -1 && TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][6])->Fill(TMath::Sqrt(nSigmaTPC * nSigmaTPC + nSigmaTOF * nSigmaTOF));
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][8])->Fill(d0BPluspion);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][9])->Fill(selectedBPlusPion->Eta());
    ((TH2F*)fDaughterHistogramArray2D[daughterType][4])->Fill(ptBPlus, pt_track);
  }

  if (isDesiredCandidate)
  {
    histType = 5;
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][0])->Fill(pt_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][1])->Fill(momentum_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][2])->Fill(numberOfITS);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][3])->Fill(numberOfTPC);

    for (Int_t j = 0; j < 10; ++j)
    {
      if (selectedBPlusPion->HasPointOnITSLayer(j)) ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(j);

    }

    if (TPCok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(nSigmaTPC);
    if (TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][5])->Fill(nSigmaTOF);
    if (TPCok != -1 && TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][6])->Fill(TMath::Sqrt(nSigmaTPC * nSigmaTPC + nSigmaTOF * nSigmaTOF));
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][8])->Fill(d0BPluspion);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][9])->Fill(selectedBPlusPion->Eta());
    ((TH2F*)fDaughterHistogramArray2D[daughterType][5])->Fill(ptBPlus, pt_track);
  }

  //we save the pdgcode of the used particle and its mother to check PID efficiency
  if (fUseMCInfo)
  {
    Float_t pdgCodeParticle = -1;
    Float_t pdgCodeParticleMother = -1;
    Int_t mcLabelParticle = -1;
    Int_t mcLabelParticleMother = -1;
    mcLabelParticle = selectedBPlusPion->GetLabel();

    if (mcLabelParticle >= 0) {

      AliAODMCParticle *mcTrackParticle = (AliAODMCParticle*)mcTrackArray->At(mcLabelParticle);
      pdgCodeParticle = TMath::Abs(mcTrackParticle->GetPdgCode());
      ((TH1F*)fDaughterHistogramArrayExtra[daughterType][2])->Fill(pdgCodeParticle);
      mcLabelParticleMother = mcTrackParticle->GetMother();

      if (mcLabelParticleMother >= 0) {
        AliAODMCParticle *mcTrackParticleMother = (AliAODMCParticle*)mcTrackArray->At(mcLabelParticleMother);
        pdgCodeParticleMother = TMath::Abs(mcTrackParticleMother->GetPdgCode());
        ((TH1F*)fDaughterHistogramArrayExtra[daughterType][3])->Fill(pdgCodeParticleMother);
      }
    }
  }

  if (!isDesiredCandidate)
  {
    ((TH1F*)(fResultsHistogramArray[11][0]))->Fill(totalNumberOfITS);
    ((TH1F*)(fResultsHistogramArray[11][2]))->Fill(totalNumberOfTPC);
    ((TH1F*)(fResultsHistogramArray[11][4]))->Fill(TMath::Sqrt(nSigmaTPCtotal + nSigmaTOFtotal));
  }
  if (isDesiredCandidate)
  {
    ((TH1F*)(fResultsHistogramArray[11][1]))->Fill(totalNumberOfITS);
    ((TH1F*)(fResultsHistogramArray[11][3]))->Fill(totalNumberOfTPC);
    ((TH1F*)(fResultsHistogramArray[11][5]))->Fill(TMath::Sqrt(nSigmaTPCtotal + nSigmaTOFtotal));
  }

  // we look at the invariant mass combinations of all daughter tracks
  AliExternalTrackParam trackBPlusPion;
  trackBPlusPion.CopyFromVTrack(selectedBPlusPion);
  AliExternalTrackParam trackD0Pion;
  trackD0Pion.CopyFromVTrack(selectedD0Pion);
  AliExternalTrackParam trackD0Kaon;
  trackD0Kaon.CopyFromVTrack(selectedD0Kaon);

  UInt_t prongs2[2] = {0};
  UInt_t prongs3[3] = {0};

  prongs2[0] = 211; prongs2[1] = 211;
  TwoTrackCombinationInfo(&trackD0Pion, &trackBPlusPion, primaryVertex, bz, isDesiredCandidate, 0, prongs2);
  prongs2[0] = 321; prongs2[1] = 211;
  TwoTrackCombinationInfo(&trackD0Kaon, &trackBPlusPion, primaryVertex, bz, isDesiredCandidate, 2, prongs2);
  prongs3[0] = 211; prongs3[1] = 321; prongs3[2] = 211;
  ThreeTrackCombinationInfo(&trackD0Pion, &trackD0Kaon, &trackBPlusPion, primaryVertex, bz, isDesiredCandidate, 4, prongs3);

  return;
}
//-------------------------------------------------------------------------------------
Double_t AliAnalysisTaskSEBPlustoD0Pi::DeltaInvMassBPlusKpipi(AliAODRecoDecayHF2Prong * BPlus) const
{
  ///
  /// 3 prong invariant mass of the D0 daughters, the soft pion, and the BPlus pion
  ///

  AliAODRecoDecayHF2Prong * D0 = (AliAODRecoDecayHF2Prong*)BPlus->GetDaughter(1);

  Double_t e[3] = {0};
  if (BPlus->Charge() == -1)
  {
    e[0] = D0->EProng(0, 211);
    e[1] = D0->EProng(1, 321);
  } else if (BPlus->Charge() == 1) {
    e[0] = D0->EProng(0, 321);
    e[1] = D0->EProng(1, 211);
  }
  e[2] = BPlus->EProng(0, 211);

  Double_t esum = e[0] + e[1] + e[2];
  Double_t invMassBPlus = TMath::Sqrt(esum * esum - BPlus->P2());

  Double_t invMassD0 = -1;

  if (BPlus->Charge() == -1) {invMassD0 = D0->InvMassD0();}
  else {invMassD0 = D0->InvMassD0bar();}
  if (invMassD0 == -1) {std::cout << "wrong invmass delta D0 BPlus" << std::endl;}

  return invMassBPlus - invMassD0;
}
//-------------------------------------------------------------------------------------
void AliAnalysisTaskSEBPlustoD0Pi::FillD0Histograms(AliAODRecoDecayHF2Prong * selectedMother, AliAODVertex *primaryVertex, Double_t bz, Int_t motherType, Int_t histType, Int_t pdgCodeMother) {

  if (histType < 0) return;

  //In this function we fill the histograms of the reconstructed mothers
  Double_t ptMother = 0.0;
  Double_t momentumMother = 0.0;
  Double_t etaMother = 0.0;
  Double_t phiMother = 0.0;
  Double_t d0Mother = 0.0;
  Double_t d0firstTrack = 0.0;
  Double_t d0secondTrack = 0.0;
  Double_t pointingAngle = 0.0;
  Double_t impactProduct = 0.0;
  Double_t impactProductXY = 0.0;
  Double_t invariantMassMother = 0.0;
  Double_t invmassDelta = 0.0;
  Double_t dcaMother = 0.0;
  AliAODVertex * vertexMother = 0x0;
  Double_t vertexDistance = 0.0;
  Double_t decayTime = 0.0;
  Double_t angleMotherFirstDaughter = 0.0;
  Double_t angleMotherSecondDaughter = 0.0;
  Double_t ptFirstDaughter = 0.0;
  Double_t ptSecondDaughter = 0.0;
  UInt_t prongs[2], prongs2[2];
  Double_t angleBetweenBothDaughters = 0;
  Double_t cosThetaStar = 0;
  Double_t normDecayLength = 0;
  Double_t pdgMassMother = TDatabasePDG::Instance()->GetParticle(421)->Mass();


  prongs[0] = 211; prongs[1] = 321;
  prongs2[0] = 321; prongs2[1] = 211;
  AliAODTrack * firstDaughter = (AliAODTrack*)selectedMother->GetDaughter(0);
  AliAODTrack * secondDaughter = (AliAODTrack*)selectedMother->GetDaughter(1);
  vertexMother = selectedMother->GetSecondaryVtx();
  ptFirstDaughter = firstDaughter->Pt();
  ptSecondDaughter = secondDaughter->Pt();

  //Topomatic
  Double_t dd0pr1 = 0.;
  Double_t dd0pr2 = 0.;
  Double_t dd0max = 0.;
  Double_t dd0min = 0.;
  for (Int_t ipr = 0; ipr < 2; ipr++)
  {
    Double_t diffIP, errdiffIP;
    selectedMother->Getd0MeasMinusExpProng(ipr, bz, diffIP, errdiffIP);
    Double_t normdd0 = 0.;
    if (errdiffIP > 0.) normdd0 = diffIP / errdiffIP;
    if (ipr == 0) dd0pr1 = normdd0;
    if (ipr == 1) dd0pr2 = normdd0;

  }
  if (TMath::Abs(dd0pr1) > TMath::Abs(dd0pr2)) {dd0max = dd0pr1; dd0min = dd0pr2;}
  else {dd0max = dd0pr2; dd0min = dd0pr1;}

  AliExternalTrackParam motherTrack;
  motherTrack.CopyFromVTrack(selectedMother);
  Double_t d0z0[2], covd0z0[3], d0[2];
  motherTrack.PropagateToDCA(primaryVertex, bz, 100., d0z0, covd0z0);
  d0[0] = d0z0[0];

  ptMother = selectedMother->Pt();
  momentumMother = selectedMother->P();
  etaMother = selectedMother->Eta();
  phiMother = selectedMother->Phi();

  d0Mother = TMath::Abs(d0[0]);
  d0firstTrack = TMath::Abs(selectedMother->Getd0Prong(0));
  d0secondTrack = TMath::Abs(selectedMother->Getd0Prong(1));
  pointingAngle = selectedMother->CosPointingAngle();
  impactProduct = selectedMother->Prodd0d0();
  impactProductXY = TMath::Abs(selectedMother->ImpParXY());
  invariantMassMother = selectedMother->InvMass(2, prongs);
  if (pdgCodeMother == -421) invariantMassMother = selectedMother->InvMass(2, prongs2);

  dcaMother = selectedMother->GetDCA();
  vertexDistance = vertexMother->DistanceToVertex(primaryVertex);
  angleMotherFirstDaughter = (selectedMother->Px() * firstDaughter->Px() + selectedMother->Py() * firstDaughter->Py() + selectedMother->Pz() * firstDaughter->Pz()) / (selectedMother->P() * firstDaughter->P());
  angleMotherSecondDaughter = (selectedMother->Px() * secondDaughter->Px() + selectedMother->Py() * secondDaughter->Py() + selectedMother->Pz() * secondDaughter->Pz()) / (selectedMother->P() * secondDaughter->P());
  cosThetaStar = selectedMother->CosThetaStar(0, 421, 211, 321);
  angleBetweenBothDaughters  = (firstDaughter->Px() * secondDaughter->Px() + firstDaughter->Py() * secondDaughter->Py() + firstDaughter->Pz() * secondDaughter->Pz()) / (firstDaughter->P() * secondDaughter->P());
  normDecayLength = selectedMother->NormalizedDecayLength();

  Double_t pseudoProperDecayLength = ((vertexMother->GetX() - primaryVertex->GetX()) * selectedMother->Px() / TMath::Abs(selectedMother->Pt())) + ((vertexMother->GetY() - primaryVertex->GetY()) * selectedMother->Py() / TMath::Abs(selectedMother->Pt()));
  Double_t pseudoProperDecayTime = pseudoProperDecayLength * pdgMassMother / ptMother;
  decayTime = vertexDistance / (299792458 * TMath::Sqrt(1 / ((pdgMassMother * pdgMassMother / (momentumMother * momentumMother)) + 1)));

  Double_t phi = selectedMother->Phi();
  Double_t theta = selectedMother->Theta();
  Double_t covMatrix[21];
  selectedMother->GetCovarianceXYZPxPyPz(covMatrix);

  Double_t cp = TMath::Cos(phi);
  Double_t sp = TMath::Sin(phi);
  Double_t ct = TMath::Cos(theta);
  Double_t st = TMath::Sin(theta);

  Double_t errorMomentum = covMatrix[9] * cp * cp * ct * ct // GetCovPxPx
                           + covMatrix[13] * 2.*cp * sp * ct * ct // GetCovPxPy
                           + covMatrix[18] * 2.*cp * ct * st // GetCovPxPz
                           + covMatrix[14] * sp * sp * ct * ct // GetCovPyPy
                           + covMatrix[19] * 2.*sp * ct * st // GetCovPyPz
                           + covMatrix[20] * st * st; // GetCovPzPz
  Double_t normalizedDecayTime = selectedMother->NormalizedDecayLength() / (299792458 * TMath::Sqrt(1 / ((pdgMassMother * pdgMassMother * errorMomentum * errorMomentum / (momentumMother * momentumMother)) + 1)));

  Double_t eKaon = selectedMother->EProng(1, 321);
  Double_t invMassKaon = TMath::Sqrt(eKaon * eKaon - secondDaughter->P() * secondDaughter->P());
  if (pdgCodeMother == -421) eKaon = selectedMother->EProng(0, 321);
  if (pdgCodeMother == -421) invMassKaon = TMath::Sqrt(eKaon * eKaon - firstDaughter->P() * firstDaughter->P());
  Double_t invMassD0 = selectedMother->InvMassD0();
  invmassDelta = invMassD0 - invMassKaon;

  Double_t vertexMotherX = vertexMother->GetX();
  Double_t vertexMotherY = vertexMother->GetY();
  Double_t vertexMotherZ = vertexMother->GetZ();

  Double_t cosPointingAngleXY = selectedMother->CosPointingAngleXY();
  Double_t distanceXYToVertex = vertexMother->DistanceXYToVertex(primaryVertex);
  Double_t normalizedDecayLengthXY = selectedMother->NormalizedDecayLengthXY();

  ((TH1F*)fMotherHistogramArray[motherType][histType][0])->Fill(ptMother);
  ((TH1F*)fMotherHistogramArray[motherType][histType][1])->Fill(ptFirstDaughter);
  ((TH1F*)fMotherHistogramArray[motherType][histType][2])->Fill(ptSecondDaughter);
  ((TH1F*)fMotherHistogramArray[motherType][histType][3])->Fill(etaMother);
  ((TH1F*)fMotherHistogramArray[motherType][histType][4])->Fill(phiMother);
  ((TH1F*)fMotherHistogramArray[motherType][histType][5])->Fill(d0Mother);
  ((TH1F*)fMotherHistogramArray[motherType][histType][6])->Fill(d0firstTrack);
  ((TH1F*)fMotherHistogramArray[motherType][histType][7])->Fill(d0secondTrack);
  ((TH1F*)fMotherHistogramArray[motherType][histType][8])->Fill(pointingAngle);
  ((TH1F*)fMotherHistogramArray[motherType][histType][9])->Fill(impactProduct);
  ((TH1F*)fMotherHistogramArray[motherType][histType][10])->Fill(impactProductXY);
  ((TH1F*)fMotherHistogramArray[motherType][histType][11])->Fill(invariantMassMother);
  ((TH1F*)fMotherHistogramArray[motherType][histType][12])->Fill(invmassDelta);
  ((TH1F*)fMotherHistogramArray[motherType][histType][13])->Fill(dcaMother);
  ((TH1F*)fMotherHistogramArray[motherType][histType][14])->Fill(vertexDistance);
  ((TH1F*)fMotherHistogramArray[motherType][histType][15])->Fill(normDecayLength);
  ((TH1F*)fMotherHistogramArray[motherType][histType][16])->Fill(pseudoProperDecayTime);
  ((TH1F*)fMotherHistogramArray[motherType][histType][17])->Fill(decayTime);
  ((TH1F*)fMotherHistogramArray[motherType][histType][18])->Fill(normalizedDecayTime);
  ((TH1F*)fMotherHistogramArray[motherType][histType][19])->Fill(angleMotherFirstDaughter);
  ((TH1F*)fMotherHistogramArray[motherType][histType][20])->Fill(angleMotherSecondDaughter);
  ((TH1F*)fMotherHistogramArray[motherType][histType][21])->Fill(angleBetweenBothDaughters);
  ((TH1F*)fMotherHistogramArray[motherType][histType][22])->Fill(cosThetaStar);

  ((TH1F*)fMotherHistogramArray[motherType][histType][23])->Fill(vertexMotherX);
  ((TH1F*)fMotherHistogramArray[motherType][histType][24])->Fill(vertexMotherY);
  ((TH1F*)fMotherHistogramArray[motherType][histType][25])->Fill(vertexMotherZ);

  ((TH1F*)fMotherHistogramArray[motherType][histType][36])->Fill(TMath::Abs(dd0pr1));
  ((TH1F*)fMotherHistogramArray[motherType][histType][37])->Fill(TMath::Abs(dd0pr2));
  ((TH1F*)fMotherHistogramArray[motherType][histType][38])->Fill(TMath::Abs(dd0max));
  ((TH1F*)fMotherHistogramArray[motherType][histType][39])->Fill(TMath::Abs(dd0min));

  ((TH1F*)fMotherHistogramArray[motherType][histType][40])->Fill(cosPointingAngleXY);
  ((TH1F*)fMotherHistogramArray[motherType][histType][41])->Fill(distanceXYToVertex);
  ((TH1F*)fMotherHistogramArray[motherType][histType][42])->Fill(normalizedDecayLengthXY);
  ((TH1F*)fMotherHistogramArray[motherType][histType][43])->Fill(vertexMother->GetChi2());
  ((TH1F*)fMotherHistogramArray[motherType][histType][44])->Fill(vertexMother->GetChi2perNDF());

  ((TH1F*)fMotherHistogramArray[motherType][histType][45])->Fill(TMath::Abs(selectedMother->Getd0Prong(0)/selectedMother->Getd0errProng(0)));
  ((TH1F*)fMotherHistogramArray[motherType][histType][46])->Fill(TMath::Abs(selectedMother->Getd0Prong(1)/selectedMother->Getd0errProng(1)));
  ((TH1F*)fMotherHistogramArray[motherType][histType][47])->Fill(TMath::Abs(d0[0]/covd0z0[0]));
  ((TH1F*)fMotherHistogramArray[motherType][histType][48])->Fill((selectedMother->Getd0Prong(0)/selectedMother->Getd0errProng(0)) * (selectedMother->Getd0Prong(1)/selectedMother->Getd0errProng(1)));

  //we fill the 2D histograms
  Int_t nFirst = 0;
  Int_t nSecond = 1;
  Int_t nVariables = 10;
  Int_t nHistograms = nVariables * (nVariables - 1) / 2;
  for (Int_t k = 0; k < nHistograms; ++k)
  {
    Double_t firstVariable = 0.0;
    Double_t secondVariable = 0.0;

    if (nFirst == 0) firstVariable = d0firstTrack;
    if (nFirst == 1) firstVariable = d0secondTrack;
    if (nFirst == 2) firstVariable = d0Mother;
    if (nFirst == 3) firstVariable = pointingAngle;
    if (nFirst == 4) firstVariable = impactProduct;
    if (nFirst == 5) firstVariable = impactProductXY;
    if (nFirst == 6) firstVariable = vertexDistance;
    if (nFirst == 7) firstVariable = normDecayLength;
    if (nFirst == 8) firstVariable = cosPointingAngleXY;
    if (nFirst == 9) firstVariable = distanceXYToVertex;
    if (nFirst == 10) firstVariable = normalizedDecayLengthXY;

    if (nSecond == 0) secondVariable = d0firstTrack;
    if (nSecond == 1) secondVariable = d0secondTrack;
    if (nSecond == 2) secondVariable = d0Mother;
    if (nSecond == 3) secondVariable = pointingAngle;
    if (nSecond == 4) secondVariable = impactProduct;
    if (nSecond == 5) secondVariable = impactProductXY;
    if (nSecond == 6) secondVariable = vertexDistance;
    if (nSecond == 7) secondVariable = normDecayLength;
    if (nSecond == 8) secondVariable = cosPointingAngleXY;
    if (nSecond == 9) secondVariable = distanceXYToVertex;
    if (nSecond == 10) secondVariable = normalizedDecayLengthXY;

    ((TH2F*)fMotherHistogramArray2D[motherType][histType][k])->Fill(firstVariable, secondVariable);

    nSecond++;
    if (nSecond > nVariables)
    {
      nFirst++;
      nSecond = nFirst + 1;
    }
  }

  return;
}
//-------------------------------------------------------------------------------------
void AliAnalysisTaskSEBPlustoD0Pi::FillBPlusHistograms(AliAODRecoDecayHF2Prong * selectedMother, AliAODVertex *primaryVertex, Double_t bz, Int_t motherType, Int_t histType) {

  //In this function we fill the histograms of the reconstructed mothers
  Double_t ptMother = 0.0;
  Double_t momentumMother = 0.0;
  Double_t etaMother = 0.0;
  Double_t phiMother = 0.0;
  Double_t d0Mother = 0.0;
  Double_t d0firstTrack = 0.0;
  Double_t d0secondTrack = 0.0;
  Double_t pointingAngle = 0.0;
  Double_t impactProduct = 0.0;
  Double_t impactProductXY = 0.0;
  Double_t invariantMassMother = 0.0;
  Double_t invmassDelta = 0.0;
  Double_t dcaMother = 0.0;
  AliAODVertex * vertexMother = 0x0;
  Double_t vertexDistance = 0.0;
  Double_t normDecayLength = 0.0;
  Double_t decayTime = 0.0;
  Double_t angleMotherFirstDaughter = 0.0;
  Double_t angleMotherSecondDaughter = 0.0;
  Double_t ptFirstDaughter = 0.0;
  Double_t ptSecondDaughter = 0.0;
  UInt_t prongs[2];
  Double_t cosThetaStar = 0;
  Double_t angleBetweenBothDaughters = 0;
  // Double_t d0MultiProduct = 0;

  Double_t pdgMassMother = 0;

  AliExternalTrackParam motherTrack;
  motherTrack.CopyFromVTrack(selectedMother);
  Double_t d0z0[2], covd0z0[3], d0[2];
  motherTrack.PropagateToDCA(primaryVertex, bz, 100., d0z0, covd0z0);
  d0[0] = d0z0[0];

  //BPlus
  prongs[1] = 421; prongs[0] = 211;
  invmassDelta = DeltaInvMassBPlusKpipi(selectedMother);
  AliAODTrack * firstDaughter = (AliAODTrack*)selectedMother->GetDaughter(0);
  AliAODRecoDecayHF2Prong * secondDaughter = (AliAODRecoDecayHF2Prong*)selectedMother->GetDaughter(1);

  ptFirstDaughter = firstDaughter->Pt();
  ptSecondDaughter = secondDaughter->Pt();
  vertexMother = selectedMother->GetSecondaryVtx();
  angleMotherFirstDaughter = (selectedMother->Px() * firstDaughter->Px() + selectedMother->Py() * firstDaughter->Py() + selectedMother->Pz() * firstDaughter->Pz()) / (selectedMother->P() * firstDaughter->P());
  angleMotherSecondDaughter = (selectedMother->Px() * secondDaughter->Px() + selectedMother->Py() * secondDaughter->Py() + selectedMother->Pz() * secondDaughter->Pz()) / (selectedMother->P() * secondDaughter->P());
  angleBetweenBothDaughters  = (firstDaughter->Px() * secondDaughter->Px() + firstDaughter->Py() * secondDaughter->Py() + firstDaughter->Pz() * secondDaughter->Pz()) / (firstDaughter->P() * secondDaughter->P());
  cosThetaStar = selectedMother->CosThetaStar(0, 521, 421, 211);
  pdgMassMother = TDatabasePDG::Instance()->GetParticle(521)->Mass();
  // d0MultiProduct = selectedMother->Getd0Prong(0) * secondDaughter->Getd0Prong(0) * secondDaughter->Getd0Prong(1);

  //Topomatic
  Double_t dd0pr1 = 0.;
  Double_t dd0pr2 = 0.;
  Double_t dd0max = 0.;
  Double_t dd0min = 0.;
  for (Int_t ipr = 0; ipr < 2; ipr++)
  {
    Double_t diffIP, errdiffIP;
    selectedMother->Getd0MeasMinusExpProng(ipr, bz, diffIP, errdiffIP);
    Double_t normdd0 = 0.;
    if (errdiffIP > 0.) normdd0 = diffIP / errdiffIP;
    if (ipr == 0) dd0pr1 = normdd0;
    if (ipr == 1) dd0pr2 = normdd0;

    // else if(TMath::Abs(normdd0)>TMath::Abs(dd0max)) dd0max=normdd0;
  }
  if (TMath::Abs(dd0pr1) > TMath::Abs(dd0pr2)) {dd0max = dd0pr1; dd0min = dd0pr2;}
  else {dd0max = dd0pr2; dd0min = dd0pr1;}



  ptMother = selectedMother->Pt();
  momentumMother = selectedMother->P();
  etaMother = selectedMother->Eta();
  phiMother = selectedMother->Phi();
  d0Mother = TMath::Abs(d0[0]);

  pointingAngle = selectedMother->CosPointingAngle();
  impactProduct = selectedMother->Prodd0d0();
  impactProductXY = TMath::Abs(selectedMother->ImpParXY());
  invariantMassMother = selectedMother->InvMass(2, prongs);
  dcaMother = selectedMother->GetDCA();
  vertexDistance = vertexMother->DistanceToVertex(primaryVertex);
  d0firstTrack = TMath::Abs(selectedMother->Getd0Prong(0));
  d0secondTrack = TMath::Abs(selectedMother->Getd0Prong(1));
  normDecayLength = selectedMother->NormalizedDecayLength();

  Double_t pseudoProperDecayLength = ((vertexMother->GetX() - primaryVertex->GetX()) * selectedMother->Px() / TMath::Abs(selectedMother->Pt())) + ((vertexMother->GetY() - primaryVertex->GetY()) * selectedMother->Py() / TMath::Abs(selectedMother->Pt()));
  Double_t pseudoProperDecayTime = pseudoProperDecayLength * pdgMassMother / ptMother;
  decayTime = vertexDistance / (299792458 * TMath::Sqrt(1 / ((pdgMassMother * pdgMassMother / (momentumMother * momentumMother)) + 1)));

  Double_t phi = selectedMother->Phi();
  Double_t theta = selectedMother->Theta();
  Double_t covMatrix[21];
  selectedMother->GetCovarianceXYZPxPyPz(covMatrix);

  Double_t cp = TMath::Cos(phi);
  Double_t sp = TMath::Sin(phi);
  Double_t ct = TMath::Cos(theta);
  Double_t st = TMath::Sin(theta);

  Double_t errorMomentum = covMatrix[9] * cp * cp * ct * ct // GetCovPxPx
                           + covMatrix[13] * 2.*cp * sp * ct * ct // GetCovPxPy
                           + covMatrix[18] * 2.*cp * ct * st // GetCovPxPz
                           + covMatrix[14] * sp * sp * ct * ct // GetCovPyPy
                           + covMatrix[19] * 2.*sp * ct * st // GetCovPyPz
                           + covMatrix[20] * st * st; // GetCovPzPz
  Double_t normalizedDecayTime = selectedMother->NormalizedDecayLength() / (299792458 * TMath::Sqrt(1 / ((pdgMassMother * pdgMassMother * errorMomentum * errorMomentum / (momentumMother * momentumMother)) + 1)));

  Double_t vertexMotherX = vertexMother->GetX();
  Double_t vertexMotherY = vertexMother->GetY();
  Double_t vertexMotherZ = vertexMother->GetZ();

  Double_t cosPointingAngleXY = selectedMother->CosPointingAngleXY();
  Double_t distanceXYToVertex = vertexMother->DistanceXYToVertex(primaryVertex);
  Double_t normalizedDecayLengthXY = selectedMother->NormalizedDecayLengthXY();

  ((TH1F*)fMotherHistogramArray[motherType][histType][0])->Fill(ptMother);
  ((TH1F*)fMotherHistogramArray[motherType][histType][1])->Fill(ptFirstDaughter);
  ((TH1F*)fMotherHistogramArray[motherType][histType][2])->Fill(ptSecondDaughter);
  ((TH1F*)fMotherHistogramArray[motherType][histType][3])->Fill(etaMother);
  ((TH1F*)fMotherHistogramArray[motherType][histType][4])->Fill(phiMother);
  ((TH1F*)fMotherHistogramArray[motherType][histType][5])->Fill(d0Mother);
  ((TH1F*)fMotherHistogramArray[motherType][histType][6])->Fill(d0firstTrack);
  ((TH1F*)fMotherHistogramArray[motherType][histType][7])->Fill(d0secondTrack);
  ((TH1F*)fMotherHistogramArray[motherType][histType][8])->Fill(pointingAngle);
  ((TH1F*)fMotherHistogramArray[motherType][histType][9])->Fill(impactProduct);
  ((TH1F*)fMotherHistogramArray[motherType][histType][10])->Fill(impactProductXY);
  ((TH1F*)fMotherHistogramArray[motherType][histType][11])->Fill(invariantMassMother);
  ((TH1F*)fMotherHistogramArray[motherType][histType][12])->Fill(invmassDelta);
  ((TH1F*)fMotherHistogramArray[motherType][histType][13])->Fill(dcaMother);
  ((TH1F*)fMotherHistogramArray[motherType][histType][14])->Fill(vertexDistance);
  ((TH1F*)fMotherHistogramArray[motherType][histType][15])->Fill(normDecayLength);
  ((TH1F*)fMotherHistogramArray[motherType][histType][16])->Fill(pseudoProperDecayTime);
  ((TH1F*)fMotherHistogramArray[motherType][histType][17])->Fill(decayTime);
  ((TH1F*)fMotherHistogramArray[motherType][histType][18])->Fill(normalizedDecayTime);
  ((TH1F*)fMotherHistogramArray[motherType][histType][19])->Fill(angleMotherFirstDaughter);
  ((TH1F*)fMotherHistogramArray[motherType][histType][20])->Fill(angleMotherSecondDaughter);
  ((TH1F*)fMotherHistogramArray[motherType][histType][21])->Fill(angleBetweenBothDaughters);
  ((TH1F*)fMotherHistogramArray[motherType][histType][22])->Fill(cosThetaStar);

  ((TH1F*)fMotherHistogramArray[motherType][histType][23])->Fill(vertexMotherX);
  ((TH1F*)fMotherHistogramArray[motherType][histType][24])->Fill(vertexMotherY);
  ((TH1F*)fMotherHistogramArray[motherType][histType][25])->Fill(vertexMotherZ);

  ((TH1F*)fMotherHistogramArray[motherType][histType][36])->Fill(TMath::Abs(dd0pr1));
  ((TH1F*)fMotherHistogramArray[motherType][histType][37])->Fill(TMath::Abs(dd0pr2));
  ((TH1F*)fMotherHistogramArray[motherType][histType][38])->Fill(TMath::Abs(dd0max));
  ((TH1F*)fMotherHistogramArray[motherType][histType][39])->Fill(TMath::Abs(dd0min));

  ((TH1F*)fMotherHistogramArray[motherType][histType][40])->Fill(cosPointingAngleXY);
  ((TH1F*)fMotherHistogramArray[motherType][histType][41])->Fill(distanceXYToVertex);
  ((TH1F*)fMotherHistogramArray[motherType][histType][42])->Fill(normalizedDecayLengthXY);
  ((TH1F*)fMotherHistogramArray[motherType][histType][43])->Fill(vertexMother->GetChi2());
  ((TH1F*)fMotherHistogramArray[motherType][histType][44])->Fill(vertexMother->GetChi2perNDF());

  ((TH1F*)fMotherHistogramArray[motherType][histType][45])->Fill(TMath::Abs(selectedMother->Getd0Prong(0)/selectedMother->Getd0errProng(0)));
  ((TH1F*)fMotherHistogramArray[motherType][histType][46])->Fill(TMath::Abs(selectedMother->Getd0Prong(1)/selectedMother->Getd0errProng(1)));
  ((TH1F*)fMotherHistogramArray[motherType][histType][47])->Fill(TMath::Abs(d0[0]/covd0z0[0]));
  ((TH1F*)fMotherHistogramArray[motherType][histType][48])->Fill((selectedMother->Getd0Prong(0)/selectedMother->Getd0errProng(0)) * (selectedMother->Getd0Prong(1)/selectedMother->Getd0errProng(1)));


  //we fill the 2D histograms
  Int_t nFirst = 0;
  Int_t nSecond = 1;
  Int_t nVariables = 10;
  Int_t nHistograms = nVariables * (nVariables - 1) / 2;
  for (Int_t k = 0; k < nHistograms; ++k)
  {
    Double_t firstVariable = 0.0;
    Double_t secondVariable = 0.0;

    if (nFirst == 0) firstVariable = d0firstTrack;
    if (nFirst == 1) firstVariable = d0secondTrack;
    if (nFirst == 2) firstVariable = d0Mother;
    if (nFirst == 3) firstVariable = pointingAngle;
    if (nFirst == 4) firstVariable = impactProduct;
    if (nFirst == 5) firstVariable = impactProductXY;
    if (nFirst == 6) firstVariable = vertexDistance;
    if (nFirst == 7) firstVariable = normDecayLength;
    if (nFirst == 8) firstVariable = cosPointingAngleXY;
    if (nFirst == 9) firstVariable = distanceXYToVertex;
    if (nFirst == 10) firstVariable = normalizedDecayLengthXY;

    if (nSecond == 0) secondVariable = d0firstTrack;
    if (nSecond == 1) secondVariable = d0secondTrack;
    if (nSecond == 2) secondVariable = d0Mother;
    if (nSecond == 3) secondVariable = pointingAngle;
    if (nSecond == 4) secondVariable = impactProduct;
    if (nSecond == 5) secondVariable = impactProductXY;
    if (nSecond == 6) secondVariable = vertexDistance;
    if (nSecond == 7) secondVariable = normDecayLength;
    if (nSecond == 8) secondVariable = cosPointingAngleXY;
    if (nSecond == 9) secondVariable = distanceXYToVertex;
    if (nSecond == 10) secondVariable = normalizedDecayLengthXY;

    ((TH2F*)fMotherHistogramArray2D[motherType][histType][k])->Fill(firstVariable, secondVariable);

    nSecond++;
    if (nSecond > nVariables)
    {
      nFirst++;
      nSecond = nFirst + 1;
    }
  }


  if (motherType == 1) {
    motherType = motherType - 1;

    AliAODRecoDecay* trackD0 = (AliAODRecoDecay*)selectedMother->GetDaughter(1);
    AliAODTrack * firstDaughterD0 = (AliAODTrack*)trackD0->GetDaughter(0);
    AliAODTrack * secondDaughterD0 = (AliAODTrack*)trackD0->GetDaughter(1);

    AliAODVertex * vertexBPlus = vertexMother;
    AliAODVertex * vertexD0 = trackD0->GetSecondaryVtx();
    vertexDistance = TMath::Abs(vertexBPlus->DistanceToVertex(vertexD0));
    pdgMassMother = TDatabasePDG::Instance()->GetParticle(421)->Mass();

    AliExternalTrackParam firstDaughterD0Track;
    AliExternalTrackParam secondDaughterD0Track;

    Double_t d0z0[2], covd0z0[3], d0[2];

    firstDaughterD0Track.CopyFromVTrack(firstDaughterD0);
    firstDaughterD0Track.PropagateToDCA(vertexBPlus, bz, 100., d0z0, covd0z0);
    d0[0] = d0z0[0];

    secondDaughterD0Track.CopyFromVTrack(secondDaughterD0);
    secondDaughterD0Track.PropagateToDCA(vertexBPlus, bz, 100., d0z0, covd0z0);
    d0[1] = d0z0[0];

    AliExternalTrackParam D0Track;
    D0Track.CopyFromVTrack(trackD0);
    Double_t d0z0D0[2], covd0z0D0[3], d0D0;
    D0Track.PropagateToDCA(vertexBPlus, bz, 100., d0z0D0, covd0z0D0);
    d0D0 = d0z0D0[0];

    Double_t impactProductToBPlus = d0[0] * d0[1];
    Double_t impactProductXYToBPlus = trackD0->ImpParXY(vertexBPlus);

    Double_t momentumMother = trackD0->P();
    Double_t pointingAngleToBPlus = trackD0->CosPointingAngle(vertexBPlus);
    Double_t d0FirstDaughterToBPlus = TMath::Abs(d0[0]);
    Double_t d0trackD0ToBPlus = TMath::Abs(d0[1]);
    Double_t normDecayLengthToBPlus = trackD0->NormalizedDecayLength(vertexBPlus);

    Double_t pseudoProperDecayLength = ((vertexD0->GetX() - vertexBPlus->GetX()) * trackD0->Px() / TMath::Abs(trackD0->Pt())) + ((vertexD0->GetY() - vertexBPlus->GetY()) * trackD0->Py() / TMath::Abs(trackD0->Pt()));
    Double_t pseudoProperDecayTimeToBPlus = pseudoProperDecayLength * pdgMassMother / ptMother;
    Double_t DecayTimeToBPlus = vertexDistance / (299792458 * TMath::Sqrt(1 / ((pdgMassMother * pdgMassMother / (momentumMother * momentumMother)) + 1)));

    Double_t phi = trackD0->Phi();
    Double_t theta = trackD0->Theta();
    Double_t covMatrix[21];
    trackD0->GetCovarianceXYZPxPyPz(covMatrix);

    Double_t cp = TMath::Cos(phi);
    Double_t sp = TMath::Sin(phi);
    Double_t ct = TMath::Cos(theta);
    Double_t st = TMath::Sin(theta);

    Double_t errorMomentum = covMatrix[9] * cp * cp * ct * ct // GetCovPxPx
                             + covMatrix[13] * 2.*cp * sp * ct * ct // GetCovPxPy
                             + covMatrix[18] * 2.*cp * ct * st // GetCovPxPz
                             + covMatrix[14] * sp * sp * ct * ct // GetCovPyPy
                             + covMatrix[19] * 2.*sp * ct * st // GetCovPyPz
                             + covMatrix[20] * st * st; // GetCovPzPz
    Double_t normDecayTimeToBPlus = trackD0->NormalizedDecayLength(vertexBPlus) / (299792458 * TMath::Sqrt(1 / ((pdgMassMother * pdgMassMother * errorMomentum * errorMomentum / (momentumMother * momentumMother)) + 1)));

    ((TH1F*)fMotherHistogramArray[motherType][histType][26])->Fill(pointingAngleToBPlus);
    ((TH1F*)fMotherHistogramArray[motherType][histType][27])->Fill(d0D0);
    ((TH1F*)fMotherHistogramArray[motherType][histType][28])->Fill(d0FirstDaughterToBPlus);
    ((TH1F*)fMotherHistogramArray[motherType][histType][29])->Fill(d0trackD0ToBPlus);
    ((TH1F*)fMotherHistogramArray[motherType][histType][30])->Fill(impactProductToBPlus);
    ((TH1F*)fMotherHistogramArray[motherType][histType][31])->Fill(impactProductXYToBPlus);
    ((TH1F*)fMotherHistogramArray[motherType][histType][32])->Fill(normDecayLengthToBPlus);
    ((TH1F*)fMotherHistogramArray[motherType][histType][33])->Fill(pseudoProperDecayTimeToBPlus);
    ((TH1F*)fMotherHistogramArray[motherType][histType][34])->Fill(DecayTimeToBPlus);
    ((TH1F*)fMotherHistogramArray[motherType][histType][35])->Fill(normDecayTimeToBPlus);



  }
  return;
}
//-------------------------------------------------------------------------------------
void AliAnalysisTaskSEBPlustoD0Pi::TwoTrackCombinationInfo(AliExternalTrackParam * firstTrack, AliExternalTrackParam * secondTrack, AliAODVertex * primaryVertex, Double_t bz, Bool_t isDesiredCandidate, Int_t histogramNumber, UInt_t prongs[2]) {

  // we calculate the vertex position
  TObjArray daughterTracks;

  daughterTracks.Add(firstTrack);
  daughterTracks.Add(secondTrack);

  Double_t dispersion = 0;
  AliAODVertex *vertex = RecalculateVertex(primaryVertex, &daughterTracks, bz, dispersion, kTRUE, kTRUE, kFALSE);
  if (!vertex) {delete vertex; vertex = NULL; return;}

  Double_t xdummy = 0., ydummy = 0., dca;
  Double_t d0z0[2], covd0z0[3], d0[2], d0err[2];

  firstTrack->PropagateToDCA(vertex, bz, 100., d0z0, covd0z0);
  secondTrack->PropagateToDCA(vertex, bz, 100., d0z0, covd0z0);
  dca = secondTrack->GetDCA(firstTrack, bz, xdummy, ydummy);

  Double_t px[2], py[2], pz[2];
  px[0] = firstTrack->Px();
  py[0] = firstTrack->Py();
  pz[0] = firstTrack->Pz();
  px[1] = secondTrack->Px();
  py[1] = secondTrack->Py();
  pz[1] = secondTrack->Pz();

  firstTrack->PropagateToDCA(primaryVertex, bz, 100., d0z0, covd0z0);
  d0[0] = d0z0[0];
  d0err[0] = TMath::Sqrt(covd0z0[0]);
  secondTrack->PropagateToDCA(primaryVertex, bz, 100., d0z0, covd0z0);
  d0[1] = d0z0[0];
  d0err[1] = TMath::Sqrt(covd0z0[0]);

  AliAODRecoDecayHF2Prong * track = new AliAODRecoDecayHF2Prong(vertex, px, py, pz, d0, d0err, dca);
  if (!track)
  {
    delete vertex; vertex = NULL;
    delete track; track = NULL;
    return;
  }

  Double_t invariantMass = track->InvMass(2, prongs);

  if (!isDesiredCandidate) ((TH1F*)(fResultsHistogramArray[12][histogramNumber]))->Fill(invariantMass);
  if (isDesiredCandidate)  ((TH1F*)(fResultsHistogramArray[12][histogramNumber + 1]))->Fill(invariantMass);

  delete vertex; vertex = NULL;
  delete track; track = NULL;
  return;
}
//-------------------------------------------------------------------------------------
void AliAnalysisTaskSEBPlustoD0Pi::ThreeTrackCombinationInfo(AliExternalTrackParam * firstTrack, AliExternalTrackParam * secondTrack, AliExternalTrackParam * thirdTrack, AliAODVertex * primaryVertex, Double_t bz, Bool_t isDesiredCandidate, Int_t histogramNumber, UInt_t prongs[3]) {

  // we calculate the vertex position
  TObjArray daughterTracks;

  daughterTracks.Add(firstTrack);
  daughterTracks.Add(secondTrack);
  daughterTracks.Add(thirdTrack);

  Double_t dispersion = 0;
  AliAODVertex *vertex = RecalculateVertex(primaryVertex, &daughterTracks, bz, dispersion, kTRUE, kTRUE, kFALSE);
  if (!vertex) {delete vertex; vertex = NULL; return;}

  Double_t xdummy = 0., ydummy = 0., dca[3];
  Double_t d0z0[2], covd0z0[3], d0[3], d0err[3];

  firstTrack->PropagateToDCA(vertex, bz, 100., d0z0, covd0z0);
  secondTrack->PropagateToDCA(vertex, bz, 100., d0z0, covd0z0);
  thirdTrack->PropagateToDCA(vertex, bz, 100., d0z0, covd0z0);

  dca[0] = firstTrack->GetDCA(secondTrack, bz, xdummy, ydummy);
  dca[1] = firstTrack->GetDCA(thirdTrack, bz, xdummy, ydummy);
  dca[2] = secondTrack->GetDCA(thirdTrack, bz, xdummy, ydummy);

  Double_t px[3], py[3], pz[3];
  px[0] = firstTrack->Px();
  py[0] = firstTrack->Py();
  pz[0] = firstTrack->Pz();
  px[1] = secondTrack->Px();
  py[1] = secondTrack->Py();
  pz[1] = secondTrack->Pz();
  px[2] = thirdTrack->Px();
  py[2] = thirdTrack->Py();
  pz[2] = thirdTrack->Pz();

  firstTrack->PropagateToDCA(primaryVertex, bz, 100., d0z0, covd0z0);
  d0[0] = d0z0[0];
  d0err[0] = TMath::Sqrt(covd0z0[0]);
  secondTrack->PropagateToDCA(primaryVertex, bz, 100., d0z0, covd0z0);
  d0[1] = d0z0[0];
  d0err[1] = TMath::Sqrt(covd0z0[0]);
  thirdTrack->PropagateToDCA(primaryVertex, bz, 100., d0z0, covd0z0);
  d0[2] = d0z0[0];
  d0err[2] = TMath::Sqrt(covd0z0[0]);

  // Double_t pos[3]; primaryVertex->GetXYZ(pos);
  // Double_t dist12=TMath::Sqrt((vertex12->GetX()-pos[0])*(vertex12->GetX()-pos[0])+(vertex12->GetY()-pos[1])*(vertex12->GetY()-pos[1])+(vertex12->GetZ()-pos[2])*(vertex12->GetZ()-pos[2]));
  // Double_t dist23=TMath::Sqrt((vertex23->GetX()-pos[0])*(vertex23->GetX()-pos[0])+(vertex23->GetY()-pos[1])*(vertex23->GetY()-pos[1])+(vertex23->GetZ()-pos[2])*(vertex23->GetZ()-pos[2]));
  Short_t charge = (Short_t)(firstTrack->Charge() + secondTrack->Charge() + thirdTrack->Charge());

  AliAODRecoDecayHF3Prong * track = new AliAODRecoDecayHF3Prong(vertex, px, py, pz, d0, d0err, dca, dispersion, 0.0, 0.0, charge); //dist12,dist23,charge);
  if (!track)
  {
    delete vertex; vertex = NULL;
    delete track; track = NULL;
    return;
  }

  Double_t invariantMass = track->InvMass(3, prongs);

  if (!isDesiredCandidate) ((TH1F*)(fResultsHistogramArray[12][histogramNumber]))->Fill(invariantMass);
  if (isDesiredCandidate)  ((TH1F*)(fResultsHistogramArray[12][histogramNumber + 1]))->Fill(invariantMass);

  delete vertex; vertex = NULL;
  delete track; track = NULL;
  return;
}
//-------------------------------------------------------------------------------------
Int_t AliAnalysisTaskSEBPlustoD0Pi::MatchCandidateToMonteCarlo(Int_t pdgabs, AliAODRecoDecayHF2Prong * candidate, TClonesArray *mcArray, TMatrix * BPlustoD0PiLabelMatrix, Bool_t bCheckLabel) const
{
  //
  // Check if this candidate is matched to a MC signal
  // If no, return -1
  // If yes, return label (>=0) of the AliAODMCParticle


  // Check number of daughters
  Int_t ndg = candidate->GetNDaughters();
  if (!ndg) { AliError("No daughters available"); return -1;}
  if (ndg != 2) return -1;

  // loop on daughters and write the labels
  Int_t dgLabels[2] = { -1};
  Int_t pdgDg[2] = {0};
  Int_t signalPosition = -1;
  if (pdgabs == 421)
  {
    AliAODTrack *trk0 = (AliAODTrack*)candidate->GetDaughter(0);
    dgLabels[0] = trk0->GetLabel();
    AliAODTrack *trk1 = (AliAODTrack*)candidate->GetDaughter(1);
    dgLabels[1] = trk1->GetLabel();
    pdgDg[0] = 211; pdgDg[1] = 321;
    signalPosition = 3;
  }
  else if (pdgabs == 521)
  {
    AliAODTrack *trk0 = (AliAODTrack*)candidate->GetDaughter(0);
    dgLabels[0] = trk0->GetLabel();
    dgLabels[1] = MatchCandidateToMonteCarlo(421, (AliAODRecoDecayHF2Prong*)candidate->GetDaughter(1), mcArray, BPlustoD0PiLabelMatrix);
    pdgDg[0] = 211; pdgDg[1] = 421;
    signalPosition = 4;
  }
  else
  {
    std::cout << "Wrong pdg supplied for function to match candidate to monte carlo signal." << std::endl;
    return -1;
  }
  if (dgLabels[0] == -1) return -1;
  if (dgLabels[1] == -1) return -1;


  Int_t labMom[2] = {0, 0};
  Int_t i, j, lab, labMother, pdgMother, pdgPart;
  AliAODMCParticle *part = 0;
  AliAODMCParticle *mother = 0;
  Double_t pxSumDgs = 0., pySumDgs = 0., pzSumDgs = 0.;
  Bool_t pdgUsed[2] = {kFALSE, kFALSE};

  // loop on daughter labels
  for (i = 0; i < ndg; i++)
  {
    labMom[i] = -1;
    lab = TMath::Abs(dgLabels[i]);
    if (lab < 0)
    {
      printf("daughter with negative label %d\n", lab);
      return -1;
    }
    part = (AliAODMCParticle*)mcArray->At(lab);
    if (!part)
    {
      printf("no MC particle\n");
      return -1;
    }

    // check the PDG of the daughter
    pdgPart = TMath::Abs(part->GetPdgCode());
    for (j = 0; j < ndg; j++)
    {
      if (!pdgUsed[j] && pdgPart == pdgDg[j])
      {
        pdgUsed[j] = kTRUE;
        break;
      }
    }


    mother = part;
    while (mother->GetMother() >= 0)
    {
      labMother = mother->GetMother();
      mother = (AliAODMCParticle*)mcArray->At(labMother);
      if (!mother)
      {
        printf("no MC mother particle\n");
        break;
      }
      pdgMother = TMath::Abs(mother->GetPdgCode());
      if (pdgMother == pdgabs)
      {
        labMom[i] = labMother;
        // keep sum of daughters' momenta, to check for mom conservation
        pxSumDgs += part->Px();
        pySumDgs += part->Py();
        pzSumDgs += part->Pz();
        break;
      }
      else break;
    }
    if (labMom[i] == -1) return -1; // mother PDG not ok for this daughter
  } // end loop on daughters

  // check if the candidate is signal
  labMother = labMom[0];
  // all labels have to be the same and !=-1
  for (i = 0; i < ndg; i++)
  {
    if (labMom[i] == -1)        return -1;
    if (labMom[i] != labMother) return -1;
  }

  // check that all daughter PDGs are matched
  for (i = 0; i < ndg; i++)
  {
    if (pdgUsed[i] == kFALSE) return -1;
  }

  // Check for mom conservation
  mother = (AliAODMCParticle*)mcArray->At(labMother);
  Double_t pxMother = mother->Px();
  Double_t pyMother = mother->Py();
  Double_t pzMother = mother->Pz();


  // check the number of daughters (we are not looking at resonant decay)
  if(mother->GetNDaughters() != 2) return -1;

  // if momentum conservation is not within 0.5%, show warning. This can be due to large propagation distance through magnetic field.
  if ((TMath::Abs(pxMother - pxSumDgs) / (TMath::Abs(pxMother) + 1.e-13)) > 0.005 ||
      (TMath::Abs(pyMother - pySumDgs) / (TMath::Abs(pyMother) + 1.e-13)) > 0.005 ||
      (TMath::Abs(pzMother - pzSumDgs) / (TMath::Abs(pzMother) + 1.e-13)) > 0.005)
  {
    std::cout << std::endl << " Momentum difference for decay pdgabs = " << pdgabs << "daughters = " << mother->GetNDaughters() << std::endl;
    std::cout << "pxMother = " << pxMother << "pyMother = " << pyMother << "pzMother = " << pzMother << std::endl;
    std::cout << "pxSumDgs = " << pxSumDgs << "pySumDgs = " << pySumDgs << "pzSumDgs = " << pzSumDgs << std::endl;
  }

  // Check if label matches a signal label
  if(bCheckLabel)
  {
    Int_t bIsSignal = kFALSE;
    TMatrix &particleMatrix = *BPlustoD0PiLabelMatrix;
    for (Int_t k = 0; k < BPlustoD0PiLabelMatrix->GetNrows(); ++k)
    {
      if(labMother == (Int_t)particleMatrix(k,signalPosition))
      {
        bIsSignal = kTRUE;
        break;
      }
    }
    if(!bIsSignal) return -1;
  }


  return labMother;
}
//-------------------------------------------------------------------------------------
Int_t AliAnalysisTaskSEBPlustoD0Pi::IsTrackInjected(AliAODTrack *part, AliAODMCHeader *header, TClonesArray *arrayMC) {

  AliVertexingHFUtils* ggg = new  AliVertexingHFUtils();

  Int_t lab = part->GetLabel();
  if (lab < 0) {delete ggg; ggg = nullptr; return 1;} //
  TString nameGen = ggg->GetGenerator(lab, header);
  TString empty = "";
  Int_t countControl = 0;
  while (nameGen.IsWhitespace()) {
    AliAODMCParticle *mcpart = (AliAODMCParticle*)arrayMC->At(lab);
    if (!mcpart) {
      printf("AliVertexingHFUtils::IsTrackInjected - BREAK: No valid AliAODMCParticle at label %i\n", lab);
      break;
    }
    Int_t mother = mcpart->GetMother();
    if (mother < 0) {
      // printf("AliVertexingHFUtils::IsTrackInjected - BREAK: Reached primary particle without valid mother\n");
      break;
    }
    lab = mother;
    nameGen = ggg->GetGenerator(mother, header);
    countControl++;
    if (countControl >= 10) { // 10 = arbitrary number; protection from infinite loops
      printf("AliVertexingHFUtils::IsTrackInjected - BREAK: Protection from infinite loop active\n");
      break;
    }
  }
  if (nameGen.IsWhitespace() || nameGen.Contains("ijing")) {delete ggg; ggg = nullptr; return 0;}

  delete ggg; ggg = nullptr;
  return 1;
}
//-------------------------------------------------------------------------------------
Bool_t AliAnalysisTaskSEBPlustoD0Pi::IsCandidateInjected(AliAODRecoDecayHF2Prong *selectedBPlus, AliAODMCHeader *header, TClonesArray *arrayMC) {

  AliAODTrack* selectedBPlusPion = (AliAODTrack*)selectedBPlus->GetDaughter(0);
  AliAODRecoDecayHF2Prong* selectedD0 = (AliAODRecoDecayHF2Prong*)selectedBPlus->GetDaughter(1);

  AliAODTrack* selectedD0FirstDaughter = (AliAODTrack*)selectedD0->GetDaughter(0);
  AliAODTrack* selectedD0SecondDaughter = (AliAODTrack*)selectedD0->GetDaughter(1);

  if (IsTrackInjected(selectedBPlusPion, header, arrayMC)) return kTRUE;
  if (IsTrackInjected(selectedD0FirstDaughter, header, arrayMC)) return kTRUE;
  if (IsTrackInjected(selectedD0SecondDaughter, header, arrayMC)) return kTRUE;

  return kFALSE;
}
//-------------------------------------------------------------------------------------
void AliAnalysisTaskSEBPlustoD0Pi::CutOptimizationLoop(Int_t variable, Int_t nVariables, Int_t nCuts, Int_t ptBin, Int_t fillNumber, Bool_t isDesiredCandidate, Int_t nSigmaBin) {

  for (Int_t iCut = 0; iCut < nCuts; ++iCut)
  {
    Int_t cutVariable = fCuts->GetCutIndexForCutOptimization(variable);
    Bool_t isUpperCut = fCuts->GetIsUpperCutForCutOptimization(variable);
    Float_t cutValue = fCuts->GetCutForCutOptimization(iCut, variable, ptBin);
    if(fCutVariableValueArray[cutVariable] > cutValue && isUpperCut == kTRUE) continue;
    if(fCutVariableValueArray[cutVariable] < cutValue && isUpperCut == kFALSE) continue;

    Int_t fill = iCut * TMath::Power(nCuts,variable) + fillNumber;
    if( (variable + 1) == nVariables) 
    {
      if(isDesiredCandidate) ((TH2F*)(fResultsHistogramArray2D[0][ptBin]))->Fill(fill,nSigmaBin);  
      if(!isDesiredCandidate) ((TH2F*)(fResultsHistogramArray2D[1][ptBin]))->Fill(fill,nSigmaBin);
    }
    else{CutOptimizationLoop(variable + 1, nVariables, nCuts, ptBin, fill, isDesiredCandidate, nSigmaBin);}
  }
  return;
}
//-------------------------------------------------------------------------------------
void AliAnalysisTaskSEBPlustoD0Pi::CutOptimizationVariableValues(AliAODRecoDecayHF2Prong * candidateBPlus, AliAODEvent*  aod) {

  if(!candidateBPlus){
    std::cout<<"candidateBPlus null"<<std::endl;
    return;
  } 

  AliAODRecoDecayHF2Prong* candidateD0 = (AliAODRecoDecayHF2Prong*)candidateBPlus->GetDaughter(1);  
  if(!candidateD0){
    std::cout<<"candidateD0 null"<<std::endl;
    return;
  }

  AliAODTrack *candidatePion = (AliAODTrack*)candidateBPlus->GetDaughter(0);
  if(!candidatePion){
    std::cout<<"candidatePion null 1"<<std::endl;
    return;
  }
  
  AliAODTrack *candidateFirstDaughter = (AliAODTrack*)candidateD0->GetDaughter(0);
  if(!candidateFirstDaughter){
    std::cout<<"candidatePion null 2"<<std::endl;
    return;
  }

  AliAODTrack *candidateSecondDaughter = (AliAODTrack*)candidateD0->GetDaughter(1);
  if(!candidateSecondDaughter){
    std::cout<<"candidateKaon null"<<std::endl;
    return;
  }

  AliAODVertex * vertexBPlus = candidateBPlus->GetSecondaryVtx();
  if(!vertexBPlus){
    std::cout<<"vertexBPlus null"<<std::endl;
    return;
  }

  AliAODVertex * vertexD0 = candidateD0->GetSecondaryVtx();
  if(!vertexD0){
    std::cout<<"vertexD0 null"<<std::endl;
    return;
  }


  AliAODVertex * primaryVertex = aod->GetPrimaryVertex();
  if(!primaryVertex){
    std::cout<<"primaryVertex null"<<std::endl;
    return;
  }

  //get the magnetic field
  Double_t bz = (Double_t)aod->GetMagneticField(); 


  // save D0 variable information
  if(kTRUE) 
  {
    // D0mass
    Double_t mD0PDG = TDatabasePDG::Instance()->GetParticle(421)->Mass();
  
    // D0 window - invariant mass
    Int_t chargeBPlus = candidateBPlus->Charge();
    UInt_t prongs[2];
    if (chargeBPlus == -1)
    {
      prongs[0] = 211;
      prongs[1] = 321;
    }
    else if (chargeBPlus == 1)
    {
      prongs[0] = 321;
      prongs[1] = 211;
    }
    else 
    {
      std::cout << "Wrong charge BPlus." << std::endl;
      return;
    }
    Double_t invMassD0 = candidateD0->InvMass(2,prongs);
    Double_t invMassDifference = TMath::Abs(mD0PDG - invMassD0);

    Double_t pointingAngle = candidateD0->CosPointingAngle();
    Double_t dcaMother = candidateD0->GetDCA();
    Double_t ptMother = candidateD0->Pt();
    Double_t momentumMother = candidateD0->P();
    Double_t ptPion = candidateFirstDaughter->Pt();
    Double_t ptKaon = candidateSecondDaughter->Pt();

    AliExternalTrackParam motherTrack;
    motherTrack.CopyFromVTrack(candidateD0);
    Double_t d0z0[2],covd0z0[3],d0[2];
    motherTrack.PropagateToDCA(primaryVertex,bz,100.,d0z0,covd0z0);
    d0[0] = d0z0[0];
    Double_t d0Mother = TMath::Abs(d0[0]);
    Double_t d0firstTrack = TMath::Abs(candidateD0->Getd0Prong(0));
    Double_t d0secondTrack = TMath::Abs(candidateD0->Getd0Prong(1));  
    
    Double_t impactProduct = candidateD0->Prodd0d0();
    Double_t impactProductXY = TMath::Abs(candidateD0->ImpParXY());  

    Double_t angleBetweenBothDaughters  = (candidateSecondDaughter->Px() * candidateFirstDaughter->Px() + candidateSecondDaughter->Py() * candidateFirstDaughter->Py() + candidateSecondDaughter->Pz() * candidateFirstDaughter->Pz()) /(candidateSecondDaughter->P() * candidateFirstDaughter->P());
    Double_t angleMotherFirstDaughter = (candidateD0->Px() * candidateFirstDaughter->Px() + candidateD0->Py() * candidateFirstDaughter->Py() + candidateD0->Pz() * candidateFirstDaughter->Pz()) /(candidateD0->P() * candidateFirstDaughter->P());
    Double_t angleMotherSecondDaughter = (candidateD0->Px() * candidateSecondDaughter->Px() + candidateD0->Py() * candidateSecondDaughter->Py() + candidateD0->Pz() * candidateSecondDaughter->Pz()) /(candidateD0->P() * candidateSecondDaughter->P());

    Double_t cosThetaStar = candidateD0->CosThetaStar(0,421,211,321);
    Double_t vertexDistance = vertexD0->DistanceToVertex(primaryVertex);
    Double_t normDecayLength = candidateD0->NormalizedDecayLength();
    Double_t pdgMassMother = TDatabasePDG::Instance()->GetParticle(421)->Mass();
    Double_t pseudoProperDecayLength = ((vertexD0->GetX() - primaryVertex->GetX()) * candidateD0->Px() / TMath::Abs(candidateD0->Pt())) + ((vertexD0->GetY() - primaryVertex->GetY()) * candidateD0->Py() / TMath::Abs(candidateD0->Pt()));
    Double_t pseudoProperDecayTime = pseudoProperDecayLength * pdgMassMother/ptMother;
    Double_t decayTime = vertexDistance / (299792458 * TMath::Sqrt(1/((pdgMassMother*pdgMassMother/(momentumMother*momentumMother)) + 1)));

    Double_t phi = candidateD0->Phi();
    Double_t theta = candidateD0->Theta();
    Double_t covMatrix[21];
    candidateD0->GetCovarianceXYZPxPyPz(covMatrix);

    Double_t cp = TMath::Cos(phi);
    Double_t sp = TMath::Sin(phi);
    Double_t ct = TMath::Cos(theta);
    Double_t st = TMath::Sin(theta);

    Double_t errorMomentum = covMatrix[9]*cp*cp*ct*ct  // GetCovPxPx
                            +covMatrix[13]*2.*cp*sp*ct*ct  // GetCovPxPy
                            +covMatrix[18]*2.*cp*ct*st  // GetCovPxPz
                            +covMatrix[14]*sp*sp*ct*ct  // GetCovPyPy
                            +covMatrix[19]*2.*sp*ct*st  // GetCovPyPz
                            +covMatrix[20]*st*st;  // GetCovPzPz
    Double_t normalizedDecayTime = candidateD0->NormalizedDecayLength() / (299792458 * TMath::Sqrt(1/((pdgMassMother*pdgMassMother*errorMomentum*errorMomentum/(momentumMother*momentumMother)) + 1)));

    Double_t cosPointingAngleXY = candidateD0->CosPointingAngleXY();
    Double_t distanceXYToVertex = vertexD0->DistanceXYToVertex(primaryVertex);
    Double_t normalizedDecayLengthXY = candidateD0->NormalizedDecayLengthXY();
    Double_t chi2Vertex = vertexD0->GetChi2perNDF();


    //Topomatic
    Double_t dd0pr1=0.;
    Double_t dd0pr2=0.;
    Double_t dd0max=0.;
    Double_t dd0min=0.;
    for(Int_t ipr=0; ipr<2; ipr++) 
    {
      Double_t diffIP, errdiffIP;
      candidateD0->Getd0MeasMinusExpProng(ipr,bz,diffIP,errdiffIP);
      Double_t normdd0=0.;
      if(errdiffIP>0.) normdd0=diffIP/errdiffIP;
      if(ipr==0) dd0pr1=normdd0;
      if(ipr==1) dd0pr2=normdd0;
    }
    if(TMath::Abs(dd0pr1)>TMath::Abs(dd0pr2)) {dd0max=dd0pr1; dd0min=dd0pr2;}
    else {dd0max=dd0pr2; dd0min=dd0pr1;}


    // We apply the cuts 
    Int_t nCutIndex = 0;

    // "inv. mass width [GeV]" --------------------------------------------
    nCutIndex = 0;
    fCutVariableValueArray[nCutIndex] = invMassDifference;
    //---------------------------------------------------------------------

    // "delta mass width [GeV]" -------------------------------------------
    nCutIndex = 1; // not used for D0
    fCutVariableValueArray[nCutIndex] = 0;
    //---------------------------------------------------------------------

    // "pointing angle [Cos(theta)]" --------------------------------------
    nCutIndex = 2;
    fCutVariableValueArray[nCutIndex] = pointingAngle;
    //---------------------------------------------------------------------

    // "dca [cm]" ---------------------------------------------------------
    nCutIndex = 3;
    fCutVariableValueArray[nCutIndex] = dcaMother;
    //---------------------------------------------------------------------

    // "Pt D0 [GeV/c]" ----------------------------------------------------
    nCutIndex = 4;
    fCutVariableValueArray[nCutIndex] = ptMother;
    //---------------------------------------------------------------------

    // "Pt Kaon [GeV/c]" -------------------------------------------------
    nCutIndex = 5;
    fCutVariableValueArray[nCutIndex] = ptKaon;
    //---------------------------------------------------------------------

    // "Pt Pion [GeV/c]" --------------------------------------------------
    nCutIndex = 6;
    fCutVariableValueArray[nCutIndex] = ptPion;
    //---------------------------------------------------------------------

    // "d0 D0 [cm]" -------------------------------------------------------
    nCutIndex = 7;
    fCutVariableValueArray[nCutIndex] = d0Mother;
    //---------------------------------------------------------------------

    // "d0 Kaon [cm]"-----------------------------------------------------
    nCutIndex = 8;
    fCutVariableValueArray[nCutIndex] = d0firstTrack;
    //---------------------------------------------------------------------

    // "d0 Pion [cm]" -----------------------------------------------------
    nCutIndex = 9;
    fCutVariableValueArray[nCutIndex] = d0secondTrack;
    //---------------------------------------------------------------------

    // "d0d0 [cm^2]" ------------------------------------------------------
    nCutIndex = 10;
    fCutVariableValueArray[nCutIndex] = impactProduct;
    //---------------------------------------------------------------------

    // "d0d0 XY [cm^2]" ---------------------------------------------------
    nCutIndex = 11;
    fCutVariableValueArray[nCutIndex] = impactProductXY;
    //---------------------------------------------------------------------

    // "angle between both daughters" -------------------------------------
    nCutIndex = 12;
    fCutVariableValueArray[nCutIndex] = angleBetweenBothDaughters;
    //---------------------------------------------------------------------

    // "angle mother with first daughter" ---------------------------------
    nCutIndex = 13;
    fCutVariableValueArray[nCutIndex] = angleMotherFirstDaughter;
    //---------------------------------------------------------------------

    // "angle mother with second daughter" --------------------------------
    nCutIndex = 14;
    fCutVariableValueArray[nCutIndex] = angleMotherSecondDaughter;
    //---------------------------------------------------------------------

    // "cosThetaStar" -----------------------------------------------------
    nCutIndex = 15;
    fCutVariableValueArray[nCutIndex] = cosThetaStar;
    //---------------------------------------------------------------------

    // "vertexDistance" ---------------------------------------------------
    nCutIndex = 16;
    fCutVariableValueArray[nCutIndex] = vertexDistance;
    //---------------------------------------------------------------------

    // "pseudoProperDecayTime" --------------------------------------------
    nCutIndex = 17;
    fCutVariableValueArray[nCutIndex] = pseudoProperDecayTime;
    //---------------------------------------------------------------------

    // "DecayTime" --------------------------------------------------------
    nCutIndex = 18;
    fCutVariableValueArray[nCutIndex] = decayTime;
    //---------------------------------------------------------------------

    // "normalizedDecayTime" ----------------------------------------------------
    nCutIndex = 19;
    fCutVariableValueArray[nCutIndex] = normalizedDecayTime;
    //---------------------------------------------------------------------

    // "normDecayLength" --------------------------------------------------
    nCutIndex = 20;
    fCutVariableValueArray[nCutIndex] = normDecayLength;
    //---------------------------------------------------------------------

    // "topomatic first daughter" -----------------------------------------
    nCutIndex = 21;
    fCutVariableValueArray[nCutIndex] = dd0pr1;
    //---------------------------------------------------------------------

    // "topomatic second daughter" ----------------------------------------
    nCutIndex = 22;
    fCutVariableValueArray[nCutIndex] = dd0pr2;
    //---------------------------------------------------------------------

    // "topomatic max" ----------------------------------------------------
    nCutIndex = 23;
    fCutVariableValueArray[nCutIndex] = dd0max;
    //---------------------------------------------------------------------

    // "topomatic min" ----------------------------------------------------
    nCutIndex = 24;
    fCutVariableValueArray[nCutIndex] = dd0min;
    //---------------------------------------------------------------------

    // "pointing angle XY" ----------------------------------------------------
    nCutIndex = 25;
    fCutVariableValueArray[nCutIndex] = cosPointingAngleXY;
    //---------------------------------------------------------------------

     // "vertex distance XY" ----------------------------------------------------
    nCutIndex = 26;
    fCutVariableValueArray[nCutIndex] = distanceXYToVertex;
    //---------------------------------------------------------------------
    
     // "normalized decay length XY" ----------------------------------------------------
    nCutIndex = 27;
    fCutVariableValueArray[nCutIndex] = normalizedDecayLengthXY;
    //---------------------------------------------------------------------
    
    // "chi squared per NDF" ----------------------------------------------------
    nCutIndex = 28;
    fCutVariableValueArray[nCutIndex] = chi2Vertex;
    //---------------------------------------------------------------------          


    
    AliAODRecoDecay* candidateD0toBPlus = (AliAODRecoDecay*)candidateD0;
    AliExternalTrackParam firstDaughterD0Track;
    AliExternalTrackParam secondDaughterD0Track;

    Double_t d0z0DSVert[2],covd0z0DSVert[3],d0DSVert[2];

    firstDaughterD0Track.CopyFromVTrack(candidateFirstDaughter);
    firstDaughterD0Track.PropagateToDCA(vertexBPlus,bz,100.,d0z0DSVert,covd0z0DSVert);
    d0DSVert[0] = d0z0DSVert[0];

    secondDaughterD0Track.CopyFromVTrack(candidateSecondDaughter);
    secondDaughterD0Track.PropagateToDCA(vertexBPlus,bz,100.,d0z0DSVert,covd0z0DSVert);
    d0DSVert[1] = d0z0DSVert[0];

    AliExternalTrackParam D0Track;
    D0Track.CopyFromVTrack(candidateD0);
    Double_t d0z0D0DSVert[2],covd0z0D0DSVert[3],d0D0DSVert;
    motherTrack.PropagateToDCA(vertexBPlus,bz,100.,d0z0D0DSVert,covd0z0D0DSVert);
    d0D0DSVert = TMath::Abs(d0z0D0DSVert[0]);

    Double_t impactProductToBPlus = d0DSVert[0]*d0DSVert[1];
    Double_t impactProductXYToBPlus = candidateD0toBPlus->ImpParXY(vertexBPlus);

    Double_t pointingAngleToBPlus = candidateD0toBPlus->CosPointingAngle(vertexBPlus);
    Double_t d0FirstDaughterToBPlus = TMath::Abs(d0DSVert[0]);
    Double_t d0SecondDaughterToBPlus = TMath::Abs(d0DSVert[1]);
    Double_t normDecayLengthToBPlus = candidateD0toBPlus->NormalizedDecayLength(vertexBPlus);

    Double_t pseudoProperDecayLengthDSVert = ((vertexD0->GetX() - vertexBPlus->GetX()) * candidateD0->Px() / TMath::Abs(candidateD0->Pt())) + ((vertexD0->GetY() - vertexBPlus->GetY()) * candidateD0->Py() / TMath::Abs(candidateD0->Pt()));
    Double_t pseudoProperDecayTimeToBPlus = pseudoProperDecayLengthDSVert * pdgMassMother/ptMother;
    Double_t DecayTimeToBPlus = vertexDistance / (299792458 * TMath::Sqrt(1/((pdgMassMother*pdgMassMother/(momentumMother*momentumMother)) + 1)));

    Double_t phiDSVert = candidateD0->Phi();
    Double_t thetaDSVert = candidateD0->Theta();
    Double_t covMatrixDSVert[21];
    candidateD0->GetCovarianceXYZPxPyPz(covMatrixDSVert);

    cp = TMath::Cos(phiDSVert);
    sp = TMath::Sin(phiDSVert);
    ct = TMath::Cos(thetaDSVert);
    st = TMath::Sin(thetaDSVert);

    errorMomentum = covMatrix[9]*cp*cp*ct*ct  // GetCovPxPx
                            +covMatrix[13]*2.*cp*sp*ct*ct  // GetCovPxPy
                            +covMatrix[18]*2.*cp*ct*st  // GetCovPxPz
                            +covMatrix[14]*sp*sp*ct*ct  // GetCovPyPy
                            +covMatrix[19]*2.*sp*ct*st  // GetCovPyPz
                            +covMatrix[20]*st*st;  // GetCovPzPz
    Double_t normalizedDecayTimeToBPlus = candidateD0toBPlus->NormalizedDecayLength(vertexBPlus) / (299792458 * TMath::Sqrt(1/((pdgMassMother*pdgMassMother*errorMomentum*errorMomentum/(momentumMother*momentumMother)) + 1)));

    // "pointingAngleToBPlus" ---------------------------------------------
    nCutIndex = 29;
    fCutVariableValueArray[nCutIndex] = pointingAngleToBPlus;
    //---------------------------------------------------------------------

    // "d0MotherToBPlus" --------------------------------------------------
    nCutIndex = 30;
    fCutVariableValueArray[nCutIndex] = d0D0DSVert;
    //---------------------------------------------------------------------

    // "d0FirstDaughterToBPlus" -------------------------------------------
    nCutIndex = 31;
    fCutVariableValueArray[nCutIndex] = d0FirstDaughterToBPlus;
    //---------------------------------------------------------------------

    // "d0SecondDaughterToBPlus" ------------------------------------------
    nCutIndex = 32;
    fCutVariableValueArray[nCutIndex] = d0SecondDaughterToBPlus;
    //---------------------------------------------------------------------

    // "impactProductToBPlus" ---------------------------------------------
    nCutIndex = 33;
    fCutVariableValueArray[nCutIndex] = impactProductToBPlus;
    //---------------------------------------------------------------------

    // "impactProductXYToBPlus" -------------------------------------------
    nCutIndex = 34;
    fCutVariableValueArray[nCutIndex] = impactProductXYToBPlus;
    //---------------------------------------------------------------------

    // "normDecayLengthToBPlus" -------------------------------------------
    nCutIndex = 35;
    fCutVariableValueArray[nCutIndex] = normDecayLengthToBPlus;
    //---------------------------------------------------------------------

    // "pseudoProperDecayTimeToBPlus" -------------------------------------
    nCutIndex = 36;
    fCutVariableValueArray[nCutIndex] = pseudoProperDecayTimeToBPlus;
    //---------------------------------------------------------------------

    // "DecayTimeToBPlus" -------------------------------------------------
    nCutIndex = 37;
    fCutVariableValueArray[nCutIndex] = DecayTimeToBPlus;
    //---------------------------------------------------------------------

    // "normalizedDecayTimeToBPlus" ---------------------------------------------
    nCutIndex = 38;
    fCutVariableValueArray[nCutIndex] = normalizedDecayTimeToBPlus;
    //---------------------------------------------------------------------
  }
 
 // save B0 variable information
  if(kTRUE) 
  {     
    // We obtain the variable values in the section below 
    // D0Mass and BPlusmass
    Double_t mD0PDG = TDatabasePDG::Instance()->GetParticle(421)->Mass();
    Double_t mBPlusPDG = TDatabasePDG::Instance()->GetParticle(521)->Mass();
    
    // delta mass PDG
    Double_t deltaPDG = mBPlusPDG - mD0PDG;
   
    // Half width BPlus mass
    UInt_t prongs[2];
    prongs[0] = 211; prongs[1] = 421;
    Double_t invMassBPlus = candidateBPlus->InvMass(2,prongs);
    Double_t invMassDifference = TMath::Abs(mBPlusPDG - invMassBPlus);
    Double_t invMassDelta = TMath::Abs(deltaPDG-(DeltaInvMassBPlusKpipi(candidateBPlus)));

    Double_t pointingAngle = candidateBPlus->CosPointingAngle();
    Double_t dcaMother = candidateBPlus->GetDCA();
    Double_t ptMother = candidateBPlus->Pt();
    Double_t momentumMother = candidateBPlus->P();
    Double_t ptD0 = candidateD0->Pt();
    Double_t ptPion = candidatePion->Pt();

    AliExternalTrackParam motherTrack;
    motherTrack.CopyFromVTrack(candidateBPlus);
    Double_t d0z0[2],covd0z0[3],d0[2];
    motherTrack.PropagateToDCA(primaryVertex,bz,100.,d0z0,covd0z0);
    d0[0] = d0z0[0];
    Double_t d0Mother = TMath::Abs(d0[0]);
    Double_t d0firstTrack = TMath::Abs(candidateBPlus->Getd0Prong(0));
    Double_t d0secondTrack = TMath::Abs(candidateBPlus->Getd0Prong(1));  
    
    Double_t impactProduct = candidateBPlus->Prodd0d0();
    Double_t impactProductXY = TMath::Abs(candidateBPlus->ImpParXY());  

    Double_t angleBetweenBothDaughters  = (candidateD0->Px() * candidatePion->Px() + candidateD0->Py() * candidatePion->Py() + candidateD0->Pz() * candidatePion->Pz()) /(candidateD0->P() * candidatePion->P());
    Double_t angleMotherFirstDaughter = (candidateBPlus->Px() * candidatePion->Px() + candidateBPlus->Py() * candidatePion->Py() + candidateBPlus->Pz() * candidatePion->Pz()) /(candidateBPlus->P() * candidatePion->P());
    Double_t angleMotherSecondDaughter = (candidateBPlus->Px() * candidateD0->Px() + candidateBPlus->Py() * candidateD0->Py() + candidateBPlus->Pz() * candidateD0->Pz()) /(candidateBPlus->P() * candidateD0->P());

    Double_t cosThetaStar = candidateBPlus->CosThetaStar(0,521,211,421);
    Double_t vertexDistance = vertexBPlus->DistanceToVertex(primaryVertex);
    Double_t normDecayLength = candidateBPlus->NormalizedDecayLength();
    Double_t pdgMassMother = TDatabasePDG::Instance()->GetParticle(521)->Mass();
    Double_t pseudoProperDecayLength = ((vertexBPlus->GetX() - primaryVertex->GetX()) * candidateBPlus->Px() / TMath::Abs(candidateBPlus->Pt())) + ((vertexBPlus->GetY() - primaryVertex->GetY()) * candidateBPlus->Py() / TMath::Abs(candidateBPlus->Pt()));
    Double_t pseudoProperDecayTime = pseudoProperDecayLength * pdgMassMother/ptMother;
    Double_t decayTime = vertexDistance / (299792458 * TMath::Sqrt(1/((pdgMassMother*pdgMassMother/(momentumMother*momentumMother)) + 1)));

    Double_t phi = candidateBPlus->Phi();
    Double_t theta = candidateBPlus->Theta();
    Double_t covMatrix[21];
    candidateBPlus->GetCovarianceXYZPxPyPz(covMatrix);

    Double_t cp = TMath::Cos(phi);
    Double_t sp = TMath::Sin(phi);
    Double_t ct = TMath::Cos(theta);
    Double_t st = TMath::Sin(theta);

    Double_t errorMomentum = covMatrix[9]*cp*cp*ct*ct  // GetCovPxPx
                            +covMatrix[13]*2.*cp*sp*ct*ct  // GetCovPxPy
                            +covMatrix[18]*2.*cp*ct*st  // GetCovPxPz
                            +covMatrix[14]*sp*sp*ct*ct  // GetCovPyPy
                            +covMatrix[19]*2.*sp*ct*st  // GetCovPyPz
                            +covMatrix[20]*st*st;  // GetCovPzPz
    Double_t normalizedDecayTime = candidateBPlus->NormalizedDecayLength() / (299792458 * TMath::Sqrt(1/((pdgMassMother*pdgMassMother*errorMomentum*errorMomentum/(momentumMother*momentumMother)) + 1)));

    Double_t cosPointingAngleXY = candidateBPlus->CosPointingAngleXY();
    Double_t distanceXYToVertex = vertexBPlus->DistanceXYToVertex(primaryVertex);
    Double_t normalizedDecayLengthXY = candidateBPlus->NormalizedDecayLengthXY();
    Double_t chi2Vertex = vertexBPlus->GetChi2perNDF();

    //Topomatic
    Double_t dd0pr1=0.;
    Double_t dd0pr2=0.;
    Double_t dd0max=0.;
    Double_t dd0min=0.;
    for(Int_t ipr=0; ipr<2; ipr++) 
    {
      Double_t diffIP, errdiffIP;
      candidateBPlus->Getd0MeasMinusExpProng(ipr,bz,diffIP,errdiffIP);
      Double_t normdd0=0.;
      if(errdiffIP>0.) normdd0=diffIP/errdiffIP;
      if(ipr==0) dd0pr1=normdd0;
      if(ipr==1) dd0pr2=normdd0;
    }
    if(TMath::Abs(dd0pr1)>TMath::Abs(dd0pr2)) {dd0max=dd0pr1; dd0min=dd0pr2;}
    else {dd0max=dd0pr2; dd0min=dd0pr1;}

    Double_t Normalizedd0D0firstdaughter = TMath::Abs(candidateD0->Getd0Prong(0)/candidateD0->Getd0errProng(0));
    Double_t Normalizedd0D0seconddaughter = TMath::Abs(candidateD0->Getd0Prong(1)/candidateD0->Getd0errProng(1));
    Double_t Normalizedd0D0 = TMath::Abs(candidateBPlus->Getd0Prong(1)/candidateBPlus->Getd0errProng(1));
    Double_t Normalizedd0BPluspion = TMath::Abs(candidateBPlus->Getd0Prong(0)/candidateBPlus->Getd0errProng(0));
    Double_t Normalizedd0BPlus = TMath::Abs(d0[0]/covd0z0[0]);
    Double_t NormalizedimpactproductD0 = (candidateD0->Getd0Prong(0)/candidateD0->Getd0errProng(0)) * (candidateD0->Getd0Prong(1)/candidateD0->Getd0errProng(1));
    Double_t NormalizedimpactproductBPlus = (candidateBPlus->Getd0Prong(0)/candidateBPlus->Getd0errProng(0)) * (candidateBPlus->Getd0Prong(1)/candidateBPlus->Getd0errProng(1));

    // We apply the cuts 
    Int_t nCutIndex = 0;

    // "inv. mass width [GeV]" --------------------------------------------
    nCutIndex = 39;
    fCutVariableValueArray[nCutIndex] = invMassDifference;
    //---------------------------------------------------------------------

    // "delta mass width [GeV]" -------------------------------------------
    nCutIndex = 40;
    fCutVariableValueArray[nCutIndex] = invMassDelta;
    //---------------------------------------------------------------------

    // "pointing angle [Cos(theta)]" --------------------------------------
    nCutIndex = 41;
    fCutVariableValueArray[nCutIndex] = pointingAngle;
    //---------------------------------------------------------------------

    // "dca [cm]" ---------------------------------------------------------
    nCutIndex = 42;
    fCutVariableValueArray[nCutIndex] = dcaMother;
    //---------------------------------------------------------------------

    // "Pt BPlus [GeV/c]" ----------------------------------------------------
    nCutIndex = 43;
    fCutVariableValueArray[nCutIndex] = ptMother;
    //---------------------------------------------------------------------

    // "Pt D0 [GeV/c]" -------------------------------------------------
    nCutIndex = 44;
    fCutVariableValueArray[nCutIndex] = ptD0;
    //---------------------------------------------------------------------

    // "Pt Pion [GeV/c]" --------------------------------------------------
    nCutIndex = 45;
    fCutVariableValueArray[nCutIndex] = ptPion;
    //---------------------------------------------------------------------

    // "d0 BPlus [cm]" -------------------------------------------------------
    nCutIndex = 46;
    fCutVariableValueArray[nCutIndex] = d0Mother;
    //---------------------------------------------------------------------

    // "d0 D0 [cm]"-----------------------------------------------------
    nCutIndex = 47;
    fCutVariableValueArray[nCutIndex] = d0firstTrack;
    //---------------------------------------------------------------------

    // "d0 Pion [cm]" -----------------------------------------------------
    nCutIndex = 48;
    fCutVariableValueArray[nCutIndex] = d0secondTrack;
    //---------------------------------------------------------------------

    // "d0d0 [cm^2]" ------------------------------------------------------
    nCutIndex = 49;
    fCutVariableValueArray[nCutIndex] = impactProduct;
    //---------------------------------------------------------------------

    // "d0d0 XY [cm^2]" ---------------------------------------------------
    nCutIndex = 50;
    fCutVariableValueArray[nCutIndex] = impactProductXY;
    //---------------------------------------------------------------------

    // "angle between both daughters" -------------------------------------
    nCutIndex = 51;
    fCutVariableValueArray[nCutIndex] = angleBetweenBothDaughters;
    //---------------------------------------------------------------------

    // "angle mother with first daughter" ---------------------------------
    nCutIndex = 52;
    fCutVariableValueArray[nCutIndex] = angleMotherFirstDaughter;
    //---------------------------------------------------------------------

    // "angle mother with second daughter" --------------------------------
    nCutIndex = 53;
    fCutVariableValueArray[nCutIndex] = angleMotherSecondDaughter;
    //---------------------------------------------------------------------

    // "cosThetaStar" -----------------------------------------------------
    nCutIndex = 54;
    fCutVariableValueArray[nCutIndex] = cosThetaStar;
    //---------------------------------------------------------------------

    // "vertexDistance" ---------------------------------------------------
    nCutIndex = 55;
    fCutVariableValueArray[nCutIndex] = vertexDistance;
    //---------------------------------------------------------------------

    // "pseudoProperDecayTime" --------------------------------------------
    nCutIndex = 56;
    fCutVariableValueArray[nCutIndex] = pseudoProperDecayTime;
    //---------------------------------------------------------------------

    // "DecayTime" --------------------------------------------------------
    nCutIndex = 57;
    fCutVariableValueArray[nCutIndex] = decayTime;
    //---------------------------------------------------------------------

    // "normalizedDecayTime" ----------------------------------------------------
    nCutIndex = 58;
    fCutVariableValueArray[nCutIndex] = normalizedDecayTime;
    //---------------------------------------------------------------------

    // "normDecayLength" --------------------------------------------------
    nCutIndex = 59;
    fCutVariableValueArray[nCutIndex] = normDecayLength;
    //---------------------------------------------------------------------

    // "topomatic first daughter" -----------------------------------------
    nCutIndex = 60;
    fCutVariableValueArray[nCutIndex] = dd0pr1;
    //---------------------------------------------------------------------

    // "topomatic second daughter" ----------------------------------------
    nCutIndex = 61;
    fCutVariableValueArray[nCutIndex] = dd0pr2;
    //---------------------------------------------------------------------

    // "topomatic max" ----------------------------------------------------
    nCutIndex = 62;
    fCutVariableValueArray[nCutIndex] = dd0max;
    //---------------------------------------------------------------------

    // "topomatic min" ----------------------------------------------------
    nCutIndex = 63;
    fCutVariableValueArray[nCutIndex] = dd0min;
    //---------------------------------------------------------------------

    // "pointing angle XY" ----------------------------------------------------
    nCutIndex = 64;
    fCutVariableValueArray[nCutIndex] = cosPointingAngleXY;
    //---------------------------------------------------------------------

     // "vertex distance XY" ----------------------------------------------------
    nCutIndex = 65;
    fCutVariableValueArray[nCutIndex] = distanceXYToVertex;
    //---------------------------------------------------------------------
    
     // "normalized decay length XY" ----------------------------------------------------
    nCutIndex = 66;
    fCutVariableValueArray[nCutIndex] = normalizedDecayLengthXY;
    //---------------------------------------------------------------------
    
    // "chi squared per NDF" ----------------------------------------------------
    nCutIndex = 67;
    fCutVariableValueArray[nCutIndex] = chi2Vertex;

    // "Normalizedd0D0firstdaughter" ----------------------------------------------------
    nCutIndex = 68;
    fCutVariableValueArray[nCutIndex] = Normalizedd0D0firstdaughter;
    //---------------------------------------------------------------------

    // "Normalizedd0D0seconddaughter" ----------------------------------------------------
    nCutIndex = 69;
    fCutVariableValueArray[nCutIndex] = Normalizedd0D0seconddaughter;
    //---------------------------------------------------------------------

    // "Normalizedd0D0" ----------------------------------------------------
    nCutIndex = 70;
    fCutVariableValueArray[nCutIndex] = Normalizedd0D0;
    //---------------------------------------------------------------------

    // "Normalizedd0BPluspion" ----------------------------------------------------
    nCutIndex = 71;
    fCutVariableValueArray[nCutIndex] = Normalizedd0BPluspion;
    //---------------------------------------------------------------------

    // "Normalizedd0BPlus" ----------------------------------------------------
    nCutIndex = 72;
    fCutVariableValueArray[nCutIndex] = Normalizedd0BPlus;
    //---------------------------------------------------------------------

    // "NormalizedimpactproductD0" ----------------------------------------------------
    nCutIndex = 73;
    fCutVariableValueArray[nCutIndex] = NormalizedimpactproductD0;
    //---------------------------------------------------------------------

    // "NormalizedimpactproductBPlus" ----------------------------------------------------
    nCutIndex = 74;
    fCutVariableValueArray[nCutIndex] = NormalizedimpactproductBPlus;
    //---------------------------------------------------------------------
  }
  return;
}

