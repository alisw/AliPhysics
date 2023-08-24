

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
//   Base class for (B0 -> DPlus pi-> K pi pi pi) Analysis
//
//
//             Cuts are centralized in AliRDHFCutsB0toDPi
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
#include "AliRDHFCutsB0toDPi.h"
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
#include "AliAnalysisTaskSEB0toDPi.h"
#include "AliAODInputHandler.h"
#include <vector>
#include <TMatrix.h>
#include <TVector3.h>
#include <TArrayI.h>
#include <bitset>
#include <TH3F.h>

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskSEB0toDPi);
/// \endcond

//__________________________________________________________________________
AliAnalysisTaskSEB0toDPi::AliAnalysisTaskSEB0toDPi():
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
  fUseSideBands(0),
  fSideBandLow(0),
  fSideBandHigh(0),
  fSaveTRHists(0),
  fOutput(0),
  fListCuts(0),
  fOutputB0Results(0),
  fOutputDPlusKaon(0),
  fOutputDPlusPions(0),
  fOutputB0Pion(0),
  fOutputDPlus(0),
  fOutputB0(0),
  fOutputDPlus_DPlusPt(0),
  fCuts(0),
  fCEvents(0),
  fCTrigger(0),
  fCRejected(0),
  fTriggerClassesCorrelated(0),
  fTriggerSubClasses(0),
  fB0PionTracks(0x0),
  fDPlusTracks(0x0),
  // fDPlusPionTracks(0x0),
  fnPtBins(0),
  fnPtBinLimits(0),
  fPtBinLimits(0x0),
  fnPtBinsDPlusforDPlusptbin(0),
  fnPtBinsDPlusforDPlusptbinLimits(0),
  fPtBinLimitsDPlusforDPlusptbin(0x0),
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
AliAnalysisTaskSEB0toDPi::AliAnalysisTaskSEB0toDPi(const Char_t* name, AliRDHFCutsB0toDPi* cuts) :
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
  fUseSideBands(0),
  fSideBandLow(0),
  fSideBandHigh(0),
  fSaveTRHists(0),
  fOutput(0),
  fListCuts(0),
  fOutputB0Results(0),
  fOutputDPlusKaon(0),
  fOutputDPlusPions(0),
  fOutputB0Pion(0),
  fOutputDPlus(0),
  fOutputB0(0),
  fOutputDPlus_DPlusPt(0),
  fCuts(0),
  fCEvents(0),
  fCTrigger(0),
  fCRejected(0),
  fTriggerClassesCorrelated(0),
  fTriggerSubClasses(0),
  fB0PionTracks(0x0),
  fDPlusTracks(0x0),
  // fDPlusPionTracks(0x0),
  fnPtBins(0),
  fnPtBinLimits(0),
  fPtBinLimits(0x0),
  fnPtBinsDPlusforDPlusptbin(0),
  fnPtBinsDPlusforDPlusptbinLimits(0),
  fPtBinLimitsDPlusforDPlusptbin(0x0),
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

  Info("AliAnalysisTaskSEB0toDPi", "Calling Constructor");

  fCuts = cuts;

  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());  // General information
  DefineOutput(2, TList::Class());  // Cuts
  DefineOutput(3, TList::Class());  // DPlus kaon output
  DefineOutput(4, TList::Class());  // DPlus pions output
  DefineOutput(5, TList::Class());  // B0 pion output
  DefineOutput(6, TList::Class());  // DPlus output
  DefineOutput(7, TList::Class());  // B0 output
  DefineOutput(8, TList::Class());  // DPlus DPlusPt output
  DefineOutput(9, TList::Class());  // B0 results output

}

//___________________________________________________________________________
AliAnalysisTaskSEB0toDPi::~AliAnalysisTaskSEB0toDPi() {
  //
  /// destructor
  //
  Info("~AliAnalysisTaskSEB0toDPi", "Calling Destructor");

  delete fOutput;
  delete fOutputDPlusKaon;
  delete fOutputDPlusPions;
  delete fOutputB0Pion;
  delete fOutputDPlus;
  delete fOutputB0;
  delete fOutputDPlus_DPlusPt;
  delete fOutputB0Results;
  delete fCuts;
  delete fCEvents;
  delete fCTrigger;
  delete fCRejected;
  delete fTriggerClassesCorrelated,
  delete fTriggerSubClasses,
  delete fDPlusTracks;
  // delete fDPlusPionTracks;
  delete fB0PionTracks;
  delete fListCuts;
  delete fPtBinLimits;
  delete fPtBinLimitsDPlusforDPlusptbin;
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
void AliAnalysisTaskSEB0toDPi::Init() {
  //
  /// Initialization
  //

  if (fDebug > 1) printf("AliAnalysisTaskSEB0toDPi::Init() \n");

  return;
}
//_________________________________________________
void AliAnalysisTaskSEB0toDPi::UserExec(Option_t *) {

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

  AliInputEventHandler *fInputHandler = (AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  //save trigger information
  for (int i = 0; i < 36; ++i)
  {
    ULong64_t compareValue = BIT(i);
    // std::cout << "Compare mask: " << std::bitset<32>(compareValue) << std::endl;
    Bool_t trigger1 = kFALSE;
    Bool_t trigger2 = kFALSE;

    if (fInputHandler->IsEventSelected() & compareValue)
    {
      // std::cout << "Match!" << std::endl;
      fCTrigger->Fill(i + 1);
      if (i == 1) trigger1 = kTRUE;
      if (i == 10) trigger2 = kTRUE;
    }
    if (trigger1 && trigger2) fCTrigger->Fill(35);
  }


  //==================================================================================
  //  EVENT INITIALIZATION - start
  //==================================================================================

  TClonesArray *array3Prong = 0;
  AliAODEvent* aodEvent = dynamic_cast<AliAODEvent*>(fInputEvent);

  fCEvents->Fill(1);

  if(!aodEvent && AODEvent() && IsStandardAOD()) {
    // In case there is an AOD handler writing a standard AOD, use the AOD 
    // event in memory rather than the input (ESD) event.    
    aodEvent = dynamic_cast<AliAODEvent*> (AODEvent());
    // in this case the braches in the deltaAOD (AliAOD.VertexingHF.root)
    // have to taken from the AOD event hold by the AliAODExtension
    AliAODHandler* aodHandler = (AliAODHandler*) 
      ((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
    if(aodHandler->GetExtensions()) {
      AliAODExtension *ext = (AliAODExtension*)aodHandler->GetExtensions()->FindObject("AliAOD.VertexingHF.root");
      AliAODEvent *aodFromExt = ext->GetAOD();
      array3Prong=(TClonesArray*)aodFromExt->GetList()->FindObject("Charm3Prong");
    }
  } else {
    array3Prong=(TClonesArray*)aodEvent->GetList()->FindObject("Charm3Prong");
  }
  

  // fix for temporary bug in ESDfilter
  // the AODs with null vertex pointer didn't pass the PhysSel
  if (!aodEvent->GetPrimaryVertex() || TMath::Abs(aodEvent->GetMagneticField()) < 0.001) return;
  fCEvents->Fill(2);

  // trigger class for PbPb C0SMH-B-NOPF-ALLNOTRD
  TString trigclass = aodEvent->GetFiredTriggerClasses();
  if (trigclass.Contains("C0SMH-B-NOPF-ALLNOTRD") || trigclass.Contains("C0SMH-B-NOPF-ALL")) fCEvents->Fill(5);

  //Save EMCAL trigger information
  // TString firedTrigClass = aodEvent->GetFiredTriggerClasses();
  
  // if (fInputHandler->IsEventSelected() & AliVEvent::kMB)fTriggerClassesCorrelated->Fill(0);
  // if (fInputHandler->IsEventSelected() & AliVEvent::kINT7)fTriggerClassesCorrelated->Fill(1);
  // if (fInputHandler->IsEventSelected() & AliVEvent::kEMC1)fTriggerClassesCorrelated->Fill(2);
  // if ((fInputHandler->IsEventSelected() & AliVEvent::kEMC7) && firedTrigClass.Contains("EMC"))fTriggerClassesCorrelated->Fill(3);
  // if (firedTrigClass.Contains("7EJE") || firedTrigClass.Contains("8EJE")) fTriggerClassesCorrelated->Fill(4);
  // if (firedTrigClass.Contains("7EJ1") || firedTrigClass.Contains("8EJ1")) fTriggerClassesCorrelated->Fill(5);
  // if (firedTrigClass.Contains("7EJ2") || firedTrigClass.Contains("8EJ2")) fTriggerClassesCorrelated->Fill(6);
  // if (firedTrigClass.Contains("7EGA") || firedTrigClass.Contains("8EGA")) fTriggerClassesCorrelated->Fill(7);
  // if (firedTrigClass.Contains("7EG1") || firedTrigClass.Contains("8EG1")) fTriggerClassesCorrelated->Fill(8);
  // if (firedTrigClass.Contains("7EG2") || firedTrigClass.Contains("8EG2")) fTriggerClassesCorrelated->Fill(9);
  // if ((fInputHandler->IsEventSelected() & AliVEvent::kEMC7) && firedTrigClass.Contains("DMC"))fTriggerClassesCorrelated->Fill(10);
  // if (firedTrigClass.Contains("7DJE") || firedTrigClass.Contains("8DJE")) fTriggerClassesCorrelated->Fill(11);
  // if (firedTrigClass.Contains("7DJ1") || firedTrigClass.Contains("8DJ1")) fTriggerClassesCorrelated->Fill(12);
  // if (firedTrigClass.Contains("7DJ2") || firedTrigClass.Contains("8DJ2")) fTriggerClassesCorrelated->Fill(13);
  // if (firedTrigClass.Contains("7DGA") || firedTrigClass.Contains("8DGA")) fTriggerClassesCorrelated->Fill(14);
  // if (firedTrigClass.Contains("7DG1") || firedTrigClass.Contains("8DG1")) fTriggerClassesCorrelated->Fill(15);
  // if (firedTrigClass.Contains("7DG2") || firedTrigClass.Contains("8DG2")) fTriggerClassesCorrelated->Fill(16);

  // if ((fInputHandler->IsEventSelected() & AliVEvent::kEMCEJE) && firedTrigClass.Contains("7EJE")) fTriggerSubClasses->Fill(0);
  // if ((fInputHandler->IsEventSelected() & AliVEvent::kEMCEJE) && firedTrigClass.Contains("7EJ1")) fTriggerSubClasses->Fill(1);
  // if ((fInputHandler->IsEventSelected() & AliVEvent::kEMCEJE) && firedTrigClass.Contains("7EJ2")) fTriggerSubClasses->Fill(2);
  // if ((fInputHandler->IsEventSelected() & AliVEvent::kEMCEGA) && firedTrigClass.Contains("7EGA")) fTriggerSubClasses->Fill(3);
  // if ((fInputHandler->IsEventSelected() & AliVEvent::kEMCEGA) && firedTrigClass.Contains("7EG1")) fTriggerSubClasses->Fill(4);
  // if ((fInputHandler->IsEventSelected() & AliVEvent::kEMCEGA) && firedTrigClass.Contains("7EG2")) fTriggerSubClasses->Fill(5);
  // if ((fInputHandler->IsEventSelected() & AliVEvent::kEMCEJE) && firedTrigClass.Contains("7DJE")) fTriggerSubClasses->Fill(6);
  // if ((fInputHandler->IsEventSelected() & AliVEvent::kEMCEJE) && firedTrigClass.Contains("7DJ1")) fTriggerSubClasses->Fill(7);
  // if ((fInputHandler->IsEventSelected() & AliVEvent::kEMCEJE) && firedTrigClass.Contains("7DJ2")) fTriggerSubClasses->Fill(8);
  // if ((fInputHandler->IsEventSelected() & AliVEvent::kEMCEGA) && firedTrigClass.Contains("7DGA")) fTriggerSubClasses->Fill(9);
  // if ((fInputHandler->IsEventSelected() & AliVEvent::kEMCEGA) && firedTrigClass.Contains("7DG1")) fTriggerSubClasses->Fill(10);
  // if ((fInputHandler->IsEventSelected() & AliVEvent::kEMCEGA) && firedTrigClass.Contains("7DG2")) fTriggerSubClasses->Fill(11);
  // if ((fInputHandler->IsEventSelected() & AliVEvent::kEMCEJE) && firedTrigClass.Contains("8EJE")) fTriggerSubClasses->Fill(12);
  // if ((fInputHandler->IsEventSelected() & AliVEvent::kEMCEJE) && firedTrigClass.Contains("8EJ1")) fTriggerSubClasses->Fill(13);
  // if ((fInputHandler->IsEventSelected() & AliVEvent::kEMCEJE) && firedTrigClass.Contains("8EJ2")) fTriggerSubClasses->Fill(14);
  // if ((fInputHandler->IsEventSelected() & AliVEvent::kEMCEGA) && firedTrigClass.Contains("8EGA")) fTriggerSubClasses->Fill(15);
  // if ((fInputHandler->IsEventSelected() & AliVEvent::kEMCEGA) && firedTrigClass.Contains("8EG1")) fTriggerSubClasses->Fill(16);
  // if ((fInputHandler->IsEventSelected() & AliVEvent::kEMCEGA) && firedTrigClass.Contains("8EG2")) fTriggerSubClasses->Fill(17);
  // if ((fInputHandler->IsEventSelected() & AliVEvent::kEMCEJE) && firedTrigClass.Contains("8DJE")) fTriggerSubClasses->Fill(18);
  // if ((fInputHandler->IsEventSelected() & AliVEvent::kEMCEJE) && firedTrigClass.Contains("8DJ1")) fTriggerSubClasses->Fill(19);
  // if ((fInputHandler->IsEventSelected() & AliVEvent::kEMCEJE) && firedTrigClass.Contains("8DJ2")) fTriggerSubClasses->Fill(20);
  // if ((fInputHandler->IsEventSelected() & AliVEvent::kEMCEGA) && firedTrigClass.Contains("8DGA")) fTriggerSubClasses->Fill(21);
  // if ((fInputHandler->IsEventSelected() & AliVEvent::kEMCEGA) && firedTrigClass.Contains("8DG1")) fTriggerSubClasses->Fill(22);
  // if ((fInputHandler->IsEventSelected() & AliVEvent::kEMCEGA) && firedTrigClass.Contains("8DG2")) fTriggerSubClasses->Fill(23);

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

  if (!array3Prong){
    AliInfo("Could not find array of Charm3Prong, skipping the event");
    return;
  }
  else AliDebug(2, Form("Found %d vertices",array3Prong->GetEntriesFast())); 
  if (array3Prong->GetEntriesFast()==0) return;

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

  //Here we create an array and vectors that we use to look how many B0s can be reconstructed after each cut
  TMatrix * B0toDPiLabelMatrix = new TMatrix(0, 6);

  //we fill the array with all B0->DPlusPi tracks
  if (fUseMCInfo) {
    B0toDPiSignalTracksInMC(mcTrackArray, aodEvent, B0toDPiLabelMatrix, fOutputB0Results);
  }

  //==================================================================================
  //  GET TRUE MC INFO - end
  //==================================================================================
  //  PARTICLE SELECTION LOOP - start
  //==================================================================================
  //
  // Here we select and reconstruct the particles for the B0 -> DPlus + Pion decay.
  //
  //==================================================================================


  // TClonesArray * DPlusTrackArray = new TClonesArray("AliAODRecoDecayHF3Prong", 5000);;
  // DaughterSelection(aodEvent, primaryVertex, bz, mcTrackArray, B0toDPiLabelMatrix, mcHeader, 0, fDPlusKaonTracks);
  // DaughterSelection(aodEvent, primaryVertex, bz, mcTrackArray, B0toDPiLabelMatrix, mcHeader, 1, fDPlusPionTracks);

  DPlusSelection(aodEvent, primaryVertex, bz, mcTrackArray, B0toDPiLabelMatrix, array3Prong, mcHeader);

  if ((Int_t)fDPlusTracks->size() > 0)
  {
    DaughterSelection(aodEvent, primaryVertex, bz, mcTrackArray, B0toDPiLabelMatrix, mcHeader, 2, fB0PionTracks);

    if ((Int_t)fB0PionTracks->size() > 0)
    {
      B0Selection(aodEvent, primaryVertex, bz, mcTrackArray, B0toDPiLabelMatrix, array3Prong, mcHeader);
    }
  }
  

  // for all entries of array get vertex and delete.
  // for (int i = 0; i < DPlusTrackArray->GetEntriesFast(); ++i) {
  //   AliAODVertex * vertex = (AliAODVertex*)((AliAODRecoDecayHF3Prong*)(*DPlusTrackArray)[i])->GetOwnSecondaryVtx();
  //   delete vertex; vertex = NULL;
  // }
  // DPlusTrackArray->Clear("C");
  // delete DPlusTrackArray;

  // Clear arrays and memory management:
  fB0PionTracks->erase(fB0PionTracks->begin(), fB0PionTracks->end());
  fDPlusTracks->erase(fDPlusTracks->begin(), fDPlusTracks->end());
  // fDPlusKaonTracks->erase(fDPlusKaonTracks->begin(), fDPlusKaonTracks->end());
  // fDPlusPionTracks->erase(fDPlusPionTracks->begin(), fDPlusPionTracks->end());



  delete B0toDPiLabelMatrix; B0toDPiLabelMatrix = NULL;

  //==================================================================================
  //  PARTICLE SELECTION LOOP - end
  //==================================================================================

  PostData(1, fOutput);
  PostData(3, fOutputDPlusKaon);
  PostData(4, fOutputDPlusPions);
  PostData(5, fOutputB0Pion);
  PostData(6, fOutputDPlus);
  PostData(7, fOutputB0);
  PostData(8, fOutputDPlus_DPlusPt);
  PostData(9, fOutputB0Results);


  //==================================================================================
  //  USER EXECUTION FUNCTION - end
  //==================================================================================

}
//___________________________________ terminate ___________________________
void AliAnalysisTaskSEB0toDPi::Terminate(Option_t*) {
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

  fCEvents                       = dynamic_cast<TH1F*>(fOutput->FindObject("fCEvents"));
  fCTrigger                      = dynamic_cast<TH1F*>(fOutput->FindObject("fCTrigger"));
  fCRejected                     = dynamic_cast<TH1F*>(fOutput->FindObject("fCRejected"));
  fTriggerClassesCorrelated      = dynamic_cast<TH1F*>(fOutput->FindObject("fTriggerClassesCorrelated"));
  fTriggerSubClasses             = dynamic_cast<TH1F*>(fOutput->FindObject("fTriggerSubClasses"));

  fListCuts = dynamic_cast<TList*> (GetOutputData(2));
  if (!fListCuts) {
    printf("ERROR: fListCuts not available\n");
    return;
  }
  fOutputDPlusKaon = dynamic_cast<TList*> (GetOutputData(3));
  if (!fOutputDPlusKaon) {
    printf("ERROR: fOutputDPlusKaon not available\n");
    return;
  }
  fOutputDPlusPions = dynamic_cast<TList*> (GetOutputData(4));
  if (!fOutputDPlusPions) {
    printf("ERROR: fOutputDPlusPions not available\n");
    return;
  }
  fOutputB0Pion = dynamic_cast<TList*> (GetOutputData(5));
  if (!fOutputB0Pion) {
    printf("ERROR: fOutputB0Pion not available\n");
    return;
  }
  fOutputDPlus = dynamic_cast<TList*> (GetOutputData(6));
  if (!fOutputDPlus) {
    printf("ERROR: fOutputDPlus not available\n");
    return;
  }
  fOutputB0 = dynamic_cast<TList*> (GetOutputData(7));
  if (!fOutputB0) {
    printf("ERROR: fOutputB0 not available\n");
    return;
  }
  fOutputDPlus_DPlusPt = dynamic_cast<TList*> (GetOutputData(8));
  if (!fOutputDPlus_DPlusPt) {
    printf("ERROR: fOutputDPlus_DPlusPt not available\n");
    return;
  }
  fOutputB0Results = dynamic_cast<TList*> (GetOutputData(9));
  if (!fOutputB0Results) {
    printf("ERROR: fOutputB0Results not available\n");
    return;
  }
  return;
}

//___________________________________________________________________________
void AliAnalysisTaskSEB0toDPi::UserCreateOutputObjects() {
  /// output
  Info("UserCreateOutputObjects", "CreateOutputObjects of task %s\n", GetName());

  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("General_information");

  fOutputDPlusKaon = new TList();
  fOutputDPlusKaon->SetOwner();
  fOutputDPlusKaon->SetName("listDPlusKaon");

  fOutputDPlusPions = new TList();
  fOutputDPlusPions->SetOwner();
  fOutputDPlusPions->SetName("listDPlusPions");

  fOutputB0Pion = new TList();
  fOutputB0Pion->SetOwner();
  fOutputB0Pion->SetName("listB0Pion");

  fOutputDPlus = new TList();
  fOutputDPlus->SetOwner();
  fOutputDPlus->SetName("listDPlus");

  fOutputB0 = new TList();
  fOutputB0->SetOwner();
  fOutputB0->SetName("listB0");

  fOutputDPlus_DPlusPt = new TList();
  fOutputDPlus_DPlusPt->SetOwner();
  fOutputDPlus_DPlusPt->SetName("listDPlus_DPlusPt");

  fOutputB0Results = new TList();
  fOutputB0Results->SetOwner();
  fOutputB0Results->SetName("listB0Results");

  // we prepare vectors that will save the positions of the daughter tracks in the track list during the reconstruction
  fB0PionTracks = new std::vector<Int_t>;
  // fDPlusPionTracks = new std::vector<Int_t>;
  fDPlusTracks = new std::vector<Int_t>;

  // we get information on the pt bins
  fnPtBins = fCuts->GetNPtBins();
  fnPtBinLimits = fnPtBins + 1;
  fPtBinLimits = fCuts->GetPtBinLimits();

  fnPtBinsDPlusforDPlusptbin = fCuts->GetNPtBinsDPlusforDPlusptbin();
  fnPtBinsDPlusforDPlusptbinLimits = fnPtBinsDPlusforDPlusptbin + 1;
  fPtBinLimitsDPlusforDPlusptbin = fCuts->GetPtBinLimitsDPlusforDPlusptbin();

  std::cout << "Nr. of B0 meson bins: " <<  fCuts->GetNPtBins() << " limits: " << std::endl;
  for (int i = 0; i < fnPtBinLimits; ++i)
  {
    std::cout << fPtBinLimits[i] << " " << std::endl;
  }
  std::cout << std::endl;
  std::cout << "Nr. of DPlus meson bins: " <<  fCuts->GetNPtBinsDPlusforDPlusptbin() << " limits: " << std::endl;
  for (int i = 0; i < fnPtBinsDPlusforDPlusptbinLimits; ++i)
  {
    std::cout << fPtBinLimitsDPlusforDPlusptbin[i] << " " << std::endl;
  }
  std::cout << std::endl;

  fListCuts = new TList();
  fListCuts->SetOwner();
  fListCuts->SetName("Cuts");
  AliRDHFCutsB0toDPi* copyfCuts = new AliRDHFCutsB0toDPi(*fCuts);
  fListCuts->Add(copyfCuts);

  // define histograms
  DefineHistograms();

  // Post the data
  PostData(1, fOutput);
  PostData(2, fListCuts);
  PostData(3, fOutputDPlusKaon);
  PostData(4, fOutputDPlusPions);
  PostData(5, fOutputB0Pion);
  PostData(6, fOutputDPlus);
  PostData(7, fOutputB0);
  PostData(8, fOutputDPlus_DPlusPt);
  PostData(9, fOutputB0Results);

  return;
}
//___________________________________ histograms _______________________________________
void  AliAnalysisTaskSEB0toDPi::DefineHistograms() {

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
  fCEvents->GetXaxis()->SetBinLabel(12, "nr. of DPlus fail to be rec");
  fOutput->Add(fCEvents);

  fCTrigger = new TH1F("fCTrigger", "Trigger counter", 36, 0, 36);
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

  fCRejected = new TH1F("fCRejected", "Rejected counter", 35, 0, 35);
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

  fTriggerClassesCorrelated = new TH1F("TriggerCorrelations", "Triggers Correlated with EMCal triggers", 17, -0.5, 16.5);
  fTriggerClassesCorrelated->GetXaxis()->SetBinLabel( 1, "kMB");
  fTriggerClassesCorrelated->GetXaxis()->SetBinLabel( 2, "kINT7");
  fTriggerClassesCorrelated->GetXaxis()->SetBinLabel( 3, "kEMC1");
  fTriggerClassesCorrelated->GetXaxis()->SetBinLabel( 4, "kEMC7");
  fTriggerClassesCorrelated->GetXaxis()->SetBinLabel( 5, "kEMCEJE");
  fTriggerClassesCorrelated->GetXaxis()->SetBinLabel( 6, "kEMCEJ1");
  fTriggerClassesCorrelated->GetXaxis()->SetBinLabel( 7, "kEMCEJ2");
  fTriggerClassesCorrelated->GetXaxis()->SetBinLabel( 8, "kEMCEGA");
  fTriggerClassesCorrelated->GetXaxis()->SetBinLabel( 9, "kEMCEG1");
  fTriggerClassesCorrelated->GetXaxis()->SetBinLabel( 10, "kEMCEG2");
  fTriggerClassesCorrelated->GetXaxis()->SetBinLabel( 11, "kDMC7");
  fTriggerClassesCorrelated->GetXaxis()->SetBinLabel( 12, "kDMCDJE");
  fTriggerClassesCorrelated->GetXaxis()->SetBinLabel( 13, "kDMCDJ1");
  fTriggerClassesCorrelated->GetXaxis()->SetBinLabel( 14, "kDMCDJ2");
  fTriggerClassesCorrelated->GetXaxis()->SetBinLabel( 15, "kDMCDGA");
  fTriggerClassesCorrelated->GetXaxis()->SetBinLabel( 16, "kDMCDG1");
  fTriggerClassesCorrelated->GetXaxis()->SetBinLabel( 17, "kDMCDG2");
  fOutput->Add(fTriggerClassesCorrelated);

  fTriggerSubClasses = new TH1F("Subtriggers", "Subtriggers Correlated with EMCal triggers", 24, -0.5, 23.5);
  fTriggerSubClasses->GetXaxis()->SetBinLabel( 1, "kEMC_7EJE");
  fTriggerSubClasses->GetXaxis()->SetBinLabel( 2, "kEMC_7EJ1");
  fTriggerSubClasses->GetXaxis()->SetBinLabel( 3, "kEMC_7EJ2");
  fTriggerSubClasses->GetXaxis()->SetBinLabel( 4, "kEMC_7EGA");
  fTriggerSubClasses->GetXaxis()->SetBinLabel( 5, "kEMC_7EG1");
  fTriggerSubClasses->GetXaxis()->SetBinLabel( 6, "kEMC_7EG2");
  fTriggerSubClasses->GetXaxis()->SetBinLabel( 7, "kDMC_7DJE");
  fTriggerSubClasses->GetXaxis()->SetBinLabel( 8, "kDMC_7DJ1");
  fTriggerSubClasses->GetXaxis()->SetBinLabel( 9, "kDMC_7DJ2");
  fTriggerSubClasses->GetXaxis()->SetBinLabel( 10, "kDMC_7DGA");
  fTriggerSubClasses->GetXaxis()->SetBinLabel( 11, "kDMC_7DG1");
  fTriggerSubClasses->GetXaxis()->SetBinLabel( 12, "kDMC_7DG2");
  fTriggerSubClasses->GetXaxis()->SetBinLabel( 13, "kEMC_8EJE");
  fTriggerSubClasses->GetXaxis()->SetBinLabel( 14, "kEMC_8EJ1");
  fTriggerSubClasses->GetXaxis()->SetBinLabel( 15, "kEMC_8EJ2");
  fTriggerSubClasses->GetXaxis()->SetBinLabel( 16, "kEMC_8EGA");
  fTriggerSubClasses->GetXaxis()->SetBinLabel( 17, "kEMC_8EG1");
  fTriggerSubClasses->GetXaxis()->SetBinLabel( 18, "kEMC_8EG2");
  fTriggerSubClasses->GetXaxis()->SetBinLabel( 19, "kDMC_8DJE");
  fTriggerSubClasses->GetXaxis()->SetBinLabel( 20, "kDMC_8DJ1");
  fTriggerSubClasses->GetXaxis()->SetBinLabel( 21, "kDMC_8DJ2");
  fTriggerSubClasses->GetXaxis()->SetBinLabel( 22, "kDMC_8DGA");
  fTriggerSubClasses->GetXaxis()->SetBinLabel( 23, "kDMC_8DG1");
  fTriggerSubClasses->GetXaxis()->SetBinLabel( 24, "kDMC_8DG2");
  fOutput->Add(fTriggerSubClasses);


  //====================================================

  TString name_mc_B0_pt = "mc_B0_pt";
  TH1F* hist_mc_B0_pt = new TH1F(name_mc_B0_pt.Data(), "Pt monte carlo B0 in B0->D*#pi; #it{p}_{T} (GeV/#it{c}); Entries", 400, 0, 20);
  hist_mc_B0_pt->Sumw2();
  hist_mc_B0_pt->SetLineColor(6);
  hist_mc_B0_pt->SetMarkerStyle(20);
  hist_mc_B0_pt->SetMarkerSize(0.6);
  hist_mc_B0_pt->SetMarkerColor(6);
  TH1F* histogram_mc_B0_pt = (TH1F*)hist_mc_B0_pt->Clone();
  fOutputB0Results->Add(histogram_mc_B0_pt);
  fResultsHistogramArray[0][0] = histogram_mc_B0_pt;

  TString name_mc_B0_pion_pt = "mc_B0_pion_pt";
  TH1F* hist_mc_B0_pion_pt = new TH1F(name_mc_B0_pion_pt.Data(), "Pt monte carlo first pion of B0 in B0->D*#pi; #it{p}_{T} (GeV/#it{c}); Entries", 400, 0, 20);
  hist_mc_B0_pion_pt->Sumw2();
  hist_mc_B0_pion_pt->SetLineColor(6);
  hist_mc_B0_pion_pt->SetMarkerStyle(20);
  hist_mc_B0_pion_pt->SetMarkerSize(0.6);
  hist_mc_B0_pion_pt->SetMarkerColor(6);
  TH1F* histogram_mc_B0_pion_pt = (TH1F*)hist_mc_B0_pion_pt->Clone();
  fOutputB0Results->Add(histogram_mc_B0_pion_pt);
  fResultsHistogramArray[0][1] = histogram_mc_B0_pion_pt;

  TString name_mc_DPlus_pt = "mc_DPlus_pt";
  TH1F* hist_mc_DPlus_pt = new TH1F(name_mc_DPlus_pt.Data(), "Pt monte carlo DPlus in B0->D*#pi; #it{p}_{T} (GeV/#it{c}); Entries", 400, 0, 20);
  hist_mc_DPlus_pt->Sumw2();
  hist_mc_DPlus_pt->SetLineColor(6);
  hist_mc_DPlus_pt->SetMarkerStyle(20);
  hist_mc_DPlus_pt->SetMarkerSize(0.6);
  hist_mc_DPlus_pt->SetMarkerColor(6);
  TH1F* histogram_mc_DPlus_pt = (TH1F*)hist_mc_DPlus_pt->Clone();
  fOutputB0Results->Add(histogram_mc_DPlus_pt);
  fResultsHistogramArray[0][2] = histogram_mc_DPlus_pt;

  TString name_mc_DPlus_first_pion_pt = "mc_DPlus_first_pion_pt";
  TH1F* hist_mc_DPlus_first_pion_pt = new TH1F(name_mc_DPlus_first_pion_pt.Data(), "Pt monte carlo first pion of DPlus in B0->D*#pi; #it{p}_{T} (GeV/#it{c}); Entries", 400, 0, 20);
  hist_mc_DPlus_first_pion_pt->Sumw2();
  hist_mc_DPlus_first_pion_pt->SetLineColor(6);
  hist_mc_DPlus_first_pion_pt->SetMarkerStyle(20);
  hist_mc_DPlus_first_pion_pt->SetMarkerSize(0.6);
  hist_mc_DPlus_first_pion_pt->SetMarkerColor(6);
  TH1F* histogram_mc_DPlus_first_pion_pt = (TH1F*)hist_mc_DPlus_first_pion_pt->Clone();
  fOutputB0Results->Add(histogram_mc_DPlus_first_pion_pt);
  fResultsHistogramArray[0][3] = histogram_mc_DPlus_first_pion_pt;

  TString name_mc_DPlus_second_pion_pt = "mc_DPlus_second_pion_pt";
  TH1F* hist_mc_DPlus_second_pion_pt = new TH1F(name_mc_DPlus_second_pion_pt.Data(), "Pt monte carlo second pion of DPlus in B0->D*#pi; #it{p}_{T} (GeV/#it{c}); Entries", 400, 0, 20);
  hist_mc_DPlus_second_pion_pt->Sumw2();
  hist_mc_DPlus_second_pion_pt->SetLineColor(6);
  hist_mc_DPlus_second_pion_pt->SetMarkerStyle(20);
  hist_mc_DPlus_second_pion_pt->SetMarkerSize(0.6);
  hist_mc_DPlus_second_pion_pt->SetMarkerColor(6);
  TH1F* histogram_mc_DPlus_second_pion_pt = (TH1F*)hist_mc_DPlus_second_pion_pt->Clone();
  fOutputB0Results->Add(histogram_mc_DPlus_second_pion_pt);
  fResultsHistogramArray[0][4] = histogram_mc_DPlus_second_pion_pt;

  TString name_mc_DPlus_kaon_pt = "mc_DPlus_kaon_pt";
  TH1F* hist_mc_DPlus_kaon_pt = new TH1F(name_mc_DPlus_kaon_pt.Data(), "Pt monte carlo kaon of DPlus in B0->D*#pi; #it{p}_{T} (GeV/#it{c}); Entries", 400, 0, 20);
  hist_mc_DPlus_kaon_pt->Sumw2();
  hist_mc_DPlus_kaon_pt->SetLineColor(6);
  hist_mc_DPlus_kaon_pt->SetMarkerStyle(20);
  hist_mc_DPlus_kaon_pt->SetMarkerSize(0.6);
  hist_mc_DPlus_kaon_pt->SetMarkerColor(6);
  TH1F* histogram_mc_DPlus_kaon_pt = (TH1F*)hist_mc_DPlus_kaon_pt->Clone();
  fOutputB0Results->Add(histogram_mc_DPlus_kaon_pt);
  fResultsHistogramArray[0][5] = histogram_mc_DPlus_kaon_pt;

  TString name_mc_B0_rapidity_true = "mc_B0_rapidity_true";
  TH1F* hist_mc_B0_rapidity_true = new TH1F(name_mc_B0_rapidity_true.Data(), "rapidity_true monte carlo B0 in B0->D*#pi; Y; Entries", 5000, -20, 20);
  hist_mc_B0_rapidity_true->Sumw2();
  hist_mc_B0_rapidity_true->SetLineColor(6);
  hist_mc_B0_rapidity_true->SetMarkerStyle(20);
  hist_mc_B0_rapidity_true->SetMarkerSize(0.6);
  hist_mc_B0_rapidity_true->SetMarkerColor(6);
  TH1F* histogram_mc_B0_rapidity_true = (TH1F*)hist_mc_B0_rapidity_true->Clone();
  fOutputB0Results->Add(histogram_mc_B0_rapidity_true);
  fResultsHistogramArray[0][7] = histogram_mc_B0_rapidity_true;

  TString name_mc_B0_pion_rapidity_true = "mc_B0_pion_rapidity_true";
  TH1F* hist_mc_B0_pion_rapidity_true = new TH1F(name_mc_B0_pion_rapidity_true.Data(), "rapidity_true monte carlo first pion of B0 in B0->D*#pi; Y; Entries", 5000, -20, 20);
  hist_mc_B0_pion_rapidity_true->Sumw2();
  hist_mc_B0_pion_rapidity_true->SetLineColor(6);
  hist_mc_B0_pion_rapidity_true->SetMarkerStyle(20);
  hist_mc_B0_pion_rapidity_true->SetMarkerSize(0.6);
  hist_mc_B0_pion_rapidity_true->SetMarkerColor(6);
  TH1F* histogram_mc_B0_pion_rapidity_true = (TH1F*)hist_mc_B0_pion_rapidity_true->Clone();
  fOutputB0Results->Add(histogram_mc_B0_pion_rapidity_true);
  fResultsHistogramArray[0][8] = histogram_mc_B0_pion_rapidity_true;

  TString name_mc_DPlus_rapidity_true = "mc_DPlus_rapidity_true";
  TH1F* hist_mc_DPlus_rapidity_true = new TH1F(name_mc_DPlus_rapidity_true.Data(), "rapidity_true monte carlo DPlus in B0->D*#pi; Y; Entries", 5000, -20, 20);
  hist_mc_DPlus_rapidity_true->Sumw2();
  hist_mc_DPlus_rapidity_true->SetLineColor(6);
  hist_mc_DPlus_rapidity_true->SetMarkerStyle(20);
  hist_mc_DPlus_rapidity_true->SetMarkerSize(0.6);
  hist_mc_DPlus_rapidity_true->SetMarkerColor(6);
  TH1F* histogram_mc_DPlus_rapidity_true = (TH1F*)hist_mc_DPlus_rapidity_true->Clone();
  fOutputB0Results->Add(histogram_mc_DPlus_rapidity_true);
  fResultsHistogramArray[0][9] = histogram_mc_DPlus_rapidity_true;

  TString name_mc_DPlus_first_pion_rapidity_true = "mc_DPlus_first_pion_rapidity_true";
  TH1F* hist_mc_DPlus_first_pion_rapidity_true = new TH1F(name_mc_DPlus_first_pion_rapidity_true.Data(), "rapidity_true monte carlo first pion of DPlus in B0->D*#pi; Y; Entries", 5000, -20, 20);
  hist_mc_DPlus_first_pion_rapidity_true->Sumw2();
  hist_mc_DPlus_first_pion_rapidity_true->SetLineColor(6);
  hist_mc_DPlus_first_pion_rapidity_true->SetMarkerStyle(20);
  hist_mc_DPlus_first_pion_rapidity_true->SetMarkerSize(0.6);
  hist_mc_DPlus_first_pion_rapidity_true->SetMarkerColor(6);
  TH1F* histogram_mc_DPlus_first_pion_rapidity_true = (TH1F*)hist_mc_DPlus_first_pion_rapidity_true->Clone();
  fOutputB0Results->Add(histogram_mc_DPlus_first_pion_rapidity_true);
  fResultsHistogramArray[0][10] = histogram_mc_DPlus_first_pion_rapidity_true;

  TString name_mc_DPlus_second_pion_rapidity_true = "mc_DPlus_second_pion_rapidity_true";
  TH1F* hist_mc_DPlus_second_pion_rapidity_true = new TH1F(name_mc_DPlus_second_pion_rapidity_true.Data(), "rapidity_true monte carlo second pion of DPlus in B0->D*#pi; Y; Entries", 5000, -20, 20);
  hist_mc_DPlus_second_pion_rapidity_true->Sumw2();
  hist_mc_DPlus_second_pion_rapidity_true->SetLineColor(6);
  hist_mc_DPlus_second_pion_rapidity_true->SetMarkerStyle(20);
  hist_mc_DPlus_second_pion_rapidity_true->SetMarkerSize(0.6);
  hist_mc_DPlus_second_pion_rapidity_true->SetMarkerColor(6);
  TH1F* histogram_mc_DPlus_second_pion_rapidity_true = (TH1F*)hist_mc_DPlus_second_pion_rapidity_true->Clone();
  fOutputB0Results->Add(histogram_mc_DPlus_second_pion_rapidity_true);
  fResultsHistogramArray[0][11] = histogram_mc_DPlus_second_pion_rapidity_true;

  TString name_mc_DPlus_kaon_rapidity_true = "mc_DPlus_kaon_rapidity_true";
  TH1F* hist_mc_DPlus_kaon_rapidity_true = new TH1F(name_mc_DPlus_kaon_rapidity_true.Data(), "rapidity_true monte carlo kaon of DPlus in B0->D*#pi; Y; Entries", 5000, -20, 20);
  hist_mc_DPlus_kaon_rapidity_true->Sumw2();
  hist_mc_DPlus_kaon_rapidity_true->SetLineColor(6);
  hist_mc_DPlus_kaon_rapidity_true->SetMarkerStyle(20);
  hist_mc_DPlus_kaon_rapidity_true->SetMarkerSize(0.6);
  hist_mc_DPlus_kaon_rapidity_true->SetMarkerColor(6);
  TH1F* histogram_mc_DPlus_kaon_rapidity_true = (TH1F*)hist_mc_DPlus_kaon_rapidity_true->Clone();
  fOutputB0Results->Add(histogram_mc_DPlus_kaon_rapidity_true);
  fResultsHistogramArray[0][12] = histogram_mc_DPlus_kaon_rapidity_true;

  TString name_mc_B0_pseudorapidity_true = "mc_B0_pseudorapidity_true";
  TH1F* hist_mc_B0_pseudorapidity_true = new TH1F(name_mc_B0_pseudorapidity_true.Data(), "pseudorapidity_true monte carlo B0 in B0->D*#pi; #eta; Entries", 5000, -20, 20);
  hist_mc_B0_pseudorapidity_true->Sumw2();
  hist_mc_B0_pseudorapidity_true->SetLineColor(6);
  hist_mc_B0_pseudorapidity_true->SetMarkerStyle(20);
  hist_mc_B0_pseudorapidity_true->SetMarkerSize(0.6);
  hist_mc_B0_pseudorapidity_true->SetMarkerColor(6);
  TH1F* histogram_mc_B0_pseudorapidity_true = (TH1F*)hist_mc_B0_pseudorapidity_true->Clone();
  fOutputB0Results->Add(histogram_mc_B0_pseudorapidity_true);
  fResultsHistogramArray[0][14] = histogram_mc_B0_pseudorapidity_true;

  TString name_mc_B0_pion_pseudorapidity_true = "mc_B0_pion_pseudorapidity_true";
  TH1F* hist_mc_B0_pion_pseudorapidity_true = new TH1F(name_mc_B0_pion_pseudorapidity_true.Data(), "pseudorapidity_true monte carlo first pion of B0 in B0->D*#pi; #eta; Entries", 5000, -20, 20);
  hist_mc_B0_pion_pseudorapidity_true->Sumw2();
  hist_mc_B0_pion_pseudorapidity_true->SetLineColor(6);
  hist_mc_B0_pion_pseudorapidity_true->SetMarkerStyle(20);
  hist_mc_B0_pion_pseudorapidity_true->SetMarkerSize(0.6);
  hist_mc_B0_pion_pseudorapidity_true->SetMarkerColor(6);
  TH1F* histogram_mc_B0_pion_pseudorapidity_true = (TH1F*)hist_mc_B0_pion_pseudorapidity_true->Clone();
  fOutputB0Results->Add(histogram_mc_B0_pion_pseudorapidity_true);
  fResultsHistogramArray[0][15] = histogram_mc_B0_pion_pseudorapidity_true;

  TString name_mc_DPlus_pseudorapidity_true = "mc_DPlus_pseudorapidity_true";
  TH1F* hist_mc_DPlus_pseudorapidity_true = new TH1F(name_mc_DPlus_pseudorapidity_true.Data(), "pseudorapidity_true monte carlo DPlus in B0->D*#pi; #eta; Entries", 5000, -20, 20);
  hist_mc_DPlus_pseudorapidity_true->Sumw2();
  hist_mc_DPlus_pseudorapidity_true->SetLineColor(6);
  hist_mc_DPlus_pseudorapidity_true->SetMarkerStyle(20);
  hist_mc_DPlus_pseudorapidity_true->SetMarkerSize(0.6);
  hist_mc_DPlus_pseudorapidity_true->SetMarkerColor(6);
  TH1F* histogram_mc_DPlus_pseudorapidity_true = (TH1F*)hist_mc_DPlus_pseudorapidity_true->Clone();
  fOutputB0Results->Add(histogram_mc_DPlus_pseudorapidity_true);
  fResultsHistogramArray[0][16] = histogram_mc_DPlus_pseudorapidity_true;

  TString name_mc_DPlus_first_pion_pseudorapidity_true = "mc_DPlus_first_pion_pseudorapidity_true";
  TH1F* hist_mc_DPlus_first_pion_pseudorapidity_true = new TH1F(name_mc_DPlus_first_pion_pseudorapidity_true.Data(), "pseudorapidity_true monte carlo first pion of DPlus in B0->D*#pi; #eta; Entries", 5000, -20, 20);
  hist_mc_DPlus_first_pion_pseudorapidity_true->Sumw2();
  hist_mc_DPlus_first_pion_pseudorapidity_true->SetLineColor(6);
  hist_mc_DPlus_first_pion_pseudorapidity_true->SetMarkerStyle(20);
  hist_mc_DPlus_first_pion_pseudorapidity_true->SetMarkerSize(0.6);
  hist_mc_DPlus_first_pion_pseudorapidity_true->SetMarkerColor(6);
  TH1F* histogram_mc_DPlus_first_pion_pseudorapidity_true = (TH1F*)hist_mc_DPlus_first_pion_pseudorapidity_true->Clone();
  fOutputB0Results->Add(histogram_mc_DPlus_first_pion_pseudorapidity_true);
  fResultsHistogramArray[0][17] = histogram_mc_DPlus_first_pion_pseudorapidity_true;

  TString name_mc_DPlus_second_pion_pseudorapidity_true = "mc_DPlus_second_pion_pseudorapidity_true";
  TH1F* hist_mc_DPlus_second_pion_pseudorapidity_true = new TH1F(name_mc_DPlus_second_pion_pseudorapidity_true.Data(), "pseudorapidity_true monte carlo second pion of DPlus in B0->D*#pi; #eta; Entries", 5000, -20, 20);
  hist_mc_DPlus_second_pion_pseudorapidity_true->Sumw2();
  hist_mc_DPlus_second_pion_pseudorapidity_true->SetLineColor(6);
  hist_mc_DPlus_second_pion_pseudorapidity_true->SetMarkerStyle(20);
  hist_mc_DPlus_second_pion_pseudorapidity_true->SetMarkerSize(0.6);
  hist_mc_DPlus_second_pion_pseudorapidity_true->SetMarkerColor(6);
  TH1F* histogram_mc_DPlus_second_pion_pseudorapidity_true = (TH1F*)hist_mc_DPlus_second_pion_pseudorapidity_true->Clone();
  fOutputB0Results->Add(histogram_mc_DPlus_second_pion_pseudorapidity_true);
  fResultsHistogramArray[0][18] = histogram_mc_DPlus_second_pion_pseudorapidity_true;

  TString name_mc_DPlus_kaon_pseudorapidity_true = "mc_DPlus_kaon_pseudorapidity_true";
  TH1F* hist_mc_DPlus_kaon_pseudorapidity_true = new TH1F(name_mc_DPlus_kaon_pseudorapidity_true.Data(), "pseudorapidity_true monte carlo kaon of DPlus in B0->D*#pi; #eta; Entries", 5000, -20, 20);
  hist_mc_DPlus_kaon_pseudorapidity_true->Sumw2();
  hist_mc_DPlus_kaon_pseudorapidity_true->SetLineColor(6);
  hist_mc_DPlus_kaon_pseudorapidity_true->SetMarkerStyle(20);
  hist_mc_DPlus_kaon_pseudorapidity_true->SetMarkerSize(0.6);
  hist_mc_DPlus_kaon_pseudorapidity_true->SetMarkerColor(6);
  TH1F* histogram_mc_DPlus_kaon_pseudorapidity_true = (TH1F*)hist_mc_DPlus_kaon_pseudorapidity_true->Clone();
  fOutputB0Results->Add(histogram_mc_DPlus_kaon_pseudorapidity_true);
  fResultsHistogramArray[0][19] = histogram_mc_DPlus_kaon_pseudorapidity_true;

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
      fOutputB0Results->Add(histogram_cut_optimization_signal);
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
      fOutputB0Results->Add(histogram_cut_optimization_background);
      fResultsHistogramArray2D[1][k] = histogram_cut_optimization_background;
    }
  }


  //==================================================

  TString name_B0s_in_analysis = "B0_in_analysis";
  TH1F* hist_B0s_in_analysis = new TH1F(name_B0s_in_analysis.Data(), "Number of B0s to kpipipi in the Analysis; Entries", 20, 0, 20);
  hist_B0s_in_analysis->Sumw2();
  hist_B0s_in_analysis->SetLineColor(6);
  hist_B0s_in_analysis->SetMarkerStyle(20);
  hist_B0s_in_analysis->SetMarkerSize(0.6);
  hist_B0s_in_analysis->SetMarkerColor(6);
  hist_B0s_in_analysis->SetStats(kTRUE);
  hist_B0s_in_analysis->GetXaxis()->SetBinLabel(1, "no. of B0");
  hist_B0s_in_analysis->GetXaxis()->SetBinLabel(2, "no. of B0 to kpipipi");
  hist_B0s_in_analysis->GetXaxis()->SetBinLabel(3, "no. with all tracks in event");
  hist_B0s_in_analysis->GetXaxis()->SetBinLabel(4, "no. ...");
  hist_B0s_in_analysis->GetXaxis()->SetBinLabel(5, "no. ...");
  hist_B0s_in_analysis->GetXaxis()->SetBinLabel(6, "no. ...");
  hist_B0s_in_analysis->GetXaxis()->SetBinLabel(7, "no. ...");
  hist_B0s_in_analysis->GetXaxis()->SetBinLabel(8, "no. ...");
  TH1F* hist_B0s_in_analysis_mc = (TH1F*)hist_B0s_in_analysis->Clone();
  fOutputB0Results->Add(hist_B0s_in_analysis_mc);
  fResultsHistogramArray[3][0] = hist_B0s_in_analysis_mc;

  TString name_B0_per_bin = "B0_per_bin";
  TH1F* hist_B0_per_bin = new TH1F(name_B0_per_bin.Data(), "Number of B0 to kpipi in the Analysis per bin; Entries", fnPtBins, 0, fnPtBins);
  for (Int_t i = 0; i < fnPtBins; ++i)
  {
    TString bin_name = "";
    bin_name += fPtBinLimits[i];
    bin_name += "-";
    bin_name += fPtBinLimits[i + 1];
    hist_B0_per_bin->GetXaxis()->SetBinLabel(i + 1, bin_name);
  }
  TH1F* hist_B0_per_bin_mc = (TH1F*)hist_B0_per_bin->Clone();
  fOutputB0Results->Add(hist_B0_per_bin_mc);
  fResultsHistogramArray[3][1] = hist_B0_per_bin_mc;

  TString name_B0_per_bin_in_Acc = "B0_per_bin_in_Acc";
  TH1F* hist_B0_per_bin_in_Acc = new TH1F(name_B0_per_bin_in_Acc.Data(), "Number of B0 to kpipi in the Analysis per bin with all daughters in acceptance; Entries", fnPtBins, 0, fnPtBins);
  for (Int_t i = 0; i < fnPtBins; ++i)
  {
    TString bin_name = "";
    bin_name += fPtBinLimits[i];
    bin_name += "-";
    bin_name += fPtBinLimits[i + 1];
    hist_B0_per_bin_in_Acc->GetXaxis()->SetBinLabel(i + 1, bin_name);
  }
  TH1F* hist_B0_per_bin_in_Acc_mc = (TH1F*)hist_B0_per_bin_in_Acc->Clone();
  fOutputB0Results->Add(hist_B0_per_bin_in_Acc_mc);
  fResultsHistogramArray[3][2] = hist_B0_per_bin_in_Acc_mc;

  TString name_Reconstructed_B0_per_bin_in_Acc = "Reconstructed_B0_per_bin_in_Acc";
  TH1F* hist_Reconstructed_B0_per_bin_in_Acc = new TH1F(name_Reconstructed_B0_per_bin_in_Acc.Data(), "Number of reconstructed B0 in the Analysis per bin with all daughters in acceptance; Entries", fnPtBins, 0, fnPtBins);
  for (Int_t i = 0; i < fnPtBins; ++i)
  {
    TString bin_name = "";
    bin_name += fPtBinLimits[i];
    bin_name += "-";
    bin_name += fPtBinLimits[i + 1];
    hist_Reconstructed_B0_per_bin_in_Acc->GetXaxis()->SetBinLabel(i + 1, bin_name);
  }
  TH1F* hist_Reconstructed_B0_per_bin_in_Acc_mc = (TH1F*)hist_Reconstructed_B0_per_bin_in_Acc->Clone();
  fOutputB0Results->Add(hist_Reconstructed_B0_per_bin_in_Acc_mc);
  fResultsHistogramArray[3][3] = hist_Reconstructed_B0_per_bin_in_Acc_mc;

  //======================================================================================================================================================

  //we make the histograms for the Pions and Kaon
  for (Int_t i = 0; i < 3; i++) {

    TString add_name = "";
    TList * listout = 0x0;
    if (i == 0) listout = fOutputDPlusKaon;
    if (i == 1) listout = fOutputDPlusPions;
    if (i == 2) listout = fOutputB0Pion;

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

      for (Int_t k = 0; k < 11; ++k)
      {
        if (k == 0) {name_Histogram = "ptTrack"; discription_Histogram = "pt track; #it{p}_{T} (GeV/#it{c}); Entries"; numberOfBins = 600; lowerBound = 0; upperBound = 30;}
        if (k == 1) {name_Histogram = "momentumTrack"; discription_Histogram = "momentum track; p (GeV/#it{c}); Entries"; numberOfBins = 600; lowerBound = 0; upperBound = 30;}
        if (k == 2) {name_Histogram = "numberOfITS"; discription_Histogram = "Number of ITS clusters track; [#]; Entries"; numberOfBins = 10; lowerBound = -0.5; upperBound = 9.5;}
        if (k == 3) {name_Histogram = "numberOfTPC"; discription_Histogram = "Number of TPC clusters track; [#]; Entries"; numberOfBins = 601; lowerBound = -0.5; upperBound = 600.5;}
        if (k == 4) {name_Histogram = "pointsOnITS"; discription_Histogram = "Number of ITS clusters track per layer; [#]; Entries"; numberOfBins = 10; lowerBound = -0.5; upperBound = 9.5;}
        if (k == 5) {name_Histogram = "nSigmaTPC"; discription_Histogram = "n sigma TPC for track PID; sigma; Entries"; numberOfBins = 500; lowerBound = -5; upperBound = 5;}
        if (k == 6) {name_Histogram = "nSigmaTOF"; discription_Histogram = "n sigma TOF for track PID; sigma; Entries"; numberOfBins = 500; lowerBound = -5; upperBound = 5;}
        if (k == 7) {name_Histogram = "nSigmaTPCandTOF"; discription_Histogram = "n sigma TPC and TOF for track PID; a.u.; Entries"; numberOfBins = 1000; lowerBound = 0; upperBound = 10;}
        if (k == 8) {name_Histogram = "impactParameter"; discription_Histogram = "Impact Parameter track;  (cm); Entries"; numberOfBins = 2000; lowerBound = 0; upperBound = 0.5;}
        if (k == 9) {name_Histogram = "NormalizedImpactParameter"; discription_Histogram = "Normalized Impact Parameter track;  (cm); Entries"; numberOfBins = 500; lowerBound = 0; upperBound = 100;}
        if (k == 10) {name_Histogram = "EtaTrack"; discription_Histogram = "Eta Track; Entries"; numberOfBins = 1000; lowerBound = -3; upperBound = 3;}

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

    TString name_ptB0_vs_ptTrack = "ptB0_vs_ptTrackBackground";
    TH2F* hist_ptB0_vs_ptTrack = new TH2F(name_ptB0_vs_ptTrack.Data(), "Pt B0 vs Pt ; #it{p}_{T} B0 (GeV/#it{c}); #it{p}_{T} track (GeV/#it{c})", 100, 0, 30, 100, 0, 30);
    hist_ptB0_vs_ptTrack->Sumw2();
    hist_ptB0_vs_ptTrack->SetLineColor(6);
    hist_ptB0_vs_ptTrack->SetMarkerStyle(20);
    hist_ptB0_vs_ptTrack->SetMarkerSize(0.6);
    hist_ptB0_vs_ptTrack->SetMarkerColor(6);
    TH2F* histogram_ptB0_vs_ptTrack = (TH2F*)hist_ptB0_vs_ptTrack->Clone();
    listout->Add(histogram_ptB0_vs_ptTrack);
    fDaughterHistogramArray2D[i][4] = histogram_ptB0_vs_ptTrack;

    TString name_ptB0_vs_ptTrackSignal = "ptB0_vs_ptTrackSignal";
    TH2F* hist_ptB0_vs_ptTrackSignal = new TH2F(name_ptB0_vs_ptTrackSignal.Data(), "Pt B0 vs Pt ; #it{p}_{T} B0 (GeV/#it{c}); #it{p}_{T} track (GeV/#it{c})", 100, 0, 30, 100, 0, 30);
    hist_ptB0_vs_ptTrackSignal->Sumw2();
    hist_ptB0_vs_ptTrackSignal->SetLineColor(4);
    hist_ptB0_vs_ptTrackSignal->SetMarkerStyle(20);
    hist_ptB0_vs_ptTrackSignal->SetMarkerSize(0.6);
    hist_ptB0_vs_ptTrackSignal->SetMarkerColor(6);
    TH2F* histogram_ptB0_vs_ptTrackSignal = (TH2F*)hist_ptB0_vs_ptTrackSignal->Clone();
    listout->Add(histogram_ptB0_vs_ptTrackSignal);
    fDaughterHistogramArray2D[i][5] = histogram_ptB0_vs_ptTrackSignal;
  }

  //we make the histograms for the reconstructed particles
  for (Int_t i = 0; i < 3; i++) {

    TString add_name = "";
    TList * listout = 0x0;
    Int_t nHistogramSets = 0;
    if (i == 0) {listout = fOutputDPlus; nHistogramSets = 6 + 2 * fnPtBins;}
    if (i == 1) {listout = fOutputB0; nHistogramSets = 6 + 2 * fnPtBins;}
    if (i == 2) {listout = fOutputDPlus_DPlusPt; nHistogramSets = 2 * fnPtBinsDPlusforDPlusptbin;}


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
        if (j % 2 == 0) {add_name = "_ptbin_"; add_name += fPtBinLimitsDPlusforDPlusptbin[j / 2]; add_name += "_to_"; add_name += fPtBinLimitsDPlusforDPlusptbin[1 + j / 2];}
        if (j % 2 == 1) {add_name = "Signal_ptbin_"; add_name += fPtBinLimitsDPlusforDPlusptbin[(j - 1) / 2]; add_name += "_to_"; add_name += fPtBinLimitsDPlusforDPlusptbin[1 + (j - 1) / 2];}
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
      // Int_t numberOfBinsTwo = 0;
      // Double_t lowerBoundTwo = 0.0;
      // Double_t upperBoundTwo = 0.0;

      if (i == 0 || i == 2)
      {
        for (Int_t k = 0; k < 49; ++k)
        {
          if (k == 0) {name_Histogram = "invariantMassMother"; discription_Histogram = "mass mother candidate; m (GeV/#it{c}^{2}); Entries"; numberOfBins = 10000; lowerBound = 0; upperBound = 10;}
          if (k == 1) {name_Histogram = "deltaMassMother"; discription_Histogram = "mass mother candidate; m (GeV/#it{c}^{2}); Entries"; numberOfBins = 10000; lowerBound = 0; upperBound = 10;}
          if (k == 2) {name_Histogram = "ptMother"; discription_Histogram = "pt mother; #it{p}_{T} (GeV/#it{c}); Entries"; numberOfBins = 300; lowerBound = 0; upperBound = 30;}
          if (k == 3) {name_Histogram = "ptFirstDaughter"; discription_Histogram = "pt first daughter; #it{p}_{T} (GeV/#it{c}); Entries"; numberOfBins = 300; lowerBound = 0; upperBound = 30;}
          if (k == 4) {name_Histogram = "ptSecondDaughter"; discription_Histogram = "pt second daughter; #it{p}_{T} (GeV/#it{c}); Entries"; numberOfBins = 300; lowerBound = 0; upperBound = 30;}
          if (k == 5) {name_Histogram = "ptThirdDaughter"; discription_Histogram = "pt third daughter; #it{p}_{T} (GeV/#it{c}); Entries"; numberOfBins = 300; lowerBound = 0; upperBound = 30;}
          if (k == 6) {name_Histogram = "d0Mother"; discription_Histogram = "d0 mother;  (cm); Entries"; numberOfBins = 500; lowerBound = 0; upperBound = 1.0;}
          if (k == 7) {name_Histogram = "d0FirstDaughter"; discription_Histogram = "d0 first daughter;  (cm); Entries"; numberOfBins = 500; lowerBound = 0; upperBound = 1;}
          if (k == 8) {name_Histogram = "d0SecondDaughter"; discription_Histogram = "d0 second daughter;  (cm); Entries"; numberOfBins = 500; lowerBound = 0; upperBound = 1;}
          if (k == 9) {name_Histogram = "d0ThirdDaughter"; discription_Histogram = "d0 third daughter;  (cm); Entries"; numberOfBins = 500; lowerBound = 0; upperBound = 1;}
          if (k == 10) {name_Histogram = "pointingAngleMother"; discription_Histogram = "pointing angle;  (Cos(#theta)); Entries"; numberOfBins = 100; lowerBound = -1; upperBound = 1;}
          if (k == 11) {name_Histogram = "pointingAngleMotherXY"; discription_Histogram = "pointing angle XY;  (Cos(#theta)); Entries"; numberOfBins = 1000; lowerBound = -1; upperBound = 1;}
          if (k == 12) {name_Histogram = "impactProduct123"; discription_Histogram = "impact product123; (cm^{3}); Entries"; numberOfBins = 500; lowerBound = -0.001; upperBound = 0.001;}
          if (k == 13) {name_Histogram = "impactProduct12"; discription_Histogram = "impact product12; (cm^{2}); Entries"; numberOfBins = 500; lowerBound = -0.01; upperBound = 0.01;}
          if (k == 14) {name_Histogram = "impactProduct13"; discription_Histogram = "impact product13; (cm^{2}); Entries"; numberOfBins = 500; lowerBound = -0.01; upperBound = 0.01;}
          if (k == 15) {name_Histogram = "impactProduct23"; discription_Histogram = "impact product23; (cm^{2}); Entries"; numberOfBins = 500; lowerBound = -0.01; upperBound = 0.01;}
          if (k == 16) {name_Histogram = "impactProductXY"; discription_Histogram = "impact product XY; (cm^{3}); Entries"; numberOfBins = 200; lowerBound = 0; upperBound = 0.5;}
          if (k == 17) {name_Histogram = "dcaMother12"; discription_Histogram = "dca mother12; distance (cm); Entries"; numberOfBins = 500; lowerBound = 0; upperBound = 0.25;}
          if (k == 18) {name_Histogram = "dcaMother13"; discription_Histogram = "dca mother13; distance (cm); Entries"; numberOfBins = 500; lowerBound = 0; upperBound = 0.25;}
          if (k == 19) {name_Histogram = "dcaMother23"; discription_Histogram = "dca mother23; distance (cm); Entries"; numberOfBins = 500; lowerBound = 0; upperBound = 0.25;}
          if (k == 20) {name_Histogram = "vertexDistance"; discription_Histogram = "vertex distance between mother and primary vertex; distance (cm); Entries"; numberOfBins = 100; lowerBound = 0; upperBound = 1;}
          if (k == 21) {name_Histogram = "vertexDistanceXY"; discription_Histogram = "vertex distance between mother and primary vertex XY; distance (cm); Entries"; numberOfBins = 1000; lowerBound = 0; upperBound = 10;}
          if (k == 22) {name_Histogram = "normDecayLength"; discription_Histogram = "Normalized decay length w.r.t primary vertex; (cm); Entries"; numberOfBins = 100; lowerBound = 0; upperBound = 50;}
          if (k == 23) {name_Histogram = "normDecayLengthXY"; discription_Histogram = "Normalized decay length w.r.t primary vertex XY; (cm); Entries"; numberOfBins = 100; lowerBound = 0; upperBound = 50;}
          if (k == 24) {name_Histogram = "pseudoProperDecayTime"; discription_Histogram = "Pseudo Proper Decay Time w.r.t primary vertex; [a.u.]; Entries"; numberOfBins = 1000; lowerBound = -10; upperBound = 10;}
          if (k == 25) {name_Histogram = "DecayTime"; discription_Histogram = "Decay Time w.r.t primary vertex; [a.u.]; Entries"; numberOfBins = 100; lowerBound = 0; upperBound = 0.00000001;}
          if (k == 26) {name_Histogram = "normDecayTime"; discription_Histogram = "Normalized Decay Time w.r.t primary vertex; [a.u.]; Entries"; numberOfBins = 100; lowerBound = 0; upperBound = 0.0000001;}
          if (k == 27) {name_Histogram = "smallestAngleMotherDaughter"; discription_Histogram = "smallest flight angle mother and daughter; (Cos(#phi)); Entries"; numberOfBins = 100; lowerBound = -1; upperBound = 1;}
          if (k == 28) {name_Histogram = "largestAngleMotherDaughter"; discription_Histogram = "largest flight angle mother and daughter; (Cos(#phi)); Entries"; numberOfBins = 100; lowerBound = -1; upperBound = 1;}
          if (k == 29) {name_Histogram = "AngleDifference"; discription_Histogram = "Angle Difference (largest - smallest); (Cos(#phi)); Entries"; numberOfBins = 100; lowerBound = -2.0; upperBound = 2;}
          if (k == 30) {name_Histogram = "Normd0FirstDaughter"; discription_Histogram = "norm d0 first daughter;  (cm); Entries"; numberOfBins = 500; lowerBound = 0; upperBound = 100;}
          if (k == 31) {name_Histogram = "Normd0SecondDaughter"; discription_Histogram = "norm d0 second daughter;  (cm); Entries"; numberOfBins = 500; lowerBound = 0; upperBound = 100;}
          if (k == 32) {name_Histogram = "Normd0ThirdDaughter"; discription_Histogram = "norm d0 third daughter;  (cm); Entries"; numberOfBins = 500; lowerBound = 0; upperBound = 100;}
          if (k == 33) {name_Histogram = "Normd0Mother"; discription_Histogram = "norm d0 mother; (cm); Entries"; numberOfBins = 500; lowerBound = 0; upperBound = 500;}
          if (k == 34) {name_Histogram = "NormimpactProduct"; discription_Histogram = "norm impact product; (cm^{3}); Entries"; numberOfBins = 500; lowerBound = -200; upperBound = 200;}
          if (k == 35) {name_Histogram = "dist12"; discription_Histogram = "dist12; (cm); Entries"; numberOfBins = 500; lowerBound = -3; upperBound = 3;}
          if (k == 36) {name_Histogram = "dist23"; discription_Histogram = "dist23; (cm); Entries"; numberOfBins = 500; lowerBound = -3; upperBound = 3;}
          if (k == 37) {name_Histogram = "sigmavertex"; discription_Histogram = "sigmavertex; []; Entries"; numberOfBins = 500; lowerBound = -200; upperBound = 200;}

          if (k == 38) {name_Histogram = "pointingAngleToB0"; discription_Histogram = "Pointing angle w.r.t. B0 decay vertex; (Cos(#theta)); Entries"; numberOfBins = 200; lowerBound = -1; upperBound = 1;}
          if (k == 39) {name_Histogram = "d0MotherToB0"; discription_Histogram = "d0 Mother w.r.t. B0 decay vertex; (cm); Entries"; numberOfBins = 500; lowerBound = 0; upperBound = 1;}
          if (k == 40) {name_Histogram = "d0FirstDaughterToB0"; discription_Histogram = "d0 first daughter w.r.t. B0 decay vertex; (cm); Entries"; numberOfBins = 500; lowerBound = 0; upperBound = 1;}
          if (k == 41) {name_Histogram = "d0SecondDaughterToB0"; discription_Histogram = "d0 second daughter w.r.t. B0 decay vertex; (cm); Entries"; numberOfBins = 500; lowerBound = 0; upperBound = 1;}
          if (k == 42) {name_Histogram = "d0ThirdDaughterToB0"; discription_Histogram = "d0 third daughter w.r.t. B0 decay vertex; (cm); Entries"; numberOfBins = 500; lowerBound = 0; upperBound = 1;}
          if (k == 43) {name_Histogram = "impactProductToB0"; discription_Histogram = "impact product w.r.t. B0 decay vertex; (cm^{3}); Entries"; numberOfBins = 400; lowerBound = -0.005; upperBound = 0.005;}
          if (k == 44) {name_Histogram = "impactProductXYToB0"; discription_Histogram = "impact product XY w.r.t. B0 decay vertex; (cm^{2}); Entries"; numberOfBins = 100; lowerBound = 0; upperBound = 0.5;}
          if (k == 45) {name_Histogram = "normDecayLengthToB0"; discription_Histogram = "Normalized decay length w.r.t. B0 decay vertex; (cm); Entries"; numberOfBins = 100; lowerBound = 0; upperBound = 50;}
          if (k == 46) {name_Histogram = "pseudoProperDecayTimeToB0"; discription_Histogram = "Pseudo Proper Decay Time w.r.t B0 vertex; [a.u.]; Entries"; numberOfBins = 1000; lowerBound = -1; upperBound = 1;}
          if (k == 47) {name_Histogram = "DecayTimeToB0"; discription_Histogram = "Decay Time w.r.t B0 vertex; [a.u.]; Entries"; numberOfBins = 100; lowerBound = 0; upperBound = 0.00000001;}
          if (k == 48) {name_Histogram = "normDecayTimeToB0"; discription_Histogram = "Normalized Decay Time w.r.t B0 vertex; [a.u.]; Entries"; numberOfBins = 100; lowerBound = 0; upperBound = 0.00000001;}

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
      }

      if (i == 1)
      {
        for (Int_t k = 0; k < 28; ++k)
        {
          if (k == 0) {name_Histogram = "invariantMassMother"; discription_Histogram = "mass mother candidate; m (GeV/#it{c}^{2}); Entries"; numberOfBins = 10000; lowerBound = 0; upperBound = 10;}
          if (k == 1) {name_Histogram = "deltaMassMother"; discription_Histogram = "mass mother candidate; m (GeV/#it{c}^{2}); Entries"; numberOfBins = 10000; lowerBound = 0; upperBound = 10;}
          if (k == 2) {name_Histogram = "ptMother"; discription_Histogram = "pt mother; #it{p}_{T} (GeV/#it{c}); Entries"; numberOfBins = 300; lowerBound = 0; upperBound = 30;}
          if (k == 3) {name_Histogram = "ptFirstDaughter"; discription_Histogram = "pt first daughter; #it{p}_{T} (GeV/#it{c}); Entries"; numberOfBins = 300; lowerBound = 0; upperBound = 30;}
          if (k == 4) {name_Histogram = "ptSecondDaughter"; discription_Histogram = "pt second daughter; #it{p}_{T} (GeV/#it{c}); Entries"; numberOfBins = 300; lowerBound = 0; upperBound = 30;}
          if (k == 5) {name_Histogram = "d0Mother"; discription_Histogram = "d0 mother;  (cm); Entries"; numberOfBins = 500; lowerBound = 0; upperBound = 1.0;}
          if (k == 6) {name_Histogram = "d0FirstDaughter"; discription_Histogram = "d0 first daughter;  (cm); Entries"; numberOfBins = 500; lowerBound = 0; upperBound = 1;}
          if (k == 7) {name_Histogram = "d0SecondDaughter"; discription_Histogram = "d0 second daughter;  (cm); Entries"; numberOfBins = 500; lowerBound = 0; upperBound = 1;}
          if (k == 8) {name_Histogram = "pointingAngleMother"; discription_Histogram = "pointing angle;  (Cos(#theta)); Entries"; numberOfBins = 100; lowerBound = -1; upperBound = 1;}
          if (k == 9) {name_Histogram = "pointingAngleMotherXY"; discription_Histogram = "pointing angle XY;  (Cos(#theta)); Entries"; numberOfBins = 1000; lowerBound = -1; upperBound = 1;}
          if (k == 10) {name_Histogram = "impactProduct"; discription_Histogram = "impact product; (cm^{2}); Entries"; numberOfBins = 500; lowerBound = -0.01; upperBound = 0.01;}
          if (k == 11) {name_Histogram = "impactProductXY"; discription_Histogram = "impact product XY; (cm^{2}); Entries"; numberOfBins = 200; lowerBound = 0; upperBound = 0.5;}
          if (k == 12) {name_Histogram = "dcaMother"; discription_Histogram = "dca mother; distance (cm); Entries"; numberOfBins = 500; lowerBound = 0; upperBound = 0.25;}
          if (k == 13) {name_Histogram = "vertexDistance"; discription_Histogram = "vertex distance between mother and primary vertex; distance (cm); Entries"; numberOfBins = 100; lowerBound = 0; upperBound = 1;}
          if (k == 14) {name_Histogram = "vertexDistanceXY"; discription_Histogram = "vertex distance between mother and primary vertex XY; distance (cm); Entries"; numberOfBins = 1000; lowerBound = 0; upperBound = 10;}
          if (k == 15) {name_Histogram = "normDecayLength"; discription_Histogram = "Normalized decay length w.r.t primary vertex; (cm); Entries"; numberOfBins = 100; lowerBound = 0; upperBound = 50;}
          if (k == 16) {name_Histogram = "normDecayLengthXY"; discription_Histogram = "Normalized decay length w.r.t primary vertex XY; (cm); Entries"; numberOfBins = 100; lowerBound = 0; upperBound = 50;}
          if (k == 17) {name_Histogram = "pseudoProperDecayTime"; discription_Histogram = "Pseudo Proper Decay Time w.r.t primary vertex; [a.u.]; Entries"; numberOfBins = 1000; lowerBound = -10; upperBound = 10;}
          if (k == 18) {name_Histogram = "DecayTime"; discription_Histogram = "Decay Time w.r.t primary vertex; [a.u.]; Entries"; numberOfBins = 100; lowerBound = 0; upperBound = 0.00000001;}
          if (k == 19) {name_Histogram = "normDecayTime"; discription_Histogram = "Normalized Decay Time w.r.t primary vertex; [a.u.]; Entries"; numberOfBins = 100; lowerBound = 0; upperBound = 0.0000001;}
          if (k == 20) {name_Histogram = "angleMotherFirstDaughter"; discription_Histogram = "flight angle mother and first daughter; [Cos(#phi)]; Entries"; numberOfBins = 100; lowerBound = -1; upperBound = 1;}
          if (k == 21) {name_Histogram = "angleMotherSecondDaughter"; discription_Histogram = "flight angle mother and second daughter; [Cos(#phi)]; Entries"; numberOfBins = 100; lowerBound = 0.5; upperBound = 1;}
          if (k == 22) {name_Histogram = "angleBetweenBothDaughters"; discription_Histogram = "angle between both daughters; [Cos(#phi)]; Entries"; numberOfBins = 100; lowerBound = -1; upperBound = 1;}
          if (k == 23) {name_Histogram = "cosThetaStar"; discription_Histogram = "cosThetaStar; [Cos(#theta*)]; Entries"; numberOfBins = 200; lowerBound = -2; upperBound = 2;}
          if (k == 24) {name_Histogram = "Normd0FirstDaughter"; discription_Histogram = "norm d0 first daughter;  (cm); Entries"; numberOfBins = 500; lowerBound = 0; upperBound = 100;}
          if (k == 25) {name_Histogram = "Normd0SecondDaughter"; discription_Histogram = "norm d0 second daughter;  (cm); Entries"; numberOfBins = 500; lowerBound = 0; upperBound = 100;}
          if (k == 26) {name_Histogram = "Normd0Mother"; discription_Histogram = "norm d0 mother; (cm); Entries"; numberOfBins = 500; lowerBound = 0; upperBound = 500;}
          if (k == 27) {name_Histogram = "NormimpactProduct"; discription_Histogram = "norm impact product; (cm^{2}); Entries"; numberOfBins = 500; lowerBound = -200; upperBound = 200;}


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
      }

      TList * listOther = new TList();
      listOther->SetOwner();
      listOther->SetName("OtherHistograms");
      listCandidates->Add(listOther);
      for (Int_t k = 49; k < 56; ++k)
      {

        if (k == 49) {name_Histogram = "yMother"; discription_Histogram = "eta mother; #eta; Entries"; numberOfBins = 100; lowerBound = -10; upperBound = 10;}
        if (k == 50) {name_Histogram = "etaMother"; discription_Histogram = "eta mother; #eta; Entries"; numberOfBins = 100; lowerBound = -2; upperBound = 2;}
        if (k == 51) {name_Histogram = "phiMother"; discription_Histogram = "phi mother; #phi; Entries"; numberOfBins = 25; lowerBound = 0; upperBound = 2 * TMath::Pi();}
        if (k == 52) {name_Histogram = "vertexX"; discription_Histogram = "Vertex position; (cm); Entries"; numberOfBins = 500; lowerBound = -5; upperBound = 5;}
        if (k == 53) {name_Histogram = "vertexY"; discription_Histogram = "Vertex position; (cm); Entries"; numberOfBins = 500; lowerBound = -5; upperBound = 5;}
        if (k == 54) {name_Histogram = "vertexZ"; discription_Histogram = "Vertex position; (cm); Entries"; numberOfBins = 500; lowerBound = -20; upperBound = 20;}
        if (k == 55) {name_Histogram = "vertexChi2NDF"; discription_Histogram = "#Chi^{2} per NDF Vertex; [#Chi^{2}/NDF]; Entries"; numberOfBins = 1000; lowerBound = 0; upperBound = 50;}

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
        listOther->Add(histogram_Clone);
        fMotherHistogramArray[i][j][k] = histogram_Clone;
      }

      // optional code for 2D histograms, to be updated for this decay
      // name_Histogram = "";
      // discription_Histogram  = "";
      // numberOfBins = 0;
      // lowerBound = 0.0;
      // upperBound = 0.0;
      // numberOfBinsTwo = 0;
      // lowerBoundTwo = 0.0;
      // upperBoundTwo = 0.0;

      // //we make the 2D histograms for the reconstructed particles
      // Int_t nFirst = 0;
      // Int_t nSecond = 1;
      // Int_t nVariables = 10;
      // Int_t nHistograms = nVariables * (nVariables - 1) / 2;

      // TList * list2D = new TList();
      // list2D->SetOwner();
      // TString name2D = "2D_Histograms";
      // name2D += add_name;
      // list2D->SetName(name2D);
      // listCandidates->Add(list2D);

      // for (Int_t k = 0; k < nHistograms; ++k)
      // {
      //   numberOfBins = 50; numberOfBinsTwo = 50;
      //   if (nFirst == 0) {name_Histogram = "dPlusFirstDaughter"; discription_Histogram = "dPlus first daughter (cm);"; lowerBound = 0; upperBound = 1;}
      //   if (nFirst == 1) {name_Histogram = "dPlusSecondDaughter"; discription_Histogram = "dPlus second daughter (cm);"; lowerBound = 0; upperBound = 1;}
      //   if (nFirst == 2) {name_Histogram = "d0Mother"; discription_Histogram = "dPlus mother (cm);"; lowerBound = 0; upperBound = 1;}
      //   if (nFirst == 3) {name_Histogram = "pointingAngleMother"; discription_Histogram = "pointing angle (Cos(#theta));"; lowerBound = -1; upperBound = 1;}
      //   if (nFirst == 4) {name_Histogram = "impactProduct"; discription_Histogram = "impact product (cm^{2});"; lowerBound = -0.01; upperBound = 0.01;}
      //   if (nFirst == 5) {name_Histogram = "impactProductXY"; discription_Histogram = "impact product XY (cm^{2});"; lowerBound = 0; upperBound = 0.5;}
      //   if (nFirst == 6) {name_Histogram = "vertexDistance"; discription_Histogram = "vertex distance between mother and primary vertex  (cm);"; lowerBound = 0; upperBound = 1;}
      //   if (nFirst == 7) {name_Histogram = "normDecayLength"; discription_Histogram = "Normalized decay length w.r.t primary vertex (cm);"; lowerBound = 0; upperBound = 50;}
      //   if (nFirst == 8) {name_Histogram = "pointingAngleMotherXY"; discription_Histogram = "pointing angle XY (Cos(#theta));"; lowerBound = -1; upperBound = 1;}
      //   if (nFirst == 9) {name_Histogram = "vertexDistanceXY"; discription_Histogram = "vertex distance between mother and primary vertex XY (cm);"; lowerBound = 0; upperBound = 1;}
      //   if (nFirst == 10) {name_Histogram = "normDecayLengthXY"; discription_Histogram = "Normalized decay length w.r.t primary vertex XY (cm);"; lowerBound = 0; upperBound = 50;}

      //   if (nSecond == 0) {name_Histogram += "dPlusFirstDaughter"; discription_Histogram += "dPlus first daughter (cm);"; lowerBoundTwo = 0; upperBoundTwo = 1;}
      //   if (nSecond == 1) {name_Histogram += "dPlusSecondDaughter"; discription_Histogram += "dPlus second daughter (cm);"; lowerBoundTwo = 0; upperBoundTwo = 1;}
      //   if (nSecond == 2) {name_Histogram += "d0Mother"; discription_Histogram += "dPlus mother (cm);"; lowerBoundTwo = 0; upperBoundTwo = 1;}
      //   if (nSecond == 3) {name_Histogram += "pointingAngleMother"; discription_Histogram += "pointing angle (Cos(#theta));"; lowerBoundTwo = -1; upperBoundTwo = 1;}
      //   if (nSecond == 4) {name_Histogram += "impactProduct"; discription_Histogram += "impact product (cm^{2});"; lowerBoundTwo = -0.01; upperBoundTwo = 0.01;}
      //   if (nSecond == 5) {name_Histogram += "impactProductXY"; discription_Histogram += "impact product XY (cm^{2});"; lowerBoundTwo = 0; upperBoundTwo = 0.5;}
      //   if (nSecond == 6) {name_Histogram += "vertexDistance"; discription_Histogram += "vertex distance between mother and primary vertex  (cm);"; lowerBoundTwo = 0; upperBoundTwo = 1;}
      //   if (nSecond == 7) {name_Histogram += "normDecayLength"; discription_Histogram += "Normalized decay length w.r.t primary vertex (cm);"; lowerBoundTwo = 0; upperBoundTwo = 50;}
      //   if (nSecond == 8) {name_Histogram += "_pointingAngleMotherXY"; discription_Histogram += "pointing angle XY (Cos(#theta));"; lowerBoundTwo = -1; upperBoundTwo = 1;}
      //   if (nSecond == 9) {name_Histogram += "_vertexDistanceXY"; discription_Histogram += "vertex distance between mother and primary vertex XY (cm);"; lowerBoundTwo = 0; upperBoundTwo = 1;}
      //   if (nSecond == 10) {name_Histogram += "_normDecayLengthXY"; discription_Histogram += "Normalized decay length w.r.t primary vertex XY (cm);"; lowerBoundTwo = 0; upperBoundTwo = 50;}

      //   name_Histogram += add_name;
      //   TH2F* histogram = new TH2F(name_Histogram.Data(), discription_Histogram.Data(), numberOfBins, lowerBound, upperBound, numberOfBinsTwo, lowerBoundTwo, upperBoundTwo);
      //   histogram->Sumw2();
      //   if (j % 2 == 0) histogram->SetLineColor(6);
      //   if (j % 2 == 1) histogram->SetLineColor(4);
      //   histogram->SetMarkerStyle(20);
      //   histogram->SetMarkerSize(0.6);
      //   histogram->SetMarkerColor(6);
      //   TH2F* histogram_Clone = (TH2F*)histogram->Clone();
      //   list2D->Add(histogram_Clone);
      //   fMotherHistogramArray2D[i][j][k] = histogram_Clone;

      //   nSecond++;
      //   if (nSecond > nVariables)
      //   {
      //     nFirst++;
      //     nSecond = nFirst + 1;
      //   }
      // }
    }

    TH1F * effectOfCuts = new TH1F("effectOfCuts", "Removal counter", 100, 0, 100);
    effectOfCuts->SetStats(kTRUE);
    effectOfCuts->GetXaxis()->SetTitle("Cut number");
    effectOfCuts->GetYaxis()->SetTitle("Particles cut");
    effectOfCuts->GetXaxis()->SetBinLabel(1, "total");
    for (Int_t i = 1; i < 100; ++i)
    {
      TString integerText = "";
      integerText += i;
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
      integerText += i;
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
    TH1F* hist_distance_vertex_from_real = new TH1F(name_distance_vertex_from_real.Data(), "Distance reconstructed vertex from real vertex; distance (cm); Entries", 500, 0, 0.5);
    hist_distance_vertex_from_real->Sumw2();
    hist_distance_vertex_from_real->SetLineColor(6);
    hist_distance_vertex_from_real->SetMarkerStyle(20);
    hist_distance_vertex_from_real->SetMarkerSize(0.6);
    hist_distance_vertex_from_real->SetMarkerColor(6);
    TH1F* histogram_distance_vertex_from_real = (TH1F*)hist_distance_vertex_from_real->Clone();
    listout->Add(histogram_distance_vertex_from_real);
    fMotherHistogramArrayExtra[i][4] = histogram_distance_vertex_from_real;

    TString name_distance_vertex_from_real_new = "distance_vertex_from_real_new";
    TH1F* hist_distance_vertex_from_real_new = new TH1F(name_distance_vertex_from_real_new.Data(), "Distance reconstructed vertex from real vertex; distance (cm); Entries", 500, 0, 0.5);
    hist_distance_vertex_from_real_new->Sumw2();
    hist_distance_vertex_from_real_new->SetLineColor(6);
    hist_distance_vertex_from_real_new->SetMarkerStyle(20);
    hist_distance_vertex_from_real_new->SetMarkerSize(0.6);
    hist_distance_vertex_from_real_new->SetMarkerColor(6);
    TH1F* histogram_distance_vertex_from_real_new = (TH1F*)hist_distance_vertex_from_real_new->Clone();
    listout->Add(histogram_distance_vertex_from_real_new);
    fMotherHistogramArrayExtra[i][5] = histogram_distance_vertex_from_real_new;

    TString name_momentum_resolution = "momentum_resolution";
    TH1F* hist_momentum_resolution = new TH1F(name_momentum_resolution.Data(), "Momentum resolution; difference between real and reconstructed momentum (GeV/#it{c}); Entries", 1000, 0, 1);
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
    fOutputB0Results->Add(listMassPlots);

    for (Int_t k = 0; k < fnPtBins + 3; ++k) {
      TString ptBinMother = "";
      if (k == 0) ptBinMother = "";
      if (k == 1) ptBinMother = "_ptbin_3_to_inf";
      if (k == 2) ptBinMother = "_ptbin_6_to_inf";
      if (k > 2) {ptBinMother += "_ptbin_"; ptBinMother += fPtBinLimits[k - 3]; ptBinMother += "_to_"; ptBinMother += fPtBinLimits[k - 2];}

      TString name_invariantMassMother = "invariantMassB0";
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

      TString name_deltainvariantMassMother = "deltainvariantMassB0";
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
  fOutputB0Results->Add(histogram_cutEffectBackground);
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
  fOutputB0Results->Add(histogram_cutEffectSignal);
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
  fOutputB0Results->Add(histogram_cutEffectUniqueBackground);
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
  fOutputB0Results->Add(histogram_cutEffectUniqueSignal);
  fResultsHistogramArray[13][1] = histogram_cutEffectUniqueSignal;

  TString name_totalITSBackground = "totalITSBackground";
  TH1F* hist_totalITSBackground = new TH1F(name_totalITSBackground.Data(), "Total nr. of ITS hits for the daughters; number [#]; Entries", 30, 0, 30);
  hist_totalITSBackground->Sumw2();
  hist_totalITSBackground->SetLineColor(6);
  hist_totalITSBackground->SetMarkerStyle(20);
  hist_totalITSBackground->SetMarkerSize(0.6);
  hist_totalITSBackground->SetMarkerColor(6);
  TH1F* histogram_totalITSBackground = (TH1F*)hist_totalITSBackground->Clone();
  fOutputB0Results->Add(histogram_totalITSBackground);
  fResultsHistogramArray[11][0] = histogram_totalITSBackground;

  TString name_totalITSSignal = "totalITSSignal";
  TH1F* hist_totalITSSignal = new TH1F(name_totalITSSignal.Data(), "Total nr. of ITS hits for the daughters; number [#]; Entries", 30, 0, 30);
  hist_totalITSSignal->Sumw2();
  hist_totalITSSignal->SetLineColor(6);
  hist_totalITSSignal->SetMarkerStyle(20);
  hist_totalITSSignal->SetMarkerSize(0.6);
  hist_totalITSSignal->SetMarkerColor(6);
  TH1F* histogram_totalITSSignal = (TH1F*)hist_totalITSSignal->Clone();
  fOutputB0Results->Add(histogram_totalITSSignal);
  fResultsHistogramArray[11][1] = histogram_totalITSSignal;

  TString name_totalTPCBackground = "totalTPCBackground";
  TH1F* hist_totalTPCBackground = new TH1F(name_totalTPCBackground.Data(), "Total nr. of TPC hits for the daughters; number [#]; Entries", 1000, 0, 1000);
  hist_totalTPCBackground->Sumw2();
  hist_totalTPCBackground->SetLineColor(6);
  hist_totalTPCBackground->SetMarkerStyle(20);
  hist_totalTPCBackground->SetMarkerSize(0.6);
  hist_totalTPCBackground->SetMarkerColor(6);
  TH1F* histogram_totalTPCBackground = (TH1F*)hist_totalTPCBackground->Clone();
  fOutputB0Results->Add(histogram_totalTPCBackground);
  fResultsHistogramArray[11][2] = histogram_totalTPCBackground;

  TString name_totalTPCSignal = "totalTPCSignal";
  TH1F* hist_totalTPCSignal = new TH1F(name_totalTPCSignal.Data(), "Total nr. of TPC hits for the daughters; number [#]; Entries", 1000, 0, 1000);
  hist_totalTPCSignal->Sumw2();
  hist_totalTPCSignal->SetLineColor(6);
  hist_totalTPCSignal->SetMarkerStyle(20);
  hist_totalTPCSignal->SetMarkerSize(0.6);
  hist_totalTPCSignal->SetMarkerColor(6);
  TH1F* histogram_totalTPCSignal = (TH1F*)hist_totalTPCSignal->Clone();
  fOutputB0Results->Add(histogram_totalTPCSignal);
  fResultsHistogramArray[11][3] = histogram_totalTPCSignal;

  TString name_totalSigmaPIDBackground = "totalSigmaPIDBackground";
  TH1F* hist_totalSigmaPIDBackground = new TH1F(name_totalSigmaPIDBackground.Data(), "Total sigma of TPC and TOF PID for the daughters; number [#]; Entries", 1000, 0, 100);
  hist_totalSigmaPIDBackground->Sumw2();
  hist_totalSigmaPIDBackground->SetLineColor(6);
  hist_totalSigmaPIDBackground->SetMarkerStyle(20);
  hist_totalSigmaPIDBackground->SetMarkerSize(0.6);
  hist_totalSigmaPIDBackground->SetMarkerColor(6);
  TH1F* histogram_totalSigmaPIDBackground = (TH1F*)hist_totalSigmaPIDBackground->Clone();
  fOutputB0Results->Add(histogram_totalSigmaPIDBackground);
  fResultsHistogramArray[11][4] = histogram_totalSigmaPIDBackground;

  TString name_totalSigmaPIDSignal = "totalSigmaPIDSignal";
  TH1F* hist_totalSigmaPIDSignal = new TH1F(name_totalSigmaPIDSignal.Data(), "Total sigma of TPC and TOF PID for the daughters; number [#]; Entries", 1000, 0, 100);
  hist_totalSigmaPIDSignal->Sumw2();
  hist_totalSigmaPIDSignal->SetLineColor(6);
  hist_totalSigmaPIDSignal->SetMarkerStyle(20);
  hist_totalSigmaPIDSignal->SetMarkerSize(0.6);
  hist_totalSigmaPIDSignal->SetMarkerColor(6);
  TH1F* histogram_totalSigmaPIDSignal = (TH1F*)hist_totalSigmaPIDSignal->Clone();
  fOutputB0Results->Add(histogram_totalSigmaPIDSignal);
  fResultsHistogramArray[11][5] = histogram_totalSigmaPIDSignal;

  // Optional code to check invariant mass of daughter combinations, to be added for this decay
  // for (int i = 0; i < 3; ++i)
  // {
  //   TString name_Histogram;
  //   TString discription_Histogram;
  //   if (i == 0) {name_Histogram = "invmassDPlusPionB0Pion"; discription_Histogram = "Invariant mass DPlus Pion and B0 Pion; inv. mass (GeV/#it{c}^{2}); Entries";}
  //   if (i == 1) {name_Histogram = "invmassDPlusKaonB0Pion"; discription_Histogram = "Invariant mass DPlus Kaon and B0 Pion; inv. mass (GeV/#it{c}^{2}); Entries";}
  //   if (i == 2) {name_Histogram = "invmassDPlusPionDPlusKaonB0Pion"; discription_Histogram = "Invariant mass DPlus Pion, DPlus Kaon, and B0 Pion; inv. mass (GeV/#it{c}^{2}); Entries";}

  //   for (int j = 0; j < 2; ++j)
  //   {
  //     TString add_name = "";
  //     if (j == 1) add_name = "_Signal";
  //     name_Histogram += add_name;
  //     TH1F* histogram = new TH1F(name_Histogram.Data(), discription_Histogram.Data(), 1000, 0, 30);
  //     histogram->Sumw2();
  //     if (j % 2 == 0) histogram->SetLineColor(6);
  //     if (j % 2 == 1) histogram->SetLineColor(4);
  //     histogram->SetMarkerStyle(20);
  //     histogram->SetMarkerSize(0.6);
  //     if (j % 2 == 0) histogram->SetMarkerColor(6);
  //     if (j % 2 == 1) histogram->SetMarkerColor(4);
  //     TH1F* histogram_Clone = (TH1F*)histogram->Clone();
  //     fOutputB0Results->Add(histogram_Clone);
  //     fResultsHistogramArray[12][2*i + j] = histogram_Clone;
  //   }
  // }

  return;
}
//-------------------------------------------------------------------------------------
AliAODVertex* AliAnalysisTaskSEB0toDPi::RecalculateVertex(const AliVVertex *primary, TObjArray *tracks, Double_t bField, Double_t dispersion, Bool_t optUseFitter, Bool_t optPropagate, Bool_t optUseDiamondConstraint, Int_t nprongs) {
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
  // Int_t nprongs = 2; //tracks->GetEntriesFast();
  vertexAOD = new AliAODVertex(pos, cov, chi2perNDF, 0x0, -1, AliAODVertex::kUndef, nprongs);

  return vertexAOD;
}
//-------------------------------------------------------------------------------------
void AliAnalysisTaskSEB0toDPi::B0toDPiSignalTracksInMC(TClonesArray * mcTrackArray, AliAODEvent*  /*aodevent*/, TMatrix * B0toDPiLabelMatrix, TList *listout) {

  TMatrix &particleMatrix = *B0toDPiLabelMatrix;
  for (Int_t i = 0; i < mcTrackArray->GetEntriesFast(); i++) {

    Int_t mcLabelKaonDPlus = 0;
    Int_t mcLabelFirstPionDPlus = 0;
    Int_t mcLabelSecondPionDPlus = 0;
    Int_t mcLabelPionB0 = 0;
    Int_t mcLabelDPlus = 0;
    Int_t mcLabelB0 = 0;

    Double_t ptMC[6] = {0.0};
    Double_t yMC[6] = {0.0};
    Double_t pseudoYMC[6] = {0.0};

    Bool_t mckaonDPlusPresent = kFALSE;
    Bool_t mcFirstPionDPlusPresent = kFALSE;
    Bool_t mcSecondPionDPlusPresent = kFALSE;
    Bool_t mcPionB0Present = kFALSE;

    AliAODMCParticle *mcTrackParticle = dynamic_cast< AliAODMCParticle*>(mcTrackArray->At(i));
    if (!mcTrackParticle) {std::cout << "no particle" << std::endl; continue;}
    Int_t pdgCodeMC = TMath::Abs(mcTrackParticle->GetPdgCode());

    if (pdgCodeMC == 511)
    { //if the track is a B0 we look at its daughters

      mcLabelB0 = i;
      Int_t nDaughterB0 = mcTrackParticle->GetNDaughters();
      ptMC[0] = mcTrackParticle->Pt();
      yMC[0] = mcTrackParticle->Y();
      pseudoYMC[0] = mcTrackParticle->Eta();

      // "B0_in_analysis";
      ((TH1F*)fResultsHistogramArray[3][0])->Fill(0);

      if (nDaughterB0 == 2)
      {
        ((TH1F*)fResultsHistogramArray[3][0])->Fill(3);
        for (Int_t iDaughterB0 = 0; iDaughterB0 < 2; iDaughterB0++)
        {
          AliAODMCParticle* daughterB0 = (AliAODMCParticle*)mcTrackArray->At(mcTrackParticle->GetDaughterLabel(0) + iDaughterB0);
          if (!daughterB0) { break;}
          Int_t pdgCodeDaughterB0 = TMath::Abs(daughterB0->GetPdgCode());
          ((TH1F*)fResultsHistogramArray[3][0])->Fill(4);
          if (pdgCodeDaughterB0 == 211) //if the track is a pion we save its monte carlo label
          {
            ((TH1F*)fResultsHistogramArray[3][0])->Fill(5);
            mcLabelPionB0 = mcTrackParticle->GetDaughterLabel(0) + iDaughterB0;
            mcPionB0Present = kTRUE;
            ptMC[1] = daughterB0->Pt();
            yMC[1] = daughterB0->Y();
            pseudoYMC[1] = daughterB0->Eta();

          } else if (pdgCodeDaughterB0 == 411) //if the track is a DPlus we look at its daughters
          {
            ((TH1F*)fResultsHistogramArray[3][0])->Fill(7);
            mcLabelDPlus = mcTrackParticle->GetDaughterLabel(0) + iDaughterB0;
            Int_t nDaughterDPlus = daughterB0->GetNDaughters();
            ptMC[2] = daughterB0->Pt();
            yMC[2] = daughterB0->Y();
            pseudoYMC[2] = daughterB0->Eta();

            if (nDaughterDPlus == 3)
            {
              Int_t numberOfDPlusPions = 0;
              ((TH1F*)fResultsHistogramArray[3][0])->Fill(8);
              for (Int_t iDaughterDPlus = 0; iDaughterDPlus < 3; iDaughterDPlus++)
              {
                AliAODMCParticle* daughterDPlus = (AliAODMCParticle*)mcTrackArray->At(daughterB0->GetDaughterLabel(0) + iDaughterDPlus);
                if (!daughterDPlus) break;
                Int_t pdgCodeDaughterDPlus = TMath::Abs(daughterDPlus->GetPdgCode());

                if (pdgCodeDaughterDPlus == 211) //if the track is a pion we save its monte carlo label
                {
                  numberOfDPlusPions++;
                  if (numberOfDPlusPions == 1)
                  {
                    ((TH1F*)fResultsHistogramArray[3][0])->Fill(9);
                    mcLabelFirstPionDPlus = daughterB0->GetDaughterLabel(0) + iDaughterDPlus;
                    ptMC[3] = daughterDPlus->Pt();
                    yMC[3] = daughterDPlus->Y();
                    pseudoYMC[3] = daughterDPlus->Eta();
                    mcFirstPionDPlusPresent = kTRUE;
                  }
                  if (numberOfDPlusPions == 2)
                  {
                    ((TH1F*)fResultsHistogramArray[3][0])->Fill(10);
                    mcLabelSecondPionDPlus = daughterB0->GetDaughterLabel(0) + iDaughterDPlus;
                    ptMC[4] = daughterDPlus->Pt();
                    yMC[4] = daughterDPlus->Y();
                    pseudoYMC[4] = daughterDPlus->Eta();
                    mcSecondPionDPlusPresent = kTRUE;
                  }
                } else if (pdgCodeDaughterDPlus == 321) //if the track is a kaon we save its monte carlo label
                {
                  ((TH1F*)fResultsHistogramArray[3][0])->Fill(11);
                  mcLabelKaonDPlus = daughterB0->GetDaughterLabel(0) + iDaughterDPlus;
                  mckaonDPlusPresent = kTRUE;
                  ptMC[5] = daughterDPlus->Pt();
                  yMC[5] = daughterDPlus->Y();
                  pseudoYMC[5] = daughterDPlus->Eta();
                } else break;
              }
            }
          } else break;
        }
      }
    }
    if (mckaonDPlusPresent && mcFirstPionDPlusPresent && mcSecondPionDPlusPresent && mcPionB0Present) {

      //  "B0_in_analysis";
      ((TH1F*)fResultsHistogramArray[3][0])->Fill(1);


      for (Int_t j = 0; j < fnPtBins; ++j)
      {
        // "B0_per_bin";
        if (fPtBinLimits[j] < ptMC[0] && ptMC[0] < fPtBinLimits[j + 1]) {((TH1F*)fResultsHistogramArray[3][1])->Fill(j); break;}
      }

      ((TH1F*)fResultsHistogramArray[0][0])->Fill(ptMC[0]);
      ((TH1F*)fResultsHistogramArray[0][1])->Fill(ptMC[1]);
      ((TH1F*)fResultsHistogramArray[0][2])->Fill(ptMC[2]);
      ((TH1F*)fResultsHistogramArray[0][3])->Fill(ptMC[3]);
      ((TH1F*)fResultsHistogramArray[0][4])->Fill(ptMC[4]);
      ((TH1F*)fResultsHistogramArray[0][5])->Fill(ptMC[5]);

      ((TH1F*)fResultsHistogramArray[0][7])->Fill(yMC[0]);
      ((TH1F*)fResultsHistogramArray[0][8])->Fill(yMC[1]);
      ((TH1F*)fResultsHistogramArray[0][9])->Fill(yMC[2]);
      ((TH1F*)fResultsHistogramArray[0][10])->Fill(yMC[3]);
      ((TH1F*)fResultsHistogramArray[0][11])->Fill(yMC[4]);
      ((TH1F*)fResultsHistogramArray[0][12])->Fill(yMC[5]);

      ((TH1F*)fResultsHistogramArray[0][14])->Fill(pseudoYMC[0]);
      ((TH1F*)fResultsHistogramArray[0][15])->Fill(pseudoYMC[1]);
      ((TH1F*)fResultsHistogramArray[0][16])->Fill(pseudoYMC[2]);
      ((TH1F*)fResultsHistogramArray[0][17])->Fill(pseudoYMC[3]);
      ((TH1F*)fResultsHistogramArray[0][18])->Fill(pseudoYMC[4]);
      ((TH1F*)fResultsHistogramArray[0][19])->Fill(pseudoYMC[5]);

      // We check if the tracks are in acceptance
      if (ptMC[1] < 0.1 || TMath::Abs(pseudoYMC[1]) > 0.9 ) continue;
      if (ptMC[3] < 0.1 || TMath::Abs(pseudoYMC[3]) > 0.9 ) continue;
      if (ptMC[4] < 0.1 || TMath::Abs(pseudoYMC[4]) > 0.9 ) continue;
      if (ptMC[5] < 0.1 || TMath::Abs(pseudoYMC[5]) > 0.9 ) continue;


      // We check if the B0 is in the fiducial region
      if (TMath::Abs(yMC[0]) > 0.8) continue;

      Int_t rows = B0toDPiLabelMatrix->GetNrows();

      B0toDPiLabelMatrix->ResizeTo(rows + 1, 6);
      particleMatrix(rows, 0) = mcLabelKaonDPlus;
      particleMatrix(rows, 1) = mcLabelFirstPionDPlus;
      particleMatrix(rows, 2) = mcLabelSecondPionDPlus;
      particleMatrix(rows, 3) = mcLabelPionB0;
      particleMatrix(rows, 4) = mcLabelDPlus;
      particleMatrix(rows, 5) = mcLabelB0;

      // "B0_in_analysis";
      ((TH1F*)fResultsHistogramArray[3][0])->Fill(2);

      for (Int_t j = 0; j < fnPtBins; ++j)
      {
        //  "B0_per_bin_in_Acc";
        if (fPtBinLimits[j] < ptMC[0] && ptMC[0] < fPtBinLimits[j + 1]) {((TH1F*)fResultsHistogramArray[3][2])->Fill(j); break;}
      }
    }
  }


  return;
}
//-------------------------------------------------------------------------------------
void AliAnalysisTaskSEB0toDPi::DaughterSelection(AliAODEvent* aodEvent, AliAODVertex *primaryVertex, Double_t bz, TClonesArray * mcTrackArray, TMatrix * B0toDPiLabelMatrix, AliAODMCHeader * header, Int_t daughterType, std::vector<Int_t> * trackVector) {

  //we keep track of the number of particles we could use and how many we actually use after cuts
  Int_t numberofparticles = 0;
  Int_t numberofparticlesused = 0;

  Int_t particleType = 2;
  if (daughterType == 0) particleType = 3;

  //we loop over all tracks in the event
  for (Int_t i = 0; i < aodEvent->GetNumberOfTracks(); i++) {
    AliAODTrack* aodTrack = dynamic_cast<AliAODTrack*>(aodEvent->GetTrack(i));
    if (!aodTrack) AliFatal("Not a standard AOD");

    //quick quality cut
    if (aodTrack->GetITSNcls() < 1) continue;
    if (aodTrack->GetTPCNcls() < 1) continue;
    if (aodTrack->GetStatus()&AliESDtrack::kITSpureSA) continue;
    if (!(aodTrack->GetStatus()&AliESDtrack::kITSin)) continue;
    if (aodTrack->GetID() < 0) continue;
    Double_t covtest[21];
    if (!aodTrack->GetCovarianceXYZPxPyPz(covtest)) continue;

    Double_t pos[3], cov[6];
    primaryVertex->GetXYZ(pos);
    primaryVertex->GetCovarianceMatrix(cov);
    const AliESDVertex vESD(pos, cov, 100., 100);
    if (!fCuts->IsDaughterSelected(aodTrack, &vESD, fCuts->GetTrackCuts(), aodEvent)) continue;

    Int_t mcLabelParticle = -1;
    mcLabelParticle = aodTrack->GetLabel();

    numberofparticles++;

    //we check if the particle is a signal track
    Bool_t isDesiredCandidate = kFALSE;
    Int_t labelPosition1 = 0;
    Int_t labelPosition2 = 0;
    if (daughterType == 1) { labelPosition1 = 1; labelPosition2 = 2;}
    if (daughterType == 2) { labelPosition1 = 3;}
    if (fUseMCInfo) {
      TMatrix &particleMatrix = *B0toDPiLabelMatrix;
      for (Int_t k = 0; k < B0toDPiLabelMatrix->GetNrows(); ++k) {
        if (mcLabelParticle == (Int_t)particleMatrix(k, labelPosition1)) {
          isDesiredCandidate = kTRUE;
          break;
        }
        if (labelPosition2 > 0)
        {
          if (mcLabelParticle == (Int_t)particleMatrix(k, labelPosition2)) {
            isDesiredCandidate = kTRUE;
            break;
          }
        }
      }
    }

    if (fUseMCInfo) {
      if (IsTrackInjected(aodTrack, header, mcTrackArray) && !isDesiredCandidate && fQuickSignalAnalysis == 2) continue;
    }

    Int_t histType = 0;
    FillDaughterHistograms(aodTrack, daughterType, histType, 0.0, primaryVertex, bz, isDesiredCandidate, mcTrackArray);

    AliExternalTrackParam particleTrack;
    particleTrack.CopyFromVTrack(aodTrack);
    Double_t d0[2], covd0[3];
    particleTrack.PropagateToDCA(primaryVertex, bz, 100., d0, covd0);

    //we apply a number of cuts on the particle
    Bool_t bCut = kFALSE;

    //we apply cut from the cutfile
    if (!(fCuts->SelectPID(aodTrack, particleType))) {
      if (isDesiredCandidate) {
        ((TH1F*)fDaughterHistogramArrayExtra[daughterType][1])->Fill(2);
      } else ((TH1F*)fDaughterHistogramArrayExtra[daughterType][0])->Fill(2);
      bCut = kTRUE;
    }

    if (aodTrack->GetITSNcls() < fCuts->GetMinITSNclsDaughterType(daughterType)) {
      if (isDesiredCandidate) {
        ((TH1F*)fDaughterHistogramArrayExtra[daughterType][1])->Fill(3);
      } else ((TH1F*)fDaughterHistogramArrayExtra[daughterType][0])->Fill(3);
      bCut = kTRUE;
    }

    if (aodTrack->GetTPCNcls() < fCuts->GetMinTPCNclsDaughterType(daughterType)) {
      if (isDesiredCandidate) {
        ((TH1F*)fDaughterHistogramArrayExtra[daughterType][1])->Fill(4);
      } else ((TH1F*)fDaughterHistogramArrayExtra[daughterType][0])->Fill(4);
      bCut = kTRUE;
    }

    if (fCuts->UseITSRefitDaughterType(daughterType) == kTRUE) {
      if (!(aodTrack->GetStatus()&AliESDtrack::kITSrefit)) {
        if (isDesiredCandidate) {
          ((TH1F*)fDaughterHistogramArrayExtra[daughterType][1])->Fill(5);
        } else ((TH1F*)fDaughterHistogramArrayExtra[daughterType][0])->Fill(5);
        bCut = kTRUE;
      }
    }

    if (fCuts->UseTPCRefitDaughterType(daughterType) == kTRUE) {
      if ((!(aodTrack->GetStatus()&AliESDtrack::kTPCrefit))) {
        if (isDesiredCandidate) {
          ((TH1F*)fDaughterHistogramArrayExtra[daughterType][1])->Fill(6);
        } else ((TH1F*)fDaughterHistogramArrayExtra[daughterType][0])->Fill(6);
        bCut = kTRUE;
      }
    }

    if (fCuts->UseFilterBitDaughterType(daughterType) == kTRUE) {
      if (!(aodTrack->TestFilterMask(BIT(fCuts->GetFilterBitDaughterType(daughterType))))) {
        if (isDesiredCandidate) {
          ((TH1F*)fDaughterHistogramArrayExtra[daughterType][1])->Fill(7);
        } else ((TH1F*)fDaughterHistogramArrayExtra[daughterType][0])->Fill(7);
        bCut = kTRUE;
      }
    }


    if (aodTrack->Pt() < fCuts->GetMinPtDaughterType(daughterType)) {
      if (isDesiredCandidate) {
        ((TH1F*)fDaughterHistogramArrayExtra[daughterType][1])->Fill(8);
      } else ((TH1F*)fDaughterHistogramArrayExtra[daughterType][0])->Fill(8);
      bCut = kTRUE;
    }


    if (TMath::Abs(d0[0]) < fCuts->GetMind0DaughterType(daughterType)) {
      if (isDesiredCandidate) {
        ((TH1F*)fDaughterHistogramArrayExtra[daughterType][1])->Fill(12);
      } else ((TH1F*)fDaughterHistogramArrayExtra[daughterType][0])->Fill(12);
      bCut = kTRUE;
    }

    if (TMath::Abs(d0[0] / TMath::Sqrt(covd0[0])) < fCuts->GetMinNormd0DaughterType(daughterType)) {
      if (isDesiredCandidate) {
        ((TH1F*)fDaughterHistogramArrayExtra[daughterType][1])->Fill(13);
      } else ((TH1F*)fDaughterHistogramArrayExtra[daughterType][0])->Fill(13);
      bCut = kTRUE;
    }

    if (TMath::Abs(aodTrack->Eta()) > fCuts->GetMaxAbsEtaDaughterType(daughterType)) {
      if (isDesiredCandidate) {
        ((TH1F*)fDaughterHistogramArrayExtra[daughterType][1])->Fill(9);
      } else ((TH1F*)fDaughterHistogramArrayExtra[daughterType][0])->Fill(9);
      bCut = kTRUE;
    }

    Bool_t bHardSelectionArrayITS[7] = {kFALSE};
    fCuts->GetHardSelectionArrayITSDaughterType(daughterType, bHardSelectionArrayITS);
    Bool_t bSoftSelectionArrayITS[7] = {kFALSE};
    fCuts->GetSoftSelectionArrayITSDaughterType(daughterType, bSoftSelectionArrayITS);

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
    if (nCounterSoftSelection < fCuts->GetNSoftITSCutDaughterType(daughterType)) bSoftITSPass = kFALSE;

    if (!bHardITSPass) {
      if (isDesiredCandidate) {
        ((TH1F*)fDaughterHistogramArrayExtra[daughterType][1])->Fill(10);
      } else ((TH1F*)fDaughterHistogramArrayExtra[daughterType][0])->Fill(10);
      bCut = kTRUE;
    }

    if (!bSoftITSPass) {
      if (isDesiredCandidate) {
        ((TH1F*)fDaughterHistogramArrayExtra[daughterType][1])->Fill(11);
      } else ((TH1F*)fDaughterHistogramArrayExtra[daughterType][0])->Fill(11);
      bCut = kTRUE;
    }

    if (!isDesiredCandidate && fQuickSignalAnalysis == 1) bCut = kTRUE;

    if (bCut) {
      if (isDesiredCandidate) {
        ((TH1F*)fDaughterHistogramArrayExtra[daughterType][1])->Fill(0);
      } else ((TH1F*)fDaughterHistogramArrayExtra[daughterType][0])->Fill(0);
      continue;
    }

    //we fill histograms with track information of the tracks that pass the cuts
    histType = 2;
    FillDaughterHistograms(aodTrack, daughterType, histType, 0.0, primaryVertex, bz, isDesiredCandidate, mcTrackArray);

    trackVector->push_back(i);
    numberofparticlesused++;
  }

  ((TH1F*)fDaughterHistogramArray[daughterType][0][10])->Fill(numberofparticles);
  ((TH1F*)fDaughterHistogramArray[daughterType][1][10])->Fill(numberofparticlesused);
  return;
}
//-------------------------------------------------------------------------------------
Bool_t AliAnalysisTaskSEB0toDPi::DPlusDaughterSelection(AliAODTrack*  aodTrack, AliAODEvent* aodEvent, AliAODVertex *primaryVertex, Double_t bz, TClonesArray * mcTrackArray, TMatrix * B0toDPiLabelMatrix, AliAODMCHeader * header, Int_t daughterType) {

  Int_t particleType = 2;
  if (daughterType == 0) particleType = 3;

  //quick quality cut
  if (aodTrack->GetITSNcls() < 1) return kFALSE;
  if (aodTrack->GetTPCNcls() < 1) return kFALSE;
  if (aodTrack->GetStatus()&AliESDtrack::kITSpureSA) return kFALSE;
  if (!(aodTrack->GetStatus()&AliESDtrack::kITSin)) return kFALSE;
  if (aodTrack->GetID() < 0) return kFALSE;
  Double_t covtest[21];
  if (!aodTrack->GetCovarianceXYZPxPyPz(covtest)) return kFALSE;

  Double_t pos[3], cov[6];
  primaryVertex->GetXYZ(pos);
  primaryVertex->GetCovarianceMatrix(cov);
  const AliESDVertex vESD(pos, cov, 100., 100);
  if (!fCuts->IsDaughterSelected(aodTrack, &vESD, fCuts->GetTrackCuts(), aodEvent)) return kFALSE;

  Int_t mcLabelParticle = -1;
  mcLabelParticle = aodTrack->GetLabel();

  //we check if the particle is a signal track
  Bool_t isDesiredCandidate = kFALSE;
  Int_t labelPosition1 = 0;
  Int_t labelPosition2 = 0;
  if (daughterType == 1) { labelPosition1 = 1; labelPosition2 = 2;}
  if (daughterType == 2) { labelPosition1 = 3;}
  if (fUseMCInfo) {
    TMatrix &particleMatrix = *B0toDPiLabelMatrix;
    for (Int_t k = 0; k < B0toDPiLabelMatrix->GetNrows(); ++k) {
      if (mcLabelParticle == (Int_t)particleMatrix(k, labelPosition1)) {
        isDesiredCandidate = kTRUE;
        break;
      }
      if (labelPosition2 > 0)
      {
        if (mcLabelParticle == (Int_t)particleMatrix(k, labelPosition2)) {
          isDesiredCandidate = kTRUE;
          break;
        }
      }
    }
  }

  if (fUseMCInfo) {
    if (IsTrackInjected(aodTrack, header, mcTrackArray) && !isDesiredCandidate && fQuickSignalAnalysis == 2) return kFALSE;
  }

  Int_t histType = 0;
  FillDaughterHistograms(aodTrack, daughterType, histType, 0.0, primaryVertex, bz, isDesiredCandidate, mcTrackArray);

  AliExternalTrackParam particleTrack;
  particleTrack.CopyFromVTrack(aodTrack);
  Double_t d0[2], covd0[3];
  particleTrack.PropagateToDCA(primaryVertex, bz, 100., d0, covd0);

  //we apply a number of cuts on the particle
  Bool_t bCut = kFALSE;

  //we apply cut from the cutfile
  if (!(fCuts->SelectPID(aodTrack, particleType))) {
    if (isDesiredCandidate) {
      ((TH1F*)fDaughterHistogramArrayExtra[daughterType][1])->Fill(2);
    } else ((TH1F*)fDaughterHistogramArrayExtra[daughterType][0])->Fill(2);
    bCut = kTRUE;
  }

  if (aodTrack->GetITSNcls() < fCuts->GetMinITSNclsDaughterType(daughterType)) {
    if (isDesiredCandidate) {
      ((TH1F*)fDaughterHistogramArrayExtra[daughterType][1])->Fill(3);
    } else ((TH1F*)fDaughterHistogramArrayExtra[daughterType][0])->Fill(3);
    bCut = kTRUE;
  }

  if (aodTrack->GetTPCNcls() < fCuts->GetMinTPCNclsDaughterType(daughterType)) {
    if (isDesiredCandidate) {
      ((TH1F*)fDaughterHistogramArrayExtra[daughterType][1])->Fill(4);
    } else ((TH1F*)fDaughterHistogramArrayExtra[daughterType][0])->Fill(4);
    bCut = kTRUE;
  }

  if (fCuts->UseITSRefitDaughterType(daughterType) == kTRUE) {
    if (!(aodTrack->GetStatus()&AliESDtrack::kITSrefit)) {
      if (isDesiredCandidate) {
        ((TH1F*)fDaughterHistogramArrayExtra[daughterType][1])->Fill(5);
      } else ((TH1F*)fDaughterHistogramArrayExtra[daughterType][0])->Fill(5);
      bCut = kTRUE;
    }
  }

  if (fCuts->UseTPCRefitDaughterType(daughterType) == kTRUE) {
    if ((!(aodTrack->GetStatus()&AliESDtrack::kTPCrefit))) {
      if (isDesiredCandidate) {
        ((TH1F*)fDaughterHistogramArrayExtra[daughterType][1])->Fill(6);
      } else ((TH1F*)fDaughterHistogramArrayExtra[daughterType][0])->Fill(6);
      bCut = kTRUE;
    }
  }

  if (fCuts->UseFilterBitDaughterType(daughterType) == kTRUE) {
    if (!(aodTrack->TestFilterMask(BIT(fCuts->GetFilterBitDaughterType(daughterType))))) {
      if (isDesiredCandidate) {
        ((TH1F*)fDaughterHistogramArrayExtra[daughterType][1])->Fill(7);
      } else ((TH1F*)fDaughterHistogramArrayExtra[daughterType][0])->Fill(7);
      bCut = kTRUE;
    }
  }


  if (aodTrack->Pt() < fCuts->GetMinPtDaughterType(daughterType)) {
    if (isDesiredCandidate) {
      ((TH1F*)fDaughterHistogramArrayExtra[daughterType][1])->Fill(8);
    } else ((TH1F*)fDaughterHistogramArrayExtra[daughterType][0])->Fill(8);
    bCut = kTRUE;
  }


  if (TMath::Abs(d0[0]) < fCuts->GetMind0DaughterType(daughterType)) {
    if (isDesiredCandidate) {
      ((TH1F*)fDaughterHistogramArrayExtra[daughterType][1])->Fill(12);
    } else ((TH1F*)fDaughterHistogramArrayExtra[daughterType][0])->Fill(12);
    bCut = kTRUE;
  }

  if (TMath::Abs(d0[0] / TMath::Sqrt(covd0[0])) < fCuts->GetMinNormd0DaughterType(daughterType)) {
    if (isDesiredCandidate) {
      ((TH1F*)fDaughterHistogramArrayExtra[daughterType][1])->Fill(13);
    } else ((TH1F*)fDaughterHistogramArrayExtra[daughterType][0])->Fill(13);
    bCut = kTRUE;
  }

  if (TMath::Abs(aodTrack->Eta()) > fCuts->GetMaxAbsEtaDaughterType(daughterType)) {
    if (isDesiredCandidate) {
      ((TH1F*)fDaughterHistogramArrayExtra[daughterType][1])->Fill(9);
    } else ((TH1F*)fDaughterHistogramArrayExtra[daughterType][0])->Fill(9);
    bCut = kTRUE;
  }

  Bool_t bHardSelectionArrayITS[7] = {kFALSE};
  fCuts->GetHardSelectionArrayITSDaughterType(daughterType, bHardSelectionArrayITS);
  Bool_t bSoftSelectionArrayITS[7] = {kFALSE};
  fCuts->GetSoftSelectionArrayITSDaughterType(daughterType, bSoftSelectionArrayITS);

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
  if (nCounterSoftSelection < fCuts->GetNSoftITSCutDaughterType(daughterType)) bSoftITSPass = kFALSE;

  if (!bHardITSPass) {
    if (isDesiredCandidate) {
      ((TH1F*)fDaughterHistogramArrayExtra[daughterType][1])->Fill(10);
    } else ((TH1F*)fDaughterHistogramArrayExtra[daughterType][0])->Fill(10);
    bCut = kTRUE;
  }

  if (!bSoftITSPass) {
    if (isDesiredCandidate) {
      ((TH1F*)fDaughterHistogramArrayExtra[daughterType][1])->Fill(11);
    } else ((TH1F*)fDaughterHistogramArrayExtra[daughterType][0])->Fill(11);
    bCut = kTRUE;
  }

  if (!isDesiredCandidate && fQuickSignalAnalysis == 1) bCut = kTRUE;

  if (bCut) {
    if (isDesiredCandidate) {
      ((TH1F*)fDaughterHistogramArrayExtra[daughterType][1])->Fill(0);
    } else ((TH1F*)fDaughterHistogramArrayExtra[daughterType][0])->Fill(0);
    return kFALSE;
  }

  //we fill histograms with track information of the tracks that pass the cuts
  histType = 2;
  FillDaughterHistograms(aodTrack, daughterType, histType, 0.0, primaryVertex, bz, isDesiredCandidate, mcTrackArray);

  // ((TH1F*)fDaughterHistogramArray[daughterType][0][10])->Fill(numberofparticles);
  // ((TH1F*)fDaughterHistogramArray[daughterType][1][10])->Fill(numberofparticlesused);
  return kTRUE;
}
//-------------------------------------------------------------------------------------
void AliAnalysisTaskSEB0toDPi::DPlusSelection(AliAODEvent* aodEvent, AliAODVertex *primaryVertex, Double_t bz, TClonesArray * mcTrackArray, TMatrix * B0toDPiLabelMatrix, TClonesArray * DPlusTrackArray, AliAODMCHeader * header) {

  AliAnalysisVertexingHF *vHF = new AliAnalysisVertexingHF();

  //next we loop over all the DPlus candidates
  for (Int_t j = 0; j < DPlusTrackArray->GetEntriesFast(); j++)
  {

    //we get the track of the DPlus
    AliAODRecoDecayHF3Prong * trackDPlus = (AliAODRecoDecayHF3Prong*)(DPlusTrackArray->At(j));
    if (!trackDPlus) {std::cout << "found none" << std::endl; continue;}
    if (trackDPlus == nullptr) {std::cout << "found nullptr" << std::endl; continue;}

    if (!(vHF->FillRecoCand(aodEvent, trackDPlus))) //Fill the data members of the candidate only if they are empty.
    {
      fCEvents->Fill(12); //monitor how often this fails
      continue;
    }

    AliAODTrack * trackDPlusFirstPion = (AliAODTrack*)(trackDPlus->GetDaughter(0));
    AliAODTrack * trackDPlusKaon = (AliAODTrack*)(trackDPlus->GetDaughter(1));
    AliAODTrack * trackDPlusSecondPion = (AliAODTrack*)(trackDPlus->GetDaughter(2));
    if (!DPlusDaughterSelection(trackDPlusFirstPion, aodEvent, primaryVertex, bz, mcTrackArray, B0toDPiLabelMatrix, header, 1)) continue;
    if (!DPlusDaughterSelection(trackDPlusKaon, aodEvent, primaryVertex, bz, mcTrackArray, B0toDPiLabelMatrix, header, 0)) continue;
    if (!DPlusDaughterSelection(trackDPlusSecondPion, aodEvent, primaryVertex, bz, mcTrackArray, B0toDPiLabelMatrix, header, 1)) continue;

    //we check if the IDs of the tracks are different
    UShort_t idProng0 = trackDPlusKaon->GetID();
    UShort_t idProng1 = trackDPlusFirstPion->GetID();
    UShort_t idProng2 = trackDPlusSecondPion->GetID();

    if (idProng0 == idProng1 || idProng0 == idProng2 || idProng1 == idProng2) continue;

    //we check if the particles have the correct charge
    if (trackDPlusKaon->Charge() == -1 && !(trackDPlusFirstPion->Charge() ==  1 && trackDPlusSecondPion->Charge() ==  1)) continue;
    if (trackDPlusKaon->Charge() ==  1 && !(trackDPlusFirstPion->Charge() == -1 && trackDPlusSecondPion->Charge() == -1)) continue;   

    AliAODVertex *vertexMother = (AliAODVertex*)trackDPlus->GetSecondaryVtx();

    //we check if the track is a desired candidate
    Int_t pdgCodeMother = -1;
    Bool_t isDesiredCandidate = kFALSE;
    Int_t motherType, histType;
    motherType = 0;
    Int_t mcLabelDPlus = -1;

    if (fUseMCInfo)
    {
      mcLabelDPlus = MatchCandidateToMonteCarlo(411, trackDPlus, mcTrackArray, B0toDPiLabelMatrix, kTRUE);

      if (mcLabelDPlus >= 0)
      {
        isDesiredCandidate = kTRUE;

        Int_t mcLabelFirstTrack = -1;
        mcLabelFirstTrack = trackDPlusKaon->GetLabel();

        if (mcLabelFirstTrack >= 0)
        {
          AliAODMCParticle *mcParticleFirstTrack = (AliAODMCParticle*)mcTrackArray->At(mcLabelFirstTrack);
          AliAODMCParticle *mcMotherParticle = (AliAODMCParticle*)mcTrackArray->At(mcLabelDPlus);

          if (mcParticleFirstTrack && mcMotherParticle)
          {
            pdgCodeMother = mcMotherParticle->GetPdgCode();

            Double_t vertex_distance = TMath::Sqrt((vertexMother->GetX() - mcParticleFirstTrack->Xv()) * (vertexMother->GetX() - mcParticleFirstTrack->Xv()) + (vertexMother->GetY() - mcParticleFirstTrack->Yv()) * (vertexMother->GetY() - mcParticleFirstTrack->Yv()) + (vertexMother->GetZ() - mcParticleFirstTrack->Zv()) * (vertexMother->GetZ() - mcParticleFirstTrack->Zv()));
            ((TH1F*)fMotherHistogramArrayExtra[motherType][4])->Fill(vertex_distance);

            // Double_t vertex_distance_new = TMath::Sqrt((vertexMotherNew->GetX() - mcParticleFirstTrack->Xv()) * (vertexMotherNew->GetX() - mcParticleFirstTrack->Xv()) + (vertexMotherNew->GetY() - mcParticleFirstTrack->Yv()) * (vertexMotherNew->GetY() - mcParticleFirstTrack->Yv()) + (vertexMotherNew->GetZ() - mcParticleFirstTrack->Zv()) * (vertexMotherNew->GetZ() - mcParticleFirstTrack->Zv()));
            // ((TH1F*)fMotherHistogramArrayExtra[motherType][5])->Fill(vertex_distance_new);

            Double_t momentum_resolution = TMath::Sqrt((trackDPlus->Px() - mcMotherParticle->Px()) * (trackDPlus->Px() - mcMotherParticle->Px()) + (trackDPlus->Py() - mcMotherParticle->Py()) * (trackDPlus->Py() - mcMotherParticle->Py()) + (trackDPlus->Pz() - mcMotherParticle->Pz()) * (trackDPlus->Pz() - mcMotherParticle->Pz()));
            ((TH1F*)fMotherHistogramArrayExtra[motherType][6])->Fill(momentum_resolution);
          }
        }
      }
    }

    // We fill the histograms
    histType = 0;
    FillDPlusHistograms(trackDPlus, primaryVertex, bz, motherType, histType);
    if (isDesiredCandidate && fUseMCInfo) {
      histType = 1;
      FillDPlusHistograms(trackDPlus, primaryVertex, bz, motherType, histType, pdgCodeMother);
    }

    // Here we apply cuts on the particle
    Bool_t cutMother = kFALSE;

    Bool_t bCutArray[39] = {0};
    Int_t cutReturnValue = fCuts->IsDPlusforDPlusptbinSelected(trackDPlus, 0, aodEvent, bCutArray);
    if (cutReturnValue == -1) cutMother = kTRUE;
    if (cutReturnValue == 0) cutMother = kTRUE;
    for (Int_t n = 0; n < 39; ++n)
    {
      if (bCutArray[n] == kTRUE) {
        if (isDesiredCandidate) {
          ((TH1F*)fMotherHistogramArrayExtra[motherType][1])->Fill(n + 1);
        } else ((TH1F*)fMotherHistogramArrayExtra[motherType][0])->Fill(n + 1);
        cutMother = kTRUE;
      }
    }

    if (!fCuts->AreDaughtersSelected(trackDPlus, aodEvent)) cutMother = kTRUE;

    if (!isDesiredCandidate && fQuickSignalAnalysis == 1) cutMother = kTRUE;

    if (cutMother) {
      if (isDesiredCandidate) {
        ((TH1F*)fMotherHistogramArrayExtra[motherType][1])->Fill(0);
      } else ((TH1F*)fMotherHistogramArrayExtra[motherType][0])->Fill(0);
      continue;
    }

    // We fill the cut histograms
    histType = 2;
    FillDPlusHistograms(trackDPlus, primaryVertex, bz, motherType, histType);
    if (isDesiredCandidate && fUseMCInfo) {
      histType = 3;
      FillDPlusHistograms(trackDPlus, primaryVertex, bz, motherType, histType, pdgCodeMother);
    }

    //we save the location of the DPlus candidate
    fDPlusTracks->push_back(j);
  }

  delete vHF; vHF = nullptr;
  return;
}
//-------------------------------------------------------------------------------------
void AliAnalysisTaskSEB0toDPi::B0Selection(AliAODEvent* aodEvent, AliAODVertex *primaryVertex, Double_t bz, TClonesArray * mcTrackArray, TMatrix * B0toDPiLabelMatrix, TClonesArray * DPlusTrackArray, AliAODMCHeader * header) {

  //we loop over all the DPlus candidates
  for (Int_t j = 0; j < (Int_t)fDPlusTracks->size(); j++)
  {

    //Save current Object count
    Int_t ObjectNumber = TProcessID::GetObjectCount();

    //we get the track of the DPlus
    AliAODRecoDecayHF3Prong * trackDPlus = (AliAODRecoDecayHF3Prong*)(DPlusTrackArray->At(fDPlusTracks->at(j)));
    if (!trackDPlus) {std::cout << "found none" << std::endl; continue;}
    if (trackDPlus == nullptr) {std::cout << "found nullptr" << std::endl; continue;}

    //we loop over all the B0 pion candidates
    for (Int_t i = 0; i < (Int_t)fB0PionTracks->size(); i++)
    {

      //we get the track of the  B0 daughter
      AliAODTrack * trackB0Pion = dynamic_cast<AliAODTrack*>(aodEvent->GetTrack(fB0PionTracks->at(i)));
      if (!trackB0Pion) continue;

      //we check if the IDs of the tracks are different
      AliAODTrack* threeProngdaughter0 = (AliAODTrack*)trackDPlus->GetDaughter(0);
      AliAODTrack* threeProngdaughter1 = (AliAODTrack*)trackDPlus->GetDaughter(1);
      AliAODTrack* threeProngdaughter2 = (AliAODTrack*)trackDPlus->GetDaughter(2);
      UShort_t idProng0 = threeProngdaughter0->GetID();
      UShort_t idProng1 = threeProngdaughter1->GetID();
      UShort_t idProng2 = threeProngdaughter2->GetID();

      UShort_t idPion1 = trackB0Pion->GetID();

      if (idPion1 == idProng0 || idPion1 == idProng1 || idPion1 == idProng2) continue;

      // UInt_t prongsDPlus[3];
      // prongsDPlus[0] = 321;
      // prongsDPlus[1] = 211;
      // prongsDPlus[2] = 211;

      // //we check if the particles have the correct charge
      Bool_t bWrongSign = kFALSE;
      if (trackDPlus->Charge() == -1 && trackB0Pion->Charge() == -1) bWrongSign = kTRUE;
      if (trackDPlus->Charge() ==  1 && trackB0Pion->Charge() ==  1) bWrongSign = kTRUE;

      if (trackDPlus->Charge() == 0 || trackB0Pion->Charge() == 0) continue;

      Int_t pdgDPlus = -1; // not relevant for DPlus

      //location B0 pion rotation around PV
      for (Int_t iRot = 0; iRot < fNumberOfRotations + 1; ++iRot)
      {
        //we create a copy of the track that we will rotate
        AliAODTrack * trackB0PionRotated = new AliAODTrack(*trackB0Pion);

        //for iRot == 0, we use the original unrotated track. For iRot > 0 we rotate the track and set the label to -1
        if (iRot != 0)
        {
          //should still check if track is already at PV
          Double_t dPhiRotated = trackB0PionRotated->Phi() + TMath::Pi() - (TMath::Pi() * fDegreePerRotation * fNumberOfRotations / (180.0 * 2.0)) + (TMath::Pi() * fDegreePerRotation * iRot / 180.0);
          trackB0PionRotated->SetPhi(dPhiRotated);
        }

        //we use the B0 pion and DPlus track to reconstruct the vertex for the B0
        AliExternalTrackParam firstTrack;
        firstTrack.CopyFromVTrack(trackDPlus);
        AliExternalTrackParam secondTrack;
        secondTrack.CopyFromVTrack(trackB0PionRotated);

        UInt_t prongs[2];
        prongs[0] = 411;
        prongs[1] = 211;

        UShort_t id[2];
        id[0] = firstTrack.GetID();
        id[1] = secondTrack.GetID();


        // we calculate the vertex of the mother candidate
        TObjArray daughterTracks;
        daughterTracks.Add(&firstTrack);
        daughterTracks.Add(&secondTrack);

        Double_t dispersion = 0;
        AliAODVertex *vertexMother = RecalculateVertex(primaryVertex, &daughterTracks, bz, dispersion, kTRUE, kTRUE, kFALSE, 2);
        if (!vertexMother)
        {
          delete vertexMother; vertexMother = nullptr;
          delete trackB0PionRotated; trackB0PionRotated = nullptr;
          continue;
        }

        //use the new vertex to create the B0 candidate
        Double_t xdummy = 0., ydummy = 0., dca;
        Double_t d0z0[2], covd0z0[3], d0[3], d0err[3];

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

        dca = firstTrack.GetDCA(&secondTrack, bz, xdummy, ydummy);

        Short_t chargeMother = trackDPlus->Charge() + trackB0PionRotated->Charge();

        AliAODRecoDecayHF2Prong trackB0(vertexMother, px, py, pz, d0, d0err, dca);

        trackB0.SetCharge(chargeMother);

        trackB0.GetSecondaryVtx()->AddDaughter(trackDPlus);
        trackB0.GetSecondaryVtx()->AddDaughter(trackB0PionRotated);
        trackB0.SetPrimaryVtxRef((AliAODVertex*)aodEvent->GetPrimaryVertex());
        trackB0.SetProngIDs(2, id);

        // Fiducial cut
        if (TMath::Abs(trackB0.Y(511)) > 0.8) {
          delete vertexMother; vertexMother = nullptr;
          delete trackB0PionRotated; trackB0PionRotated = nullptr;
          continue;
        }

        // We check if the signal is injected, optionally we can reject injected signals
        Bool_t bIsInjected = kFALSE;
        if (fUseMCInfo) {
          bIsInjected = IsCandidateInjected(&trackB0, header, mcTrackArray);
          if (fRemoveInjected && bIsInjected) {
            delete vertexMother; vertexMother = nullptr;
            delete trackB0PionRotated; trackB0PionRotated = nullptr;
            continue;
          }
        }

        // We check if the B0 candidate is a true signal in Monte Carlo
        Bool_t isDesiredCandidate = kFALSE;
        Int_t mcLabelB0 = -1;

        Int_t motherType, histType;
        motherType = 1;

        if (fUseMCInfo)
        {
          mcLabelB0 = MatchCandidateToMonteCarlo(511, &trackB0, mcTrackArray, B0toDPiLabelMatrix);
          if (mcLabelB0 >= 0 && trackB0PionRotated->GetLabel() >= 0 && iRot == 0)
          {
            isDesiredCandidate = kTRUE;
          }
        }

        histType = 0;

        if (isDesiredCandidate)
        {
          AliAODMCParticle *mcTrackB0Pion = (AliAODMCParticle*)mcTrackArray->At(trackB0Pion->GetLabel());
          AliAODMCParticle *mcTrackB0 = (AliAODMCParticle*)mcTrackArray->At(mcLabelB0);

          // Double_t vertex_distance = TMath::Sqrt((vertexMotherTemp->GetX() - mcTrackB0Pion->Xv()) * (vertexMotherTemp->GetX() - mcTrackB0Pion->Xv()) + (vertexMotherTemp->GetY() - mcTrackB0Pion->Yv()) * (vertexMotherTemp->GetY() - mcTrackB0Pion->Yv()) + (vertexMotherTemp->GetZ() - mcTrackB0Pion->Zv()) * (vertexMotherTemp->GetZ() - mcTrackB0Pion->Zv()));
          // ((TH1F*)fMotherHistogramArrayExtra[motherType][4])->Fill(vertex_distance);

          Double_t vertex_distance_new = TMath::Sqrt((vertexMother->GetX() - mcTrackB0Pion->Xv()) * (vertexMother->GetX() - mcTrackB0Pion->Xv()) + (vertexMother->GetY() - mcTrackB0Pion->Yv()) * (vertexMother->GetY() - mcTrackB0Pion->Yv()) + (vertexMother->GetZ() - mcTrackB0Pion->Zv()) * (vertexMother->GetZ() - mcTrackB0Pion->Zv()));
          ((TH1F*)fMotherHistogramArrayExtra[motherType][5])->Fill(vertex_distance_new);

          Double_t momentum_resolution = TMath::Sqrt((trackB0.Px() - mcTrackB0->Px()) * (trackB0.Px() - mcTrackB0->Px()) + (trackB0.Py() - mcTrackB0->Py()) * (trackB0.Py() - mcTrackB0->Py()) + (trackB0.Pz() - mcTrackB0->Pz()) * (trackB0.Pz() - mcTrackB0->Pz()));
          ((TH1F*)fMotherHistogramArrayExtra[motherType][6])->Fill(momentum_resolution);
        }

        Bool_t bFillUncutHistograms = kTRUE;
        if (iRot != 0 && !fSaveTRHists) bFillUncutHistograms = kFALSE;
        if (iRot == 0 && fSaveTRHists) bFillUncutHistograms = kFALSE;

        if (!bWrongSign && bFillUncutHistograms)
        {
          // We fill the histograms
          histType = 0;
          FillB0Histograms(&trackB0, primaryVertex, bz, motherType, histType);

          if (isDesiredCandidate)
          {
            histType = 1;
            FillB0Histograms(&trackB0, primaryVertex, bz, motherType, histType);
          }
        }


        // We apply cuts
        Bool_t cutMother = kFALSE;

        Bool_t bCutArray[78] = {0};
        Int_t numberOfCuts = 78;
        Int_t cutReturnValue = fCuts->IsSelected(&trackB0, 0, aodEvent, bCutArray);
        if (cutReturnValue == -1) cutMother = kTRUE;
        if (cutReturnValue == 0) cutMother = kTRUE;


        // We save information about the cuts
        TString histName = "";
        Double_t invariantMassMother = trackB0.InvMass(2, prongs);
        Double_t pdgMassMother = TDatabasePDG::Instance()->GetParticle(511)->Mass();
        Double_t massWindow = fHistMassWindow; //GeV/c^2

        for (Int_t n = 0; n < 78; ++n)
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
            for (Int_t m = 0; m < numberOfCuts; ++m)
            {
              if (bCutArray[m] == kFALSE) continue;
              if (isDesiredCandidate == kFALSE) ((TH2F*)(fResultsHistogramArray2D[2][0]))->Fill(l, m);
              if (isDesiredCandidate == kTRUE) ((TH2F*)(fResultsHistogramArray2D[2][1]))->Fill(l, m);
            }
          }

          for (Int_t l = 0; l < numberOfCuts; ++l) //unique
          {
            if (bCutArray[l] == kFALSE) continue;
            Bool_t bFill = kTRUE;
            for (Int_t m = 0; m < numberOfCuts; ++m)
            {
              if (l == m) continue;
              if (bCutArray[m] == kTRUE)
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
          delete trackB0PionRotated; trackB0PionRotated = nullptr;
          continue;
        }

        if (!bWrongSign)
        {
          histType = 2;
          FillB0Histograms(&trackB0, primaryVertex, bz, motherType, histType);
          if (fUseMCInfo && isDesiredCandidate)
          {
            //fill mc histograms
            histType = 3;
            FillB0Histograms(&trackB0, primaryVertex, bz, motherType, histType);
          }
        }

        Bool_t bFillhistograms = kFALSE;
        Bool_t bPassMassWindow = (TMath::Abs(invariantMassMother - pdgMassMother) < massWindow);
        Bool_t bPassSideBand = (TMath::Abs(invariantMassMother - pdgMassMother) > fSideBandLow) && (TMath::Abs(invariantMassMother - pdgMassMother) < fSideBandHigh);
        if (isDesiredCandidate && bPassMassWindow) bFillhistograms = kTRUE;
        if (!isDesiredCandidate && bPassMassWindow && !fUseSideBands) bFillhistograms = kTRUE;
        if (!isDesiredCandidate && bPassSideBand && fUseSideBands) bFillhistograms = kTRUE;
        if (iRot != 0 && !fSaveTRHists) bFillhistograms = kFALSE;
        if (iRot == 0 && fSaveTRHists) bFillhistograms = kFALSE;

        if (bFillhistograms)
        {
          if (!bWrongSign)
          {
            FillFinalTrackHistograms(&trackB0, primaryVertex, bz, isDesiredCandidate, mcTrackArray);
            if (!isDesiredCandidate)
            {
              motherType = 0; histType = 4; FillDPlusHistograms(trackDPlus, primaryVertex, bz, motherType, histType, pdgDPlus);
              motherType = 1; histType = 4; FillB0Histograms(&trackB0, primaryVertex, bz, motherType, histType);
            }
            if (isDesiredCandidate)
            {
              motherType = 0; histType = 5; FillDPlusHistograms(trackDPlus, primaryVertex, bz, motherType, histType, pdgDPlus);
              motherType = 1; histType = 5; FillB0Histograms(&trackB0, primaryVertex, bz, motherType, histType);
            }
          }
        }


        // Here we fill the histograms per pt bin and apply the same sign method
        TString ptBinMother = "";
        Int_t ptBin = fCuts->PtBin(trackB0.Pt());
        ptBinMother += "_ptbin_"; ptBinMother += fPtBinLimits[ptBin]; ptBinMother += "_to_"; ptBinMother += fPtBinLimits[ptBin + 1];
        histType = 6 + 2 * ptBin;

        Int_t dPlusPtBin = fCuts->PtBinDPlusforDPlusptbin(trackDPlus->Pt());
        Int_t histTypeDPlus = 2 * dPlusPtBin;


        if (bFillhistograms)
        {
          if (!bWrongSign && histType > 5)
          {
            if (!isDesiredCandidate)
            {
              motherType = 0; FillDPlusHistograms(trackDPlus, primaryVertex, bz, motherType, histType, pdgDPlus);
              motherType = 1; FillB0Histograms(&trackB0, primaryVertex, bz, motherType, histType);
              motherType = 2; FillDPlusHistograms(trackDPlus, primaryVertex, bz, motherType, histTypeDPlus, pdgDPlus);
            }

            if (isDesiredCandidate)
            {
              histType += 1;
              motherType = 0; FillDPlusHistograms(trackDPlus, primaryVertex, bz, motherType, histType, pdgDPlus);
              motherType = 1; FillB0Histograms(&trackB0, primaryVertex, bz, motherType, histType);
              motherType = 2; FillDPlusHistograms(trackDPlus, primaryVertex, bz, motherType, histTypeDPlus + 1, pdgDPlus);
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
            CutOptimizationVariableValues(&trackB0, aodEvent);
            CutOptimizationLoop(nStartVariable, nVariables, nCuts, ptBin, nStartFillNumber, isDesiredCandidate, nSigmaBin);
          }
        }

        if (isDesiredCandidate) ((TH1F*)fResultsHistogramArray[3][3])->Fill(ptBin);

        Double_t invmassDelta = DeltaInvMassB0Kpipipi(&trackB0);

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

          ((TH1F*)fResultsHistogramArray[4 + m][0])->Fill(invariantMassMother, fillfactor);
          ((TH1F*)fResultsHistogramArray[4 + m][0 + fnPtBins + 3])->Fill(invmassDelta, fillfactor);

          if (ptBin > 0)
          {
            ((TH1F*)fResultsHistogramArray[4 + m][1])->Fill(invariantMassMother, fillfactor);
            ((TH1F*)fResultsHistogramArray[4 + m][1 + fnPtBins + 3])->Fill(invmassDelta, fillfactor);
          }
          if (ptBin > 1)
          {
            ((TH1F*)fResultsHistogramArray[4 + m][2])->Fill(invariantMassMother, fillfactor);
            ((TH1F*)fResultsHistogramArray[4 + m][2 + fnPtBins + 3])->Fill(invmassDelta, fillfactor);
          }

          ((TH1F*)fResultsHistogramArray[4 + m][ptBin + 3])->Fill(invariantMassMother, fillfactor);
          ((TH1F*)fResultsHistogramArray[4 + m][ptBin + 3 + fnPtBins + 3])->Fill(invmassDelta, fillfactor);
        }
        delete vertexMother; vertexMother = nullptr;
        delete trackB0PionRotated; trackB0PionRotated = nullptr;
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
void AliAnalysisTaskSEB0toDPi::FillDaughterHistograms(AliAODTrack* daughterTrack, Int_t daughterType, Int_t histType, Double_t ptB0, AliAODVertex *primaryVertex, Double_t bz, Bool_t isDesiredCandidate, TClonesArray * mcTrackArray) {

  //In this function we fill histograms with the properties of all the daughters of our selected signal candidate

  Int_t particleType = 2;
  if (daughterType == 0) particleType = 3;

  Double_t pt_track = 0;
  Double_t momentum_track = 0;
  Int_t numberOfITS = 0;
  Int_t numberOfTPC = 0;
  Int_t totalNumberOfITS = 0;
  Int_t totalNumberOfTPC = 0;
  Double_t nSigmaTPC = 0;
  Double_t nSigmaTOF = 0;
  Double_t nSigmaTPCtotal = 0;
  Double_t nSigmaTOFtotal = 0;
  Int_t TPCok = 0;
  Int_t TOFok = 0;

  if (isDesiredCandidate) histType++;

  AliExternalTrackParam particleTrack;
  particleTrack.CopyFromVTrack(daughterTrack);
  Double_t d0[2], covd0[3];
  particleTrack.PropagateToDCA(primaryVertex, bz, 100., d0, covd0);

  AliAODPidHF* trackPIDHF = (AliAODPidHF*)fCuts->GetPidHF();

  //fill the DPlus pion info
  pt_track = daughterTrack->Pt();
  momentum_track = daughterTrack->P();
  numberOfITS = daughterTrack->GetITSNcls();
  numberOfTPC = daughterTrack->GetTPCNcls();
  totalNumberOfITS += numberOfITS;
  totalNumberOfTPC += numberOfTPC;
  if (trackPIDHF) TPCok = trackPIDHF->GetnSigmaTPC(daughterTrack, particleType, nSigmaTPC);
  if (trackPIDHF) TOFok = trackPIDHF->GetnSigmaTOF(daughterTrack, particleType, nSigmaTOF);
  if (TPCok != -1) nSigmaTPCtotal += nSigmaTPC * nSigmaTPC;
  if (TOFok != -1) nSigmaTOFtotal += nSigmaTOF * nSigmaTOF;

  ((TH1F*)fDaughterHistogramArray[daughterType][histType][0])->Fill(pt_track);
  ((TH1F*)fDaughterHistogramArray[daughterType][histType][1])->Fill(momentum_track);
  ((TH1F*)fDaughterHistogramArray[daughterType][histType][2])->Fill(numberOfITS);
  ((TH1F*)fDaughterHistogramArray[daughterType][histType][3])->Fill(numberOfTPC);

  for (Int_t j = 0; j < 10; ++j)
  {
    if (daughterTrack->HasPointOnITSLayer(j)) ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(j);

  }

  if (TPCok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][5])->Fill(nSigmaTPC);
  if (TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][6])->Fill(nSigmaTOF);
  if (TPCok != -1 && TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][7])->Fill(TMath::Sqrt(nSigmaTPC * nSigmaTPC + nSigmaTOF * nSigmaTOF));
  ((TH1F*)fDaughterHistogramArray[daughterType][histType][8])->Fill(TMath::Abs(d0[0]));
  ((TH1F*)fDaughterHistogramArray[daughterType][histType][9])->Fill(TMath::Abs(d0[0] / TMath::Sqrt(covd0[0])));
  ((TH1F*)fDaughterHistogramArray[daughterType][histType][10])->Fill(daughterTrack->Eta());

  if (ptB0 > 0.0) ((TH2F*)fDaughterHistogramArray2D[daughterType][4])->Fill(ptB0, pt_track);


  //we save the pdgcode of the used particle and its mother to check PID efficiency
  if (fUseMCInfo)
  {
    Float_t pdgCodeParticle = -1;
    Float_t pdgCodeParticleMother = -1;
    Int_t mcLabelParticle = -1;
    Int_t mcLabelParticleMother = -1;
    mcLabelParticle = daughterTrack->GetLabel();

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
  return;
}
//-------------------------------------------------------------------------------------
void AliAnalysisTaskSEB0toDPi::FillFinalTrackHistograms(AliAODRecoDecayHF2Prong * selectedB0, AliAODVertex *primaryVertex, Double_t bz, Bool_t isDesiredCandidate, TClonesArray * mcTrackArray) {

  //In this function we fill histograms with the properties of all the daughters of our selected signal candidate

  AliAODRecoDecayHF3Prong* selectedDPlus = (AliAODRecoDecayHF3Prong*)selectedB0->GetDaughter(0);
  AliAODTrack* selectedB0Pion = (AliAODTrack*)selectedB0->GetDaughter(1);

  AliAODTrack* selectedDPlusKaon = (AliAODTrack*)selectedDPlus->GetDaughter(1);
  AliAODTrack* selectedDPlusFirstPion = (AliAODTrack*)selectedDPlus->GetDaughter(0);
  AliAODTrack* selectedDPlusSecondPion = (AliAODTrack*)selectedDPlus->GetDaughter(2);


  FillDaughterHistograms(selectedB0Pion, 2, 4, selectedB0->Pt(), primaryVertex, bz, isDesiredCandidate, mcTrackArray);
  FillDaughterHistograms(selectedDPlusKaon, 0, 4, selectedB0->Pt(), primaryVertex, bz, isDesiredCandidate, mcTrackArray);
  FillDaughterHistograms(selectedDPlusFirstPion, 1, 4, selectedB0->Pt(), primaryVertex, bz, isDesiredCandidate, mcTrackArray);
  FillDaughterHistograms(selectedDPlusSecondPion, 1, 4, selectedB0->Pt(), primaryVertex, bz, isDesiredCandidate, mcTrackArray);

  // not yet implemented for this decay
  // if (!isDesiredCandidate)
  // {
  //   ((TH1F*)(fResultsHistogramArray[11][0]))->Fill(totalNumberOfITS);
  //   ((TH1F*)(fResultsHistogramArray[11][2]))->Fill(totalNumberOfTPC);
  //   ((TH1F*)(fResultsHistogramArray[11][4]))->Fill(TMath::Sqrt(nSigmaTPCtotal + nSigmaTOFtotal));
  // }
  // if (isDesiredCandidate)
  // {
  //   ((TH1F*)(fResultsHistogramArray[11][1]))->Fill(totalNumberOfITS);
  //   ((TH1F*)(fResultsHistogramArray[11][3]))->Fill(totalNumberOfTPC);
  //   ((TH1F*)(fResultsHistogramArray[11][5]))->Fill(TMath::Sqrt(nSigmaTPCtotal + nSigmaTOFtotal));
  // }

  // we look at the invariant mass combinations of all daughter tracks
  // AliExternalTrackParam trackB0Pion;
  // trackB0Pion.CopyFromVTrack(selectedB0Pion);
  // AliExternalTrackParam trackDPlusPion;
  // trackDPlusPion.CopyFromVTrack(selectedDPlusPion);
  // AliExternalTrackParam trackDPlusKaon;
  // trackDPlusKaon.CopyFromVTrack(selectedDPlusKaon);

  // UInt_t prongs2[2] = {0};
  // UInt_t prongs3[3] = {0};

  // prongs2[0] = 211; prongs2[1] = 211;
  // TwoTrackCombinationInfo(&trackDPlusPion, &trackB0Pion, primaryVertex, bz, isDesiredCandidate, 0, prongs2);
  // prongs2[0] = 321; prongs2[1] = 211;
  // TwoTrackCombinationInfo(&trackDPlusKaon, &trackB0Pion, primaryVertex, bz, isDesiredCandidate, 2, prongs2);
  // prongs3[0] = 211; prongs3[1] = 321; prongs3[2] = 211;
  // ThreeTrackCombinationInfo(&trackDPlusPion, &trackDPlusKaon, &trackB0Pion, primaryVertex, bz, isDesiredCandidate, 4, prongs3);

  return;
}
//-------------------------------------------------------------------------------------
Double_t AliAnalysisTaskSEB0toDPi::DeltaInvMassB0Kpipipi(AliAODRecoDecayHF2Prong *Bzero) const
{
  ///
  /// delta invariant mass
  ///

  AliAODRecoDecayHF3Prong * DPlus = (AliAODRecoDecayHF3Prong*)Bzero->GetDaughter(0);

  Double_t e[4] = {0};
  e[0] = DPlus->EProng(0, 211); // to be done: Check effect of taking energies of DPlus daughters at B vertex by propagating D to B first ////////////////////////////////////////////////////////////////////////////////////////////////////////
  e[1] = DPlus->EProng(1, 321);
  e[2] = DPlus->EProng(2, 211);
  e[3] = Bzero->EProng(1, 211);

  Double_t esum = e[0] + e[1] + e[2] + e[3];
  Double_t invMassB0 = TMath::Sqrt(esum * esum - Bzero->P2());

  Double_t eD[3] = {0};
  eD[0] = DPlus->EProng(0, 211);
  eD[1] = DPlus->EProng(1, 321);
  eD[2] = DPlus->EProng(2, 211);

  Double_t esumD = eD[0] + eD[1] + eD[2];
  Double_t invMassDPlus = TMath::Sqrt(esumD * esumD - DPlus->P2());

  return invMassB0 - invMassDPlus;
}
//-------------------------------------------------------------------------------------
void AliAnalysisTaskSEB0toDPi::FillDPlusHistograms(AliAODRecoDecayHF3Prong * selectedMother, AliAODVertex *primaryVertex, Double_t bz, Int_t motherType, Int_t histType, Int_t pdgCodeMother) {

  if (histType < 0) return;

  //In this function we fill the histograms of the reconstructed mothers
  Double_t ptMother = 0.0;
  Double_t momentumMother = 0.0;
  Double_t etaMother = 0.0;
  Double_t phiMother = 0.0;
  Double_t d0Mother = 0.0;
  Double_t d0firstTrack = 0.0;
  Double_t d0secondTrack = 0.0;
  Double_t d0thirdTrack = 0.0;
  Double_t pointingAngle = 0.0;
  Double_t impactProduct123 = 0.0;
  Double_t impactProduct12 = 0.0;
  Double_t impactProduct13 = 0.0;
  Double_t impactProduct23 = 0.0;
  Double_t impactProductXY = 0.0;
  Double_t invariantMassMother = 0.0;
  Double_t invmassDelta = 0.0;
  Double_t dcaMother12 = 0.0;
  Double_t dcaMother13 = 0.0;
  Double_t dcaMother23 = 0.0;
  AliAODVertex * vertexMother = 0x0;
  Double_t vertexDistance = 0.0;
  Double_t decayTime = 0.0;
  Double_t ptFirstDaughter = 0.0;
  Double_t ptSecondDaughter = 0.0;
  Double_t ptThirdDaughter = 0.0;
  UInt_t prongs[3];

  Double_t angleMotherFirstDaughter = 0.0;
  Double_t angleMotherSecondDaughter = 0.0;
  Double_t angleMotherThirdDaughter = 0.0;

  Double_t smallestAngleMotherDaughter = 0.0;
  Double_t largestAngleMotherDaughter = 0.0;

  Double_t normDecayLength = 0;
  Double_t pdgMassMother = TDatabasePDG::Instance()->GetParticle(411)->Mass();

  prongs[0] = 211;
  prongs[1] = 321;
  prongs[2] = 211;

  AliAODTrack * firstDaughter = (AliAODTrack*)selectedMother->GetDaughter(0);
  AliAODTrack * secondDaughter = (AliAODTrack*)selectedMother->GetDaughter(1);
  AliAODTrack * thirdDaughter = (AliAODTrack*)selectedMother->GetDaughter(2);

  vertexMother = selectedMother->GetSecondaryVtx();
  ptFirstDaughter = firstDaughter->Pt();
  ptSecondDaughter = secondDaughter->Pt();
  ptThirdDaughter = thirdDaughter->Pt();

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
  d0thirdTrack = TMath::Abs(selectedMother->Getd0Prong(2));
  pointingAngle = selectedMother->CosPointingAngle();
  impactProduct123 = selectedMother->Getd0Prong(0) * selectedMother->Getd0Prong(1) * selectedMother->Getd0Prong(2);
  impactProduct12 = selectedMother->Getd0Prong(0) * selectedMother->Getd0Prong(1);
  impactProduct13 = selectedMother->Getd0Prong(0) * selectedMother->Getd0Prong(2);
  impactProduct23 = selectedMother->Getd0Prong(1) * selectedMother->Getd0Prong(2);

  impactProductXY = TMath::Abs(selectedMother->ImpParXY());
  invariantMassMother = selectedMother->InvMass(3, prongs);

  dcaMother12 = selectedMother->GetDCA(0);
  dcaMother13 = selectedMother->GetDCA(1);
  dcaMother23 = selectedMother->GetDCA(2);

  vertexDistance = vertexMother->DistanceToVertex(primaryVertex);

  angleMotherFirstDaughter = (selectedMother->Px() * firstDaughter->Px() + selectedMother->Py() * firstDaughter->Py() + selectedMother->Pz() * firstDaughter->Pz()) / (selectedMother->P() * firstDaughter->P());
  angleMotherSecondDaughter = (selectedMother->Px() * secondDaughter->Px() + selectedMother->Py() * secondDaughter->Py() + selectedMother->Pz() * secondDaughter->Pz()) / (selectedMother->P() * secondDaughter->P());
  angleMotherThirdDaughter = (selectedMother->Px() * secondDaughter->Px() + selectedMother->Py() * secondDaughter->Py() + selectedMother->Pz() * secondDaughter->Pz()) / (selectedMother->P() * secondDaughter->P());

  if (angleMotherFirstDaughter > angleMotherSecondDaughter) {smallestAngleMotherDaughter = angleMotherSecondDaughter;} else {smallestAngleMotherDaughter = angleMotherFirstDaughter;}
  if (angleMotherThirdDaughter < smallestAngleMotherDaughter) smallestAngleMotherDaughter = angleMotherThirdDaughter;
  if (angleMotherFirstDaughter > angleMotherSecondDaughter) {largestAngleMotherDaughter = angleMotherFirstDaughter;} else {largestAngleMotherDaughter = angleMotherSecondDaughter;}
  if (angleMotherThirdDaughter > largestAngleMotherDaughter) largestAngleMotherDaughter = angleMotherThirdDaughter;

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

  Double_t eKaon1 = selectedMother->EProng(1, 321);
  Double_t ePion2 = selectedMother->EProng(2, 211);
  Double_t eSum = eKaon1 + ePion2;
  Double_t pxPions = (selectedMother->PxProng(1) + selectedMother->PxProng(2)) * (selectedMother->PxProng(1) + selectedMother->PxProng(2));
  Double_t pyPions = (selectedMother->PyProng(1) + selectedMother->PyProng(2)) * (selectedMother->PyProng(1) + selectedMother->PyProng(2));
  Double_t pzPions = (selectedMother->PzProng(1) + selectedMother->PzProng(2)) * (selectedMother->PzProng(1) + selectedMother->PzProng(2));
  Double_t pSum = pxPions + pyPions + pzPions;
  Double_t invMassPions = TMath::Sqrt(eSum * eSum - pSum);
  Double_t invMassDPlus = selectedMother->InvMass(3, prongs);;
  invmassDelta = invMassDPlus - invMassPions;

  Double_t vertexMotherX = vertexMother->GetX();
  Double_t vertexMotherY = vertexMother->GetY();
  Double_t vertexMotherZ = vertexMother->GetZ();

  Double_t cosPointingAngleXY = selectedMother->CosPointingAngleXY();
  Double_t distanceXYToVertex = vertexMother->DistanceXYToVertex(primaryVertex);
  Double_t normalizedDecayLengthXY = selectedMother->NormalizedDecayLengthXY();

  ((TH1F*)fMotherHistogramArray[motherType][histType][0])->Fill(invariantMassMother);
  ((TH1F*)fMotherHistogramArray[motherType][histType][1])->Fill(invmassDelta);
  ((TH1F*)fMotherHistogramArray[motherType][histType][2])->Fill(ptMother);
  ((TH1F*)fMotherHistogramArray[motherType][histType][3])->Fill(ptFirstDaughter);
  ((TH1F*)fMotherHistogramArray[motherType][histType][4])->Fill(ptSecondDaughter);
  ((TH1F*)fMotherHistogramArray[motherType][histType][5])->Fill(ptThirdDaughter);
  ((TH1F*)fMotherHistogramArray[motherType][histType][6])->Fill(d0Mother);
  ((TH1F*)fMotherHistogramArray[motherType][histType][7])->Fill(d0firstTrack);
  ((TH1F*)fMotherHistogramArray[motherType][histType][8])->Fill(d0secondTrack);
  ((TH1F*)fMotherHistogramArray[motherType][histType][9])->Fill(d0thirdTrack);
  ((TH1F*)fMotherHistogramArray[motherType][histType][10])->Fill(pointingAngle);
  ((TH1F*)fMotherHistogramArray[motherType][histType][11])->Fill(cosPointingAngleXY);
  ((TH1F*)fMotherHistogramArray[motherType][histType][12])->Fill(impactProduct123);
  ((TH1F*)fMotherHistogramArray[motherType][histType][13])->Fill(impactProduct12);
  ((TH1F*)fMotherHistogramArray[motherType][histType][14])->Fill(impactProduct13);
  ((TH1F*)fMotherHistogramArray[motherType][histType][15])->Fill(impactProduct23);
  ((TH1F*)fMotherHistogramArray[motherType][histType][16])->Fill(impactProductXY);
  ((TH1F*)fMotherHistogramArray[motherType][histType][17])->Fill(dcaMother12);
  ((TH1F*)fMotherHistogramArray[motherType][histType][18])->Fill(dcaMother13);
  ((TH1F*)fMotherHistogramArray[motherType][histType][19])->Fill(dcaMother23);
  ((TH1F*)fMotherHistogramArray[motherType][histType][20])->Fill(vertexDistance);
  ((TH1F*)fMotherHistogramArray[motherType][histType][21])->Fill(distanceXYToVertex);
  ((TH1F*)fMotherHistogramArray[motherType][histType][22])->Fill(normDecayLength);
  ((TH1F*)fMotherHistogramArray[motherType][histType][23])->Fill(normalizedDecayLengthXY);
  ((TH1F*)fMotherHistogramArray[motherType][histType][24])->Fill(pseudoProperDecayTime);
  ((TH1F*)fMotherHistogramArray[motherType][histType][25])->Fill(decayTime);
  ((TH1F*)fMotherHistogramArray[motherType][histType][26])->Fill(normalizedDecayTime);

  ((TH1F*)fMotherHistogramArray[motherType][histType][27])->Fill(smallestAngleMotherDaughter);
  ((TH1F*)fMotherHistogramArray[motherType][histType][28])->Fill(largestAngleMotherDaughter);
  ((TH1F*)fMotherHistogramArray[motherType][histType][29])->Fill(largestAngleMotherDaughter - smallestAngleMotherDaughter);

  ((TH1F*)fMotherHistogramArray[motherType][histType][30])->Fill(TMath::Abs(selectedMother->Getd0Prong(0) / selectedMother->Getd0errProng(0)));
  ((TH1F*)fMotherHistogramArray[motherType][histType][31])->Fill(TMath::Abs(selectedMother->Getd0Prong(1) / selectedMother->Getd0errProng(1)));
  ((TH1F*)fMotherHistogramArray[motherType][histType][32])->Fill(TMath::Abs(selectedMother->Getd0Prong(2) / selectedMother->Getd0errProng(2)));
  ((TH1F*)fMotherHistogramArray[motherType][histType][33])->Fill(TMath::Abs(d0[0] / TMath::Sqrt(covd0z0[0])));
  ((TH1F*)fMotherHistogramArray[motherType][histType][34])->Fill((selectedMother->Getd0Prong(0) / selectedMother->Getd0errProng(0)) * (selectedMother->Getd0Prong(1) / selectedMother->Getd0errProng(1)) * (selectedMother->Getd0Prong(2) / selectedMother->Getd0errProng(2)));

  ((TH1F*)fMotherHistogramArray[motherType][histType][35])->Fill(selectedMother->GetDist12toPrim());
  ((TH1F*)fMotherHistogramArray[motherType][histType][36])->Fill(selectedMother->GetDist23toPrim());
  ((TH1F*)fMotherHistogramArray[motherType][histType][37])->Fill(selectedMother->GetSigmaVert());

  ((TH1F*)fMotherHistogramArray[motherType][histType][49])->Fill(selectedMother->Y(411));
  ((TH1F*)fMotherHistogramArray[motherType][histType][50])->Fill(etaMother);
  ((TH1F*)fMotherHistogramArray[motherType][histType][51])->Fill(phiMother);
  ((TH1F*)fMotherHistogramArray[motherType][histType][52])->Fill(vertexMotherX);
  ((TH1F*)fMotherHistogramArray[motherType][histType][53])->Fill(vertexMotherY);
  ((TH1F*)fMotherHistogramArray[motherType][histType][54])->Fill(vertexMotherZ);
  ((TH1F*)fMotherHistogramArray[motherType][histType][55])->Fill(vertexMother->GetChi2perNDF());

  //we fill the 2D histograms // to be done
  // Int_t nFirst = 0;
  // Int_t nSecond = 1;
  // Int_t nVariables = 10;
  // Int_t nHistograms = nVariables * (nVariables - 1) / 2;
  // for (Int_t k = 0; k < nHistograms; ++k)
  // {
  //   Double_t firstVariable = 0.0;
  //   Double_t secondVariable = 0.0;

  //   if (nFirst == 0) firstVariable = d0firstTrack;
  //   if (nFirst == 1) firstVariable = d0secondTrack;
  //   if (nFirst == 2) firstVariable = d0Mother;
  //   if (nFirst == 3) firstVariable = pointingAngle;
  //   if (nFirst == 4) firstVariable = impactProduct;
  //   if (nFirst == 5) firstVariable = impactProductXY;
  //   if (nFirst == 6) firstVariable = vertexDistance;
  //   if (nFirst == 7) firstVariable = normDecayLength;
  //   if (nFirst == 8) firstVariable = cosPointingAngleXY;
  //   if (nFirst == 9) firstVariable = distanceXYToVertex;
  //   if (nFirst == 10) firstVariable = normalizedDecayLengthXY;

  //   if (nSecond == 0) secondVariable = d0firstTrack;
  //   if (nSecond == 1) secondVariable = d0secondTrack;
  //   if (nSecond == 2) secondVariable = d0Mother;
  //   if (nSecond == 3) secondVariable = pointingAngle;
  //   if (nSecond == 4) secondVariable = impactProduct;
  //   if (nSecond == 5) secondVariable = impactProductXY;
  //   if (nSecond == 6) secondVariable = vertexDistance;
  //   if (nSecond == 7) secondVariable = normDecayLength;
  //   if (nSecond == 8) secondVariable = cosPointingAngleXY;
  //   if (nSecond == 9) secondVariable = distanceXYToVertex;
  //   if (nSecond == 10) secondVariable = normalizedDecayLengthXY;

  //   ((TH2F*)fMotherHistogramArray2D[motherType][histType][k])->Fill(firstVariable, secondVariable);

  //   nSecond++;
  //   if (nSecond > nVariables)
  //   {
  //     nFirst++;
  //     nSecond = nFirst + 1;
  //   }
  // }

  return;
}
//-------------------------------------------------------------------------------------
void AliAnalysisTaskSEB0toDPi::FillB0Histograms(AliAODRecoDecayHF2Prong * selectedMother, AliAODVertex *primaryVertex, Double_t bz, Int_t motherType, Int_t histType) {

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
  Double_t ptFirstDaughter = 0.0;
  Double_t ptSecondDaughter = 0.0;
  UInt_t prongs[2];

  Double_t cosThetaStar = 0;
  Double_t angleMotherFirstDaughter = 0.0;
  Double_t angleMotherSecondDaughter = 0.0;
  Double_t angleBetweenBothDaughters = 0;

  Double_t normDecayLength = 0;
  Double_t pdgMassMother = TDatabasePDG::Instance()->GetParticle(511)->Mass();

  prongs[0] = 411;
  prongs[1] = 211;

  AliAODRecoDecayHF3Prong * firstDaughter = (AliAODRecoDecayHF3Prong*)selectedMother->GetDaughter(0);
  AliAODTrack * secondDaughter = (AliAODTrack*)selectedMother->GetDaughter(1);

  vertexMother = selectedMother->GetSecondaryVtx();
  ptFirstDaughter = firstDaughter->Pt();
  ptSecondDaughter = secondDaughter->Pt();

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
  impactProduct = selectedMother->Getd0Prong(0) * selectedMother->Getd0Prong(1);

  impactProductXY = TMath::Abs(selectedMother->ImpParXY());
  invariantMassMother = selectedMother->InvMass(2, prongs);

  dcaMother = selectedMother->GetDCA();

  vertexDistance = vertexMother->DistanceToVertex(primaryVertex);

  angleMotherFirstDaughter = (selectedMother->Px() * firstDaughter->Px() + selectedMother->Py() * firstDaughter->Py() + selectedMother->Pz() * firstDaughter->Pz()) / (selectedMother->P() * firstDaughter->P());
  angleMotherSecondDaughter = (selectedMother->Px() * secondDaughter->Px() + selectedMother->Py() * secondDaughter->Py() + selectedMother->Pz() * secondDaughter->Pz()) / (selectedMother->P() * secondDaughter->P());
  angleBetweenBothDaughters  = (firstDaughter->Px() * secondDaughter->Px() + firstDaughter->Py() * secondDaughter->Py() + firstDaughter->Pz() * secondDaughter->Pz()) / (firstDaughter->P() * secondDaughter->P());

  cosThetaStar = selectedMother->CosThetaStar(0, 511, 411, 211);
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

  invmassDelta = DeltaInvMassB0Kpipipi(selectedMother);

  Double_t vertexMotherX = vertexMother->GetX();
  Double_t vertexMotherY = vertexMother->GetY();
  Double_t vertexMotherZ = vertexMother->GetZ();

  Double_t cosPointingAngleXY = selectedMother->CosPointingAngleXY();
  Double_t distanceXYToVertex = vertexMother->DistanceXYToVertex(primaryVertex);
  Double_t normalizedDecayLengthXY = selectedMother->NormalizedDecayLengthXY();

  ((TH1F*)fMotherHistogramArray[motherType][histType][0])->Fill(invariantMassMother);
  ((TH1F*)fMotherHistogramArray[motherType][histType][1])->Fill(invmassDelta);
  ((TH1F*)fMotherHistogramArray[motherType][histType][2])->Fill(ptMother);
  ((TH1F*)fMotherHistogramArray[motherType][histType][3])->Fill(ptFirstDaughter);
  ((TH1F*)fMotherHistogramArray[motherType][histType][4])->Fill(ptSecondDaughter);
  ((TH1F*)fMotherHistogramArray[motherType][histType][5])->Fill(d0Mother);
  ((TH1F*)fMotherHistogramArray[motherType][histType][6])->Fill(d0firstTrack);
  ((TH1F*)fMotherHistogramArray[motherType][histType][7])->Fill(d0secondTrack);
  ((TH1F*)fMotherHistogramArray[motherType][histType][8])->Fill(pointingAngle);
  ((TH1F*)fMotherHistogramArray[motherType][histType][9])->Fill(cosPointingAngleXY);
  ((TH1F*)fMotherHistogramArray[motherType][histType][10])->Fill(impactProduct);
  ((TH1F*)fMotherHistogramArray[motherType][histType][11])->Fill(impactProductXY);
  ((TH1F*)fMotherHistogramArray[motherType][histType][12])->Fill(dcaMother);
  ((TH1F*)fMotherHistogramArray[motherType][histType][13])->Fill(vertexDistance);
  ((TH1F*)fMotherHistogramArray[motherType][histType][14])->Fill(distanceXYToVertex);
  ((TH1F*)fMotherHistogramArray[motherType][histType][15])->Fill(normDecayLength);
  ((TH1F*)fMotherHistogramArray[motherType][histType][16])->Fill(normalizedDecayLengthXY);
  ((TH1F*)fMotherHistogramArray[motherType][histType][17])->Fill(pseudoProperDecayTime);
  ((TH1F*)fMotherHistogramArray[motherType][histType][18])->Fill(decayTime);
  ((TH1F*)fMotherHistogramArray[motherType][histType][19])->Fill(normalizedDecayTime);

  ((TH1F*)fMotherHistogramArray[motherType][histType][20])->Fill(angleMotherFirstDaughter);
  ((TH1F*)fMotherHistogramArray[motherType][histType][21])->Fill(angleMotherSecondDaughter);
  ((TH1F*)fMotherHistogramArray[motherType][histType][22])->Fill(angleBetweenBothDaughters);
  ((TH1F*)fMotherHistogramArray[motherType][histType][23])->Fill(cosThetaStar);

  ((TH1F*)fMotherHistogramArray[motherType][histType][24])->Fill(TMath::Abs(selectedMother->Getd0Prong(0) / selectedMother->Getd0errProng(0)));
  ((TH1F*)fMotherHistogramArray[motherType][histType][25])->Fill(TMath::Abs(selectedMother->Getd0Prong(1) / selectedMother->Getd0errProng(1)));
  ((TH1F*)fMotherHistogramArray[motherType][histType][26])->Fill(TMath::Abs(d0[0] / TMath::Sqrt(covd0z0[0])));
  ((TH1F*)fMotherHistogramArray[motherType][histType][27])->Fill((selectedMother->Getd0Prong(0) / selectedMother->Getd0errProng(0)) * (selectedMother->Getd0Prong(1) / selectedMother->Getd0errProng(1)));


  ((TH1F*)fMotherHistogramArray[motherType][histType][49])->Fill(selectedMother->Y(511));
  ((TH1F*)fMotherHistogramArray[motherType][histType][50])->Fill(etaMother);
  ((TH1F*)fMotherHistogramArray[motherType][histType][51])->Fill(phiMother);
  ((TH1F*)fMotherHistogramArray[motherType][histType][52])->Fill(vertexMotherX);
  ((TH1F*)fMotherHistogramArray[motherType][histType][53])->Fill(vertexMotherY);
  ((TH1F*)fMotherHistogramArray[motherType][histType][54])->Fill(vertexMotherZ);
  ((TH1F*)fMotherHistogramArray[motherType][histType][55])->Fill(vertexMother->GetChi2perNDF());

  //we fill the 2D histograms // to be done
  // Int_t nFirst = 0;
  // Int_t nSecond = 1;
  // Int_t nVariables = 10;
  // Int_t nHistograms = nVariables * (nVariables - 1) / 2;
  // for (Int_t k = 0; k < nHistograms; ++k)
  // {
  //   Double_t firstVariable = 0.0;
  //   Double_t secondVariable = 0.0;

  //   if (nFirst == 0) firstVariable = d0firstTrack;
  //   if (nFirst == 1) firstVariable = d0secondTrack;
  //   if (nFirst == 2) firstVariable = d0Mother;
  //   if (nFirst == 3) firstVariable = pointingAngle;
  //   if (nFirst == 4) firstVariable = impactProduct;
  //   if (nFirst == 5) firstVariable = impactProductXY;
  //   if (nFirst == 6) firstVariable = vertexDistance;
  //   if (nFirst == 7) firstVariable = normDecayLength;
  //   if (nFirst == 8) firstVariable = cosPointingAngleXY;
  //   if (nFirst == 9) firstVariable = distanceXYToVertex;
  //   if (nFirst == 10) firstVariable = normalizedDecayLengthXY;

  //   if (nSecond == 0) secondVariable = d0firstTrack;
  //   if (nSecond == 1) secondVariable = d0secondTrack;
  //   if (nSecond == 2) secondVariable = d0Mother;
  //   if (nSecond == 3) secondVariable = pointingAngle;
  //   if (nSecond == 4) secondVariable = impactProduct;
  //   if (nSecond == 5) secondVariable = impactProductXY;
  //   if (nSecond == 6) secondVariable = vertexDistance;
  //   if (nSecond == 7) secondVariable = normDecayLength;
  //   if (nSecond == 8) secondVariable = cosPointingAngleXY;
  //   if (nSecond == 9) secondVariable = distanceXYToVertex;
  //   if (nSecond == 10) secondVariable = normalizedDecayLengthXY;

  //   ((TH2F*)fMotherHistogramArray2D[motherType][histType][k])->Fill(firstVariable, secondVariable);

  //   nSecond++;
  //   if (nSecond > nVariables)
  //   {
  //     nFirst++;
  //     nSecond = nFirst + 1;
  //   }
  // }


  if (motherType == 1) {
    motherType = motherType - 1;

    AliAODRecoDecay* trackDPlus = (AliAODRecoDecay*)selectedMother->GetDaughter(0);
    AliAODTrack * firstDaughterDPlus = (AliAODTrack*)trackDPlus->GetDaughter(0);
    AliAODTrack * secondDaughterDPlus = (AliAODTrack*)trackDPlus->GetDaughter(1);
    AliAODTrack * thirdDaughterDPlus = (AliAODTrack*)trackDPlus->GetDaughter(2);

    AliAODVertex * vertexB0 = vertexMother;
    AliAODVertex * vertexDPlus = trackDPlus->GetSecondaryVtx();
    vertexDistance = TMath::Abs(vertexB0->DistanceToVertex(vertexDPlus));
    pdgMassMother = TDatabasePDG::Instance()->GetParticle(411)->Mass();

    AliExternalTrackParam firstDaughterDPlusTrack;
    AliExternalTrackParam secondDaughterDPlusTrack;
    AliExternalTrackParam thirdDaughterDPlusTrack;

    Double_t d0z0[2], covd0z0[3], d0[3];

    firstDaughterDPlusTrack.CopyFromVTrack(firstDaughterDPlus);
    firstDaughterDPlusTrack.PropagateToDCA(vertexB0, bz, 100., d0z0, covd0z0);
    d0[0] = d0z0[0];

    secondDaughterDPlusTrack.CopyFromVTrack(secondDaughterDPlus);
    secondDaughterDPlusTrack.PropagateToDCA(vertexB0, bz, 100., d0z0, covd0z0);
    d0[1] = d0z0[0];

    thirdDaughterDPlusTrack.CopyFromVTrack(thirdDaughterDPlus);
    thirdDaughterDPlusTrack.PropagateToDCA(vertexB0, bz, 100., d0z0, covd0z0);
    d0[2] = d0z0[0];

    AliExternalTrackParam DPlusTrack;
    DPlusTrack.CopyFromVTrack(trackDPlus);
    Double_t d0z0DPlus[2], covd0z0DPlus[3], d0DPlus;
    DPlusTrack.PropagateToDCA(vertexB0, bz, 100., d0z0DPlus, covd0z0DPlus);
    d0DPlus = d0z0DPlus[0];

    Double_t impactProductToB0 = d0[0] * d0[1] * d0[2];
    Double_t impactProductXYToB0 = trackDPlus->ImpParXY(vertexB0);

    Double_t momentumMother = trackDPlus->P();
    Double_t pointingAngleToB0 = trackDPlus->CosPointingAngle(vertexB0);
    Double_t d0FirstDaughterToB0 = TMath::Abs(d0[0]);
    Double_t d0SecondDaughterToB0 = TMath::Abs(d0[1]);
    Double_t d0ThirdDaughterToB0 = TMath::Abs(d0[2]);

    Double_t normDecayLengthToB0 = trackDPlus->NormalizedDecayLength(vertexB0);

    Double_t pseudoProperDecayLength = ((vertexDPlus->GetX() - vertexB0->GetX()) * trackDPlus->Px() / TMath::Abs(trackDPlus->Pt())) + ((vertexDPlus->GetY() - vertexB0->GetY()) * trackDPlus->Py() / TMath::Abs(trackDPlus->Pt()));
    Double_t pseudoProperDecayTimeToB0 = pseudoProperDecayLength * pdgMassMother / ptMother;
    Double_t DecayTimeToB0 = vertexDistance / (299792458 * TMath::Sqrt(1 / ((pdgMassMother * pdgMassMother / (momentumMother * momentumMother)) + 1)));

    Double_t phi = trackDPlus->Phi();
    Double_t theta = trackDPlus->Theta();
    Double_t covMatrix[21];
    trackDPlus->GetCovarianceXYZPxPyPz(covMatrix);

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
    Double_t normDecayTimeToB0 = trackDPlus->NormalizedDecayLength(vertexB0) / (299792458 * TMath::Sqrt(1 / ((pdgMassMother * pdgMassMother * errorMomentum * errorMomentum / (momentumMother * momentumMother)) + 1)));

    ((TH1F*)fMotherHistogramArray[motherType][histType][38])->Fill(pointingAngleToB0);
    ((TH1F*)fMotherHistogramArray[motherType][histType][39])->Fill(d0DPlus);
    ((TH1F*)fMotherHistogramArray[motherType][histType][40])->Fill(d0FirstDaughterToB0);
    ((TH1F*)fMotherHistogramArray[motherType][histType][41])->Fill(d0SecondDaughterToB0);
    ((TH1F*)fMotherHistogramArray[motherType][histType][42])->Fill(d0ThirdDaughterToB0);
    ((TH1F*)fMotherHistogramArray[motherType][histType][43])->Fill(impactProductToB0);
    ((TH1F*)fMotherHistogramArray[motherType][histType][44])->Fill(impactProductXYToB0);
    ((TH1F*)fMotherHistogramArray[motherType][histType][45])->Fill(normDecayLengthToB0);
    ((TH1F*)fMotherHistogramArray[motherType][histType][46])->Fill(pseudoProperDecayTimeToB0);
    ((TH1F*)fMotherHistogramArray[motherType][histType][47])->Fill(DecayTimeToB0);
    ((TH1F*)fMotherHistogramArray[motherType][histType][48])->Fill(normDecayTimeToB0);



  }
  return;
}
//-------------------------------------------------------------------------------------
void AliAnalysisTaskSEB0toDPi::TwoTrackCombinationInfo(AliExternalTrackParam * firstTrack, AliExternalTrackParam * secondTrack, AliAODVertex * primaryVertex, Double_t bz, Bool_t isDesiredCandidate, Int_t histogramNumber, UInt_t prongs[2]) {

  // we calculate the vertex position
  TObjArray daughterTracks;

  daughterTracks.Add(firstTrack);
  daughterTracks.Add(secondTrack);

  Double_t dispersion = 0;
  AliAODVertex *vertex = RecalculateVertex(primaryVertex, &daughterTracks, bz, dispersion, kTRUE, kTRUE, kFALSE, 3);
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
void AliAnalysisTaskSEB0toDPi::ThreeTrackCombinationInfo(AliExternalTrackParam * firstTrack, AliExternalTrackParam * secondTrack, AliExternalTrackParam * thirdTrack, AliAODVertex * primaryVertex, Double_t bz, Bool_t isDesiredCandidate, Int_t histogramNumber, UInt_t prongs[3]) {

  // we calculate the vertex position
  TObjArray daughterTracks;

  daughterTracks.Add(firstTrack);
  daughterTracks.Add(secondTrack);
  daughterTracks.Add(thirdTrack);

  Double_t dispersion = 0;
  AliAODVertex *vertex = RecalculateVertex(primaryVertex, &daughterTracks, bz, dispersion, kTRUE, kTRUE, kFALSE, 3);
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
Int_t AliAnalysisTaskSEB0toDPi::MatchCandidateToMonteCarlo(Int_t pdgabs, AliAODRecoDecayHF * candidate, TClonesArray *mcArray, TMatrix * B0toDPiLabelMatrix, Bool_t bCheckLabel) const
{
  //
  // Check if this candidate is matched to a MC signal
  // If no, return -1
  // If yes, return label (>=0) of the AliAODMCParticle

  // Check number of daughters
  Int_t ndg = candidate->GetNDaughters();
  if (!ndg) { AliError("No daughters available"); return -1;}

  // loop on daughters and write the labels
  Int_t dgLabels[3] = { -1};
  Int_t pdgDg[3] = {0};
  Int_t signalPosition = -1;
  if (pdgabs == 411)
  {
    if (ndg != 3) return -1;
    AliAODTrack *trk0 = (AliAODTrack*)candidate->GetDaughter(0);
    dgLabels[0] = trk0->GetLabel();
    AliAODTrack *trk1 = (AliAODTrack*)candidate->GetDaughter(1);
    dgLabels[1] = trk1->GetLabel();
    AliAODTrack *trk2 = (AliAODTrack*)candidate->GetDaughter(2);
    dgLabels[2] = trk2->GetLabel();
    pdgDg[0] = 211; pdgDg[1] = 321; pdgDg[2] = 211;
    signalPosition = 4;
  }
  else if (pdgabs == 511)
  {
    if (ndg != 2) return -1;
    dgLabels[0] = MatchCandidateToMonteCarlo(411, (AliAODRecoDecayHF3Prong*)candidate->GetDaughter(0), mcArray, B0toDPiLabelMatrix);
    AliAODTrack *trk1 = (AliAODTrack*)candidate->GetDaughter(1);
    dgLabels[1] = trk1->GetLabel();
    pdgDg[0] = 411; pdgDg[1] = 211; pdgDg[2] = -1;
    signalPosition = 5;
  }
  else
  {
    std::cout << "Wrong pdg supplied for function to match candidate to monte carlo signal." << std::endl;
    return -1;
  }
  if (dgLabels[0] == -1) return -1;
  if (dgLabels[1] == -1) return -1;

  Int_t labMom[3] = {0, 0, 0};
  Int_t i, j, lab, labMother, pdgMother, pdgPart;
  AliAODMCParticle *part = 0;
  AliAODMCParticle *mother = 0;
  Double_t pxSumDgs = 0., pySumDgs = 0., pzSumDgs = 0.;
  Bool_t pdgUsed[3] = {kFALSE, kFALSE, kFALSE};

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
  if (pdgabs == 411  && mother->GetNDaughters() != 3) return -1;
  if (pdgabs == 511  && mother->GetNDaughters() != 2) return -1;

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
  if (bCheckLabel)
  {
    Int_t bIsSignal = kFALSE;
    TMatrix &particleMatrix = *B0toDPiLabelMatrix;
    for (Int_t k = 0; k < B0toDPiLabelMatrix->GetNrows(); ++k)
    {
      if (labMother == (Int_t)particleMatrix(k, signalPosition))
      {
        bIsSignal = kTRUE;
        break;
      }
    }
    if (!bIsSignal) return -1;
  }

  return labMother;
}
//-------------------------------------------------------------------------------------
Int_t AliAnalysisTaskSEB0toDPi::IsTrackInjected(AliAODTrack *part, AliAODMCHeader *header, TClonesArray *arrayMC) {

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
Bool_t AliAnalysisTaskSEB0toDPi::IsCandidateInjected(AliAODRecoDecayHF2Prong *selectedB0, AliAODMCHeader *header, TClonesArray *arrayMC) {


  AliAODRecoDecayHF3Prong* selectedDPlus = (AliAODRecoDecayHF3Prong*)selectedB0->GetDaughter(0);
  AliAODTrack* selectedB0Pion = (AliAODTrack*)selectedB0->GetDaughter(1);

  AliAODTrack* selectedDPlusFirstDaughter = (AliAODTrack*)selectedDPlus->GetDaughter(0);
  AliAODTrack* selectedDPlusSecondDaughter = (AliAODTrack*)selectedDPlus->GetDaughter(1);
  AliAODTrack* selectedDPlusThirdDaughter = (AliAODTrack*)selectedDPlus->GetDaughter(2);

  if (IsTrackInjected(selectedB0Pion, header, arrayMC)) return kTRUE;
  if (IsTrackInjected(selectedDPlusFirstDaughter, header, arrayMC)) return kTRUE;
  if (IsTrackInjected(selectedDPlusSecondDaughter, header, arrayMC)) return kTRUE;
  if (IsTrackInjected(selectedDPlusThirdDaughter, header, arrayMC)) return kTRUE;

  return kFALSE;
}
//-------------------------------------------------------------------------------------
void AliAnalysisTaskSEB0toDPi::CutOptimizationLoop(Int_t variable, Int_t nVariables, Int_t nCuts, Int_t ptBin, Int_t fillNumber, Bool_t isDesiredCandidate, Int_t nSigmaBin) {

  for (Int_t iCut = 0; iCut < nCuts; ++iCut)
  {
    Int_t cutVariable = fCuts->GetCutIndexForCutOptimization(variable);
    Bool_t isUpperCut = fCuts->GetIsUpperCutForCutOptimization(variable);
    Float_t cutValue = fCuts->GetCutForCutOptimization(iCut, variable, ptBin);
    if (fCutVariableValueArray[cutVariable] > cutValue && isUpperCut == kTRUE) continue;
    if (fCutVariableValueArray[cutVariable] < cutValue && isUpperCut == kFALSE) continue;

    Int_t fill = iCut * TMath::Power(nCuts, variable) + fillNumber;
    if ( (variable + 1) == nVariables)
    {
      if (isDesiredCandidate) ((TH2F*)(fResultsHistogramArray2D[0][ptBin]))->Fill(fill, nSigmaBin);
      if (!isDesiredCandidate) ((TH2F*)(fResultsHistogramArray2D[1][ptBin]))->Fill(fill, nSigmaBin);
    }
    else {CutOptimizationLoop(variable + 1, nVariables, nCuts, ptBin, fill, isDesiredCandidate, nSigmaBin);}
  }
  return;
}
//-------------------------------------------------------------------------------------
void AliAnalysisTaskSEB0toDPi::CutOptimizationVariableValues(AliAODRecoDecayHF2Prong * candidateB0, AliAODEvent*  aod) {

  if (!candidateB0) {
    std::cout << "candidateB0 null" << std::endl;
    return;
  }

  AliAODRecoDecayHF3Prong* candidateDPlus = (AliAODRecoDecayHF3Prong*)candidateB0->GetDaughter(0);
  if (!candidateDPlus) {
    std::cout << "candidateDPlus null" << std::endl;
    return;
  }

  AliAODTrack *candidateB0Pion = (AliAODTrack*)candidateB0->GetDaughter(1);
  if (!candidateB0Pion) {
    std::cout << "candidateB0Pion null" << std::endl;
    return;
  }

  AliAODTrack *candidateFirstDaughter = (AliAODTrack*)candidateDPlus->GetDaughter(0);
  if (!candidateFirstDaughter) {
    std::cout << "candidatePion null" << std::endl;
    return;
  }

  AliAODTrack *candidateSecondDaughter = (AliAODTrack*)candidateDPlus->GetDaughter(1);
  if (!candidateSecondDaughter) {
    std::cout << "candidateKaon null" << std::endl;
    return;
  }

  AliAODTrack *candidateThirdDaughter = (AliAODTrack*)candidateDPlus->GetDaughter(2);
  if (!candidateSecondDaughter) {
    std::cout << "candidatePion null" << std::endl;
    return;
  }

  AliAODVertex * vertexB0 = candidateB0->GetSecondaryVtx();
  if (!vertexB0) {
    std::cout << "vertexB0 null" << std::endl;
    return;
  }

  AliAODVertex * vertexDPlus = candidateDPlus->GetSecondaryVtx();
  if (!vertexDPlus) {
    std::cout << "vertexDPlus null" << std::endl;
    return;
  }

  AliAODVertex * primaryVertex = aod->GetPrimaryVertex();
  if (!primaryVertex) {
    std::cout << "primaryVertex null" << std::endl;
    return;
  }

  //get the magnetic field
  Double_t bz = (Double_t)aod->GetMagneticField();


  // save DPlus variable information
  if (kTRUE)
  {
    // DPlusmass
    Double_t mDPlusPDG = TDatabasePDG::Instance()->GetParticle(411)->Mass();
    Double_t mPionPDG = TDatabasePDG::Instance()->GetParticle(211)->Mass();

    // DPlus window - invariant mass
    UInt_t prongs[3];
    prongs[0] = 211;
    prongs[1] = 321;
    prongs[2] = 211;

    // delta mass PDG
    Double_t deltaPDG = mDPlusPDG - 2 * mPionPDG;

    Double_t invMassDPlus = candidateDPlus->InvMass(3, prongs);
    Double_t invMassDifference = TMath::Abs(mDPlusPDG - invMassDPlus);
    Double_t eKaon1 = candidateDPlus->EProng(1, 321);
    Double_t ePion2 = candidateDPlus->EProng(2, 211);
    Double_t eSum = eKaon1 + ePion2;
    Double_t pxPions = (candidateDPlus->PxProng(1) + candidateDPlus->PxProng(2)) * (candidateDPlus->PxProng(1) + candidateDPlus->PxProng(2));
    Double_t pyPions = (candidateDPlus->PyProng(1) + candidateDPlus->PyProng(2)) * (candidateDPlus->PyProng(1) + candidateDPlus->PyProng(2));
    Double_t pzPions = (candidateDPlus->PzProng(1) + candidateDPlus->PzProng(2)) * (candidateDPlus->PzProng(1) + candidateDPlus->PzProng(2));
    Double_t pSum = pxPions + pyPions + pzPions;
    Double_t invMassPions = TMath::Sqrt(eSum * eSum - pSum);
    Double_t invMassDelta = invMassDPlus - invMassPions;
    Double_t invMassDeltaDifference = TMath::Abs(deltaPDG - invMassDelta);

    Double_t pointingAngle = candidateDPlus->CosPointingAngle();
    Double_t dcaMother12 = candidateDPlus->GetDCA(0);
    Double_t dcaMother13 = candidateDPlus->GetDCA(1);
    Double_t dcaMother23 = candidateDPlus->GetDCA(2);
    Double_t ptMother = candidateDPlus->Pt();
    Double_t momentumMother = candidateDPlus->P();
    Double_t ptKaon = candidateSecondDaughter->Pt();
    Double_t ptFirstPion = candidateFirstDaughter->Pt();
    Double_t ptSecondPion = candidateThirdDaughter->Pt();

    AliExternalTrackParam motherTrack;
    motherTrack.CopyFromVTrack(candidateDPlus);
    Double_t d0z0[2], covd0z0[3], d0[2];
    motherTrack.PropagateToDCA(primaryVertex, bz, 100., d0z0, covd0z0);
    d0[0] = d0z0[0];
    Double_t d0Mother = TMath::Abs(d0[0]);
    Double_t d0firstTrack = TMath::Abs(candidateDPlus->Getd0Prong(0));
    Double_t d0secondTrack = TMath::Abs(candidateDPlus->Getd0Prong(1));
    Double_t d0thirdTrack = TMath::Abs(candidateDPlus->Getd0Prong(2));

    Double_t impactProduct123 = candidateDPlus->Getd0Prong(0) * candidateDPlus->Getd0Prong(1) * candidateDPlus->Getd0Prong(2);
    Double_t impactProduct12 = candidateDPlus->Getd0Prong(0) * candidateDPlus->Getd0Prong(1);
    Double_t impactProduct13 = candidateDPlus->Getd0Prong(0) * candidateDPlus->Getd0Prong(2);
    Double_t impactProduct23 = candidateDPlus->Getd0Prong(1) * candidateDPlus->Getd0Prong(2);
    Double_t impactProductXY = TMath::Abs(candidateDPlus->ImpParXY());

    Double_t angleMotherFirstDaughter = (candidateDPlus->Px() * candidateFirstDaughter->Px() + candidateDPlus->Py() * candidateFirstDaughter->Py() + candidateDPlus->Pz() * candidateFirstDaughter->Pz()) / (candidateDPlus->P() * candidateFirstDaughter->P());
    Double_t angleMotherSecondDaughter = (candidateDPlus->Px() * candidateSecondDaughter->Px() + candidateDPlus->Py() * candidateSecondDaughter->Py() + candidateDPlus->Pz() * candidateSecondDaughter->Pz()) / (candidateDPlus->P() * candidateSecondDaughter->P());
    Double_t angleMotherThirdDaughter = (candidateDPlus->Px() * candidateThirdDaughter->Px() + candidateDPlus->Py() * candidateThirdDaughter->Py() + candidateDPlus->Pz() * candidateThirdDaughter->Pz()) / (candidateDPlus->P() * candidateThirdDaughter->P());

    Double_t smallestAngleMotherDaughter = 0.0;
    Double_t largestAngleMotherDaughter = 0.0;

    if (angleMotherFirstDaughter > angleMotherSecondDaughter) {smallestAngleMotherDaughter = angleMotherSecondDaughter;} else {smallestAngleMotherDaughter = angleMotherFirstDaughter;}
    if (angleMotherThirdDaughter < smallestAngleMotherDaughter) smallestAngleMotherDaughter = angleMotherThirdDaughter;
    if (angleMotherFirstDaughter > angleMotherSecondDaughter) {largestAngleMotherDaughter = angleMotherFirstDaughter;} else {largestAngleMotherDaughter = angleMotherSecondDaughter;}
    if (angleMotherThirdDaughter > largestAngleMotherDaughter) largestAngleMotherDaughter = angleMotherThirdDaughter;

    Double_t vertexDistance = vertexDPlus->DistanceToVertex(primaryVertex);
    Double_t normDecayLength = candidateDPlus->NormalizedDecayLength();
    Double_t pdgMassMother = TDatabasePDG::Instance()->GetParticle(411)->Mass();
    Double_t pseudoProperDecayLength = ((vertexDPlus->GetX() - primaryVertex->GetX()) * candidateDPlus->Px() / TMath::Abs(candidateDPlus->Pt())) + ((vertexDPlus->GetY() - primaryVertex->GetY()) * candidateDPlus->Py() / TMath::Abs(candidateDPlus->Pt()));
    Double_t pseudoProperDecayTime = pseudoProperDecayLength * pdgMassMother / ptMother;
    Double_t decayTime = vertexDistance / (299792458 * TMath::Sqrt(1 / ((pdgMassMother * pdgMassMother / (momentumMother * momentumMother)) + 1)));

    Double_t phi = candidateDPlus->Phi();
    Double_t theta = candidateDPlus->Theta();
    Double_t covMatrix[21];
    candidateDPlus->GetCovarianceXYZPxPyPz(covMatrix);

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
    Double_t normalizedDecayTime = candidateDPlus->NormalizedDecayLength() / (299792458 * TMath::Sqrt(1 / ((pdgMassMother * pdgMassMother * errorMomentum * errorMomentum / (momentumMother * momentumMother)) + 1)));

    Double_t cosPointingAngleXY = candidateDPlus->CosPointingAngleXY();
    Double_t distanceXYToVertex = vertexDPlus->DistanceXYToVertex(primaryVertex);
    Double_t normalizedDecayLengthXY = candidateDPlus->NormalizedDecayLengthXY();
    Double_t chi2Vertex = vertexDPlus->GetChi2perNDF();


    Double_t Normalizedd0DPlusfirstdaughter = TMath::Abs(candidateDPlus->Getd0Prong(0) / candidateDPlus->Getd0errProng(0));
    Double_t Normalizedd0DPlusseconddaughter = TMath::Abs(candidateDPlus->Getd0Prong(1) / candidateDPlus->Getd0errProng(1));
    Double_t Normalizedd0DPlusthirddaughter = TMath::Abs(candidateDPlus->Getd0Prong(2) / candidateDPlus->Getd0errProng(2));
    Double_t Normalizedd0DPlus = TMath::Abs(candidateB0->Getd0Prong(1) / candidateB0->Getd0errProng(1));
    Double_t NormalizedimpactproductDPlus = (candidateDPlus->Getd0Prong(0) / candidateDPlus->Getd0errProng(0)) * (candidateDPlus->Getd0Prong(1) / candidateDPlus->Getd0errProng(1)) * (candidateDPlus->Getd0Prong(2) / candidateDPlus->Getd0errProng(2));

    Double_t Dist12 = candidateDPlus->GetDist12toPrim();
    Double_t Dist23 = candidateDPlus->GetDist23toPrim();
    Double_t SigmaVertex = candidateDPlus->GetSigmaVert();

    // We apply the cuts
    Int_t nCutIndex = 0;

    // "inv. mass width [GeV]" --------------------------------------------
    nCutIndex = 0;
    fCutVariableValueArray[nCutIndex] = invMassDifference;
    //---------------------------------------------------------------------

    // "delta mass width [GeV]" -------------------------------------------
    nCutIndex = 1;
    fCutVariableValueArray[nCutIndex] = invMassDeltaDifference;
    //---------------------------------------------------------------------

    // "pointing angle [Cos(theta)]" --------------------------------------
    nCutIndex = 2;
    fCutVariableValueArray[nCutIndex] = pointingAngle;
    //---------------------------------------------------------------------

    // "pointing angle XY" ----------------------------------------------------
    nCutIndex = 3;
    fCutVariableValueArray[nCutIndex] = cosPointingAngleXY;
    //---------------------------------------------------------------------

    // "dca12 (cm)" ---------------------------------------------------------
    nCutIndex = 4;
    fCutVariableValueArray[nCutIndex] = dcaMother12;
    //---------------------------------------------------------------------

    // "dca13 (cm)" ---------------------------------------------------------
    nCutIndex = 5;
    fCutVariableValueArray[nCutIndex] = dcaMother13;
    //---------------------------------------------------------------------

    // "dca23 (cm)" ---------------------------------------------------------
    nCutIndex = 6;
    fCutVariableValueArray[nCutIndex] = dcaMother23;
    //---------------------------------------------------------------------

    // "Pt DPlus [GeV/c]" ----------------------------------------------------
    nCutIndex = 7;
    fCutVariableValueArray[nCutIndex] = ptMother;
    //---------------------------------------------------------------------

    // "Pt Kaon [GeV/c]" -------------------------------------------------
    nCutIndex = 8;
    fCutVariableValueArray[nCutIndex] = ptKaon;
    //---------------------------------------------------------------------

    // "Pt First Pion [GeV/c]" --------------------------------------------------
    nCutIndex = 9;
    fCutVariableValueArray[nCutIndex] = ptFirstPion;
    //---------------------------------------------------------------------

    // "Pt Second Pion [GeV/c]" --------------------------------------------------
    nCutIndex = 10;
    fCutVariableValueArray[nCutIndex] = ptSecondPion;
    //---------------------------------------------------------------------

    // "d0 DPlus (cm)" -------------------------------------------------------
    nCutIndex = 11;
    fCutVariableValueArray[nCutIndex] = d0Mother;
    //---------------------------------------------------------------------

    // "d0 Kaon (cm)"-----------------------------------------------------
    nCutIndex = 12;
    fCutVariableValueArray[nCutIndex] = d0firstTrack;
    //---------------------------------------------------------------------

    // "d0 First Pion (cm)" -----------------------------------------------------
    nCutIndex = 13;
    fCutVariableValueArray[nCutIndex] = d0secondTrack;
    //---------------------------------------------------------------------

    // "d0 Second Pion (cm)" -----------------------------------------------------
    nCutIndex = 14;
    fCutVariableValueArray[nCutIndex] = d0thirdTrack;
    //---------------------------------------------------------------------

    // "d0d0d0 [cm^3]" ------------------------------------------------------
    nCutIndex = 15;
    fCutVariableValueArray[nCutIndex] = impactProduct123;
    //---------------------------------------------------------------------

    // "d0d0 [cm^2]" ------------------------------------------------------
    nCutIndex = 16;
    fCutVariableValueArray[nCutIndex] = impactProduct12;
    //---------------------------------------------------------------------

    // "d0d0 [cm^2]" ------------------------------------------------------
    nCutIndex = 17;
    fCutVariableValueArray[nCutIndex] = impactProduct13;
    //---------------------------------------------------------------------

    // "d0d0 [cm^2]" ------------------------------------------------------
    nCutIndex = 18;
    fCutVariableValueArray[nCutIndex] = impactProduct23;
    //---------------------------------------------------------------------

    // "d0 XY [cm^2]" ---------------------------------------------------
    nCutIndex = 19;
    fCutVariableValueArray[nCutIndex] = impactProductXY;
    //---------------------------------------------------------------------

    // "smallestAngleMotherDaughter" -------------------------------------
    nCutIndex = 20;
    fCutVariableValueArray[nCutIndex] = smallestAngleMotherDaughter;
    //---------------------------------------------------------------------

    // "largestAngleMotherDaughter" ---------------------------------
    nCutIndex = 21;
    fCutVariableValueArray[nCutIndex] = largestAngleMotherDaughter;
    //---------------------------------------------------------------------

    // "angle difference" --------------------------------
    nCutIndex = 22;
    fCutVariableValueArray[nCutIndex] = largestAngleMotherDaughter - smallestAngleMotherDaughter;
    //---------------------------------------------------------------------

    // "vertexDistance" ---------------------------------------------------
    nCutIndex = 23;
    fCutVariableValueArray[nCutIndex] = vertexDistance;
    //---------------------------------------------------------------------

    // "vertex distance XY" ----------------------------------------------------
    nCutIndex = 24;
    fCutVariableValueArray[nCutIndex] = distanceXYToVertex;
    //---------------------------------------------------------------------

    // "pseudoProperDecayTime" --------------------------------------------
    nCutIndex = 25;
    fCutVariableValueArray[nCutIndex] = pseudoProperDecayTime;
    //---------------------------------------------------------------------

    // "DecayTime" --------------------------------------------------------
    nCutIndex = 26;
    fCutVariableValueArray[nCutIndex] = decayTime;
    //---------------------------------------------------------------------

    // "normalizedDecayTime" ----------------------------------------------------
    nCutIndex = 27;
    fCutVariableValueArray[nCutIndex] = normalizedDecayTime;
    //---------------------------------------------------------------------

    // "normDecayLength" --------------------------------------------------
    nCutIndex = 28;
    fCutVariableValueArray[nCutIndex] = normDecayLength;
    //---------------------------------------------------------------------

    // "normalized decay length XY" ----------------------------------------------------
    nCutIndex = 29;
    fCutVariableValueArray[nCutIndex] = normalizedDecayLengthXY;
    //---------------------------------------------------------------------

    // "chi squared per NDF" ----------------------------------------------------
    nCutIndex = 30;
    fCutVariableValueArray[nCutIndex] = chi2Vertex;
    //---------------------------------------------------------------------

    // "Normalizedd0DPlusfirstdaughter" ----------------------------------------------------
    nCutIndex = 31;
    fCutVariableValueArray[nCutIndex] = Normalizedd0DPlusfirstdaughter;
    //---------------------------------------------------------------------

    // "Normalizedd0DPlusseconddaughter" ----------------------------------------------------
    nCutIndex = 32;
    fCutVariableValueArray[nCutIndex] = Normalizedd0DPlusseconddaughter;
    //---------------------------------------------------------------------

    // "Normalizedd0DPlusthirddaughter" ----------------------------------------------------
    nCutIndex = 33;
    fCutVariableValueArray[nCutIndex] = Normalizedd0DPlusthirddaughter;
    //---------------------------------------------------------------------

    // "Normalizedd0DPlus" ----------------------------------------------------
    nCutIndex = 34;
    fCutVariableValueArray[nCutIndex] = Normalizedd0DPlus;
    //---------------------------------------------------------------------

    // "NormalizedimpactproductDPlus" ----------------------------------------------------
    nCutIndex = 35;
    fCutVariableValueArray[nCutIndex] = NormalizedimpactproductDPlus;
    //---------------------------------------------------------------------

    // "Dist12" ----------------------------------------------------
    nCutIndex = 36;
    fCutVariableValueArray[nCutIndex] = Dist12;
    //---------------------------------------------------------------------

    // "Dist23" ----------------------------------------------------
    nCutIndex = 37;
    fCutVariableValueArray[nCutIndex] = Dist23;
    //---------------------------------------------------------------------

    // "Sigma vertex" ----------------------------------------------------
    nCutIndex = 38;
    fCutVariableValueArray[nCutIndex] = SigmaVertex;
    //---------------------------------------------------------------------



    AliAODRecoDecay* candidateDPlustoB0 = (AliAODRecoDecay*)candidateDPlus;
    AliExternalTrackParam firstDaughterDPlusTrack;
    AliExternalTrackParam secondDaughterDPlusTrack;
    AliExternalTrackParam thirdDaughterDPlusTrack;

    Double_t d0z0DSVert[2], covd0z0DSVert[3], d0DSVert[3];

    firstDaughterDPlusTrack.CopyFromVTrack(candidateFirstDaughter);
    firstDaughterDPlusTrack.PropagateToDCA(vertexB0, bz, 100., d0z0DSVert, covd0z0DSVert);
    d0DSVert[0] = d0z0DSVert[0];

    secondDaughterDPlusTrack.CopyFromVTrack(candidateSecondDaughter);
    secondDaughterDPlusTrack.PropagateToDCA(vertexB0, bz, 100., d0z0DSVert, covd0z0DSVert);
    d0DSVert[1] = d0z0DSVert[0];

    thirdDaughterDPlusTrack.CopyFromVTrack(candidateThirdDaughter);
    thirdDaughterDPlusTrack.PropagateToDCA(vertexB0, bz, 100., d0z0DSVert, covd0z0DSVert);
    d0DSVert[2] = d0z0DSVert[0];

    AliExternalTrackParam DPlusTrack;
    DPlusTrack.CopyFromVTrack(candidateDPlus);
    Double_t d0z0D0DSVert[2], covd0z0D0DSVert[3];
    motherTrack.PropagateToDCA(vertexB0, bz, 100., d0z0D0DSVert, covd0z0D0DSVert);
    Double_t d0d0DSVert = TMath::Abs(d0z0D0DSVert[0]);

    Double_t impactProductToB0 = d0DSVert[0] * d0DSVert[1] * d0DSVert[2];
    Double_t impactProductXYToB0 = candidateDPlustoB0->ImpParXY(vertexB0);

    Double_t pointingAngleToB0 = candidateDPlustoB0->CosPointingAngle(vertexB0);
    Double_t d0FirstDaughterToB0 = TMath::Abs(d0DSVert[0]);
    Double_t d0SecondDaughterToB0 = TMath::Abs(d0DSVert[1]);
    Double_t d0ThirdDaughterToB0 = TMath::Abs(d0DSVert[2]);
    Double_t normDecayLengthToB0 = candidateDPlustoB0->NormalizedDecayLength(vertexB0);

    Double_t pseudoProperDecayLengthDSVert = ((vertexDPlus->GetX() - vertexB0->GetX()) * candidateDPlus->Px() / TMath::Abs(candidateDPlus->Pt())) + ((vertexDPlus->GetY() - vertexB0->GetY()) * candidateDPlus->Py() / TMath::Abs(candidateDPlus->Pt()));
    Double_t pseudoProperDecayTimeToB0 = pseudoProperDecayLengthDSVert * pdgMassMother / ptMother;
    Double_t DecayTimeToB0 = vertexDistance / (299792458 * TMath::Sqrt(1 / ((pdgMassMother * pdgMassMother / (momentumMother * momentumMother)) + 1)));

    Double_t phiDSVert = candidateDPlus->Phi();
    Double_t thetaDSVert = candidateDPlus->Theta();
    Double_t covMatrixDSVert[21];
    candidateDPlus->GetCovarianceXYZPxPyPz(covMatrixDSVert);

    cp = TMath::Cos(phiDSVert);
    sp = TMath::Sin(phiDSVert);
    ct = TMath::Cos(thetaDSVert);
    st = TMath::Sin(thetaDSVert);

    errorMomentum = covMatrix[9] * cp * cp * ct * ct // GetCovPxPx
                    + covMatrix[13] * 2.*cp * sp * ct * ct // GetCovPxPy
                    + covMatrix[18] * 2.*cp * ct * st // GetCovPxPz
                    + covMatrix[14] * sp * sp * ct * ct // GetCovPyPy
                    + covMatrix[19] * 2.*sp * ct * st // GetCovPyPz
                    + covMatrix[20] * st * st; // GetCovPzPz
    Double_t normalizedDecayTimeToB0 = candidateDPlustoB0->NormalizedDecayLength(vertexB0) / (299792458 * TMath::Sqrt(1 / ((pdgMassMother * pdgMassMother * errorMomentum * errorMomentum / (momentumMother * momentumMother)) + 1)));

    // "pointingAngleToB0" ---------------------------------------------
    nCutIndex = 39;
    fCutVariableValueArray[nCutIndex] = pointingAngleToB0;
    //---------------------------------------------------------------------

    // "d0MotherToB0" --------------------------------------------------
    nCutIndex = 40;
    fCutVariableValueArray[nCutIndex] = d0d0DSVert;
    //---------------------------------------------------------------------

    // "d0FirstDaughterToB0" -------------------------------------------
    nCutIndex = 41;
    fCutVariableValueArray[nCutIndex] = d0FirstDaughterToB0;
    //---------------------------------------------------------------------

    // "d0SecondDaughterToB0" ------------------------------------------
    nCutIndex = 42;
    fCutVariableValueArray[nCutIndex] = d0SecondDaughterToB0;
    //---------------------------------------------------------------------

    // "d0ThirdDaughterToB0" ------------------------------------------
    nCutIndex = 43;
    fCutVariableValueArray[nCutIndex] = d0ThirdDaughterToB0;
    //---------------------------------------------------------------------

    // "impactProductToB0" ---------------------------------------------
    nCutIndex = 44;
    fCutVariableValueArray[nCutIndex] = impactProductToB0;
    //---------------------------------------------------------------------

    // "impactProductXYToB0" -------------------------------------------
    nCutIndex = 45;
    fCutVariableValueArray[nCutIndex] = impactProductXYToB0;
    //---------------------------------------------------------------------

    // "normDecayLengthToB0" -------------------------------------------
    nCutIndex = 46;
    fCutVariableValueArray[nCutIndex] = normDecayLengthToB0;
    //---------------------------------------------------------------------

    // "pseudoProperDecayTimeToB0" -------------------------------------
    nCutIndex = 47;
    fCutVariableValueArray[nCutIndex] = pseudoProperDecayTimeToB0;
    //---------------------------------------------------------------------

    // "DecayTimeToB0" -------------------------------------------------
    nCutIndex = 48;
    fCutVariableValueArray[nCutIndex] = DecayTimeToB0;
    //---------------------------------------------------------------------

    // "normalizedDecayTimeToB0" ---------------------------------------------
    nCutIndex = 49;
    fCutVariableValueArray[nCutIndex] = normalizedDecayTimeToB0;
    //---------------------------------------------------------------------
  }

  // save B0 variable information
  if (kTRUE)
  {
    // We obtain the variable values in the section below
    // DPlusMass and B0mass
    Double_t mDPlusPDG = TDatabasePDG::Instance()->GetParticle(411)->Mass();
    Double_t mB0PDG = TDatabasePDG::Instance()->GetParticle(511)->Mass();

    // delta mass PDG
    Double_t deltaPDG = mB0PDG - mDPlusPDG;

    // Half width B0 mass
    UInt_t prongs[2];
    prongs[0] = 411; prongs[1] = 211;
    Double_t invMassB0 = candidateB0->InvMass(2, prongs);
    Double_t invMassDifference = TMath::Abs(mB0PDG - invMassB0);
    Double_t invMassDelta = TMath::Abs(deltaPDG - (DeltaInvMassB0Kpipipi(candidateB0)));

    Double_t pointingAngle = candidateB0->CosPointingAngle();
    Double_t dcaMother = candidateB0->GetDCA();
    Double_t ptMother = candidateB0->Pt();
    Double_t momentumMother = candidateB0->P();
    Double_t ptDPlus = candidateDPlus->Pt();
    Double_t ptPion = candidateB0Pion->Pt();

    AliExternalTrackParam motherTrack;
    motherTrack.CopyFromVTrack(candidateB0);
    Double_t d0z0[2], covd0z0[3], d0[2];
    motherTrack.PropagateToDCA(primaryVertex, bz, 100., d0z0, covd0z0);
    d0[0] = d0z0[0];
    Double_t d0Mother = TMath::Abs(d0[0]);
    Double_t d0firstTrack = TMath::Abs(candidateB0->Getd0Prong(0));
    Double_t d0secondTrack = TMath::Abs(candidateB0->Getd0Prong(1));

    Double_t impactProduct = candidateB0->Getd0Prong(0) * candidateB0->Getd0Prong(1);
    Double_t impactProductXY = TMath::Abs(candidateB0->ImpParXY());

    Double_t angleMotherFirstDaughter = (candidateB0->Px() * candidateDPlus->Px() + candidateB0->Py() * candidateDPlus->Py() + candidateB0->Pz() * candidateDPlus->Pz()) / (candidateB0->P() * candidateDPlus->P());
    Double_t angleMotherSecondDaughter = (candidateB0->Px() * candidateB0Pion->Px() + candidateB0->Py() * candidateB0Pion->Py() + candidateB0->Pz() * candidateB0Pion->Pz()) / (candidateB0->P() * candidateB0Pion->P());

    Double_t angleBetweenBothDaughters  = (candidateDPlus->Px() * candidateB0Pion->Px() + candidateDPlus->Py() * candidateB0Pion->Py() + candidateDPlus->Pz() * candidateB0Pion->Pz()) / (candidateDPlus->P() * candidateB0Pion->P());

    Double_t cosThetaStar = candidateB0->CosThetaStar(0, 511, 411, 211);

    Double_t vertexDistance = vertexB0->DistanceToVertex(primaryVertex);
    Double_t normDecayLength = candidateB0->NormalizedDecayLength();
    Double_t pdgMassMother = TDatabasePDG::Instance()->GetParticle(511)->Mass();
    Double_t pseudoProperDecayLength = ((vertexB0->GetX() - primaryVertex->GetX()) * candidateB0->Px() / TMath::Abs(candidateB0->Pt())) + ((vertexB0->GetY() - primaryVertex->GetY()) * candidateB0->Py() / TMath::Abs(candidateB0->Pt()));
    Double_t pseudoProperDecayTime = pseudoProperDecayLength * pdgMassMother / ptMother;
    Double_t decayTime = vertexDistance / (299792458 * TMath::Sqrt(1 / ((pdgMassMother * pdgMassMother / (momentumMother * momentumMother)) + 1)));

    Double_t phi = candidateB0->Phi();
    Double_t theta = candidateB0->Theta();
    Double_t covMatrix[21];
    candidateB0->GetCovarianceXYZPxPyPz(covMatrix);

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
    Double_t normalizedDecayTime = candidateB0->NormalizedDecayLength() / (299792458 * TMath::Sqrt(1 / ((pdgMassMother * pdgMassMother * errorMomentum * errorMomentum / (momentumMother * momentumMother)) + 1)));

    Double_t cosPointingAngleXY = candidateB0->CosPointingAngleXY();
    Double_t distanceXYToVertex = vertexB0->DistanceXYToVertex(primaryVertex);
    Double_t normalizedDecayLengthXY = candidateB0->NormalizedDecayLengthXY();
    Double_t chi2Vertex = vertexB0->GetChi2perNDF();

    Double_t Normalizedd0B0pion = TMath::Abs(candidateB0->Getd0Prong(1) / candidateB0->Getd0errProng(1));
    Double_t Normalizedd0B0 = TMath::Abs(d0[0] / TMath::Sqrt(covd0z0[0]));
    Double_t NormalizedimpactproductB0 = (candidateB0->Getd0Prong(0) / candidateB0->Getd0errProng(0)) * (candidateB0->Getd0Prong(1) / candidateB0->Getd0errProng(1));


    // We apply the cuts
    Int_t nCutIndex = 0;

    // "inv. mass width [GeV]" --------------------------------------------
    nCutIndex = 50;
    fCutVariableValueArray[nCutIndex] = invMassDifference;
    //---------------------------------------------------------------------

    // "delta mass width [GeV]" -------------------------------------------
    nCutIndex = 51;
    fCutVariableValueArray[nCutIndex] = invMassDelta;
    //---------------------------------------------------------------------

    // "pointing angle [Cos(theta)]" --------------------------------------
    nCutIndex = 52;
    fCutVariableValueArray[nCutIndex] = pointingAngle;
    //---------------------------------------------------------------------

    // "pointing angle XY" ----------------------------------------------------
    nCutIndex = 53;
    fCutVariableValueArray[nCutIndex] = cosPointingAngleXY;
    //---------------------------------------------------------------------

    // "dca12 (cm)" ---------------------------------------------------------
    nCutIndex = 54;
    fCutVariableValueArray[nCutIndex] = dcaMother;
    //---------------------------------------------------------------------

    // "Pt B0 [GeV/c]" ----------------------------------------------------
    nCutIndex = 55;
    fCutVariableValueArray[nCutIndex] = ptMother;
    //---------------------------------------------------------------------

    // "Pt DPlus [GeV/c]" -------------------------------------------------
    nCutIndex = 56;
    fCutVariableValueArray[nCutIndex] = ptDPlus;
    //---------------------------------------------------------------------

    // "Pt Pion [GeV/c]" --------------------------------------------------
    nCutIndex = 57;
    fCutVariableValueArray[nCutIndex] = ptPion;
    //---------------------------------------------------------------------

    // "d0 B0 (cm)" -------------------------------------------------------
    nCutIndex = 58;
    fCutVariableValueArray[nCutIndex] = d0Mother;
    //---------------------------------------------------------------------

    // "d0 DPlus (cm)"-----------------------------------------------------
    nCutIndex = 59;
    fCutVariableValueArray[nCutIndex] = d0firstTrack;
    //---------------------------------------------------------------------

    // "d0 First Pion (cm)" -----------------------------------------------------
    nCutIndex = 60;
    fCutVariableValueArray[nCutIndex] = d0secondTrack;
    //---------------------------------------------------------------------

    // "d0d0d0 [cm^2]" ------------------------------------------------------
    nCutIndex = 61;
    fCutVariableValueArray[nCutIndex] = impactProduct;
    //---------------------------------------------------------------------

    // "d0 XY [cm^2]" ---------------------------------------------------
    nCutIndex = 62;
    fCutVariableValueArray[nCutIndex] = impactProductXY;
    //---------------------------------------------------------------------

    // "angleMotherFirstDaughter" -------------------------------------
    nCutIndex = 63;
    fCutVariableValueArray[nCutIndex] = angleMotherFirstDaughter;
    //---------------------------------------------------------------------

    // "angleMotherSecondDaughter" ---------------------------------
    nCutIndex = 64;
    fCutVariableValueArray[nCutIndex] = angleMotherSecondDaughter;
    //---------------------------------------------------------------------

    // "angleBetweenBothDaughters" --------------------------------
    nCutIndex = 65;
    fCutVariableValueArray[nCutIndex] = angleBetweenBothDaughters;
    //---------------------------------------------------------------------

    // "cosThetaStar" --------------------------------
    nCutIndex = 66;
    fCutVariableValueArray[nCutIndex] = cosThetaStar;
    //---------------------------------------------------------------------

    // "vertexDistance" ---------------------------------------------------
    nCutIndex = 67;
    fCutVariableValueArray[nCutIndex] = vertexDistance;
    //---------------------------------------------------------------------

    // "vertex distance XY" ----------------------------------------------------
    nCutIndex = 68;
    fCutVariableValueArray[nCutIndex] = distanceXYToVertex;
    //---------------------------------------------------------------------

    // "pseudoProperDecayTime" --------------------------------------------
    nCutIndex = 69;
    fCutVariableValueArray[nCutIndex] = pseudoProperDecayTime;
    //---------------------------------------------------------------------

    // "DecayTime" --------------------------------------------------------
    nCutIndex = 70;
    fCutVariableValueArray[nCutIndex] = decayTime;
    //---------------------------------------------------------------------

    // "normalizedDecayTime" ----------------------------------------------------
    nCutIndex = 71;
    fCutVariableValueArray[nCutIndex] = normalizedDecayTime;
    //---------------------------------------------------------------------

    // "normDecayLength" --------------------------------------------------
    nCutIndex = 72;
    fCutVariableValueArray[nCutIndex] = normDecayLength;
    //---------------------------------------------------------------------

    // "normalized decay length XY" ----------------------------------------------------
    nCutIndex = 73;
    fCutVariableValueArray[nCutIndex] = normalizedDecayLengthXY;
    //---------------------------------------------------------------------

    // "chi squared per NDF" ----------------------------------------------------
    nCutIndex = 74;
    fCutVariableValueArray[nCutIndex] = chi2Vertex;

    // "Normalizedd0B0Pion" ----------------------------------------------------
    nCutIndex = 75;
    fCutVariableValueArray[nCutIndex] = Normalizedd0B0pion;
    //---------------------------------------------------------------------

    // "Normalizedd0B0" ----------------------------------------------------
    nCutIndex = 76;
    fCutVariableValueArray[nCutIndex] = Normalizedd0B0;
    //---------------------------------------------------------------------

    // "NormalizedimpactproductB0" ----------------------------------------------------
    nCutIndex = 77;
    fCutVariableValueArray[nCutIndex] = NormalizedimpactproductB0;
    //---------------------------------------------------------------------

  }
  return;
}

