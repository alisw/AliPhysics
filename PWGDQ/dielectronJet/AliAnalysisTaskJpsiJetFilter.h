/* Copyright(c) 1998-2019, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/**
 * J/psi filter for Jpsi in jets analysis
 * 
 * Impelement <merge> interface for AliAODExtension 
 * Based on AliAnalysisTaskDielectronFilter
 * Refer: AliAnalysisTaskUpcFilter/Replicator 
 * 
 * Filter event by dielectron candidates
 * Fill all central barrel branch without filter
 *
 * By: Yìtāo WÚ <yitao@cern.ch> @ USTC, CN
*/

#ifndef ALIANALYSISTASK_JPSIJET_FILTER_H
#define ALIANALYSISTASK_JPSIJET_FILTER_H

#include "TH1D.h"

#include "AliLog.h"
#include "AliAODHandler.h"
#include "AliAODExtension.h"
#include "AliESDInputHandler.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisManager.h"
#include "AliTriggerAnalysis.h"

#include "AliDielectron.h"
class AliDielectron;
#include "AliDielectronVarManager.h"
class AliDielectronVarManager;

#include "AliEmcalJet.h"
#include "AliJetContainer.h"

#include "AliAnalysisTaskSE.h"

class AliESDInputHandler;
class AliAODInputHandler;
class AliAnalysisManager;
class AliTriggerAnalysis;

class AliAnalysisTaskJpsiJetFilter : public AliAnalysisTaskSE {

public:
  AliAnalysisTaskJpsiJetFilter();
  AliAnalysisTaskJpsiJetFilter(const char* name);
  virtual ~AliAnalysisTaskJpsiJetFilter();

  virtual void Init();
  virtual void UserExec(Option_t* option);

  Bool_t IsToMerge() { return fIsToMerge;}
  void SetToMerge(Bool_t isToMerge = kTRUE){ fIsToMerge = isToMerge;}
  Bool_t IsToReplace() { return fIsToReplace;}
  void SetToReplace(Bool_t isToReplace = kTRUE){ fIsToReplace = isToReplace;}

   Bool_t RejectPileup(){return fRejectPileup;}
   void SetRejectPileup(Bool_t p = kTRUE){fRejectPileup = p;}

  void InitDielectron();
  void InitHistogramsForDielectron(const char* histMgrName);

private:
  Bool_t fIsToMerge; // Option for AliAODExtension
  Bool_t fIsToReplace; // Option for e+e- Pair-Track replacing
  TString fOutputFileName;  // File name of filtered AOD file
  AliAODExtension* fExtAOD; // Filtered nano AOD

  AliAODTracklets* fSPD;
  AliAODCaloTrigger* fEMCALTrigger;
  AliAODCaloTrigger* fPHOSTrigger;
  AliAODCaloCells* fEMCalCells;
  AliAODCaloCells* fPHOSCells;
  AliAODZDC* fAODZDC;
  AliAODAD* fAODAD;
  AliAODTZERO* fAODTZERO;
  
  TClonesArray* fPairs;  // J/psi candidates - AliDielectronPair
  TClonesArray* fDaughters; // Daughters of e+e- pairs
  TClonesArray* fJets; // Jet found with R=0.4 - AliEmcalJet

private:
  void FillJets(AliAODEvent* aodEv, TClonesArray* jetArray, TString jetName);
  Bool_t FindDaughters(AliVTrack* trk);
  void SetTrackFromPair(AliDielectronPair* pair, AliAODTrack* tmp);

/*Copy from AliAnalysisTaskDielectronFilter*/
public:
  enum ETriggerLogig {kAny, kExact};
  virtual void UserCreateOutputObjects();
  virtual void LocalInit() {Init();}
  void NotifyRun(){AliDielectronPID::SetCorrVal((Double_t)fCurrentRunNumber);}

  void UsePhysicsSelection(Bool_t phy=kTRUE) {fSelectPhysics=phy;}
  void SetTriggerMask(UInt_t mask) {fTriggerMask=mask;}
  UInt_t GetTriggerMask() const { return fTriggerMask; }
  void SetExcludeTriggerMask(ULong64_t mask) {fExcludeTriggerMask=mask;}
  UInt_t GetExcludeTriggerMask() const { return fExcludeTriggerMask; }
  void SetTriggerLogic(ETriggerLogig log) {fTriggerLogic=log;}
  ETriggerLogig GetTriggerLogic() const {return fTriggerLogic;}


  void SetDielectron(AliDielectron * const die) { fDielectron = die; }

  void SetStoreLikeSignCandidates(Bool_t storeLS) { fStoreLikeSign = storeLS; }
  void SetStoreRotatedPairs(Bool_t storeTR) { fStoreRotatedPairs = storeTR; }
  void SetStoreEventsWithSingleTracks(Bool_t storeSingleTrk) { fStoreEventsWithSingleTracks = storeSingleTrk; }
  void SetCreateNanoAODs(Bool_t storeTrackRef) { fCreateNanoAOD = storeTrackRef; }

  void SetEventFilter(AliAnalysisCuts * const filter) {fEventFilter=filter;}

private:
  enum {kAllEvents=0, kSelectedEvents, kV0andEvents, kFilteredEvents, kPileupEvents, kNbinsEvent};

  AliDielectron *fDielectron;             // J/psi framework object

  Bool_t fSelectPhysics;                  // Whether to use physics selection
  UInt_t fTriggerMask;               // Event trigger mask
  UInt_t fExcludeTriggerMask;        // Triggers to exclude from the analysis
  Bool_t fTriggerOnV0AND;            // if to trigger on V0and
  Bool_t fRejectPileup;              // pileup rejection wanted

  TH1D *fEventStat;                  //! Histogram with event statistics

  ETriggerLogig fTriggerLogic;       // trigger logic: any or all bits need to be matching

  AliTriggerAnalysis *fTriggerAnalysis; //! trigger analysis class

  Bool_t fStoreLikeSign;        // flag to store like-sign candidates
  Bool_t fStoreRotatedPairs;    // flag to store rotation
  Bool_t fStoreEventsWithSingleTracks;    // flag to store events with a least one reconstructed track
  Bool_t fCreateNanoAOD;        // flag to create nanoAODs

  AliAnalysisCuts *fEventFilter;     // event filter

  ClassDef(AliAnalysisTaskJpsiJetFilter, 2);

};
#endif // #ifndef ALIANALYSISTASK_JPSIJET_FILTER_H
