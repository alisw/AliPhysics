#ifndef ALIANALYSISTASKDIELECTRONFILTER_H
#define ALIANALYSISTASKDIELECTRONFILTER_H
/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//#####################################################
//#                                                   # 
//#        Dielectron even filter task                #
//#                                                   #
//#                                                   #
//#  by WooJin J. Park, GSI / W.J.Park@gsi.de         #
//#     Ionut C. Arsene, GSI / I.C.Arsene@gsi.de      #
//#     Magnus Mager, CERN / Magnus.Mager@cern.ch     #
//#     Jens Wiechula, Uni HD / Jens.Wiechula@cern.ch #
//#                                                   #
//#####################################################
/*
Filter Event based on cuts provided in the AliDielectron class.

Write an AOD file containing events with Dielectron candidates.
Add a sattelite AOD with the array of candidates.
*/



#include "AliAnalysisTaskSE.h"

#include "AliDielectronPID.h"

class AliDielectron;
class TH1D;
class AliTriggerAnalysis;
class AliAODCaloCluster;

class AliAnalysisTaskDielectronFilter : public AliAnalysisTaskSE {
  
public:
  AliAnalysisTaskDielectronFilter();
  AliAnalysisTaskDielectronFilter(const char *name);
  virtual ~AliAnalysisTaskDielectronFilter(){}

  enum ETriggerLogig {kAny, kExact};

  virtual void UserExec(Option_t *option);
  virtual void Init();
  virtual void UserCreateOutputObjects();
  virtual void LocalInit() {Init();}
  //temporary
  virtual void NotifyRun(){AliDielectronPID::SetCorrVal((Double_t)fCurrentRunNumber);}
  
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
  void SetStoreHeader(Bool_t storeHeader) { fStoreHeader = storeHeader; }

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
  Bool_t fStoreHeader;          // flag to store header for all events 

  AliAnalysisCuts *fEventFilter;     // event filter
  
  AliAnalysisTaskDielectronFilter(const AliAnalysisTaskDielectronFilter &c);
  AliAnalysisTaskDielectronFilter& operator= (const AliAnalysisTaskDielectronFilter &c);
  
  ClassDef(AliAnalysisTaskDielectronFilter, 1);
};
#endif
