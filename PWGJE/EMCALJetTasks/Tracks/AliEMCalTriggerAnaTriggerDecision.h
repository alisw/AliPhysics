#ifndef ALIEMCALTRIGGERANATRIGGERDECISION_H
#define ALIEMCALTRIGGERANATRIGGERDECISION_H
/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Author: Markus Fasel
#include <TObject.h>

class TClonesArray;
class TString;
class AliEmcalTriggerPatchInfo;

namespace EMCalTriggerPtAnalysis {

class AliEMCalTriggerEventData;

class AliEMCalTriggerAnaTriggerDecision : public TObject {
public:
  enum ETATriggerType{
    kTAEMCJHigh       = 0,
    kTAEMCJLow        = 1,
    kTAEMCGHigh       = 2,
    kTAEMCGLow        = 3
  };
  AliEMCalTriggerAnaTriggerDecision();
  virtual ~AliEMCalTriggerAnaTriggerDecision(){}

  void Create(const AliEMCalTriggerEventData * const data);
  Bool_t IsTriggered(ETATriggerType trigger, Bool_t fromPatches = kFALSE) const {
    if(fromPatches) return fDecisionFromPatches[trigger];
    return fDecisionFromString[trigger];
  }
  void SetSwapThresholds(Bool_t doSwap = kTRUE) { fSwapThresholds = doSwap; }
  void SetIsMinBias(Bool_t isMB = kTRUE) { fIsMinBias = isMB; }
  void SetUseOfflinePatches(bool useOffline = kTRUE) { fUseOfflinePatches = useOffline; }
  void SetOfflineEnergyThreshold(ETATriggerType triggerClass, double threshold){
	  fEnergyThresholds[triggerClass] = threshold;
  }
  void Reset();

  Bool_t IsMinBias() const { return fIsMinBias; }
  void Print(Option_t * opt = NULL) const;

  void SetDebugMode(Bool_t doDebug = true) { fDoDebug = doDebug; }

  bool CheckConsistency() const;

protected:
  void MakeDecisionFromString(const TString &triggerstring);
  void MakeDecisionFromPatches(const TClonesArray &listOfPatches);

  Bool_t SelectTriggerPatch(ETATriggerType trigger, const AliEmcalTriggerPatchInfo * const recpatch) const;

  Bool_t fSwapThresholds;
  Bool_t fIsMinBias;
  Bool_t fUseOfflinePatches;
  Bool_t fDecisionFromPatches[4];
  Bool_t fDecisionFromString[4];

  Double_t fEnergyThresholds[4];

  Bool_t fDoDebug;

  ClassDef(AliEMCalTriggerAnaTriggerDecision, 1);     // EMCal trigger decision
};

} /* namespace EMCalTriggerPtAnalysis */

#endif /* ALIEMCALTRIGGERANATRIGGERDECISION_H */
