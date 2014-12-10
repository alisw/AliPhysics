#ifndef ALIEMCALTRIGGERANATRIGGERDECISION_H
#define ALIEMCALTRIGGERANATRIGGERDECISION_H
/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Author: Markus Fasel
#include <TObject.h>

class TClonesArray;
class TString;

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
  void Reset();

  Bool_t IsMinBias() const { return fIsMinBias; }

protected:
  void MakeDecisionFromString(const TString &triggerstring);
  void MakeDecisionFromPatches(const TClonesArray &listOfPatches);

  Bool_t fSwapThresholds;
  Bool_t fIsMinBias;
  Bool_t fDecisionFromPatches[4];
  Bool_t fDecisionFromString[4];

  ClassDef(AliEMCalTriggerAnaTriggerDecision, 1);     // EMCal trigger decision
};

} /* namespace EMCalTriggerPtAnalysis */

#endif /* ALIEMCALTRIGGERANATRIGGERDECISION_H */
