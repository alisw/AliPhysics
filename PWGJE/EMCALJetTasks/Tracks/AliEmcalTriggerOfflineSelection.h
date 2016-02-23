/*
 * AliEmcalTriggerOfflineSelection.h
 *
 *  Created on: Feb 23, 2016
 *      Author: markus
 */

#ifndef ALIEMCALTRIGGEROFFLINESELECTION_H
#define ALIEMCALTRIGGEROFFLINESELECTION_H

#include <TObject.h>

class TClonesArray;

namespace EMCalTriggerPtAnalysis {

class AliEmcalTriggerOfflineSelection {
public:
  enum EmcalTriggerClass{
    kTrgEL0 = 0,
    kTrgEG1,
    kTrgEG2,
    kTrgEJ1,
    kTrgEJ2,
    kTrgDL0,
    kTrgDG1,
    kTrgDG2,
    kTrgDJ1,
    kTrgDJ2,
    kTrgn
  };

  AliEmcalTriggerOfflineSelection();
  virtual ~AliEmcalTriggerOfflineSelection() {}

  void                        SetOfflineEnergyThreshold(EmcalTriggerClass trgcls, double threshold) { fOfflineEnergyThreshold[trgcls] = threshold; }
  Bool_t                      IsOfflineSelected(EmcalTriggerClass trgcls, const TClonesArray * const triggerpatches) const;
  Double_t                    GetThresholdForTrigger(EmcalTriggerClass trgcls) const {return fOfflineEnergyThreshold[trgcls]; }

  static Bool_t IsSingleShower(EmcalTriggerClass cls);
  static Bool_t IsDCAL(EmcalTriggerClass cls);

protected:
  Double_t                    fOfflineEnergyThreshold[kTrgn];

  ClassDef(AliEmcalTriggerOfflineSelection, 1);
};

} /* namespace EMCalTriggerPtAnalysis */

#endif /* ALIEMCALTRIGGEROFFLINESELECTION_H */
