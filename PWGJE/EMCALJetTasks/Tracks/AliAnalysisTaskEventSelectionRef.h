#ifndef ALIANALYSISTASKEVENTSELECTIONREF_H
#define ALIANALYSISTASKEVENTSELECTIONREF_H
/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliAnalysisTaskSE.h"

class TClonesArray;

class AliAnalysisUtils;

namespace EMCalTriggerPtAnalysis {

class AliEMCalHistoContainer;

class AliAnalysisTaskEventSelectionRef : public AliAnalysisTaskSE {
public:
  enum EmcalTriggerClass{
    kCPREL0 = 0,
    kCPREG1,
    kCPREG2,
    kCPREJ1,
    kCPREJ2,
    kCPRntrig
  };
  AliAnalysisTaskEventSelectionRef();
  AliAnalysisTaskEventSelectionRef(const char *name);
  virtual ~AliAnalysisTaskEventSelectionRef();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *);

  void SetOfflineEnergyThreshold(EmcalTriggerClass trgcls, double threshold) { fOfflineEnergyThreshold[trgcls] = threshold; }

protected:
  void FillEventCounterHists(const char *triggerclass, double vtxz, bool isSelected, bool isOfflineSelected);
  Bool_t IsOfflineSelected(EmcalTriggerClass trgcls, const TClonesArray * const triggerpatches) const;

  AliAnalysisUtils              *fAnalysisUtils;
  AliEMCalHistoContainer        *fHistos;
  Double_t                      fOfflineEnergyThreshold[kCPRntrig];

  ClassDef(AliAnalysisTaskEventSelectionRef, 1);
};

} /* namespace EMCalTriggerPtAnalysis */

#endif
