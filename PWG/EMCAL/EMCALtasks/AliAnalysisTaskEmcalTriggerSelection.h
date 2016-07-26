#ifndef ALIANALYSISTASKEMCALTRIGGERSELECTION_H
#define ALIANALYSISTASKEMCALTRIGGERSELECTION_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Author: Markus Fasel

#include <TList.h>
#include <TString.h>
#include "AliAnalysisTaskEmcal.h"

class AliEmcalTriggerDecisionContainer;
class AliEmcalTriggerSelection;

class AliAnalysisTaskEmcalTriggerSelection: public AliAnalysisTaskEmcal {
public:
  AliAnalysisTaskEmcalTriggerSelection();
  AliAnalysisTaskEmcalTriggerSelection(const char *name);
  virtual ~AliAnalysisTaskEmcalTriggerSelection() {}

  void AddTriggerSelection(AliEmcalTriggerSelection * const selection);
  void SetGlobalDecisionContainerName(const char *name) { fGlobalDecisionContainerName = name; }

  virtual Bool_t Run();

protected:
  AliEmcalTriggerDecisionContainer *GetGlobalTriggerDecisionContainer();

  TString fGlobalDecisionContainerName;     // Name of the global trigger selection
  TList fTriggerSelections;                 // List of trigger selections

  ClassDef(AliAnalysisTaskEmcalTriggerSelection, 1);    // Task running different EMCAL trigger selections
};

#endif /* ALIANALYSISTASKEMCALTRIGGERSELECTION_H */
