#ifndef ALIEMCALTRIGGERSELECTION_H
#define ALIEMCALTRIGGERSELECTION_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Author: Markus Fasel

#include <TNamed.h>
#include <TString.h>

class AliEmcalTriggerDecision;
class AliEMCALTriggerPatchInfo;
class AliEmcalTriggerSelectionCuts;
class TClonesArray;

class AliEmcalTriggerSelection: public TNamed {
public:
  AliEmcalTriggerSelection();
  AliEmcalTriggerSelection(const char *name, const AliEmcalTriggerSelectionCuts * const cuts);
  virtual ~AliEmcalTriggerSelection() {}

  const AliEmcalTriggerSelectionCuts *GetSelectionCuts() const { return fSelectionCuts; }

  void SetOutputName(const char *name) { fOutputName = name; }
  void SetSelectionCuts(const AliEmcalTriggerSelectionCuts * const cuts) { fSelectionCuts = cuts; }

  AliEmcalTriggerDecision * MakeDecison(const TClonesArray * const reconstructedPatches) const;
protected:
  const AliEmcalTriggerSelectionCuts  *fSelectionCuts;    // Cuts used for the trigger patch selection
  TString                              fOutputName;       // Name of the output object (AliEmcalTriggerDecision)

  ClassDef(AliEmcalTriggerSelection, 1);    // EMCAL trigger selection component
private:
  AliEmcalTriggerSelection(const AliEmcalTriggerSelection &ref);
  AliEmcalTriggerSelection &operator=(const AliEmcalTriggerSelection &ref);
};

#endif /* ALIEMCALTRIGGERSELECTION_H */
