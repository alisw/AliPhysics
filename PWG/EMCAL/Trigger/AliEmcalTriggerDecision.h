#ifndef ALIEMCALTRIGGERDECISION_H
#define ALIEMCALTRIGGERDECISION_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Author: Markus Fasel

#include <TList.h>
#include <TNamed.h>

class AliEmcalTriggerPatchInfo;
class AliEmcalTriggerSelectionCuts;

class AliEmcalTriggerDecision: public TNamed {
public:
  AliEmcalTriggerDecision();
  AliEmcalTriggerDecision(const char *name, const char *title = "");
  virtual ~AliEmcalTriggerDecision() {}

  const AliEmcalTriggerPatchInfo *GetMainPatch() const { return fMainPatch; }
  const AliEmcalTriggerSelectionCuts *GetSelectionCuts() const { return fSelectionCuts; }
  const TList *GetAcceptedPatches() const { return &fAcceptedPatches; }
  Bool_t IsSelected() const { return fMainPatch != NULL; }

  void SetSelectionCuts(const AliEmcalTriggerSelectionCuts * const cuts) { fSelectionCuts = cuts; }
  void SetMainPatch(const AliEmcalTriggerPatchInfo * const mainpatch) { fMainPatch = mainpatch; }
  void AddAcceptedPatch(AliEmcalTriggerPatchInfo * const acceptedPatch);

protected:
  const AliEmcalTriggerPatchInfo          *fMainPatch;         // Main trigger patch which fires the decision
  const AliEmcalTriggerSelectionCuts      *fSelectionCuts;     // Pointer to the cuts used for the trigger selection
  TList                                    fAcceptedPatches;   // All trigger patches which are accepted as well

  ClassDef(AliEmcalTriggerDecision, 1);               // Container of the trigger decision information
private:
  AliEmcalTriggerDecision(const AliEmcalTriggerDecision &ref);
  AliEmcalTriggerDecision &operator=(const AliEmcalTriggerDecision &ref);
};

#endif /* ALIEMCALTRIGGERDECISION_H */
