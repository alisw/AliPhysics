#ifndef ALIEMCALTRIGGERSELECTIONCUTS_H
#define ALIEMCALTRIGGERSELECTIONCUTS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Author: Markus Fasel

#include <TObject.h>

class AliEMCALTriggerPatchInfo;

class AliEmcalTriggerSelectionCuts: public TObject {
public:
  enum SelectionMethod_t {
    kADC = 0,
    kEnergyRough = 1,
    kEnergyOffline = 2
  };
  enum PatchType_t {
    kAnyPatch = 0,
    kL1JetPatch = 1,
    kL1GammaPatch = 2,
    kL0Patch = 3,
    kL1JetLowPatch = 4,
    kL1JetHighPatch = 5,
    kL1GammaLowPatch = 6,
    kL1GammaHighPatch = 7
  };

  AliEmcalTriggerSelectionCuts();
  virtual ~AliEmcalTriggerSelectionCuts() {}

  PatchType_t GetPatchType() const { return fPatchType; }
  SelectionMethod_t GetSelectionMethod() const { return fSelectionMethod; }
  Double_t GetThreshold() const { return fThreshold; }
  Bool_t IsRequestingSimpleOfflinePatches() const { return fUseSimpleOffline; }

  void SetPatchType(PatchType_t patchType) { fPatchType = patchType; }
  void SetSelectionMethod(SelectionMethod_t selectionMethod) { fSelectionMethod = selectionMethod; }
  void SetThreshold(Double_t threshold) { fThreshold = threshold; }
  void SetUseSimpleOfflinePatches(Bool_t doUse = kTRUE) { fUseSimpleOffline = doUse; }

  Bool_t IsSelected(const AliEMCALTriggerPatchInfo * const patch) const;
  Int_t CompareTriggerPatches(const AliEMCALTriggerPatchInfo *first, const AliEMCALTriggerPatchInfo *second) const;

protected:
  Double_t GetCutPrimitive(const AliEMCALTriggerPatchInfo * const patch) const;
  Bool_t SelectPatchType(const AliEMCALTriggerPatchInfo * const patch) const;

  SelectionMethod_t     fSelectionMethod;           // Variable to cut on
  PatchType_t           fPatchType;                 // Type of the patch to be selected
  Double_t              fThreshold;                 // Threshold used
  Bool_t                fUseSimpleOffline;          // Request simple offline patches

  ClassDef(AliEmcalTriggerSelectionCuts, 1);         // Cuts for the EMCAL Trigger selection
};

#endif /* ALIEMCALTRIGGERSELECTIONCUTS_H */
