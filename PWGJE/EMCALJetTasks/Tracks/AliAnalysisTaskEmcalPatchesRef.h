#ifndef ALIANALYSISTASKEMCALPATCHESREF_H
#define ALIANALYSISTASKEMCALPATCHESREF_H
/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliAnalysisTaskSE.h"

class AliAnalysisUtils;
class TArrayD;
class TString;

namespace EMCalTriggerPtAnalysis {

class AliAnalysisTaskEmcalPatchesRef : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskEmcalPatchesRef();
  AliAnalysisTaskEmcalPatchesRef(const char *name);
  virtual ~AliAnalysisTaskEmcalPatchesRef();

  void UserCreateOutputObjects();
  void UserExec(Option_t *);
  void Terminate(Option_t *) {}

  void SetCreateTriggerStringFromPatches(Bool_t doUsePatches) { fTriggerStringFromPatches = doUsePatches; }

protected:
  void CreateEnergyBinning(TArrayD& binning) const;
  void FillPatchHistograms(TString triggerclass, TString patchname, double energy, double eta, double phi);
  TString GetFiredTriggerClassesFromPatches(const TClonesArray* triggerpatches) const;

  AliAnalysisUtils                    *fAnalysisUtil;
  AliEMCalHistoContainer              *fHistos;

  Bool_t                              fTriggerStringFromPatches;

private:
  AliAnalysisTaskEmcalPatchesRef(const AliAnalysisTaskEmcalPatchesRef &);
  AliAnalysisTaskEmcalPatchesRef &operator=(const AliAnalysisTaskEmcalPatchesRef &);

  ClassDef(AliAnalysisTaskEmcalPatchesRef, 1);

};

} /* namespace EMCalTriggerPtAnalysis */

#endif /* ALIANALYSISTASKEMCALPATCHESREF_H */
