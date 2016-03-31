#ifndef ALIANALYSISTASKEMCALPATCHESREF_H
#define ALIANALYSISTASKEMCALPATCHESREF_H
/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliAnalysisTaskSE.h"
#include "AliCutValueRange.h"

class AliAnalysisUtils;
class THistManager;
class TArrayD;
class TString;

namespace EMCalTriggerPtAnalysis {

class AliEmcalTriggerOfflineSelection;

class AliAnalysisTaskEmcalPatchesRef : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskEmcalPatchesRef();
  AliAnalysisTaskEmcalPatchesRef(const char *name);
  virtual ~AliAnalysisTaskEmcalPatchesRef();

  void UserCreateOutputObjects();
  void UserExec(Option_t *);
  void Terminate(Option_t *) {}

  void SetOfflineTriggerSelection(AliEmcalTriggerOfflineSelection *sel) { fTriggerSelection = sel; }
  void SetCreateTriggerStringFromPatches(Bool_t doUsePatches) { fTriggerStringFromPatches = doUsePatches; }
  void SetRequestAnalysisUtil(bool doUse) { fRequestAnalysisUtil = doUse; }
  void SetCentralityRange(double min, double max) { fCentralityRange.SetLimits(min,max); }
  void SetVertexRange(double min, double max) { fVertexRange.SetLimits(min, max); }

protected:

  void GetPatchBoundaries(TObject *o, Double_t *boundaries) const;
  bool IsOfflineSimplePatch(TObject *o) const;
  bool SelectDCALPatch(TObject *o) const;
  bool SelectSingleShowerPatch(TObject *o) const;
  bool SelectJetPatch(TObject *o) const;
  double GetPatchEnergy(TObject *o) const;

  void CreateEnergyBinning(TArrayD& binning) const;
  void CreateLinearBinning(TArrayD& binning, int nbins, double min, double max) const;
  void FillPatchHistograms(TString triggerclass, TString patchname, double energy, double transverseenergy, double eta, double phi);
  TString GetFiredTriggerClassesFromPatches(const TClonesArray* triggerpatches) const;

  AliAnalysisUtils                    *fAnalysisUtil;
  AliEmcalTriggerOfflineSelection     *fTriggerSelection;
  THistManager                        *fHistos;

  Bool_t                              fRequestAnalysisUtil;
  Bool_t                              fTriggerStringFromPatches;

  AliCutValueRange<double>            fCentralityRange;
  AliCutValueRange<double>            fVertexRange;

private:
  AliAnalysisTaskEmcalPatchesRef(const AliAnalysisTaskEmcalPatchesRef &);
  AliAnalysisTaskEmcalPatchesRef &operator=(const AliAnalysisTaskEmcalPatchesRef &);

  ClassDef(AliAnalysisTaskEmcalPatchesRef, 1);

};

} /* namespace EMCalTriggerPtAnalysis */

#endif /* ALIANALYSISTASKEMCALPATCHESREF_H */
