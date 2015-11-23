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
  enum EmcalTriggerClass{
    kEPREL0 = 0,
    kEPREG1,
    kEPREG2,
    kEPREJ1,
    kEPREJ2,
    kEPRDL0,
    kEPRDG1,
    kEPRDG2,
    kEPRDJ1,
    kEPRDJ2,
    kEPRntrig
  };
  AliAnalysisTaskEmcalPatchesRef();
  AliAnalysisTaskEmcalPatchesRef(const char *name);
  virtual ~AliAnalysisTaskEmcalPatchesRef();

  void UserCreateOutputObjects();
  void UserExec(Option_t *);
  void Terminate(Option_t *) {}

  void SetCreateTriggerStringFromPatches(Bool_t doUsePatches) { fTriggerStringFromPatches = doUsePatches; }
  void SetOfflineEnergyThreshold(EmcalTriggerClass trgcls, double threshold) { fOfflineEnergyThreshold[trgcls] = threshold; }

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
  Bool_t IsOfflineSelected(EmcalTriggerClass trgcls, const TClonesArray * const triggerpatches) const;

  AliAnalysisUtils                    *fAnalysisUtil;
  AliEMCalHistoContainer              *fHistos;

  Bool_t                              fRequestAnalysisUtil;
  Bool_t                              fTriggerStringFromPatches;
  Double_t                            fOfflineEnergyThreshold[kEPRntrig];

private:
  AliAnalysisTaskEmcalPatchesRef(const AliAnalysisTaskEmcalPatchesRef &);
  AliAnalysisTaskEmcalPatchesRef &operator=(const AliAnalysisTaskEmcalPatchesRef &);

  ClassDef(AliAnalysisTaskEmcalPatchesRef, 1);

};

} /* namespace EMCalTriggerPtAnalysis */

#endif /* ALIANALYSISTASKEMCALPATCHESREF_H */
