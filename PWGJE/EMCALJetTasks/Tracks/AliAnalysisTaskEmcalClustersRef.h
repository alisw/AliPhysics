#ifndef ALIANALYSISTASKEMCALCLUSTERSREF_H
#define ALIANALYSISTASKEMCALCLUSTERSREF_H
/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliAnalysisTaskSE.h"
#include "AliCutValueRange.h"
#include <TCustomBinning.h>

class AliAnalysisUtils;
class AliEMCALGeometry;

class TArrayD;
class TClonesArray;
class THistManager;
class TList;
class TString;

namespace EMCalTriggerPtAnalysis {

class AliEmcalTriggerOfflineSelection;

class AliAnalysisTaskEmcalClustersRef : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskEmcalClustersRef();
  AliAnalysisTaskEmcalClustersRef(const char *name);
  virtual ~AliAnalysisTaskEmcalClustersRef();

  void UserCreateOutputObjects();
  void UserExec(Option_t *);
  void Terminate(Option_t *) {}

  void SetOfflineTriggerSelection(AliEmcalTriggerOfflineSelection *sel) { fTriggerSelection = sel; }
  void SetClusterContainer(TString clustercontname) { fClusterContainer = clustercontname; }
  void SetCreateTriggerStringFromPatches(Bool_t doUsePatches) { fTriggerStringFromPatches = doUsePatches; }

  void SetRequestAnalysisUtil(Bool_t doRequest) { fRequestAnalysisUtil = doRequest; }
  void SetCentralityRange(double min, double max) { fCentralityRange.SetLimits(min, max); fRequestCentrality = true; }
  void SetVertexRange(double min, double max) { fVertexRange.SetLimits(min, max); }

protected:

  void GetPatchBoundaries(TObject *o, Double_t *boundaries) const;
  bool IsOfflineSimplePatch(TObject *o) const;
  bool SelectDCALPatch(TObject *o) const;
  bool SelectSingleShowerPatch(TObject *o) const;
  bool SelectJetPatch(TObject *o) const;
  double GetPatchEnergy(TObject *o) const;

  void FillClusterHistograms(TString triggerclass, double energy, double transversenergy, double eta, double phi, TList *triggerpatches);
  TString GetFiredTriggerClassesFromPatches(const TClonesArray* triggerpatches) const;
  void FindPatchesForTrigger(TString triggerclass, const TClonesArray * triggerpatches, TList &foundpatches) const;
  Bool_t CorrelateToTrigger(Double_t etaclust, Double_t phiclust, TList *triggerpatches) const;

  AliAnalysisUtils                    *fAnalysisUtil;
  THistManager                        *fHistos;
  AliEmcalTriggerOfflineSelection     *fTriggerSelection;
  AliEMCALGeometry                    *fGeometry;
  TString                             fClusterContainer;

  Bool_t                              fRequestAnalysisUtil;
  Bool_t                              fTriggerStringFromPatches;
  AliCutValueRange<double>            fCentralityRange;
  AliCutValueRange<double>            fVertexRange;
  Bool_t                              fRequestCentrality;

private:

  class EnergyBinning : public TCustomBinning {
  public:
    EnergyBinning();
    virtual ~EnergyBinning() {}
  };

  AliAnalysisTaskEmcalClustersRef(const AliAnalysisTaskEmcalClustersRef &);
  AliAnalysisTaskEmcalClustersRef &operator=(const AliAnalysisTaskEmcalClustersRef &);

  ClassDef(AliAnalysisTaskEmcalClustersRef, 1);
};

} /* namespace EMCalTriggerPtAnalysis */

#endif /* ALIANALYSISTASKEMCALCLUSTERSREF_H */
