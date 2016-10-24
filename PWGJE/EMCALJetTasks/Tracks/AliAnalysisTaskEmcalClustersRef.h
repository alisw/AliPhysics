#ifndef ALIANALYSISTASKEMCALCLUSTERSREF_H
#define ALIANALYSISTASKEMCALCLUSTERSREF_H
/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliAnalysisTaskEmcalTriggerBase.h"
#include "AliCutValueRange.h"
#include <TCustomBinning.h>
#include <TString.h>

class TClonesArray;

namespace EMCalTriggerPtAnalysis {

class AliAnalysisTaskEmcalClustersRef : public AliAnalysisTaskEmcalTriggerBase {
public:
  AliAnalysisTaskEmcalClustersRef();
  AliAnalysisTaskEmcalClustersRef(const char *name);
  virtual ~AliAnalysisTaskEmcalClustersRef();

  void SetClusterContainer(TString clustercontname) { fNameClusterContainer = clustercontname; }
  void SetCreateTriggerStringFromPatches(Bool_t doUsePatches) { fTriggerStringFromPatches = doUsePatches; }

  void SetCentralityRange(double min, double max) { fCentralityRange.SetLimits(min, max); fRequestCentrality = true; }

protected:
  virtual void CreateUserObjects() {}
  virtual void CreateUserHistos();
  virtual bool IsUserEventSelected();
  virtual bool Run();
  virtual void UserFillHistosAfterEventSelection();

  void GetPatchBoundaries(TObject *o, Double_t *boundaries) const;
  bool IsOfflineSimplePatch(TObject *o) const;
  bool SelectDCALPatch(TObject *o) const;
  bool SelectSingleShowerPatch(TObject *o) const;
  bool SelectJetPatch(TObject *o) const;
  double GetPatchEnergy(TObject *o) const;

  void FillClusterHistograms(const TString &triggerclass, double energy, double transversenergy, double eta, double phi, TList *triggerpatches);
  void FindPatchesForTrigger(TString triggerclass, const TClonesArray * triggerpatches, TList &foundpatches) const;
  Bool_t CorrelateToTrigger(Double_t etaclust, Double_t phiclust, TList *triggerpatches) const;

  TString                             fNameClusterContainer;      ///< Name of the cluster container in the event

  AliCutValueRange<double>            fCentralityRange;           ///< Selected centrality range
  Bool_t                              fRequestCentrality;         ///< Switch on request for centrality range
  Double_t                            fEventCentrality;           //!<! Current event centrality

private:

  class EnergyBinning : public TCustomBinning {
  public:
    EnergyBinning();
    virtual ~EnergyBinning() {}
  };

  AliAnalysisTaskEmcalClustersRef(const AliAnalysisTaskEmcalClustersRef &);
  AliAnalysisTaskEmcalClustersRef &operator=(const AliAnalysisTaskEmcalClustersRef &);

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskEmcalClustersRef, 1);
  /// \endcond
};

} /* namespace EMCalTriggerPtAnalysis */

#endif /* ALIANALYSISTASKEMCALCLUSTERSREF_H */
