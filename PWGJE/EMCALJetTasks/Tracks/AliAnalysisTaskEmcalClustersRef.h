#ifndef ALIANALYSISTASKEMCALCLUSTERSREF_H
#define ALIANALYSISTASKEMCALCLUSTERSREF_H
/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliAnalysisTaskSE.h"

class AliAnalysisUtils;
class AliEMCALGeometry;

class TArrayD;
class TClonesArray;
class TList;
class TString;

namespace EMCalTriggerPtAnalysis {

class AliEMCalHistoContainer;

class AliAnalysisTaskEmcalClustersRef : public AliAnalysisTaskSE {
public:
  enum EmcalTriggerClass{
    kECREL0 = 0,
    kECREG1,
    kECREG2,
    kECREJ1,
    kECREJ2,
    kECRntrig
  };
  AliAnalysisTaskEmcalClustersRef();
  AliAnalysisTaskEmcalClustersRef(const char *name);
  virtual ~AliAnalysisTaskEmcalClustersRef();

  void UserCreateOutputObjects();
  void UserExec(Option_t *);
  void Terminate(Option_t *) {}

  void SetClusterContainer(TString clustercontname) { fClusterContainer = clustercontname; }
  void SetCreateTriggerStringFromPatches(Bool_t doUsePatches) { fTriggerStringFromPatches = doUsePatches; }

  void SetOfflineEnergyThreshold(EmcalTriggerClass trgcls, double threshold) { fOfflineEnergyThreshold[trgcls] = threshold; }

protected:

  void CreateEnergyBinning(TArrayD& binning) const;
  void CreateLinearBinning(TArrayD& binning, int nbins, double min, double max) const;
  void FillClusterHistograms(TString triggerclass, double energy, double eta, double phi, TList *triggerpatches);
  TString GetFiredTriggerClassesFromPatches(const TClonesArray* triggerpatches) const;
  void FindPatchesForTrigger(TString triggerclass, const TClonesArray * triggerpatches, TList &foundpatches) const;
  Bool_t CorrelateToTrigger(Double_t etaclust, Double_t phiclust, TList *triggerpatches) const;
  Bool_t IsOfflineSelected(EmcalTriggerClass trgcls, const TClonesArray * const triggerpatches) const;

  AliAnalysisUtils                    *fAnalysisUtil;
  AliEMCalHistoContainer              *fHistos;
  AliEMCALGeometry                    *fGeometry;
  TString                             fClusterContainer;

  Bool_t                              fTriggerStringFromPatches;
  Double_t                            fOfflineEnergyThreshold[kECRntrig];

private:
  AliAnalysisTaskEmcalClustersRef(const AliAnalysisTaskEmcalClustersRef &);
  AliAnalysisTaskEmcalClustersRef &operator=(const AliAnalysisTaskEmcalClustersRef &);

  ClassDef(AliAnalysisTaskEmcalClustersRef, 1);
};

} /* namespace EMCalTriggerPtAnalysis */

#endif /* ALIANALYSISTASKEMCALCLUSTERSREF_H */
