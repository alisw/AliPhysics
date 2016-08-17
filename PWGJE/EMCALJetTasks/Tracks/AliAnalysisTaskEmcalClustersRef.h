#ifndef ALIANALYSISTASKEMCALCLUSTERSREF_H
#define ALIANALYSISTASKEMCALCLUSTERSREF_H
/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliAnalysisTaskSE.h"
#include "AliCutValueRange.h"
#include <TCustomBinning.h>
#include <TString.h>

class AliAnalysisUtils;
class AliEMCALGeometry;
class AliOADBContainer;

class TArrayD;
class TClonesArray;
class THistManager;
class TList;
class TObjArray;
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
  void SetDownscaleOADB(TString oadbname) { fNameDownscaleOADB = oadbname; }

protected:

  virtual void ExecOnce();
  virtual void RunChanged(Int_t runnumber);

  Double_t GetTriggerWeight(const TString &triggerclass) const;
  void GetPatchBoundaries(TObject *o, Double_t *boundaries) const;
  bool IsOfflineSimplePatch(TObject *o) const;
  bool SelectDCALPatch(TObject *o) const;
  bool SelectSingleShowerPatch(TObject *o) const;
  bool SelectJetPatch(TObject *o) const;
  double GetPatchEnergy(TObject *o) const;

  void FillClusterHistograms(const TString &triggerclass, double energy, double transversenergy, double eta, double phi, TList *triggerpatches);
  void FillEventHistograms(const TString &triggerclass, double centrality, double vertexz);
  TString GetFiredTriggerClassesFromPatches(const TClonesArray* triggerpatches) const;
  void FindPatchesForTrigger(TString triggerclass, const TClonesArray * triggerpatches, TList &foundpatches) const;
  Bool_t CorrelateToTrigger(Double_t etaclust, Double_t phiclust, TList *triggerpatches) const;

  AliAnalysisUtils                    *fAnalysisUtil;             ///< Analysis utils for additional event selection / pileup rejection
  THistManager                        *fHistos;                   //!<! Histogram handler
  AliEmcalTriggerOfflineSelection     *fTriggerSelection;         ///< EMCAL offline trigger selection tool
  AliEMCALGeometry                    *fGeometry;                 //!<! EMCAL geometry
  TClonesArray                        *fTriggerPatches;           //!<! Container with trigger patches
  TString                             fClusterContainer;          ///< Name of the cluster container in the event

  Bool_t                              fRequestAnalysisUtil;       ///< Switch on request for event selection using analysis utils
  Bool_t                              fTriggerStringFromPatches;  ///< Build trigger string from trigger patches
  AliCutValueRange<double>            fCentralityRange;           ///< Selected centrality range
  AliCutValueRange<double>            fVertexRange;               ///< Selected vertex range
  Bool_t                              fRequestCentrality;         ///< Swich on request for centrality range

  TString                             fNameDownscaleOADB;         ///< Name of the downscale OADB container
  AliOADBContainer                    *fDownscaleOADB;            //!<! Container with downscale factors for different triggers
  TObjArray                           *fDownscaleFactors;         //!<! Downscalfactors for given run

  Int_t                               fCurrentRun;                ///< Current run number (for RunChange method)
  Bool_t                              fInitialized;               ///< Check for initialized

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
