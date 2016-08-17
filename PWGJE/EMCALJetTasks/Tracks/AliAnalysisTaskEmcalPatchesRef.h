#ifndef ALIANALYSISTASKEMCALPATCHESREF_H
#define ALIANALYSISTASKEMCALPATCHESREF_H
/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliAnalysisTaskSE.h"
#include "AliCutValueRange.h"
#include <TCustomBinning.h>
#include <TString.h>

class AliAnalysisUtils;
class AliOADBContainer;
class THistManager;
class TObjArray;

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
  void SetCentralityRange(double min, double max) { fCentralityRange.SetLimits(min,max); fRequestCentrality = true; }
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

  void FillPatchHistograms(TString triggerclass, TString patchname, double energy, double transverseenergy, double eta, double phi, int col, int row);
  void FillEventHistograms(const TString &triggerclass, double centrality, double vertexz);
  TString GetFiredTriggerClassesFromPatches(const TClonesArray* triggerpatches) const;

  AliAnalysisUtils                    *fAnalysisUtil;             ///< Analysis utils for event selection
  AliEmcalTriggerOfflineSelection     *fTriggerSelection;         ///< Offline trigger selection tool
  THistManager                        *fHistos;                   //!<! Histogram handler
  TClonesArray                        *fTriggerPatches;           //!<! Container with trigger patches

  Bool_t                              fRequestAnalysisUtil;       ///< Switch on request for analysis util
  Bool_t                              fTriggerStringFromPatches;  ///< Switch on building a trigger string based on available trigger patches

  AliCutValueRange<double>            fCentralityRange;           ///< Range of accepted event centralities
  AliCutValueRange<double>            fVertexRange;               ///< Range of the z-position in order to accept a primary vertex
  Bool_t                              fRequestCentrality;         ///< Switch for request of centrality selection

  TString                             fNameDownscaleOADB;         ///< Name of the downscale OADB container
  AliOADBContainer                    *fDownscaleOADB;            //!<! Container with downscale factors for different triggers
  TObjArray                           *fDownscaleFactors;         //!<! Downscalefactors for given run

  Int_t                               fCurrentRun;                ///< Current run number (for RunChange method)
  Bool_t                              fInitialized;               ///< Check for initialized

private:

  class EnergyBinning : public TCustomBinning {
  public:
    EnergyBinning();
    virtual ~EnergyBinning() {}
  };

  AliAnalysisTaskEmcalPatchesRef(const AliAnalysisTaskEmcalPatchesRef &);
  AliAnalysisTaskEmcalPatchesRef &operator=(const AliAnalysisTaskEmcalPatchesRef &);

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskEmcalPatchesRef, 1);
  /// \endcond
};

} /* namespace EMCalTriggerPtAnalysis */

#endif /* ALIANALYSISTASKEMCALPATCHESREF_H */
