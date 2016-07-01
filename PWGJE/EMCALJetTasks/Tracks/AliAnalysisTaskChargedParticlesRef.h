#ifndef ALIANALYSISTASKCHARGEDPARTICLESREF_H
#define ALIANALYSISTASKCHARGEDPARTICLESREF_H
/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliAnalysisTaskSE.h"
#include <string>
#include <TCustomBinning.h>

class TArrayD;
class THistManager;
class TString;
class AliAnalysisUtils;
class AliEmcalTrackSelection;
class AliEMCALGeometry;

namespace EMCalTriggerPtAnalysis {

class AliEmcalTriggerOfflineSelection;

/**
 * \class AliAnalysisTaskChargedParticlesRef
 * \brief Unit test class for charged particle distributions
 *
 * As generator of reference distributions this task is written as "must work",
 * meaning as simple as possible and as independent as possible. Using only
 * well-tested components. Objects responsible for further problems, i.e. the
 * usage of THnSparse due to memory problems in several places, are forbidden.
 */
class AliAnalysisTaskChargedParticlesRef : public AliAnalysisTaskSE {
public:
  enum BeamDirection_t{
    kpPb = 1,
    kPbp = -1
  };
  AliAnalysisTaskChargedParticlesRef();
  AliAnalysisTaskChargedParticlesRef(const char *name);
  virtual ~AliAnalysisTaskChargedParticlesRef();

  void UserCreateOutputObjects();
  void UserExec(Option_t *);
  void Terminate(Option_t *) {}

  void SetRapidityShift(Double_t yshift) { fYshift = yshift; }
  void SetBeamDirection(BeamDirection_t beamdir) { fEtaSign = static_cast<Double_t>(beamdir); }
  void UseTriggerPatches(Bool_t doUse) { fTriggerStringFromPatches = doUse; }

  void SetTrackSelection(AliEmcalTrackSelection *sel) { fTrackCuts = sel; }
  void InitializeTrackCuts(TString cutname, bool isAOD);
  void SetOfflineTriggerSelection(AliEmcalTriggerOfflineSelection *sel) { fTriggerSelection = sel; }
  void SetAnalysisUtil(AliAnalysisUtils *util) { fAnalysisUtil = util; }
  void SetEtaLabCut(double etamin, double etamax) { fEtaLabCut[0] = etamin; fEtaLabCut[1] = etamax; }
  void SetEtaCMSCut(double etamin, double etamax) { fEtaCmsCut[0] = etamin; fEtaCmsCut[1] = etamax; }

protected:
  void FillEventCounterHists(const std::string &triggerclass, double vtxz, bool isSelected);
  void FillTrackHistos(const std::string &eventclass, Double_t pt, Double_t eta, Double_t etacent, Double_t phi, Bool_t etacut, Bool_t inEmcal, Bool_t hasTRD);
  void FillPIDHistos(const std::string &eventclass, const AliVTrack &track);
  TString GetFiredTriggerClassesFromPatches(const TClonesArray* triggerpatches) const;

  AliEmcalTrackSelection          *fTrackCuts;                ///< Standard track selection
  AliAnalysisUtils                *fAnalysisUtil;             ///< Event selection
  AliEmcalTriggerOfflineSelection *fTriggerSelection;         ///< Offline trigger selection
  THistManager                    *fHistos;                   ///< Histogram manager
  AliEMCALGeometry                *fGeometry;                 ///< EMCAL geometry methods

  Bool_t                          fTriggerStringFromPatches;  ///< Do rebuild the trigger string from trigger patches
  Double_t                        fYshift;                    ///< Rapidity shift
  Double_t                        fEtaSign;                   ///< Sign of the eta distribution (swaps when beam directions swap): p-Pb: +1, Pb-p: -1

  Double_t                        fEtaLabCut[2];              ///< Cut applied in Eta Lab frame
  Double_t                        fEtaCmsCut[2];              ///< Cut applied in Eta centre-of-mass frame

private:

  class OldPtBinning : public TCustomBinning{
  public:
    OldPtBinning();
    virtual ~OldPtBinning() {}
  };
  class NewPtBinning : public TCustomBinning{
  public:
    NewPtBinning();
    virtual ~NewPtBinning() {}
  };

  AliAnalysisTaskChargedParticlesRef(const AliAnalysisTaskChargedParticlesRef &);
  AliAnalysisTaskChargedParticlesRef &operator=(const AliAnalysisTaskChargedParticlesRef &);

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskChargedParticlesRef, 1);
  /// \endcond
};

} /* namespace EMCalTriggerPtAnalysis */

#endif /* ALIANALYSISTASKCHARGEDPARTICLESREF_H */
