#ifndef ALIANALYSISTASKCHARGEDPARTICLESREF_H
#define ALIANALYSISTASKCHARGEDPARTICLESREF_H
/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliAnalysisTaskSE.h"

class TArrayD;
class TString;
class AliAnalysisUtils;
class AliAODTrack;
class AliEMCALGeometry;
class AliESDtrackCuts;
class AliESDtrack;

namespace EMCalTriggerPtAnalysis {

class AliEMCalHistoContainer;

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
  enum EmcalTriggerClass{
    kCPREL0 = 0,
    kCPREG1,
    kCPREG2,
    kCPREJ1,
    kCPREJ2,
    kCPRntrig
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

  void SetEtaLabCut(double etamin, double etamax) { fEtaLabCut[0] = etamin; fEtaLabCut[1] = etamax; }
  void SetEtaCMSCut(double etamin, double etamax) { fEtaCmsCut[0] = etamin; fEtaCmsCut[1] = etamax; }

  void SetOfflineEnergyThreshold(EmcalTriggerClass trgcls, double threshold) { fOfflineEnergyThreshold[trgcls] = threshold; }

protected:
  void CreateOldPtBinning(TArrayD &binning) const;
  void CreateNewPtBinning(TArrayD &binning) const;

  void FillEventCounterHists(const char *triggerclass, double vtxz, bool isSelected);
  void FillTrackHistos(const char *eventclass, Double_t pt, Double_t eta, Double_t etacent, Double_t phi, Bool_t etacut, Bool_t inEmcal, Bool_t hasTRD);
  TString GetFiredTriggerClassesFromPatches(const TClonesArray* triggerpatches) const;
  Bool_t IsOfflineSelected(EmcalTriggerClass trgcls, const TClonesArray * const triggerpatches) const;

  Bool_t TrackSelectionESD(AliESDtrack *track);
  Bool_t TrackSelectionAOD(AliAODTrack *track);

  AliESDtrackCuts                 *fTrackCuts;                ///< Standard track selection
  AliAnalysisUtils                *fAnalysisUtil;             ///< Event selection
  AliEMCalHistoContainer          *fHistos;                   ///< Histogram manager
  AliEMCALGeometry                *fGeometry;                 ///< EMCAL geometry methods

  Bool_t                          fTriggerStringFromPatches;  ///< Do rebuild the trigger string from trigger patches
  Double_t                        fYshift;                    ///< Rapidity shift
  Double_t                        fEtaSign;                   ///< Sign of the eta distribution (swaps when beam directions swap): p-Pb: +1, Pb-p: -1

  Double_t                        fEtaLabCut[2];              ///< Cut applied in Eta Lab frame
  Double_t                        fEtaCmsCut[2];              ///< Cut applied in Eta centre-of-mass frame

  Double_t                        fOfflineEnergyThreshold[kCPRntrig];    ///< Threhold applied on offline patches

private:
  AliAnalysisTaskChargedParticlesRef(const AliAnalysisTaskChargedParticlesRef &);
  AliAnalysisTaskChargedParticlesRef &operator=(const AliAnalysisTaskChargedParticlesRef &);

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskChargedParticlesRef, 1);
  /// \endcond
};

} /* namespace EMCalTriggerPtAnalysis */

#endif /* ALIANALYSISTASKCHARGEDPARTICLESREF_H */
