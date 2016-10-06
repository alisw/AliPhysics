#ifndef ALIANALYSISTASKCHARGEDPARTICLESREF_H
#define ALIANALYSISTASKCHARGEDPARTICLESREF_H
/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliAnalysisTaskEmcalTriggerBase.h"
#include "AliCutValueRange.h"
#include <vector>
#include <TCustomBinning.h>
#include <TString.h>

class THistManager;
class AliEmcalTrackSelection;

namespace EMCalTriggerPtAnalysis {

/**
 * \class AliAnalysisTaskChargedParticlesRef
 * \brief Unit test class for charged particle distributions
 *
 * As generator of reference distributions this task is written as "must work",
 * meaning as simple as possible and as independent as possible. Using only
 * well-tested components. Objects responsible for further problems, i.e. the
 * usage of THnSparse due to memory problems in several places, are forbidden.
 */
class AliAnalysisTaskChargedParticlesRef : public AliAnalysisTaskEmcalTriggerBase {
public:
  enum BeamDirection_t{
    kpPb = 1,
    kPbp = -1
  };
  AliAnalysisTaskChargedParticlesRef();
  AliAnalysisTaskChargedParticlesRef(const char *name);
  virtual ~AliAnalysisTaskChargedParticlesRef();

  void SetRapidityShift(Double_t yshift) { fYshift = yshift; }
  void SetBeamDirection(BeamDirection_t beamdir) { fEtaSign = static_cast<Double_t>(beamdir); }

  void SetKineCorrelation(Bool_t useKineCorrelation) { fKineCorrelation = useKineCorrelation; }
  void SetPlotPID(Bool_t plotPID) { fStudyPID = plotPID; }
  void SetEMCALTrackSelection(AliEmcalTrackSelection *sel) { fTrackCuts = sel; }
  void InitializeTrackCuts(TString cutname, bool isAOD);
  void SetEtaLabCut(double etamin, double etamax) { fEtaLabCut.SetLimits(etamin, etamax); }
  void SetEtaCMSCut(double etamin, double etamax) { fEtaCmsCut.SetLimits(etamin, etamax); }
  void SetTrackPhiCut(double phimin, double phimax) { fPhiCut.SetLimits(phimin, phimax); }

protected:
  virtual void CreateUserHistos();
  virtual void CreateUserObjects();
  virtual bool Run();

  virtual void UserFillHistosBeforeEventSelection();
  virtual void UserFillHistosAfterEventSelection();

  void FillEventCounterHists(const TString &triggerclass, double vtxz, bool isSelected);
  void FillTrackHistos(const TString &eventclass, Double_t pt, Double_t eta, Double_t etacent, Double_t phi, Bool_t etacut, Bool_t inEmcal, Bool_t hasTRD);
  void FillPIDHistos(const TString &eventclass, const AliVTrack &track);

  AliEmcalTrackSelection          *fTrackCuts;                ///< Standard track selection

  Double_t                        fYshift;                    ///< Rapidity shift
  Double_t                        fEtaSign;                   ///< Sign of the eta distribution (swaps when beam directions swap): p-Pb: +1, Pb-p: -1

  AliCutValueRange<double>        fEtaLabCut;                 ///< Cut applied in \f$\eta_{lab}\f$ frame
  AliCutValueRange<double>        fEtaCmsCut;                 ///< Cut applied in \f$\eta_{centre-of-mass}\f$ frame
  AliCutValueRange<double>        fPhiCut;                    ///< Cut applied in \f$\phi_{track}\f$

  Bool_t                          fKineCorrelation;           ///< Use kinematics correlation histograms
  Bool_t                          fStudyPID;                  ///< Use kinematics correlation histograms

private:

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
