#ifndef ALIANALYSISTASKCHARGEDPARTICLESREFMC_H
#define ALIANALYSISTASKCHARGEDPARTICLESREFMC_H
/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliAnalysisTaskEmcal.h"
#include "AliCutValueRange.h"
#include <TString.h>
#include <TCustomBinning.h>

#include <string>
#include <vector>

class TArrayD;
class TClonesArray;
class THistManager;

class AliAnalysisUtils;
class AliAODTrack;
class AliEMCALGeometry;
class AliESDtrack;
class AliEmcalTrackSelection;
class AliGenPythiaEventHeader;
class AliVParticle;
class AliMCEvent;

namespace EMCalTriggerPtAnalysis {

class AliEMCalTriggerWeightHandler;

/**
 * \class AliAnalysisTaskChargedParticlesRefMC
 * \brief Unit test class for charged particle distributions (MC case)
 *
 * As generator of reference distributions this task is written as "must work",
 * meaning as simple as possible and as independent as possible. Using only
 * well-tested components. Objects responsible for further problems, i.e. the
 * usage of THnSparse due to memory problems in several places, are forbidden.
 */
class AliAnalysisTaskChargedParticlesRefMC: public AliAnalysisTaskEmcal {
public:
  enum BeamDirection_t{
    kpPb = 1,
    kPbp = -1
  };

  AliAnalysisTaskChargedParticlesRefMC();
  AliAnalysisTaskChargedParticlesRefMC(const char *name);
  virtual ~AliAnalysisTaskChargedParticlesRefMC();

  void                        SetRapidityShift(Double_t yshift) { fYshift = yshift; }
  void                        SetBeamDirection(BeamDirection_t beamdir) { fEtaSign = static_cast<Double_t>(beamdir); }

  void                        SetAnalysisUtil(AliAnalysisUtils *util) { fAliAnalysisUtils = util; }
  void                        SetTrackSelection(AliEmcalTrackSelection * sel) { fTrackCuts = sel; }
  void                        SetEtaLabCut(double etamin, double etamax) { fEtaLabCut.SetLimits(etamin, etamax); }
  void                        SetEtaCMSCut(double etamin, double etamax) { fEtaCmsCut.SetLimits(etamin, etamax); }
  void                        SetOfflineTriggerSelection(AliEmcalTriggerOfflineSelection *sel) { fTriggerSelection = sel; }
  void                        InitializeTrackCuts(TString cutname, bool isAOD);
  void                        SetWeightHandler(const AliEMCalTriggerWeightHandler * wh) { fWeightHandler = wh; }

protected:
  virtual void                UserCreateOutputObjects();
  virtual bool                Run();
  virtual bool                IsEventSelected();

  void                        FillTrackHistos(const char *eventclass, Double_t weight, Double_t pt, Double_t eta, Double_t etacent, Double_t phi, Bool_t etacut, Bool_t inEmcal, Bool_t hasTRD, const char *pid);
  void                        FillTriggerJetHistograms(Bool_t aftercut, AliGenPythiaEventHeader *const header);

  TString                     GetFiredTriggerClasses(const TClonesArray * triggerpatches);
  Bool_t                      IsPhysicalPrimary(const AliVParticle *const part, AliMCEvent *const mcevent);

private:

  class PtBinning : public TCustomBinning{
  public:
    PtBinning();
    virtual ~PtBinning() {}
  };

  AliAnalysisTaskChargedParticlesRefMC(const AliAnalysisTaskChargedParticlesRefMC &);
  AliAnalysisTaskChargedParticlesRefMC &operator=(const AliAnalysisTaskChargedParticlesRefMC &);

  AliEmcalTrackSelection                *fTrackCuts;                ///< Standard track selection
  AliEmcalTriggerOfflineSelection       *fTriggerSelection;         ///< Offline trigger selection
  THistManager                          *fHistos;                   ///< Histogram manager
  const AliEMCalTriggerWeightHandler    *fWeightHandler;            ///< Weight handler (optional)

  std::vector<std::string>              fEventTriggers;             //!<! Temporary container for selected triggers
  Double_t                              fEventWeight;               ///< Event weight

  Double_t                              fYshift;                    ///< Rapidity shift
  Double_t                              fEtaSign;                   ///< Sign of the eta distribution (swaps when beam directions swap): p-Pb: +1, Pb-p: -1
  AliCutValueRange<double>              fEtaLabCut;                 ///< Cut applied in Eta Lab frame
  AliCutValueRange<double>              fEtaCmsCut;                 ///< Cut applied in Eta centre-of-mass frame
  Double_t                              fFracPtHard;                ///< Cut on the maximum fraction of pt hard of any trigger jet

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskChargedParticlesRefMC, 1);
  /// \endcond
};

} /* namespace EMCalTriggerPtAnalysis */

#endif /* ALIANALYSISTASKCHARGEDPARTICLESREFMC_H */
