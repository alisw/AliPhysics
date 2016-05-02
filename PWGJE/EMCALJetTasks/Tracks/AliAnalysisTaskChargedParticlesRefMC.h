#ifndef ALIANALYSISTASKCHARGEDPARTICLESREFMC_H
#define ALIANALYSISTASKCHARGEDPARTICLESREFMC_H
/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliAnalysisTaskSE.h"
#include <TString.h>

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
class AliAnalysisTaskChargedParticlesRefMC: public AliAnalysisTaskSE {
public:
  enum BeamDirection_t{
    kpPb = 1,
    kPbp = -1
  };

  AliAnalysisTaskChargedParticlesRefMC();
  AliAnalysisTaskChargedParticlesRefMC(const char *name);
  virtual ~AliAnalysisTaskChargedParticlesRefMC();

  void                        UserCreateOutputObjects();
  void                        UserExec(Option_t *);
  Bool_t                      UserNotify();
  void                        Terminate(Option_t *) {}

  void                        SetRapidityShift(Double_t yshift) { fYshift = yshift; }
  void                        SetBeamDirection(BeamDirection_t beamdir) { fEtaSign = static_cast<Double_t>(beamdir); }

  void                        SetAnalysisUtil(AliAnalysisUtils *util) { fAnalysisUtil = util; }
  void                        SetTrackSelection(AliEmcalTrackSelection * sel) { fTrackCuts = sel; }
  void                        SetEtaLabCut(double etamin, double etamax) { fEtaLabCut[0] = etamin; fEtaLabCut[1] = etamax; }
  void                        SetEtaCMSCut(double etamin, double etamax) { fEtaCmsCut[0] = etamin; fEtaCmsCut[1] = etamax; }
  void                        SetOutlierCut(double fracpthard = 1.2) { fFracPtHard = fracpthard; }
  void                        SetOfflineTriggerSelection(AliEmcalTriggerOfflineSelection *sel) { fTriggerSelection = sel; }
  void                        InitializeTrackCuts(TString cutname, bool isAOD);
  void                        SetWeightHandler(const AliEMCalTriggerWeightHandler * wh) { fWeightHandler = wh; }

protected:
  void                        CreateOldPtBinning(TArrayD &binning) const;
  void                        CreateNewPtBinning(TArrayD &binning) const;

  void                        FillTrackHistos(const char *eventclass, Double_t weight, Double_t pt, Double_t eta, Double_t etacent, Double_t phi, Bool_t etacut, Bool_t inEmcal, Bool_t hasTRD, const char *pid);
  void                        FillTriggerJetHistograms(Bool_t aftercut, AliGenPythiaEventHeader *const header);

  TString                     GetFiredTriggerClasses(const TClonesArray * triggerpatches);
  Bool_t                      PythiaInfoFromFile(const char* currFile, Float_t &fXsec, Float_t &fTrials, Int_t &pthard) const;
  AliGenPythiaEventHeader    *GetPythiaHeader() const;
  Bool_t                      IsPhysicalPrimary(const AliVParticle *const part, AliMCEvent *const mcevent);
  Bool_t                      IsOutlier(AliGenPythiaEventHeader * const header) const;


  AliEmcalTrackSelection                *fTrackCuts;                ///< Standard track selection
  AliAnalysisUtils                      *fAnalysisUtil;             ///< Event selection
  AliEmcalTriggerOfflineSelection       *fTriggerSelection;         ///< Offline trigger selection
  THistManager                          *fHistos;                   ///< Histogram manager
  AliEMCALGeometry                      *fGeometry;                 ///< EMCAL geometry methods
  const AliEMCalTriggerWeightHandler    *fWeightHandler;            ///< Weight handler (optional)

  // Monte-Carlo specific information
  Double_t                        fPtHard;                    ///< event pt hard
  Int_t                           fPtHardBin;                 ///< event pt hard bin
  Int_t                           fNTrials;                   ///< event trials
  Float_t                         fXsection;                  ///< x-section from pythia header

  Double_t                        fYshift;                    ///< Rapidity shift
  Double_t                        fEtaSign;                   ///< Sign of the eta distribution (swaps when beam directions swap): p-Pb: +1, Pb-p: -1
  Double_t                        fEtaLabCut[2];              ///< Cut applied in Eta Lab frame
  Double_t                        fEtaCmsCut[2];              ///< Cut applied in Eta centre-of-mass frame
  Double_t                        fFracPtHard;                ///< Cut on the maximum fraction of pt hard of any trigger jet

private:
  AliAnalysisTaskChargedParticlesRefMC(const AliAnalysisTaskChargedParticlesRefMC &);
  AliAnalysisTaskChargedParticlesRefMC &operator=(const AliAnalysisTaskChargedParticlesRefMC &);

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskChargedParticlesRefMC, 1);
  /// \endcond
};

} /* namespace EMCalTriggerPtAnalysis */

#endif /* ALIANALYSISTASKCHARGEDPARTICLESREFMC_H */
