#ifndef ALIANALYSISTASKCHARGEDPARTICLESMCTRIGGERMIMIC_H
#define ALIANALYSISTASKCHARGEDPARTICLESMCTRIGGERMIMIC_H
/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliAnalysisTaskEmcal.h"
#include "AliCutValueRange.h"
#include <TCustomBinning.h>
#include <TString.h>

class TClonesArray;
class THistManager;

class AliAnalysisUtils;
class AliEmcalTrackSelection;
class AliGenPythiaEventHeader;
class AliMCEvent;
class AliVParticle;

namespace EMCalTriggerPtAnalysis {

class AliEMCalTriggerWeightHandler;

/**
 * @class AliAnalysisTaskChargedParticlesMCTriggerMimic
 * @brief Lightweight task, mimicing EMCAL triggers in MC
 * @author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 *
 * Similar to AliAnalysisTaskChargedParticlesRefMC, however:
 * - No split in trigger classes - this has to be don on subwagon level
 */
class AliAnalysisTaskChargedParticlesMCTriggerMimic : public AliAnalysisTaskEmcal {
public:

  /**
   * @enum BeamDirection_t
   * @brief Direction of the beam, defined as Pb-going
   */
  enum BeamDirection_t{
    kpPb = 1,   //!< p-Pb (\f$ \eta \f$-sign positive)
    kPbp = -1   //!< Pb-p (\f$ \eta \f$-sign negative)
  };

  /**
   * @enum PatchType_t
   * @brief Type of the patch used in the EMCAL trigger selection
   */
  enum PatchType_t {
    kEMCEGA = 0,    //!< EMCAL gamma patch
    kEMCEJE = 1,    //!< EMCAL jet patch
    kUndef = -1     //!< Patch type not defined
  };

  AliAnalysisTaskChargedParticlesMCTriggerMimic();
  AliAnalysisTaskChargedParticlesMCTriggerMimic(const char *name);
  virtual ~AliAnalysisTaskChargedParticlesMCTriggerMimic();

  void                        SetRapidityShift(Double_t yshift) { fYshift = yshift; }
  void                        SetBeamDirection(BeamDirection_t beamdir) { fEtaSign = static_cast<Double_t>(beamdir); }

  void                        SetAnalysisUtil(AliAnalysisUtils *util) { fAliAnalysisUtils = util; }
  void                        SetTrackSelection(AliEmcalTrackSelection * sel) { fTrackCuts = sel; }
  void                        SetEtaLabCut(double etamin, double etamax) { fEtaLabCut.SetLimits(etamin, etamax); }
  void                        SetEtaCMSCut(double etamin, double etamax) { fEtaCmsCut.SetLimits(etamin, etamax); }
  void                        InitializeTrackCuts(TString cutname, bool isAOD);
  void                        SetWeightHandler(const AliEMCalTriggerWeightHandler * wh) { fWeightHandler = wh; }

  void                        SetEMCALTrigger(PatchType_t patchtype, Double_t threshold) { fPatchType = patchtype; fEnergyThreshold = threshold; }

protected:

  virtual void                UserCreateOutputObjects();
  virtual Bool_t              IsEventSelected();
  virtual Bool_t              Run();

  void                        FillTrackHistos(const char *eventclass, Double_t weight, Double_t pt, Double_t eta, Double_t etacent, Double_t phi, Bool_t etacut, Bool_t inEmcal, const char *pid);
  Bool_t                      IsPhysicalPrimary(const AliVParticle* const part, AliMCEvent* const mcevent);

  Bool_t                      SelectEmcalTrigger(const TClonesArray * triggerpatches);

private:

  class PtBinning : public TCustomBinning{
  public:
    PtBinning();
    virtual ~PtBinning() {}
  };

  AliAnalysisTaskChargedParticlesMCTriggerMimic(const AliAnalysisTaskChargedParticlesMCTriggerMimic &);
  AliAnalysisTaskChargedParticlesMCTriggerMimic &operator=(const AliAnalysisTaskChargedParticlesMCTriggerMimic &);


  AliEmcalTrackSelection                *fTrackCuts;                ///< Standard track selection
  THistManager                          *fHistos;                   ///< Histogram manager
  const AliEMCalTriggerWeightHandler    *fWeightHandler;            ///< Weight handler (optional)

  // Kinematical selection
  Double_t                              fYshift;                    ///< Rapidity shift
  Double_t                              fEtaSign;                   ///< Sign of the eta distribution (swaps when beam directions swap): p-Pb: +1, Pb-p: -1
  AliCutValueRange<double>              fEtaLabCut;                 ///< Cut applied in Eta Lab frame
  AliCutValueRange<double>              fEtaCmsCut;                 ///< Cut applied in Eta centre-of-mass frame

  // Trigger specific information
  PatchType_t                           fPatchType;                 ///< Type of the trigger patch (default: kUndef)
  Double_t                              fEnergyThreshold;           ///< Energy threshold in patch


  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskChargedParticlesMCTriggerMimic, 1);
  /// \endcond
};

} /* namespace EMCalTriggerPtAnalysis */

#endif /* ALIANALYSISTASKCHARGEDPARTICLEMCTRIGGERMIMIC_H */
