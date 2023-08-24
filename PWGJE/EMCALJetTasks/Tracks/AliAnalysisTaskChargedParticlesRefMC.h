#ifndef ALIANALYSISTASKCHARGEDPARTICLESREFMC_H
#define ALIANALYSISTASKCHARGEDPARTICLESREFMC_H
/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliAnalysisTaskEmcal.h"
#include "AliCutValueRange.h"
#include <TString.h>
#include <TCustomBinning.h>

#include <vector>

class TArrayD;
class TClonesArray;
class THistManager;

class AliAnalysisUtils;
class AliAODTrack;
class AliEMCALGeometry;
class AliESDtrack;
class AliGenPythiaEventHeader;
class AliVParticle;
class AliMCEvent;
class AliEmcalTrackSelection;

namespace PWGJE{
  
namespace EMCALJetTasks {

class AliEMCalTriggerWeightHandler;
class AliEmcalTriggerOfflineSelection;

/**
 * @class AliAnalysisTaskChargedParticlesRefMC
 * @brief Test class for charged particle distributions (MC case)
 * @ingroup PWGJETASKS
 * @author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 * @since July 13, 2015
 *
 * As generator of reference distributions this task is written as "must work",
 * meaning as simple as possible and as independent as possible. Using only
 * well-tested components. Objects responsible for further problems, i.e. the
 * usage of THnSparse due to memory problems in several places, are forbidden.
 */
class AliAnalysisTaskChargedParticlesRefMC: public AliAnalysisTaskEmcal {
public:
  /**
   * @enum BeamDirection_t
   * @brief Direction of the beams
   *
   * In asymmetric collision systems the beam direction matters when calculating
   * the rapidity in the cms-frame. The sign will be used to correct the rapidity
   * shift.
   */
  enum BeamDirection_t{
    kpPb = 1,       ///< p-Pb
    kPbp = -1       ///< Pb-p
  };

  /**
   * @brief Dummy constructor
   */
  AliAnalysisTaskChargedParticlesRefMC();

  /**
   * @brief Main constructor
   * @param[in] name Name of the task
   */
  AliAnalysisTaskChargedParticlesRefMC(const char *name);

  /**
   * @brief Destuctor
   */
  virtual ~AliAnalysisTaskChargedParticlesRefMC();

  /**
   * @brief Enable Sumw2 when creating the histograms.
   *
   * Attention: Enabling Sumw2 will increase memory consumption significantly. Option should only be
   * used in case histograms are filled with a weight.
   * @param[in] doEnable If true Sumw2 is enabled for all histograms
   */
  void EnableSumw2(Bool_t doEnable) { fEnableSumw2 = doEnable; }

  /**
   * @brief Set rapidity shift originating from the asymmetric collision system.
   *
   * The rapidity shift is applied to tracks in order to get the rapidity
   * in the cms system.
   *
   * @param[in] yshift Rapidity shift applied to tracks
   */
  void SetRapidityShift(Double_t yshift) { fYshift = yshift; }

  /**
   * Define direction of the beams in the asymmetric collision system (p-Pb or Pb-p)
   * @param[in] beamdir Direction of the beams
   */
  void SetBeamDirection(BeamDirection_t beamdir) { fEtaSign = static_cast<Double_t>(beamdir); }

  /**
   * Set AliAnalysisUtils object. It will be used in the event selection.
   * @param[in] util Analysis Utils object used in the event selection
   */
  void  SetAnalysisUtil(AliAnalysisUtils *util) { fAliAnalysisUtils = util; }

  /**
   * @brief Set the virtual track selection.
   *
   * Assumes that the track selection object is already fully configured
   * @param[in] sel Virtual track selection
   */
  void SetEMCALTrackSelection(AliEmcalTrackSelection *sel) { fTrackCuts = sel; }

  /**
   * Define kinematic cut to particles and tracks on \f$ \eta \f$ in the lab frame. By default
   * the accepted range is defined to be within [-0.8, 0.8].
   * @param[in] etamin Min. allowed \f$ \eta \f$
   * @param[in] etamax Max. allowed \f$ \eta \f$
   */
  void SetEtaLabCut(double etamin, double etamax) { fEtaLabCut.SetLimits(etamin, etamax); }

  /**
   * Define kinematic cut to particles and tracks on \f$ \eta \f$ in the cms frame. By default
   * no cut is applied.
   * @param[in] etamin Min. allowed \f$ \eta \f$
   * @param[in] etamax Max. allowed \f$ \eta \f$
   */
  void SetEtaCMSCut(double etamin, double etamax) { fEtaCmsCut.SetLimits(etamin, etamax); }

  /**
   * Define kinematic cut the particles and tracks on \f$ \phi \f$ . By default no cut is applied.
   * @param[in] phimin Min. allowed \f$ \phi \f$
   * @param[in] phimax Max. allowed \f$ \phi \f$
   */
  void SetTrackPhiCut(double phimin, double phimax) { fPhiCut.SetLimits(phimin, phimax); }

  /**
   * @brief Set minimum \f$ p_{t}\f$ used to select track candidate
   * @param[in] minpt Minimum \f$ p_{t}\f$
   */
  void SetMinPtTracks(Double_t minpt) { fMinPt = minpt; }

  /**
   * @brief Set offline trigger selection.
   *
   * The trigger selection will be used to select events for the given trigger class given the presence
   * of a patch with energy above threshold.
   *
   * @param[in] sel Trigger offline selection
   */
  void SetOfflineTriggerSelection(AliEmcalTriggerOfflineSelection *sel) { fTriggerSelection = sel; }

  /**
   * @brief Require bunch crossing information of track obtained from TOF (if available)
   * matches the bunch crossing ID of the event
   *
   * Pileup is not simulated in ALICE. This means on simulation level only the TOF matching is required
   * in order to check the efficiency of the TOF matching
   *
   * @param[in] doRequire If true the track is only selected if it has a TOF match
   */
  void SetRequireTOFBunchCrossing(Bool_t doRequire) { fRequireTOFBunchCrossing = doRequire; }

  /**
   * @brief Get the trigger offline selection.
   *
   * Providing access to the trigger offline selection. Note that
   * it must be initialized from outside the task.
   *
   * @return Trigger offline selection of this instance (nullptr if not defined)
   */
  AliEmcalTriggerOfflineSelection *GetOfflineTriggerSelection() const { return fTriggerSelection; }

  /**
   * @brief Enable PID-related plots.
   *
   * Monitoring track counts for different particle species.
   * Attention: PID plots are stored as THnSparses, therefore
   * they might be big.
   *
   * @param[in] plotPID If true then PID-related plots are enabled
   */
  void SetPlotPID(Bool_t plotPID) { fStudyPID = plotPID; }

  /**
   * Initializing pre-defined track cuts based on a track cuts name. The virtual track selection
   * is handled for ESDs or AODs according to whether isAOD is false or true
   * @param[in] cutname Name of the track cut configuration
   * @param[in] isAOD Switch distinguishing between ESD and AOD
   */
  void InitializeTrackCuts(TString cutname, bool isAOD);

  /**
   * Set external weight handler. The weight handler is used to handle cross section weights
   * for different \f$ p_{t} \f$-hard bins
   * @param[in] wh Weight handler for \f$ p_{t} \f$-hard bins
   */
  void SetWeightHandler(const AliEMCalTriggerWeightHandler * wh) { fWeightHandler = wh; }

  /**
   * Set the name of the OADB container with the trigger acceptance maps. It will be used
   * in the trigger offline selection to mimic the acceptance profile observed in data.
   * @param[in] name Name of the OADB container with the trigger acceptance map
   */
  void SetTriggerAcceptanceOADB(const TString &name) { fNameAcceptanceOADB = name; }

  /**
   * @brief Add histograms for tracks pointing to EMCAL supermodules
   * @param[in] doStudy If true histograms are added
   */
  void SetStudyEMCALgeo(Bool_t doStudy) { fStudyEMCALgeo = doStudy; }

  /**
   * @brief Switch for whether the analysis runs only in min. bias mode
   * 
   * In min. bias mode EMCAL triggers are not evaluated. Also the connection
   * to the EMCAL trigger patch container is disabled.
   * 
   * @param[in] exclusiveMinBias Switch for whether the analysis runs in min. bias or full mode
   */
  void SetExclusiveMinBias(Bool_t exclusiveMinBias) { fExclusiveMinBias = exclusiveMinBias; this->SetCaloTriggerPatchInfoName(fExclusiveMinBias ? "" : "EmcalTriggers"); }

  /**
   * Preconfigure task. Intended for subwagon configuration inside the train.
   * @param[in] suffix Suffix of the subwagon
   * @return Preconfigured task
   */
  static AliAnalysisTaskChargedParticlesRefMC *AddTaskChargedParticlesRefMC(const TString &suffix);

  /**
   * Configure task fully with a default configuration. Not intended for subwagons
   * @param[in] cutname Name of the cut configuration
   * @return Fully-configured task
   */
  static AliAnalysisTaskChargedParticlesRefMC *AddTaskChargedParticlesRefMCDefault(const TString &cutname = "standard");

protected:

  /**
   * @brief Create the output histograms.
   *
   * Histograms are available for three different categories
   * -# Events
   * -# Tracks
   * -# Trigger jets
   * Histograms are internally handled by a THistManager, and appended to the
   * common output container.
   */
  virtual void UserCreateOutputObjects();

  /**
   * @brief Main event loop

   * Fill histograms for MC-true particles and reconstructed
   * true particles (reconstructed particles with assigned MC-truth information),
   * and fill several kinematic histograms.
   * -# track-\f$ p_{t} \f$
   * -# track-\$f \eta \f$, both cms and lab frame
   * -# track-\f$ \phi \f$
   * In addition restrictions are done to the EMCAL acceptance and the presence
   * of TRD material.
   *
   * Histograms are filled only for events accepted in the member function
   * IsEventSelected().
   *
   * @return Always true
   */
  virtual bool  Run();

  /**
   * @brief Apply standard event selection.
   *
   * The standard event selection is consisting of
   * - Trigger selection
   * - Vertex selection
   * - Pileup rejection
   *
   * In addition also a cut on the MC outliers based on the PYTHIA trigger jet
   * information is applied.
   *
   * Trigger classes connected to the event are stored in the event and can be
   * used in the run method.
   *
   * @return True if the event is selected
   */
  virtual bool IsEventSelected();

  /**
   * @brief Loading initial settings at the beginning of the job.
   *
   * Implentation of framework function ExceOnce: Loading acceptance map from the OADB
   * in case a valid OADB container is provided, and forwarding it to the EMCAL trigger
   * offline selection.
   */
  virtual void ExecOnce();

  /**
   * @brief Fill track histograms
   * @param[in] eventclass Trigger class fired
   * @param[in] weight \f$ p_{t} \f$-hard dependent weight
   * @param[in] posCharge Switch whether track charge is positive or negative
   * @param[in] pt track \f$ p_{t} \f$
   * @param[in] etalab Track \f$ \eta \f$ in lab frame
   * @param[in] etacent Track \f$ \eta \f$ in cms frame
   * @param[in] phi Track \f$ \eta \f$ in lab frame
   * @param[in] inEmcal Track in EMCAL \f$ \phi \f$ acceptance
   * @param[in] isPrimary If true then the particle is a primary particle
   * @param[in] pid True particle species
   */
  void FillTrackHistos(const TString &eventclass, Double_t weight, Bool_t posCharge, Double_t pt, Double_t eta, Double_t etacent, Double_t phi, Bool_t inEmcal, Bool_t isPrimary, const TString &pid);

  /**
   * Apply trigger selection using offline patches and trigger thresholds based on offline ADC Amplitude
   * @param[in] triggerpatches Trigger patches found by the trigger maker
   * @return String with EMCAL trigger decision
   */
  TString GetFiredTriggerClasses(const TClonesArray * triggerpatches);

  /**
   * Check in a transparent way for ESDs and AODs whether the particle is physical primary or not
   * -# AOD: Information stored in the AliAODMCParticle
   * -# ESD: Information needs to be retrieved from the stack via the label of the MC particle
   * @param[in] part The particle to check
   * @param[in] mcevent The MC event containing the stack (ESD only)
   * @return True if particle is a physical primary particle, false otherwise
   */
  Bool_t IsPhysicalPrimary(const AliVParticle *const part, AliMCEvent *const mcevent);

private:

  /**
   * @class PtBinning
   * @brief Default \f$ p_{t} \f$ binning for the charged particle \f$ p_{t} \f$ analysis
   *
   * Object oriented approach for \f$ p_{t} \f$ binning: Handling in a class,
   * inheriting from TCustomBinning. Can be used for all histograms handled by
   * the hist manager which rely on the given binning.
   */
  class PtBinning : public TCustomBinning{
  public:

    /**
     * Constructor
     */
    PtBinning();

    /**
     * Destructor
     */
    virtual ~PtBinning() {}
  };

  AliAnalysisTaskChargedParticlesRefMC(const AliAnalysisTaskChargedParticlesRefMC &);
  AliAnalysisTaskChargedParticlesRefMC &operator=(const AliAnalysisTaskChargedParticlesRefMC &);

  AliEmcalTrackSelection    *fTrackCuts;                ///< Standard track selection
  AliEmcalTriggerOfflineSelection       *fTriggerSelection;         ///< Offline trigger selection
  THistManager                          *fHistos;                   ///< Histogram manager
  const AliEMCalTriggerWeightHandler    *fWeightHandler;            ///< Weight handler (optional)

  std::vector<TString>                  fEventTriggers;             //!<! Temporary container for selected triggers
  Double_t                              fEventWeight;               ///< Event weight

  Double_t                              fYshift;                    ///< Rapidity shift
  Double_t                              fEtaSign;                   ///< Sign of the \f$\eta\f$-distribution (swaps when beam directions swap): p-Pb: +1, Pb-p: -1
  Double_t                              fMinPt;                     ///< Minimum \f$ p_{t}\f$-used to select tracks
  AliCutValueRange<double>              fEtaLabCut;                 ///< Cut applied in \f$\eta_{Lab}\f$-frame
  AliCutValueRange<double>              fEtaCmsCut;                 ///< Cut applied in \f$\eta_{centre-of-mass}\f$ frame
  AliCutValueRange<double>              fPhiCut;                    ///< Track \f$\phi\f$ cut
  Double_t                              fFracPtHard;                ///< Cut on the maximum fraction of pt hard of any trigger jet
  Bool_t                                fEnableSumw2;               ///< Enable sumw2 when filling histograms (by default off)
  Bool_t                                fStudyPID;                  ///< Fill kinematics histograms with information of true particle species (default: off)
  Bool_t                                fStudyEMCALgeo;             ///< Add histograms for tracks pointing to the EMCAL acceptance
  Bool_t                                fExclusiveMinBias;          ///< Study only min. bias;
  Bool_t                                fRequireTOFBunchCrossing;   ///< Require TOF bunch crossing information (reducing to simple TOF matching in simulation)

  TString                               fNameAcceptanceOADB;        ///< Name of the OADB container with trigger acceptance maps

  ClassDef(AliAnalysisTaskChargedParticlesRefMC, 1);
};

} /* namespace EMCALJetTasks */

} /* namespace PWGJE */

#endif /* ALIANALYSISTASKCHARGEDPARTICLESREFMC_H */
