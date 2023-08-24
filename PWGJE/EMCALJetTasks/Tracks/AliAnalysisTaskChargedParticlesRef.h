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

namespace PWGJE {

namespace EMCALJetTasks {

/**
 * @class AliAnalysisTaskChargedParticlesRef
 * @brief Test class for charged particle distributions
 * @ingroup PWGJETASKS
 * @author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 * @since Jul 9, 2015
 *
 * As generator of reference distributions this task is written as "must work",
 * meaning as simple as possible and as independent as possible. Using only
 * well-tested components. Objects responsible for further problems, i.e. the
 * usage of THnSparse due to memory problems in several places, are forbidden.
 */
class AliAnalysisTaskChargedParticlesRef : public AliAnalysisTaskEmcalTriggerBase {
public:

  /**
   * @enum BeamDirection_t
   * @brief Direction of the beams in an asymmetric collision system
   *
   * In case of an asymmetric collision system, a rapidity shift is
   * applied to convert from the lab frame to the centre-of-mass frame.
   * This enum is used to determine the sign of the rapidity shift.
   */
  enum BeamDirection_t{
    kpPb = 1,       ///< p-Pb
    kPbp = -1       ///< Pb-p
  };

  /**
   * @brief Dummy constructor
   */
  AliAnalysisTaskChargedParticlesRef();

  /**
   * @brief Main constructor
   * @param[in] name Name of the task
   */
  AliAnalysisTaskChargedParticlesRef(const char *name);

  /**
   * @brief Destuctor
   */
  virtual ~AliAnalysisTaskChargedParticlesRef();

  /**
   * @brief Enable Sumw2 when creating the histograms.
   *
   * Attention: Enabling Sumw2 will increase memory consumption
   * significantly. Option should only be used in case histograms
   * are filled with a weight.
   * @param[in] doEnable If true Sumw2 is enabled for all histograms
   */
  void EnableSumw2(Bool_t doEnable) { fEnableSumw2 = doEnable; }

  /**
   * @brief Set the track selection.
   *
   * The track cuts are generated according to predefined
   * recipes in the AliEmcalAnalysisFactory track cuts factory
   * @param[in] cutname Name of the track cuts
   * @param[in] isAOD check whether we run on ESDs or AODs
   */
  void InitializeTrackCuts(TString cutname, bool isAOD);

  /**
   * @brief Define rapidity shift applied to convert \f$ \eta_{lab} \f$
   * to \f$ \eta_{cms} \f$
   * @param[in] yshift Rapidity shift
   */
  void SetRapidityShift(Double_t yshift) { fYshift = yshift; }

  /**
   * @brief Specify directions of the incoming beams (p-Pb or Pb-p)
   * @param[in] beamdir Directions of the incoming beams
   */
  void SetBeamDirection(BeamDirection_t beamdir) { fEtaSign = static_cast<Double_t>(beamdir); }

  /**
   * @brief Enable PID-related plots (THnSparses, might be big)
   * @param[in] plotPID If true then PID-related plots are enabled
   */
  void SetPlotPID(Bool_t plotPID) { fStudyPID = plotPID; }

  /**
   * @brief Switch on-/off- study of exotic triggers (default: off)
   * @param[in] doStudy If true study of exotic triggers is enabled
   */
  void SetStudyExoticTriggers(Bool_t doStudy) { fStudyExoticTriggers = doStudy; SetCaloTriggerPatchInfoName("EmcalTriggers"); }

  /**
   * @brief Set the virtual track selection.
   *
   * Assumes that the track selection object is already fully configured
   * @param[in] sel Virtual track selection
   */
  void SetEMCALTrackSelection(AliEmcalTrackSelection *sel) { fTrackCuts = sel; }

  /**
   * @brief Define cut on \f$ \eta_{lab} \f$ used to select tracks
   * @param[in] etamin Minimum allowed \f$ \eta \f$
   * @param[in] etamax Maximum allowed \f$ \eta \f$
   */
  void SetEtaLabCut(double etamin, double etamax) { fEtaLabCut.SetLimits(etamin, etamax); }

  /**
   * @brief Define cut on \f$ \eta_{cms} \f$ used to select tracks
   * @param etamin[in] Minimum allowed \f$ \eta \f$
   * @param etamax[in] Maximum allowed \f$ \eta \f$
   */
  void SetEtaCMSCut(double etamin, double etamax) { fEtaCmsCut.SetLimits(etamin, etamax); }

  /**
   * @brief Define cut on \f$ \phi \f$ used to select tracks
   * @param[in] phimin Minimum allowed \f$ \phi \f$
   * @param[in] phimax Maximum allowed \f$ \phi \f$
   */
  void SetTrackPhiCut(double phimin, double phimax) { fPhiCut.SetLimits(phimin, phimax); }

  /**
   * @brief Set minimum \f$ p_{t}\f$ used to select track candidate
   * @param[in] minpt Minimum \f$ p_{t}\f$
   */
  void SetMinPtTracks(Double_t minpt) { fMinPt = minpt; }

  /**
   * @brief Add histograms for tracks pointing to EMCAL supermodules
   * @param[in] doStudy If true histograms are added
   */
  void SetStudyEMCALgeo(Bool_t doStudy) { fStudyEMCALgeo = doStudy; }

  /**
   * @brief Require bunch crossing information of track obtained from TOF (if available)
   * matches the bunch crossing ID of the event
   *
   * @param[in] doRequire If true the track is only selected if it is from the right bunch crossing according to TOF information
   */
  void SetRequireTOFBunchCrossing(Bool_t doRequire) { fRequireTOFBunchCrossing = doRequire; }

  /**
   * @brief Pre-configure task so that it can be used in subwagons
   * @param[in] suffix Suffix of the subwagon
   * @return Pre-configured task
   */
  static AliAnalysisTaskChargedParticlesRef *AddTaskChargedParticlesRef(const TString &suffix);

  /**
   * @brief Fully configure task, not intended for subwagons. Using a default cut configuration.
   * @param[in] cutname Name of the track cuts
   * @return Fully configured task
   */
  static AliAnalysisTaskChargedParticlesRef *AddTaskChargedParticlesRefDefault(const TString &cutname = "standard");

protected:

  /**
   * @brief Create the output histograms
   */
  virtual void CreateUserHistos();

  /**
   * @brief Creating user objects.
   *
   * In case no track selection is applied from outside,
   * the default track selection (RAA standard cuts is
   * created in this class)
   */
  virtual void CreateUserObjects();

  /**
   * @brief Simple task testing particle spectra in triggered events
   * - Select event using AliAnalysisUtil (in IsEventSelected)
   * - Assigning trigger type (Request INT7, EJ*, EG*)
   * - Loop over tracks, select particles
   * - Fill distributions
   * @param option Not used
   */
  virtual bool Run();

  /**
   * @brief Implementation of framework function UserFillHistosBeforeEventSelection.
   *
   * Filling event counter histograms before vertex selection
   */
  virtual void UserFillHistosBeforeEventSelection();

  /**
   * @brief Implementation of framework function UserFillHistosAfterEventSelection.
   *
   * Filling event counter histograms after vertex selection
   */
  virtual void UserFillHistosAfterEventSelection();

  /**
   * @brief Fill track (kinematic) histograms
   * @param[in] eventclass Trigger class fired
   * @param[in] histtag Tag in the histogram
   * @param[in] posCharge True if charge is positive, false if charge is negative
   * @param[in] pt track \f$ p_{t} \f$
   * @param[in] etalab Track \f$ \eta \f$ in lab frame
   * @param[in] etacent Track \f$ \eta \f$ in cms frame
   * @param[in] phi Track \f$ \eta \f$ in lab frame
   * @param[in] inEmcal Track in EMCAL \f$ \phi \f$ acceptance
   */
  void FillTrackHistos(const TString &eventclass, const TString &histtag, Bool_t posCharge, Double_t pt, Double_t eta, Double_t etacent, Double_t phi, Bool_t inEmcal);

  /**
   * @brief Fill PID-related histograms
   * @param[in] eventclass Trigger class of the event
   * @param[in] track Track containing PID information
   */
  void FillPIDHistos(const TString &eventclass, const AliVTrack &track);

  /**
   * Converted from AliESDtrack::GetTOFBunchCrossing, allowing for different values of bunch spacing
   * @param trk Track to check
   * @param b
   * @param spacing Bunch spacing (default: 25 ns - corresponding to pp)
   * @param pidTPConly
   * @return
   */
  Int_t GetTOFBunchCrossing(const AliVTrack *trk, Double_t b = 0, Double_t spacing = 25, Bool_t pidTPConly = kTRUE) const;

  /**
   * @brief Check if the trigger is fired by an exotic cluster.
   *
   * An event is considered as fired by exotics if there are only patches with exotic clusters
   * above threshold. This means every patch not overlapping with an exotic cluster is a valid trigger
   * patch. For good events one should therefor select events which have *false* as return value.
   *
   * @param[in] trg Trigger class to be checked
   * @return True if the event is purely fired by an exotic cluster, false otherwise
   */
  bool IsExoticsTrigger(const TString &trg);

  AliEmcalTrackSelection *fTrackCuts;                ///< Standard track selection

  Double_t                        fYshift;                    ///< Rapidity shift
  Double_t                        fEtaSign;                   ///< Sign of the eta distribution (swaps when beam directions swap): p-Pb: +1, Pb-p: -1

  Double_t                        fMinPt;                     ///< Minimum \f$ p_{t}\f$-used to select tracks
  AliCutValueRange<double>        fEtaLabCut;                 ///< Cut applied in \f$\eta_{lab}\f$ frame
  AliCutValueRange<double>        fEtaCmsCut;                 ///< Cut applied in \f$\eta_{centre-of-mass}\f$ frame
  AliCutValueRange<double>        fPhiCut;                    ///< Cut applied in \f$\phi_{track}\f$

  Bool_t                          fStudyPID;                  ///< Use kinematics correlation histograms
  Bool_t                          fStudyEMCALgeo;             ///< Add histograms for tracks pointing to the EMCAL acceptance
  Bool_t                          fEnableSumw2;               ///< Enable sumw2 during histogram creation
  Bool_t                          fRequireTOFBunchCrossing;   ///< Require that the bunch crossing ID determined by TOF matches the bunch crossing ID of the event
  Bool_t                          fStudyExoticTriggers;       ///< Switch on whether to study exotic triggers or not

private:

  /**
   * @class PtBinning
   * @brief Definition of the \f$ p_{t} \f$-binning used in the histograms
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

  AliAnalysisTaskChargedParticlesRef(const AliAnalysisTaskChargedParticlesRef &);
  AliAnalysisTaskChargedParticlesRef &operator=(const AliAnalysisTaskChargedParticlesRef &);

  ClassDef(AliAnalysisTaskChargedParticlesRef, 1);
};

} /* namespace EMCALJetTasks */

} /* namespace PWGJE */ 

#endif /* ALIANALYSISTASKCHARGEDPARTICLESREF_H */
