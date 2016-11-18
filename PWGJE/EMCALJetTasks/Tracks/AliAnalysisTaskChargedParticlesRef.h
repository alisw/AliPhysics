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
   * Dummy constructor
   */
  AliAnalysisTaskChargedParticlesRef();

  /**
   * Main constructor
   * @param[in] name Name of the task
   */
  AliAnalysisTaskChargedParticlesRef(const char *name);

  /**
   * Destuctor
   */
  virtual ~AliAnalysisTaskChargedParticlesRef();

  /**
   * Enable Sumw2 when creating the histograms. Attention: Enabling Sumw2
   * will increase memory consumption significantly. Option should only be
   * used in case histograms are filled with a weight.
   * @param[in] doEnable If true Sumw2 is enabled for all histograms
   */
  void EnableSumw2(Bool_t doEnable) { fEnableSumw2 = doEnable; }

  /**
   * Set the track selection: The track cuts are generated according to predefined
   * recipes in the AliEmcalAnalysisFactory track cuts factory
   * @param[in] cutname Name of the track cuts
   * @param[in] isAOD check whether we run on ESDs or AODs
   */
  void InitializeTrackCuts(TString cutname, bool isAOD);

  /**
   * Define rapidity shift applied to convert \f$ \eta_{lab} \f$
   * to \f$ \eta_{cms} \f$
   * @param[in] yshift Rapidity shift
   */
  void SetRapidityShift(Double_t yshift) { fYshift = yshift; }

  /**
   * Specify directions of the incoming beams (p-Pb or Pb-p)
   * @param[in] beamdir Directions of the incoming beams
   */
  void SetBeamDirection(BeamDirection_t beamdir) { fEtaSign = static_cast<Double_t>(beamdir); }

  /**
   * Enable PID-related plots (THnSparses, might be big)
   * @param[in] plotPID If true then PID-related plots are enabled
   */
  void SetPlotPID(Bool_t plotPID) { fStudyPID = plotPID; }

  /**
   * Set the virtual track selection. Assumes that the track selection object
   * is already fully configured
   * @param[in] sel Virtual track selection
   */
  void SetEMCALTrackSelection(AliEmcalTrackSelection *sel) { fTrackCuts = sel; }

  /**
   * Define cut on \f$ \eta_{lab} \f$ used to select tracks
   * @param[in] etamin Minimum allowed \f$ \eta \f$
   * @param[in] etamax Maximum allowed \f$ \eta \f$
   */
  void SetEtaLabCut(double etamin, double etamax) { fEtaLabCut.SetLimits(etamin, etamax); }

  /**
   * Define cut on \f$ \eta_{cms} \f$ used to select tracks
   * @param etamin[in] Minimum allowed \f$ \eta \f$
   * @param etamax[in] Maximum allowed \f$ \eta \f$
   */
  void SetEtaCMSCut(double etamin, double etamax) { fEtaCmsCut.SetLimits(etamin, etamax); }

  /**
   * Define cut on \f$ \phi \f$ used to select tracks
   * @param[in] phimin Minimum allowed \f$ \phi \f$
   * @param[in] phimax Maximum allowed \f$ \phi \f$
   */
  void SetTrackPhiCut(double phimin, double phimax) { fPhiCut.SetLimits(phimin, phimax); }

  /**
   * Preconfigure task so that it can be used in subwagons
   * @param[in] suffix Suffix of the subwagon
   * @return Preconfigured task
   */
  static AliAnalysisTaskChargedParticlesRef *AddTaskChargedParticlesRef(const TString &suffix);

  /**
   * Fully configure task, not intended for subwagons. Using a default cut configuration.
   * @param[in] cutname Name of the track cuts
   * @return Fully configured task
   */
  static AliAnalysisTaskChargedParticlesRef *AddTaskChargedParticlesRefDefault(const TString &cutname = "standard");

protected:

  /**
   * Create the output histograms
   */
  virtual void CreateUserHistos();

  /**
   * Creating user objects: In case no track selection
   * is applied from outside, the default track selection
   * (RAA standard cuts is created in this class)
   */
  virtual void CreateUserObjects();

  /**
   * Simple task testing particle spectra in triggered events
   * - Select event using AliAnalysisUtil (in IsEventSelected)
   * - Assing trigger type (Request INT7, EJ*, EG*)
   * - Loop over tracks, select particles
   * - Fill distributions
   * @param option Not used
   */
  virtual bool Run();

  /**
   * Implementation of framework function UserFillHistosBeforeEventSelection.
   * Filling event counter histograms before vertex selection
   */
  virtual void UserFillHistosBeforeEventSelection();

  /**
   * Implementation of framework function UserFillHistosAfterEventSelection.
   * Filling event counter histograms after vertex selection
   */
  virtual void UserFillHistosAfterEventSelection();

  /**
   * Fill track (kinematic) histograms
   * @param[in] eventclass Trigger class fired
   * @param[in] pt track \f$ p_{t} \f$
   * @param[in] etalab Track \f$ \eta \f$ in lab frame
   * @param[in] etacent Track \f$ \eta \f$ in cms frame
   * @param[in] phi Track \f$ \eta \f$ in lab frame
   * @param[in] inEmcal Track in EMCAL \f$ \phi \f$ acceptance
   */
  void FillTrackHistos(const TString &eventclass, Double_t pt, Double_t eta, Double_t etacent, Double_t phi, Bool_t inEmcal);

  /**
   * Fill PID-related histograms
   * @param[in] eventclass Trigger class of the event
   * @param[in] track Track containing PID information
   */
  void FillPIDHistos(const TString &eventclass, const AliVTrack &track);

  AliEmcalTrackSelection          *fTrackCuts;                ///< Standard track selection

  Double_t                        fYshift;                    ///< Rapidity shift
  Double_t                        fEtaSign;                   ///< Sign of the eta distribution (swaps when beam directions swap): p-Pb: +1, Pb-p: -1

  AliCutValueRange<double>        fEtaLabCut;                 ///< Cut applied in \f$\eta_{lab}\f$ frame
  AliCutValueRange<double>        fEtaCmsCut;                 ///< Cut applied in \f$\eta_{centre-of-mass}\f$ frame
  AliCutValueRange<double>        fPhiCut;                    ///< Cut applied in \f$\phi_{track}\f$

  Bool_t                          fStudyPID;                  ///< Use kinematics correlation histograms
  Bool_t                          fEnableSumw2;               ///< Enable sumw2 during histogram creation

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

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskChargedParticlesRef, 1);
  /// \endcond
};

} /* namespace EMCalTriggerPtAnalysis */

#endif /* ALIANALYSISTASKCHARGEDPARTICLESREF_H */
