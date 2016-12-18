#ifndef ALIANALYSISTASKEMCALPATCHESREF_H
#define ALIANALYSISTASKEMCALPATCHESREF_H
/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <string>
#include <vector>
#include "AliAnalysisTaskEmcalTriggerBase.h"
#include "AliCutValueRange.h"
#include <TCustomBinning.h>
#include <TString.h>

class AliOADBContainer;
class AliEMCALTriggerPatchInfo;
class THistManager;
class TObjArray;

namespace EMCalTriggerPtAnalysis {

class AliEmcalTriggerOfflineSelection;

/**
 * @class AliAnalysisTaskEmcalPatchesRef
 * @brief Analysis of trigger patches in min. bias and triggered events
 * @author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 * @ingroup PWGJETASKS
 * @since Aug 31. 2015
 */
class AliAnalysisTaskEmcalPatchesRef : public AliAnalysisTaskEmcalTriggerBase {
public:
  /**
   * Dummy (I/O) onstructor
   */
  AliAnalysisTaskEmcalPatchesRef();

  /**
   * Named constructor
   * @param[in] name Name of the task
   */
  AliAnalysisTaskEmcalPatchesRef(const char *name);

  /**
   * Destructor
   */
  virtual ~AliAnalysisTaskEmcalPatchesRef(){}

  /**
   * Enable Sumw2 when creating the histograms. Attention: Enabling Sumw2
   * will increase memory consumption significantly. Option should only be
   * used in case histograms are filled with a weight.
   * @param[in] doEnable If true Sumw2 is enabled for all histograms
   */
  void EnableSumw2(Bool_t doEnable) { fEnableSumw2 = doEnable; }

  /**
   * Set centrality selection.
   * Note: Needs multiplicity task to run in front
   * @param[in] min Min. value of the centrality interval
   * @param[in] max Max. value of the centrality interval
   */
  void SetCentralityRange(double min, double max) { fCentralityRange.SetLimits(min,max); fRequestCentrality = true; }

protected:

  /**
   * Task does not need user objects, so the function is not implemented
   */
  virtual void CreateUserObjects() {}

  /**
   * Creating output histograms:
   * - Patch (calibrated) energy spectrum - separated by patch type - for different trigger classes
   * - Patch eta-phi map - separated by patch type - for different trigger classes and different min. energies
   */
  virtual void CreateUserHistos();

  /**
   * Handle centrality selection (if requested, and multiplicity information available)
   * @return True if event is in proper centrality class
   */
  virtual bool IsUserEventSelected();

  /**
   * Run event loop: Loop over all patches, assign them to trigger classes, and
   * fill relevant histograms
   * @param Not used
   */
  virtual bool Run();

  /**
   * Fill event-based histograms. Monitored are
   * - Number of events
   * - Centrality percentile (if available)
   * - z-position of the primary vertex
   * In case a downscaling correction is avaiable it is applied to all
   * histograms as a weight.
   */
  virtual void UserFillHistosAfterEventSelection();

  /**
   * Helper function reading the patch boundaries in \f$\eta\f$ and \f$\phi\f$
   * from a trigger patch
   * @param[in] patch Patch to check
   * @param[out] boundaries Patch boundaries
   */
  void GetPatchBoundaries(const AliEMCALTriggerPatchInfo *patch, Double_t *boundaries) const;

  /**
   * Helper function, selecting offline simple single shower (gamma or L0) patches
   * @param[in] patch Patch to check
   * @return True if the patch is selected
   */
  bool SelectSingleShowerPatch(const AliEMCALTriggerPatchInfo *patch) const;

  /**
   * Helper function, selecting offline simple jet patches
   * @param[in] patch Patch to check
   * @return True if the patch is selected
   */
  bool SelectJetPatch(const AliEMCALTriggerPatchInfo  *patch) const;

  /**
   * Filling patch related histogram. In case a downscaling correction is
   * available it is applied to the histograms as weight
   * @param[in] triggerclass Name of the trigger class firing the event
   * @param[in] patchname Name of the patchtype
   * @param[in] energy Calibrated energy of the patch
   * @param[in] eta Patch eta at the geometrical center
   * @param[in] phi Patch phi at the geometrical center
   */
  void FillPatchHistograms(TString triggerclass, TString patchname, double energy, double transverseenergy, double eta, double phi, int col, int row);

  AliCutValueRange<double>            fCentralityRange;           ///< Range of accepted event centralities
  Bool_t                              fEnableSumw2;               ///< Enable sumw2 during histogram creation
  Bool_t                              fRequestCentrality;         ///< Switch for request of centrality selection
  Double_t                            fEventCentrality;           //!<1! Event centrality

private:

  class EnergyBinning : public TCustomBinning {
  public:
    /**
     * Constructor
     */
    EnergyBinning();

    /**
     * Destructor
     */
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
