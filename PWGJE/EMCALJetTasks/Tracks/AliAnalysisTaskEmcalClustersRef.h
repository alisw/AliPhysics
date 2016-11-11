#ifndef ALIANALYSISTASKEMCALCLUSTERSREF_H
#define ALIANALYSISTASKEMCALCLUSTERSREF_H
/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliAnalysisTaskEmcalTriggerBase.h"
#include "AliCutValueRange.h"
#include <TCustomBinning.h>
#include <TString.h>

class TClonesArray;

namespace EMCalTriggerPtAnalysis {

/**
 * @class AliAnalysisTaskEmcalClustersRef
 * @brief Simple monitoring task for cluster-related quantities in EMCAL-triggered events
 * @ingroup PWGJETASKS
 * @author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 * @since Aug 31, 2015
 *
 * Simple monitoring class monitoring
 * - energy / transverse energy spectra
 * - differential energy/transverse energy as function of
 *   - \f$ \eta \f$
 *   - tracking sector
 *   - supermodule
 * Clusters are selected from a provided cluster container (argument "usedefault" should
 * be used). The energy can be the uncorrected, corrected energy for non-linearity or the
 * energy corrected for hadronic contribution.
 */
class AliAnalysisTaskEmcalClustersRef : public AliAnalysisTaskEmcalTriggerBase {
public:

  /**
   * @enum EnergyDefinition_t
   * @brief Type of the energy used in the monitoring histograms
   */
  enum EnergyDefinition_t {
    kDefaultEnergy,       ///< Uncorrected energy measurement
    kNonLinCorrEnergy,    ///< Energy corrected for non-linearity
    kHadCorrEnergy        ///< Energy corrected for the hadronic contribution
  };

  /**
   * Dummy (I/O) constructor
   */
  AliAnalysisTaskEmcalClustersRef();

  /**
   * Named constructor
   * @param[in] name Name of the task
   */
  AliAnalysisTaskEmcalClustersRef(const char *name);

  /**
   * Destructor
   */
  virtual ~AliAnalysisTaskEmcalClustersRef();

  /**
   * Define type of energy used in the monitoring histograms
   * @param[in] edef type of energy definition
   */
  void SetEnergyDefinition(EnergyDefinition_t edef) { fEnergyDefinition = edef; }

  /**
   * Define name of the cluster container used to read EMCAL cluster information
   * from
   * @param[in] clustercontname Name of the cluster container
   */
  void SetClusterContainer(TString clustercontname) { fNameClusterContainer = clustercontname; }

  /**
   * Define centrality range used to select
   * @param[in] min Min. allowed centrality
   * @param[in] max Max. allowed centrality
   */
  void SetCentralityRange(double min, double max) { fCentralityRange.SetLimits(min, max); fRequestCentrality = true; }

  /**
   * Preconfigure task so that it can be used in subwagons
   * @param[in] nClusters Name of the input cluster container
   * @param[in] suffix Suffix of the subwagon
   * @return Preconfigure task
   */
  static AliAnalysisTaskEmcalClustersRef *AddTaskEmcalClustersRef(const TString &nClusters = "usedefault", const TString &suffix = "");

  /**
   *
   * @param[in] nClusters Name of the input cluster container
   * @return Fully configured task
   */
  static AliAnalysisTaskEmcalClustersRef *AddTaskEmcalClustersRefDefault(const TString &nClusters = "usedefault");

protected:

  /**
   * Task has no user-defined objects
   */
  virtual void CreateUserObjects() {}

  /**
   * Creating histograms for the distributions monitored by the task
   */
  virtual void CreateUserHistos();

  /**
   * User event selection: Select event in maching centrality range (if requested)
   * @return True if the event is in the matching centrality range
   */
  virtual bool IsUserEventSelected();

  /**
   * Main function: If the event was selected previously loop
   * over all clusters in the cluster container and fill distributions
   * for all different selected triggers
   * @return Always true
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
   * Get the boundaries of the trigger patch
   * @param[in] o Patch o check
   * @param[out] boundaries Patch boundaries [etamin, phimin, etamax, phimax]
   */
  void GetPatchBoundaries(TObject *o, Double_t *boundaries) const;

  /**
   * Select patch as offline simple (FEE) patch
   * @param[in] o Patch to check
   * @return True if the patch is an offline simple patch
   */
  bool IsOfflineSimplePatch(TObject *o) const;

  /**
   * Selection of patches in the DCAL acceptance
   * @param[in] o Patch to check
   * @return True if the patch is in the DCAL acceptance
   */
  bool SelectDCALPatch(TObject *o) const;

  /**
   * Selection of patches as single shower (gamma/Level0) patches
   * @param[in] o Patch to check
   * @return True if the patch is a single shower patch
   */
  bool SelectSingleShowerPatch(TObject *o) const;

  /**
   * Selection of jet patches
   * @param[in] o Patch to check
   * @return True if the patch is a jet patch
   */
  bool SelectJetPatch(TObject *o) const;

  /**
   * Get the energy of a trigger patch
   * @param[in] o Patch to check
   * @return Energy of the trigger patch
   */
  double GetPatchEnergy(TObject *o) const;

  void FillClusterHistograms(const TString &triggerclass, double energy, double transversenergy, double eta, double phi, TList *triggerpatches);

  /**
   * Find all patches in an event which could have fired the trigger
   * Attention: This task groups into single shower triggers (L0, EG1, EG2) and jet triggers (EJ1 and EJ2).
   * Per convention the low threshold patch is selected. No energy cut should be applied in the trigger maker
   * @param[in] triggerclass EMCAL trigger class firing
   * @param[in] fTriggerPatches Trigger patches found in the event
   * @return List of patches which could have fired the trigger
   */
  void FindPatchesForTrigger(TString triggerclass, const TClonesArray * triggerpatches, TList &foundpatches) const;

  /**
   * Check whether cluster is inside a trigger patch which has fired the trigger
   * @param[in] etaclust \f$ \eta \f$ of the cluster at center
   * @param[in] phiclust \f$ \phi \f$ of the cluster at center
   * @param[in] fTriggerPatches List of trigger patches which have fired the trigger
   * @return[in] True if the cluster can be correlated to a triggerpatch fired the trigger, false otherwise
   */
  Bool_t CorrelateToTrigger(Double_t etaclust, Double_t phiclust, TList *triggerpatches) const;

  TString                             fNameClusterContainer;      ///< Name of the cluster container in the event

  AliCutValueRange<double>            fCentralityRange;           ///< Selected centrality range
  Bool_t                              fRequestCentrality;         ///< Switch on request for centrality range
  Double_t                            fEventCentrality;           //!<! Current event centrality

  EnergyDefinition_t                  fEnergyDefinition;          ///< Energy definition used for a given cluster

private:

  class EnergyBinning : public TCustomBinning {
  public:
    EnergyBinning();
    virtual ~EnergyBinning() {}
  };

  AliAnalysisTaskEmcalClustersRef(const AliAnalysisTaskEmcalClustersRef &);
  AliAnalysisTaskEmcalClustersRef &operator=(const AliAnalysisTaskEmcalClustersRef &);

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskEmcalClustersRef, 1);
  /// \endcond
};

} /* namespace EMCalTriggerPtAnalysis */

#endif /* ALIANALYSISTASKEMCALCLUSTERSREF_H */
