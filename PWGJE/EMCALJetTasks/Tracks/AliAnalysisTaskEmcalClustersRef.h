/************************************************************************************
 * Copyright (C) 2017, Copyright Holders of the ALICE Collaboration                 *
 * All rights reserved.                                                             *
 *                                                                                  *
 * Redistribution and use in source and binary forms, with or without               *
 * modification, are permitted provided that the following conditions are met:      *
 *     * Redistributions of source code must retain the above copyright             *
 *       notice, this list of conditions and the following disclaimer.              *
 *     * Redistributions in binary form must reproduce the above copyright          *
 *       notice, this list of conditions and the following disclaimer in the        *
 *       documentation and/or other materials provided with the distribution.       *
 *     * Neither the name of the <organization> nor the                             *
 *       names of its contributors may be used to endorse or promote products       *
 *       derived from this software without specific prior written permission.      *
 *                                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND  *
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED    *
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE           *
 * DISCLAIMED. IN NO EVENT SHALL ALICE COLLABORATION BE LIABLE FOR ANY              *
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES       *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;     *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND      *
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT       *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS    *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                     *
 ************************************************************************************/
#ifndef ALIANALYSISTASKEMCALCLUSTERSREF_H
#define ALIANALYSISTASKEMCALCLUSTERSREF_H

#include "AliAnalysisEmcalTriggerSelectionHelper.h"
#include "AliAnalysisTaskEmcalTriggerBase.h"
#include "AliCutValueRange.h"
#include <TCustomBinning.h>
#include <TString.h>
#include <vector>

class TClonesArray;

namespace PWGJE {

namespace EMCALJetTasks {

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
class AliAnalysisTaskEmcalClustersRef : public AliAnalysisTaskEmcalTriggerBase, public AliAnalysisEmcalTriggerSelectionHelperImpl {
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
   * @brief Dummy (I/O) constructor
   */
  AliAnalysisTaskEmcalClustersRef();

  /**
   * @brief Named constructor
   * @param[in] name Name of the task
   */
  AliAnalysisTaskEmcalClustersRef(const char *name);

  /**
   * @brief Destructor
   */
  virtual ~AliAnalysisTaskEmcalClustersRef();

  /**
   * @brief Enable Sumw2 when creating the histograms. 
   * 
   * Attention: Enabling Sumw2
   * will increase memory consumption significantly. Option should only be
   * used in case histograms are filled with a weight.
   * 
   * @param[in] doEnable If true Sumw2 is enabled for all histograms
   */
  void EnableSumw2(Bool_t doEnable) { fEnableSumw2 = doEnable; }

  /**
   * @brief Define type of energy used in the monitoring histograms
   * @param[in] edef type of energy definition
   */
  void SetEnergyDefinition(EnergyDefinition_t edef) { fEnergyDefinition = edef; }

  /**
   * @brief Add dimensions for Eta-phi to the cluster THnSparse
   * @param[in] doMonitor If true dimensions for eta-phi are added
   */
  void SetMonitorEtaPhi(bool doMonitor) { fMonitorEtaPhi = doMonitor; }

  /**
   * @brief Define cut on the time of the leading cell in the cluster
   * @param[in] mintime Minimum selected time for cluster
   * @param[in] maxtime Maximum selected time for cluster
   */
  void SetClusterTimeRange(double mintime, double maxtime) { fClusterTimeRange.SetLimits(mintime, maxtime); }

  /**
   * @brief Define centrality range used to select
   * @param[in] min Min. allowed centrality
   * @param[in] max Max. allowed centrality
   */
  void SetCentralityRange(double min, double max) { fCentralityRange.SetLimits(min, max); fRequestCentrality = true; }

  /**
   * @brief Set the centrality estimator used to select centrality ranges. 
   * 
   * By default V0M will be used
   * 
   * @param[in] centest Name of the centrality estimator
   */
  void SetUserCentralityEstimator(TString centest) { fCentralityEstimator = centest; }

  /**
   * @brief Select events only from certain bunch crossings
   * @param[in] bunchCrossingIndex index of then bunch crossing selected
   */
  void SetBunchCrossingIndex(Int_t bunchCrossingIndex) { fBunchCrossingIndex = bunchCrossingIndex; };

  /**
   * @brief Switch for filling multiplicity correlation histograms
   * @param[in] doFill if true multiplicity histograms are filled
   */
  void SetFillMultiplicityHistograms(Bool_t doFill) { fDoFillMultiplicityHistograms = doFill; }

  /**
   * @brief Switch on histograms for fired clusters
   * @param[in] doUse If true histograms are switched on
   */
  void SetUsedFiredClusters(Bool_t doUse) { fUseFiredTriggers = doUse; }

  /**
   * @brief Set filling trigger cluster dimension
   * @param[in] doFill If true trigger cluster axis is filled (only relevant in pp for livetime correction) 
   */
  void SetFillTriggerClusters(Bool_t doFill) { fFillTriggerClusters = doFill; }

  /**
  * @brief Switch on/off exclusive triggers
  * 
  * Exclusive triggers do not contain lower threshold triggers. Switching
  * off is for memory considerations
  * 
  * @param doUse If true exclusive triggers are used 
  */
  void SetUseExclusiveTriggers(Bool_t doUse) { fUseExclusiveTriggers = doUse; }

  /**
   * @brief Add trigger for which overlap is required
   * 
   * Event is only selected if the trigger fires as well
   * 
   * @param[in] trigger Trigger for which overlap is required 
   */
  void AddRequiredTriggerOverlap(const char *trigger);

  /**
   * @brief Add trigger for which overlap is excluded
   * 
   * Event is only selected if the trigger does not fire the event in addition
   * (histos of this trigger will of course be empty)
   * 
   * @param[in] trigger Trigger for which overlap is excluded 
   */
  void AddExcludedTriggerOverlap(const char *trigger);

  /**
   * @brief Preconfigure task so that it can be used in subwagons
   * @param[in] nClusters Name of the input cluster container
   * @param[in] suffix Suffix of the subwagon
   * @return Preconfigure task
   */
  static AliAnalysisTaskEmcalClustersRef *AddTaskEmcalClustersRef(const TString &nClusters = "usedefault", const TString &suffix = "");

  /**
   * @brief Preconfigure task and add it to the analysis manager
   * @param[in] nClusters Name of the input cluster container
   * @return Fully configured task
   */
  static AliAnalysisTaskEmcalClustersRef *AddTaskEmcalClustersRefDefault(const TString &nClusters = "usedefault");

protected:

  /**
   * @brief Task has no user-defined objects
   */
  virtual void CreateUserObjects() {}

  /**
   * @brief Creating histograms for the distributions monitored by the task
   */
  virtual void CreateUserHistos();

  /**
   * @brief User event selection: Select event in maching centrality range (if requested)
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
   * @brief Get the boundaries of the trigger patch
   * @param[in] o Patch o check
   * @param[out] boundaries Patch boundaries [etamin, phimin, etamax, phimax]
   */
  void GetPatchBoundaries(AliEMCALTriggerPatchInfo &o, Double_t *boundaries) const;

  void FillClusterHistograms(const TString &triggerclass, double energy, double eta, double phi, double clustertime, int ncell, int trgcluster, const TList *triggerpatches, int energycomp);

  /**
   * @brief Check whether cluster is inside a trigger patch which has fired the trigger
   * @param[in] etaclust \f$ \eta \f$ of the cluster at center
   * @param[in] phiclust \f$ \phi \f$ of the cluster at center
   * @param[in] fTriggerPatches List of trigger patches which have fired the trigger
   * @return[in] List with all patches with geometric overlap with the cluster
   */
  std::vector<AliEMCALTriggerPatchInfo *> CorrelateToTrigger(Double_t etaclust, Double_t phiclust, const TList &triggerpatches) const;


  AliCutValueRange<double>            fCentralityRange;           ///< Selected centrality range
  Bool_t                              fRequestCentrality;         ///< Switch on request for centrality range
  Double_t                            fEventCentrality;           //!<! Current event centrality
  TString                             fCentralityEstimator;       ///< Centrality estimator (default: V0M for PbPb)
  Char_t                              fBunchCrossingIndex;        ///< Bunch Crossing index

  EnergyDefinition_t                  fEnergyDefinition;          ///< Energy definition used for a given cluster
  Bool_t                              fEnableSumw2;               ///< Enable sumw2 when creating histograms
  Bool_t                              fDoFillMultiplicityHistograms;    ///< Swich for multiplcity histograms
  Bool_t                              fUseFiredTriggers;          ///< Study clusters connected with patches
  Bool_t                              fUseExclusiveTriggers;      ///< Include exclusive triggers (without lower threshold triggers)
  Bool_t                              fFillTriggerClusters;       ///< Fill trigger cluster histograms
  Bool_t                              fMonitorEtaPhi;             ///< Add dimensions for eta-phi in the THnSparses
  AliCutValueRange<double>            fClusterTimeRange;          ///< Selected range on cluster time
  std::vector<TriggerCluster_t>       fTriggerClusters;           //!<! Detected trigger clusters for event
  TObjArray                           fRequiredOverlaps;          ///< Add option to require overlap with certain triggers
  TObjArray                           fExcludedOverlaps;          ///< Add option to exclude overlap with certain triggers

private:

  class EnergyBinning : public TCustomBinning {
  public:
    EnergyBinning();
    virtual ~EnergyBinning() {}
  };

  int CountTracklets(double etamin, double etamax, double phimin, double phimax);
  int CountEmcalClusters(double ecut);
  int GetEMCALCellOccupancy(double ecut);

  AliAnalysisTaskEmcalClustersRef(const AliAnalysisTaskEmcalClustersRef &);
  AliAnalysisTaskEmcalClustersRef &operator=(const AliAnalysisTaskEmcalClustersRef &);

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskEmcalClustersRef, 1);
  /// \endcond
};

} /* namespace EMCALJetTasks */

} /* namespace PWGJE */

#endif /* ALIANALYSISTASKEMCALCLUSTERSREF_H */
