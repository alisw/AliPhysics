#ifndef ALIANALYSISTASKEMCALCLUSTERMATCHED_H
#define ALIANALYSISTASKEMCALCLUSTERMATCHED_H
/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliAnalysisTaskEmcalTriggerBase.h"
#include "AliCutValueRange.h"
#include <TCustomBinning.h>
#include <TString.h>

class AliEmcalTrackSelection;

namespace PWGJE { 

namespace EMCALJetTasks {

/**
 * @class AliAnalysisTaskEmcalClusterMatched
 * @brief Simple task checking energy spectra of clusters matched to tracks
 * @ingroup PWGJETASKS
 * @author Markus Fasel <markus.fasel@cern.ch>, Oak Ridge National Laboratory
 * @since Dec 16, 2016
 *
 * Simple task comparing energy spectra of EMCAL clusters
 * - all
 * - matched to global tracks
 * - matched to TPC-only tracks
 */
class AliAnalysisTaskEmcalClusterMatched : public AliAnalysisTaskEmcalTriggerBase {
public:

  /**
   * ROOT I/O constructor
   */
  AliAnalysisTaskEmcalClusterMatched();

  /**
   * Named constructor - performs initial task configuration
   * @param[in] name
   */
  AliAnalysisTaskEmcalClusterMatched(const char *name);

  /**
   * Destructor, deleting track selections
   */
  virtual ~AliAnalysisTaskEmcalClusterMatched();

  /**
   * Set time range used to accept clusters
   * @param[in] mintime Minimum allowed cluster time
   * @param[in] maxtime Maximum allowed cluster time
   */
  void SetClusterTimeRange(double mintime, double maxtime) { fTimeCut.SetLimits(mintime, maxtime); }

  /**
   * Enable sumw2 option
   * @param[in] doenable If true sumw2 is enabled for all histograms
   */
  void EnableSumw2(bool doenable) { fEnableSumw2 = doenable; }

  /**
   * Create and configure task skeleton
   * @param[in] suffix Suffix in task and output container name
   * @return Pre-configured task
   */
  static AliAnalysisTaskEmcalClusterMatched *AddTaskEmcalClusterMatched(const char *suffix);

protected:

  /**
   * Initializes the two track selections (global and TPConly) and
   * sets the default track container
   */
  virtual void CreateUserObjects();

  /**
   * Create histograms for the cluster energy per trigger and supermodule
   */
  virtual void CreateUserHistos();

  /**
   * Perform analysis
   * @return always true
   */
  virtual bool Run();

  /**
   * Fill event counter histograms
   */
  virtual void UserFillHistosAfterEventSelection();

  /**
   * Initialize the selection of global and TPC-only tracks
   * @param[in] isAOD Switch between ESD and AOD tracks
   */
  void InitializeTrackSelections(bool isAOD);

private:

  /**
   * @class EnergyBinning
   * @brief Binning for cluster energy
   */
  class EnergyBinning : public TCustomBinning{
  public:

    /**
     * Constructor, define the binning
     */
    EnergyBinning();

    /**
     * Destructor
     */
    virtual ~EnergyBinning() {}
  };

  AliAnalysisTaskEmcalClusterMatched(const AliAnalysisTaskEmcalClusterMatched &);
  AliAnalysisTaskEmcalClusterMatched &operator=(const AliAnalysisTaskEmcalClusterMatched &);

  AliEmcalTrackSelection       *fTrackSelectionGlobal;                  ///< Global track cuts (strong condition)
  AliEmcalTrackSelection       *fTrackSelectionTPConly;                 ///< TPC-only track cuts (weak condition)
  AliCutValueRange<double>                  fTimeCut;                               ///< Cut on cluster time
  Bool_t                                    fEnableSumw2;                           ///< Enable Sumw2

  ClassDef(AliAnalysisTaskEmcalClusterMatched, 1)
};

} /* namespace EMCALJetTasks */

} /* namespace PWGJE */

#endif /* ALIANALYSISTASKEMCALCLUSTERMATCHED_H */
