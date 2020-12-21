/**************************************************************************************
 * Copyright (C) 2017, Copyright Holders of the ALICE Collaboration                   *
 * All rights reserved.                                                               *
 *                                                                                    *
 * Redistribution and use in source and binary forms, with or without                 *
 * modification, are permitted provided that the following conditions are met:        *
 *     * Redistributions of source code must retain the above copyright               *
 *       notice, this list of conditions and the following disclaimer.                *
 *     * Redistributions in binary form must reproduce the above copyright            *
 *       notice, this list of conditions and the following disclaimer in the          *
 *       documentation and/or other materials provided with the distribution.         *
 *     * Neither the name of the <organization> nor the                               *
 *       names of its contributors may be used to endorse or promote products         *
 *       derived from this software without specific prior written permission.        *
 *                                                                                    *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND    *
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED      *
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE             *
 * DISCLAIMED. IN NO EVENT SHALL ALICE COLLABORATION BE LIABLE FOR ANY                *
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES         *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;       *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND        *
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT         *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS      *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                       *
 **************************************************************************************/
#ifndef ALIEMCALCELLMONITOR_H_
#define ALIEMCALCELLMONITOR_H_

#include "AliAnalysisTaskSE.h"
#include <TCustomBinning.h>
#include <TString.h>
#include <vector>

class TArrayD;
class THistManager;
class AliEMCALGeometry;

/**
 * @namespace PWG
 * @brief Namespace for PWG framework classes
 */
namespace PWG {

/**
 * @namespace EMCAL
 * @ingroup EMCALFW
 * @brief Namespace for EMCAL framework classes and task
 */
namespace EMCAL {

/**
 * @class AliEmcalCellMonitorTask
 * @brief Simple monitoring task for cell related quantities
 * @author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 * @since May 19, 2016
 * @ingroup EMCALFWTASKS
 *
 * This task monitors simple cell quantities like
 * - Amplitude distribution
 * - Time distribution
 * - Count rate in col-row space
 * - Integrated amplitude in col-row space
 * In addition, in case the task has access to reconstructed clusters,
 * it plots
 * - Frequency of a cell contributing to clusters
 * - Summed cell amplitude fraction in a cluster
 * In order to switch on cluster handling a cluster input container needs
 * to be provided:
 *
 * ~~~{.cxx}
 * ...
 * task->SetNameClusterContainer("caloClusters");
 * ~~~
 *
 * Note: When working with the default branches the names are:
 * - caloClusters for ESDs
 * - CaloClusters for AODs
 *
 * The task should run separately for different trigger
 * classes. This is handled in the event selection provided
 * from outside:
 *
 * ~~~{.cxx}
 * ...
 * // Selecting EMCAL gamma triggered events, high threshold
 * task->SetRequestTrigger(AliVEvent::kEGA, "EG1");
 * ~~~
 *
 * It can be added to the train using the add macro
 * ~~~
 * $ALICE_PHYSICS/PWG/EMCAL/AddEmcalCellMonitorTask.C
 * ~~~
 */
class AliEmcalCellMonitorTask : public AliAnalysisTaskSE {
public:

  /**
   * Dummy constructor, for ROOT I/O only
   */
  AliEmcalCellMonitorTask();

  /**
   * Default constructor, setting also the name and defining the output type
   * @param name Name of the task
   */
  AliEmcalCellMonitorTask(const char *name);

  /**
   * Destructor
   */
  virtual ~AliEmcalCellMonitorTask();

  /**
   * Set the minimum amplitude requested for a cell in order to use the
   * cell time in the cell time plot
   * @param minAmp Min. cell amplitude
   */
  void SetMinCellAmplitudeForCellTime(Double_t minAmp) { fMinCellAmplitude = minAmp; }

  /**
   * Define trigger selection. It can consist of trigger bits and a trigger string
   * @param triggerbits Trigger bit selection
   * @param triggerstring Trigger string (i.e. EG1, EG2, DG1, ...)
   */
  void SetRequestTrigger(ULong_t triggerbits, TString triggerstring = "") {
    fRequestTrigger = triggerbits;
    fTriggerString = triggerstring;
  }

  /**
   * Set number of Cells in order to match run2 geometry
   */
  void SetRun2() { fNumberOfCells = 17664; }

  /**
   * Mark cell as bad (i.e) exclude from running
   * @param cellId Cell to be masked
   */
  void SetBadCell(Int_t cellId);

  /**
   * Set the name of the cluster container. If provided the following
   * additional histograms will be provided:
   * - Occurrency of a cell inside clusters:  cellClusterOccurrency
   * - Summed energy distribution of a cell in clusters: cellSumAmplitudeCluster
   * In case of the summed energy the value is calculated summing up the
   * energy of all occurrences within one cluster
   * @param[in] nameclusters Name of the cluster container
   */
  void SetNameClusterContainer(const TString &nameclusters) { fNameClusters = nameclusters; }

  /**
   * Read bad channels from OADB container and set the cell with the ID
   * as bad.
   * @param[in] containername Name of the OADB container.
   */
  void InitBadChannelsFromContainer(const TString &containername) { fBadChannelContainer = containername; }

protected:

  /**
   * @class AliEmcalCellMonitorAmplitudeBinning
   * @brief Defining binning in amplitude direction
   *
   * Create binning for the amplitude:
   * 0 - 2 GeV: 100 MeV bins
   * 2 - 5 GeV: 200 MeV bins
   * 5 - 10 GeV: 500 MeV bins
   * 10 - 20 GeV: 1 GeV bins
   * 20 - 50 GeV: 2 GeV bins
   * 50 - 100 GeV: 5 GeV bins
   * 100 - 200 GeV: 10 GeV bins
   */
  class AliEmcalCellMonitorAmplitudeBinning : public TCustomBinning {
  public:

    /**
     * Constructor, defining minimum and ranges
     */
    AliEmcalCellMonitorAmplitudeBinning();

    /**
     * Destructor
     */
    virtual ~AliEmcalCellMonitorAmplitudeBinning() {}
  };

  /**
   * Prepare histogram manager for later initialization.
   * As the histograms depend on the number of cells which is
   * available only for after the first event is initialized,
   * the histograms are not created here but in ExecOnce.
   */
  virtual void UserCreateOutputObjects();

  /**
   * Event loop.
   *
   * Running over all cells and filling the cell-related
   * quantities specified in the class documentation.
   *
   * @param[in]
   */
  virtual void UserExec(Option_t *);

  /**
   * Perform initializations of the task which require a
   * run number (only available as soon as the first event
   * is available). This contains:
   * - Loading EMCAL geometry
   * - Creating histograms
   * As the dimension for most of the histograms is based on
   * the amount of cells, which is obtained from the EMCAL
   * geometry, the histograms are not initialized in
   * UserCreateOutputObjects but in ExecOnce.
   */
  virtual void ExecOnce();

  /**
   * Perform tasks necessary when the run number changes. In this
   * case a new masked cell map is loaded and masked cells are marked
   * in a historgam.
   */
  virtual void RunChanged();

  /**
   * Create the output histograms
   *
   * For all supermodules the following histograms will be created:
   * - cellAmplitude with the cell amplitude distribution
   * - cellTime with the cell time distribution
   *
   * For each supermodule the following histograms will be created:
   * - cellAmpSM with the integrated amplitude for cells within a supermodule in col and row
   * - cellCountSM with the count rate for cells within a supermodule in col and row
   *
   * Function is called in ExecOnce
   */
  void CreateHistograms();

  /**
   * Check whether cell with a given ID is masked by the user
   * @param cellId Cell ID to be checked
   * @return True if the cell is masked, false otherwise
   */
  bool IsCellMasked(Int_t cellId) const;

  /**
   * Load the cell masking from the OADB container into the task
   */
  void LoadCellMasking();

private:
  Bool_t                              fLocalInitialized;    ///< Check whether task is initialized (for ExecOnce)
  THistManager                        *fHistManager;        //!<! Histogram handler
  AliEMCALGeometry                    *fGeometry;           //!<! EMCAL geometry

  Double_t                            fMinCellAmplitude;    ///< Min. cell amplitude requested for cell time and frequency
  ULong_t                             fRequestTrigger;      ///< Trigger selection
  TString                             fTriggerString;       ///< Trigger string in addition to trigger selection
  TString                             fBadChannelContainer; ///< Bad channel container name
  TString                             fNameClusters;        ///< Name of the cluster container (as TClonesArray)
  Int_t                               fNumberOfCells;       ///< Number of cells
  Int_t                               fOldRun;              //!<! Old Run number (for run change check)

  std::vector<Int_t>                  fMaskedCells;         ///< Vector of masked cells

  AliEmcalCellMonitorTask(const AliEmcalCellMonitorTask &ref);
  AliEmcalCellMonitorTask &operator=(const AliEmcalCellMonitorTask &ref);

  /// \cond CLASSIMP
  ClassDef(AliEmcalCellMonitorTask, 1);
  /// \endcond
};

} /* namespace EMCAL */

} /* namespace PWG */

#endif /* ALIEMCALCELLMONITORTASK_H_ */
