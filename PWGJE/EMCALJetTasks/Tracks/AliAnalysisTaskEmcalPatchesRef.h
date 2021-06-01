/************************************************************************************
 * Copyright (C) 2015, Copyright Holders of the ALICE Collaboration                 *
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
#ifndef ALIANALYSISTASKEMCALPATCHESREF_H
#define ALIANALYSISTASKEMCALPATCHESREF_H

#include <string>
#include <vector>
#include "AliAnalysisTaskEmcalTriggerBase.h"
#include "AliEMCALTriggerDataGrid.h"
#include "AliCutValueRange.h"
#include <TCustomBinning.h>
#include <TString.h>

class AliOADBContainer;
class AliEMCALTriggerPatchInfo;
class THistManager;
class TObjArray;

namespace PWGJE {
  
namespace EMCALJetTasks {

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
   * @brief Dummy (I/O) onstructor
   */
  AliAnalysisTaskEmcalPatchesRef();

  /**
   * @brief Named constructor
   * @param[in] name Name of the task
   */
  AliAnalysisTaskEmcalPatchesRef(const char *name);

  /**
   * @brief Destructor
   */
  virtual ~AliAnalysisTaskEmcalPatchesRef(){}

  /**
   * @brief Enable Sumw2 when creating the histograms.
   *
   * Attention: Enabling Sumw2 will increase memory consumption
   * significantly. Option should only be used in case histograms
   * are filled with a weight.
   *
   * @param[in] doEnable If true Sumw2 is enabled for all histograms
   */
  void EnableSumw2(Bool_t doEnable) { fEnableSumw2 = doEnable; }

  /**
   * @brief Set centrality selection.
   *
   * Note: Needs multiplicity task to run in front
   *
   * @param[in] min Min. value of the centrality interval
   * @param[in] max Max. value of the centrality interval
   */
  void SetCentralityRange(double min, double max) { fCentralityRange.SetLimits(min,max); fRequestCentrality = true; }

  /**
   * @brief Set time range for cell amplitudes for filling the offline trigger data grid
   * @param min Min. cell time
   * @param max Max. cell time
   */
  void SetCellTimeRange(double min, double max) { fCellTimeRange.SetLimits(min, max); }

  /**
   * @brief Set minimum number of FastORs used to accept the patch
   * @param min Min. number of FastORs
   */
  void SetMinNumberOfFastorsPatch(int min) { fMinNumberFastors = min; }

  /**
   * @brief Switch for recalc patches
   *
   * By default the task analyzes offline patches and offline energy (meaning
   * energy calculated from FEE cells).  Switching on recalc mode the task
   * is using recalc patches (calculated from L1 time sums, but ignoring STU
   * decision).
   *
   * @param[in] doUse If true then recalc patches (and recalc energy) is used
   */
  void SetUseRecalcPatches(Bool_t doUse = kTRUE) { fUseRecalcPatches = doUse; }

protected:

  /**
   * Task does not need user objects, so the function is not implemented
   */
  virtual void CreateUserObjects() {}

  /**
   * @brief Creating output histograms.
   *
   * Supported user histograms are:
   * - Patch (calibrated) energy spectrum - separated by patch type - for different trigger classes
   * - Patch eta-phi map - separated by patch type - for different trigger classes and different min. energies
   */
  virtual void CreateUserHistos();

  /**
   * @brief Apply additional event selection requested by the user
   *
   * Handle only centrality selection (if requested, and multiplicity
   * information available). Supported is only new centrality framework.
   *
   * @return True if event is in proper centrality class
   */
  virtual bool IsUserEventSelected();

  /**
   * @brief Main event loop.
   *
   * Loop over all patches, assign them to trigger classes, and
   * fill relevant histograms. Event loop is only called once the event
   * selection (common and user) is passed.
   *
   * @param Not used
   */
  virtual bool Run();

  /**
   * @brief Fill event-based histograms.
   *
   * Monitored are:
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
   * @brief Selecting offline/recalc simple single shower (gamma or L0) patches
   * @param[in] patch Patch to check
   * @return True if the patch is selected
   */
  bool SelectSingleShowerPatch(const AliEMCALTriggerPatchInfo *patch) const;

  /**
   * @brief Selecting offline/recalc simple jet patches.
   * @param[in] patch Patch to check
   * @return True if the patch is selected
   */
  bool SelectJetPatch(const AliEMCALTriggerPatchInfo  *patch) const;

  /**
   * @brief Filling patch related histogram.
   *
   * In case a downscaling correction is available it is applied to the
   * histograms as weight.
   *
   * @param[in] triggerclass Name of the trigger class firing the event
   * @param[in] patchname Name of the patchtype
   * @param[in] energy Calibrated energy (calculated from cells) of the patch
   * @param[in] transverseenergy Calibrated transverse energy (calculated from cells) of the patch
   * @param[in] smearedenergy Calibrated energy (from cells) after gaussian smearing at fastor level
   * @param[in] eta Patch \f$ eta \f$ at the geometrical center
   * @param[in] phi Patch \f$ varphi \f$ at the geometrical center
   * @param[in] col Column (\f$ \eta \f$-direction) of the patch
   * @param[in] row Row (\f$ \varphi \f$-direction) of the patch
   */
  void FillPatchHistograms(TString triggerclass, TString patchname, double energy, double transverseenergy, double smearedenergy, double eta, double phi, int col, int row);

  /**
   * @brief Fill trigger data grid
   * 
   * Fill offline FEE grid in FastOR dimenision in order to 
   * provide access single FastOR energies included in the
   * patch
   */
  void FillDataGrid();

  /**
   * @brief Get the number of non-0 FastORs contributing to the trigger patch
   * @param patch Trigger patch to analyse
   * @return Number of non-0 FastORs
   */
  Int_t GetNumberOfFastORs(const AliEMCALTriggerPatchInfo *patch) const;

  AliCutValueRange<double>            fCentralityRange;           ///< Range of accepted event centralities
  AliCutValueRange<double>            fCellTimeRange;             ///< Range of accepted cell time (for data grid)
  AliEMCALTriggerDataGrid<double>     fOfflineTriggerData;        ///< Trigger energy data grid       
  Double_t                            fMinNumberFastors;          ///< Min. number of non-0 fastors contributing to the patch
  Bool_t                              fEnableSumw2;               ///< Enable sumw2 during histogram creation
  Bool_t                              fUseRecalcPatches;          ///< Switch between offline (FEE) and recalc (L1) patches
  Bool_t                              fRequestCentrality;         ///< Switch for request of centrality selection
  Double_t                            fEventCentrality;           //!<! Event centrality

private:

  /**
   * @class EnergyBinning
   * @brief Common energy binning used in the patch QA task
   */
  class EnergyBinning : public TCustomBinning {
  public:
    /**
     * @brief Constructor
     */
    EnergyBinning();

    /**
     * @brief Destructor
     */
    virtual ~EnergyBinning() {}
  };

  AliAnalysisTaskEmcalPatchesRef(const AliAnalysisTaskEmcalPatchesRef &);
  AliAnalysisTaskEmcalPatchesRef &operator=(const AliAnalysisTaskEmcalPatchesRef &);

  ClassDef(AliAnalysisTaskEmcalPatchesRef, 1);
};

} /* namespace EMCALJetTasks */

} /* namespace PWGJE */

#endif /* ALIANALYSISTASKEMCALPATCHESREF_H */
