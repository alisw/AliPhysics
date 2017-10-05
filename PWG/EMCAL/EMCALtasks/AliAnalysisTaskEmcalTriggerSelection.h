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
#ifndef ALIANALYSISTASKEMCALTRIGGERSELECTION_H
#define ALIANALYSISTASKEMCALTRIGGERSELECTION_H

#include <TList.h>
#include <TNamed.h>
#include <TString.h>
#include "AliAnalysisTaskEmcal.h"

namespace PWG{
namespace EMCAL {

class AliEmcalTriggerDecision;
class AliEmcalTriggerDecisionContainer;
class AliEmcalTriggerSelection;

/**
 * @class AliAnalysisTaskEmcalTriggerSelection
 * @brief Task providing an event selection for EMCAL-triggered events based on the
 * reconstructed EMCAL trigger patches
 * @ingroup EMCALFWTASKS
 * @author Markus Fasel <markus.fasel@cern.ch> Oak Ridge National Laboratory
 * @since Dec 17, 2014
 */
class AliAnalysisTaskEmcalTriggerSelection: public AliAnalysisTaskEmcal {
public:
  /**
   * @brief Dummy constructor
   *
   * Only for I/O, not to be used by the user
   */
  AliAnalysisTaskEmcalTriggerSelection();

  /**
   * @brief Main constructor
   *
   * To be called by the users
   */
  AliAnalysisTaskEmcalTriggerSelection(const char *name);

  /**
   * @brief Destructor
   */
  virtual ~AliAnalysisTaskEmcalTriggerSelection() {}

  /**
   * Add trigger selection to the trigger selection task
   *
   * @param[in] selection the trigger selection to be added
   */
  void AddTriggerSelection(AliEmcalTriggerSelection * const selection);

  /**
   * @brief Set the name of the global trigger decision container
   *
   * Other tasks have to connet to the container via this name.
   *
   * @param[in] name Name of the trigger decision container
   */
  void SetGlobalDecisionContainerName(const char *name) { fGlobalDecisionContainerName = name; }

  /**
   * @brief Automatically configure trigger decision handler for different periods
   */
  void AutoConfigure(const char *period);

  /**
   * @brief Configure the trigger selection task for pp anchored to 2016
   *
   * Using recalc patches (recalculated from FASTOR, no STU trigger decision) and
   * ADC cut, where the cut values are set to the nominal ADC thresholds
   */
  void ConfigurePP2016();

  /**
   * @brief Configure the trigger selection task for MC anchored to pp 2016
   *
   * Using simple offline patches (from cells) and energy cut, where the cut
   * values are set to the expected energy cuts
   */
  void ConfigureMCPP2016();

protected:

  /**
   * @class AliEmcalTriggerSelectionQA
   * @brief Helper class for the trigger selection
   * @author Markus Fasel <markus.fasel@cern.ch>, Oak Ridge National Laboratory
   * @since Oct 3, 2017
   */
  class AliEmcalTriggerSelectionQA : public TNamed {
  public:

    /**
     * @brief Dummy constructor
     */
    AliEmcalTriggerSelectionQA();

    /**
     * @brief Constructor
     *
     * Initializing component matching to the trigger selection
     *
     * @param[in] sel Corresponding trigger selection to be monitored by this QA component
     */
    AliEmcalTriggerSelectionQA(const AliEmcalTriggerSelection * const sel);

    /**
     * @brief Copy constructor
     *
     * Histograms will be shared between the two QA components
     * @param[in] ref Reference for copy
     */
    AliEmcalTriggerSelectionQA(const AliEmcalTriggerSelectionQA &ref);

    /**
     * @brief Assignment operator
     *
     * Histograms will be shared between the two QA components
     * @param[in] ref Reference for assignment
     * @return This object after assignment
     */
    AliEmcalTriggerSelectionQA &operator=(const AliEmcalTriggerSelectionQA &ref);

    /**
     * @brief Destructor
     */
    virtual ~AliEmcalTriggerSelectionQA() {}

    /**
     * Fill histograms for the main patch for the given trigger decision
     * @param[in] decision
     */
    void Fill(const AliEmcalTriggerDecision * const decision);

    /**
     * @brief Fill histograms of this QA component into the targetlist
     *
     * The outputlist will take ownership over the histograms
     *
     * @param[out] targetlist list toe
     */
    void GetHistos(TList *targetlist) const;

  private:
     TH1 *fMaxPatchADC;                 ///< Histogram with patch ADC of the max patch
     TH1 *fMaxPatchEnergy;              ///< Histogram with patch energy of the max patch
     TH1 *fMaxPatchEnergySmeared;       ///< Histogram with smeared patch energy of the max patch
  };

  /**
   * @brief Initialize QA histograms
   */
  virtual void UserCreateOutputObjects();

  /**
   * @brief Initializing common output container for trigger decision
   */
  virtual void UserExecOnce();

  /**
   * @brief Run over all trigger selections, and append the selection to the global trigger selection container
   * @return Always true
   */
  virtual Bool_t Run();

  /**
   * @brief Filling basic QA Histograms of the trigger selection task
   *
   * The QA histograms are connected to the main patch and monitor
   * - ADC Amplitude
   * - Energy
   * - Smeared energy
   * @return Always true
   */
  virtual Bool_t FillHistograms();

  /**
   * @brief Find the main trigger container in the input event.
   *
   * If not available, create it and add it to the input event
   */
  AliEmcalTriggerDecisionContainer *GetGlobalTriggerDecisionContainer() const { return fTriggerDecisionContainer; }

  /**
   * @brief Fill QA histograms for the event
   * @param cont
   */
  void MakeQA(const AliEmcalTriggerDecisionContainer *cont);

  /**
   * @brief Initialize QA histograms for trigger selection
   * @param sel
   */
  void InitQA(const AliEmcalTriggerSelection *const sel);

  AliEmcalTriggerDecisionContainer          *fTriggerDecisionContainer;        ///<
  TString                                    fGlobalDecisionContainerName;     ///< Name of the global trigger selection
  TList                                      fTriggerSelections;               ///< List of trigger selections
  TList                                      fSelectionQA;                     ///< Trigger selection QA

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskEmcalTriggerSelection, 1);    // Task running different EMCAL trigger selections
  /// \endcond
};

}
}

#endif /* ALIANALYSISTASKEMCALTRIGGERSELECTION_H */
