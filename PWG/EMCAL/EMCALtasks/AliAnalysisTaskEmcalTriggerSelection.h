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
#include <TString.h>
#include "AliAnalysisTaskEmcal.h"

namespace PWG{
namespace EMCAL {

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
  virtual ~AliAnalysisTaskEmcalTriggerSelection() {}

  /**
   * Add trigger selection to the trigger selection task
   *
   * @param[in] selection the trigger selection to be added
   */
  void AddTriggerSelection(AliEmcalTriggerSelection * const selection);
  void SetGlobalDecisionContainerName(const char *name) { fGlobalDecisionContainerName = name; }

  /**
   * @brief Run over all trigger selections, and append the selection to the global trigger selection container
   */
  virtual Bool_t Run();

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
   * @brief Find the main trigger container in the input event.
   *
   * If not available, create it and add it to the input event
   */
  AliEmcalTriggerDecisionContainer *GetGlobalTriggerDecisionContainer();

  TString fGlobalDecisionContainerName;     ///< Name of the global trigger selection
  TList fTriggerSelections;                 ///< List of trigger selections

  ClassDef(AliAnalysisTaskEmcalTriggerSelection, 1);    // Task running different EMCAL trigger selections
};

}
}

#endif /* ALIANALYSISTASKEMCALTRIGGERSELECTION_H */
