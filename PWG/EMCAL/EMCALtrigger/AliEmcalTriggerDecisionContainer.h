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
#ifndef ALIEMCALTRIGGERDECISIONCONTAINER_H
#define ALIEMCALTRIGGERDECISIONCONTAINER_H

#include <TList.h>
#include <TNamed.h>

namespace PWG{

namespace EMCAL {

class AliEmcalTriggerDecision;

/**
 * @class AliEmcalTriggerDecisionContainer
 * @brief Container for trigger decision object
 * @ingroup EMCALTRGFW
 * @author Markus Fasel <markus.fasel@cern.ch> Oak Ridge National Laboratory
 * @since Dec 17, 2014
 * 
 * # Event selection status of EMCAL level 1 triggers
 * 
 * In case the event selection of EMCAL Level1 triggers is done from trigger patches,
 * the selection results for all supported triggers are stored in the AliEmcalTriggerDecision
 * container for easy event selection. In case the event is selected as triggered event
 * the primitives firing the trigger are also linked so they can be accessed in user analyses.
 * 
 * ## 1) Accessing the trigger decision container in user task
 * 
 * In order to get the trigger decision the AliAnalysisTaskEmcalTriggerSelection needs to
 * run together with the AliEmcalTriggerMakerTask before the actual user task. In this 
 * case the trigger selection task creates an object of type 
 * PWG::EMCAL::AliEmcalTriggerDecisionContainer, if not specified differently with name
 * *EmcalTriggerDecision*, and adds it to the input event as list object. Users can then
 * access it in the task via
 * 
 * ~~~{.cxx}
 * auto trgsel = static_cast<PWG::EMCAL::AliEmcalTriggerDecisionContainer *>(fInputEvent->FindListObject("EmcalTriggerDecision"));
 * ~~~ 
 * 
 * ## 2) Selection of Level1 triggers from the trigger decision container
 * 
 * Once the trigger selection container is available, EMCAL Level1 triggers can be queried
 * using the function IsEventSelected(trigger) with a valid trigger name. Trigger names
 * are expected to be valid Level1 trigger names according to the EMCAL convention
 * 
 * - pp, 2012: EGA (gamma trigger), EJE (jet trigger)
 * - p-Pb, 2013: EG1, EG2 (gamma triggers), EJ1, EJ2 (jet triggers)
 * - run2 pp/p-Pb: EG1, EG2 (gamma triggers, EMCAL), DG1, DG2 (gamma triggers, DCAL)
 *                 EJ1, EJ2 (jet trigger, EMCAL), DJ1, DJ2 (jet triggers, DCAL)
 * 
 * The following example selects events as EMCAL EG1 events:
 * 
 * ~~~{.cxx}
 * bool SelectEG1(const AliVEvent *input){
 *   auto trgsel = static_cast<PWG::EMCAL::AliEmcalTriggerDecisionContainer *>(input->FindListObject("EmcalTriggerDecision"))
 *   if(!trgsel->IsEventSelected("EG1")) return false;
 *   return true;
 * }
 * ~~~
 * 
 * ## 3) Retrieving trigger primitives for fired triggers
 * 
 * More information about the trigger selection for a given triggers are attached as
 * objects of type PWG::EMCAL::AliEmcalTriggerDecision to the AliEmcalTriggerDecisionContainer.
 * They can be queried from the trigger decision container via the function FindTriggerDecision(triggername).
 * 
 * The following example extracts the maximum trigger patch firing the trigger for the EG1 trigger:
 * 
 * ~~~{.cxx}
 * auto trgdec = trgsel->FindTriggerDecision("EG1");
 * if(trgsel){
 *   auto maxpatch = trgdec->GetMainPatch();
 *   if(maxpatch){
 *     std::cout << "Max patch found - event was triggered" << std::endl;
 *   }
 * }
 * ~~~
 */
class AliEmcalTriggerDecisionContainer: public TNamed {
public:
  /**
   * @brief Dummy constructor
   *
   * For I/O, not to be called by the user
   */
  AliEmcalTriggerDecisionContainer();

  /**
   * @brief Main constructor
   *
   * Called by the user. The name of the trigger decision container
   * is used in user analyses for the query of the container from
   * the event.
   * 
   * @param[in] name Name of the trigger decision container
   */
  AliEmcalTriggerDecisionContainer(const char *name);

  /**
   * @brief Destructor
   */
  virtual ~AliEmcalTriggerDecisionContainer() {}

  /**
   * @brief Clear container with trigger decisions
   * 
   * Removing all trigger decision objects from the container. The container
   * will be empty at this point.
   */
  void Reset();

  /**
   * @brief Clear function of the trigger decision container
   * 
   * Overwriting the Clear function preventing to reset the name.
   * Only clearing the trigger decision objects content.
   * @param option Not used
   */
  virtual void Clear(Option_t *option="") { Reset(); }

  /**
   * @brief Add trigger decision to the container
   * 
   * The name of the trigger decision corresponds to the name of the
   * trigger class according to the EMCAL naming convention and is used
   * by the user to query the event trigger selection status.
   *
   * @param[in] decision Trigger decision, created by the trigger selection task
   */
  void AddTriggerDecision(AliEmcalTriggerDecision * const decision);

  /**
   * @brief Find a trigger decision with a given name in the trigger decision container
   * 
   * Objects of type PWG::EMCAL::AliEmcalTriggerDecision contain additional information
   * beside the event selection status
   * 
   * - Main patch (highest-energetic patch)
   * - All patches firing the trigger condition
   * 
   * The name of the trigger decision object must be a valid trigger class name according
   * to the EMCAL convention, where the trigger class is supported for the dataset to
   * be analyzed.
   *
   * @param[in] decname Name of the EMCAL Level1 trigger class
   * @return the trigger decision (NULL if not found)
   */
  const AliEmcalTriggerDecision *FindTriggerDecision(const char *name) const;

  /**
   * @brief Get container with trigger decision results
   * 
   * Objects of type PWG::EMCAL::AliEmcalTriggerDecision contain additional information
   * beside the event selection status
   * 
   * - Main patch (highest-energetic patch)
   * - All patches firing the trigger condition
   * 
   * The list concerns the trigger decision objects for all Level1 triggers supported in
   * the dataset, no matter whether the trigger was fired or not. For a specific dataset 
   * the trigger decision object can be found in the output list by the name of the trigger
   * class, following the EMCAL nameing convention.
   * 
   * @return List of trigger decision objects
   */
  const TList *GetListOfTriggerDecisions() const { return &fContainer; }

  /**
   * @brief Checks whether the events is selected for a given trigger type
   * 
   * Selecting the trigger selection status for the given Level1 trigger.
   * The trigger decision is done base on the existence of a trigger patch
   * of a given type in the event above threshold, where patch type and
   * threshold are defined in the configuration of task 
   * PWG::EMCAL::AliAnalysisTaskEmcalTriggerSelection as in the trigger
   * selection cuts. The name of the trigger must be a valid Level1 trigger
   * name according to the EMCAL naming convention.
   * 
   * @param[in] name Name of the EMCAL trigger
   * @return True if the event is selected, false otherwise
   */
  bool IsEventSelected(const char *name)  const;

protected:
  TList     fContainer;         ///< List of trigger decisions

  /// \cond CLASSIMP
  ClassDef(AliEmcalTriggerDecisionContainer, 1);    // Container for trigger decisions
  /// \endcond
};

}
}

#endif /* ALIEMCALTRIGGERDECISIONCONTAINER_H */
