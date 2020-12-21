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
#ifndef ALIEMCALTRIGGERDECISION_H
#define ALIEMCALTRIGGERDECISION_H

#include <TList.h>
#include <TNamed.h>

class AliEMCALTriggerPatchInfo;

namespace PWG {

namespace EMCAL {

class AliEmcalTriggerAlias;
class AliEmcalTriggerSelectionCuts;

/**
 * @class AliEmcalTriggerDecision
 * @brief Container for trigger decision
 * @ingroup EMCALTRGFW
 * @author Markus Fasel <markus.fasel@cern.ch>
 * @since Dec 17, 2014
 * 
 * # Status of the trigger selection process for a given Level1 trigger
 * 
 * AliEmcalTriggerDecision object collect all relevant information for 
 * a given Level1 trigger:
 * 
 * - Maximum patch firing the trigger
 * - All other patches firing the trigger
 * - A corresponding trigger selection cuts object defining the selection process
 * 
 * An AliEmcalTriggerDecision object handles only one Level1 trigger class,
 * each trigger class supported by the dataset has its own AliEmcalTriggerDecision
 * object.
 * 
 * ## Checking whether an event was triggered for the given Level1 trigger
 * 
 * The presence of a maximum patch marks an event as triggered. Consequently it
 * is sufficient to check for the presence of this. The following example selects
 * events triggered by the EG1 trigger:
 * 
 * ~~~{.cxx}
 * auto trgcont = static_cast<PWG::EMCAL::AliEmcalTriggerDecisionContainer*>(fInputEvent->FindListObject("EmcalTriggerDecision"));
 * auto eg1 = trgcont->FindTriggerDecision("EG1");
 * if(eg1->GetMainPatch()) {
 *   std::cout << "Event is an EG1 event";
 * }
 * ~~~
 * 
 * ## Getting informatinon about all patches firing the trigger
 *
 * The trigger is fired if at least one patch above threshold according to the definition in the associated
 * trigger selection cuts is found. Consequently all patches above  threshold are valid trigger patches. 
 * They are attached to the event and can be queried via  GetAcceptedPatches(). The following example draws 
 * the energy spectrum of all accepted trigger patches:
 * 
 * ~~~{.cxx}
 * auto hPatchEnergy = new TH1F("hPatchEnergy", "PatchEnergy", 200, 0., 200.);
 * auto trgcont = static_cast<PWG::EMCAL::AliEmcalTriggerDecisionContainer*>(fInputEvent->FindListObject("EmcalTriggerDecision"));
 * auto eg1 = trgcont->FindTriggerDecision("EG1");
 * for(auto p : eg1) {
 *   auto patch = static_cast<AliEMCALTriggerPatchInfo *>(p);
 *   hPatchEnergy->Fill(patch->GetPatchE());  
 * } 
 * ~~~
 * 
 * The AliEmcalTriggerDecision object is not owner of the patches. Delete calls would probably lead to
 * double delets. Users must not delete patches attached.
 * 
 * ## Specification of the level1 trigger
 * 
 * The level1 trigger is specified in the corresponding AliEmcalTriggerSelectionCuts object. The specifications
 * consist of 
 * 
 * - Patch type (Gamma or Jet patch)
 * - Energy definition (FastOR ADC amplitude or FEE energy, smeared or non-smeared)
 * - Trigger threshold
 * 
 * The associated trigger selection cut object is used in the trigger selection process and linked to 
 * this trigger selection object. It can be queried via GetSelectionCuts(). The following example indicates 
 * how to query the trigger threshold:
 * 
 * ~~~{.cxx}
 * auto trgcont = static_cast<PWG::EMCAL::AliEmcalTriggerDecisionContainer*>(fInputEvent->FindListObject("EmcalTriggerDecision"));
 * auto eg1 = trgcont->FindTriggerDecision("EG1"); 
 * auto cuts = eg1->GetSelectionCuts();
 * std::cout << "Threshold for EG1: " << cuts->GetThreshold() << std::endl;
 * ~~~
 * 
 * The trigger decision object is not owner of the trigger selection cuts. Users must not delete them.
 */
class AliEmcalTriggerDecision: public TNamed {
public:

  /**
   * @brief Dummy constructor
   *
   * Needed for I/O, not to be used by the user
   */
  AliEmcalTriggerDecision();

  /**
   * @brief The main (named) constructor.
   *
   * The decision object can be read out later by the consumer
   * task according to the name. The name has to be a valid Level1
   * trigger name of the corresponding Level1 trigger supported
   * in the data set, according to the EMCAL naming convention.
   *
   * @param[in] name Name of the decision object
   * @param[in] title Title of the decision object
   */
  AliEmcalTriggerDecision(const char *name, const char *title = "");

  /**
   * @brief Destructor
   * 
   * As this class is not owner of the trigger patches the destructor will
   * not delete them. Deleting trigger patches has to be the responsibility
   * of the owner (typically the TClonesArray created by the trigger maker
   * attached to the input event).
   */
  virtual ~AliEmcalTriggerDecision();


  /**
   * @brief Get the trigger alias
   * @return corresponding trigger alias (nullptr if not set)
   */
  const AliEmcalTriggerAlias *GetTriggerAlias() const { return fTriggerAlias; }

  /**
   * @brief Get the highest energetic trigger patch of the event firing the trigger
   * 
   * The highest energic patch is defined according to the energy definition (FastOR
   * ADC amplitude or FEE energy) and matching patch type above threshold. The presence
   * of a main patch indicates that the event was fired for the Level1 trigger handled
   * by this trigger selection. 
   * 
   * @return Highest energetic trigger patch of the event. In case the event was not triggered
   *         returning null.
   */
  const AliEMCALTriggerPatchInfo *GetMainPatch() const { return fMainPatch; }

  /**
   * @brief Get the selection cuts used in the trigger selection
   * 
   * Selection cuts specify the trigger patch selection configuration for the
   * Level1 trigger producing this trigger selection result (patch type, energy
   * definition, trigger threshold) and process the trigger patch selection.
   * 
   * @return Selection cuts used for the corresponding trigger class
   */
  const AliEmcalTriggerSelectionCuts *GetSelectionCuts() const { return fSelectionCuts; }

  /**
   * @brief Get the list of all patches in the event satisfying the trigger condition
   * 
   * The trigger patch selection is specified in the PWG::EMCAL::AliEmcalTriggerSelectionCuts 
   * object which can be obtained via GetSelectionCuts() and is configured in the task
   * PWG::EMCAL::AliAnalysisTaskEmcalTrigger selection. Selected trigger patches must
   * be above the trigger threshold for the given trigger class according to the energy
   * definition specified in the corresponding AliEmcalTriggerSelectionCuts (ADC amplitude
   * or offline patch energy)
   * 
   * @return List of acccepted patches: Containing all patches which would fire the trigger.
   *         Returning an empty list in case the event was not selected.
   */
  const TList *GetAcceptedPatches() const { return &fAcceptedPatches; }

  /**
   * @brief Check whether event is selected under the given trigger
   *
   * An event is selected if a main (highest energy) patch was found,
   * indicating at least one patch was above nominal threshold according
   * to the energy definition specified in the corresponding AliEmcalTriggerSelectionCuts
   * object.
   * 
   * @return True if the event was selected, false otherwise
   */
  Bool_t IsSelected() const { return fMainPatch != NULL; }

  /**
   * @brief Set trigger alias (names of trigger classes corresponding to the trigger condintion)
   * @param alias Trigger alias object
   */
  void SetTriggerAlias(const AliEmcalTriggerAlias *alias) { fTriggerAlias = alias; }

  /**
   * @brief Set the selection cuts used in the trigger selection
   * 
   * Selection cuts specify the trigger patch selection configuration for the
   * Level1 trigger producing this trigger selection result (patch type, energy
   * definition, trigger threshold) and process the trigger patch selection.
   * 
   * Users can access the trigger selection cuts in thier task via the corresponding 
   * getter (GetSelectionCuts)
   * 
   * @param[in] cuts Selection cuts for the given trigger class
   */
  void SetSelectionCuts(const AliEmcalTriggerSelectionCuts * const cuts) { fSelectionCuts = cuts; }

  /**
   * @brief Set the main (highest-energetic) trigger patch
   * 
   * The main patch is defined as the highest energetic trigger patch
   * according to patch type and energy definition selected by the 
   * associated trigger selection cuts.
   * 
   * Setting the main patch marks the event as triggered for the 
   * given Level1 trigger producing this result.
   * 
   * @param[in] mainpatch Highest energetic trigger patch of the event firing the trigger
   */
  void SetMainPatch(const AliEMCALTriggerPatchInfo * const mainpatch) { fMainPatch = mainpatch; }

  /**
   * Add accepted patch to the trigger decision
   * 
   * Patches added to the trigger selection object must comply with
   * the trigger definition specified in the corresponding 
   * AliEmcalTriggerDecision object. All accepted patches would have 
   * fired the corresponding Level1 trigger.
   *
   * @param[in] patch the accepted patch
   */
  void AddAcceptedPatch(AliEMCALTriggerPatchInfo * const acceptedPatch);

protected:
  const AliEMCALTriggerPatchInfo          *fMainPatch;         ///< Main trigger patch which fires the decision
  const AliEmcalTriggerSelectionCuts      *fSelectionCuts;     ///< Pointer to the cuts used for the trigger selection
  const AliEmcalTriggerAlias              *fTriggerAlias;      ///< Pointer to trigger alias object with names of corresponding classes
  TList                                    fAcceptedPatches;   ///< All trigger patches which are accepted as well

  /// \cond CLASSIMP
  ClassDef(AliEmcalTriggerDecision, 1);               // Container of the trigger decision information
  /// \endcond
private:
  AliEmcalTriggerDecision(const AliEmcalTriggerDecision &ref);
  AliEmcalTriggerDecision &operator=(const AliEmcalTriggerDecision &ref);
};

}
}

#endif /* ALIEMCALTRIGGERDECISION_H */
