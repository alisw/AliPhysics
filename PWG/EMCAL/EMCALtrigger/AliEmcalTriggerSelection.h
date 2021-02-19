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
#ifndef ALIEMCALTRIGGERSELECTION_H
#define ALIEMCALTRIGGERSELECTION_H

#include <iosfwd>
#include <TNamed.h>
#include <TString.h>
#include "AliEmcalTriggerSelectionCuts.h"

class AliEMCALTriggerPatchInfo;
class TClonesArray;

// operator<< has to be forward declared carefully to stay in the global namespace so that it works with CINT.
// For generally how to keep the operator in the global namespace, See: https://stackoverflow.com/a/38801633
namespace PWG { namespace EMCAL { class AliEmcalTriggerSelection; } }
std::ostream & operator<< (std::ostream &in, const PWG::EMCAL::AliEmcalTriggerSelection &sel);

namespace PWG {

namespace EMCAL {

class AliEmcalTriggerAlias;
class AliEmcalTriggerDecision;

/**
 * @class AliEmcalTriggerSelection
 * @brief Object performing offline EMCAL trigger selection
 * @ingroup EMCALTRGFW
 * @author Markus Fasel <markus.fasel@cern.ch> Oak Ridge National Laboratory
 * @since Dec 17, 2014
 *
 * Object performing an offline EMCAL trigger decision based on user defined criterions
 * (trigger patch type, energy threshold,...). The main method MakeTriggerDecision performs
 * an event selection and creates a trigger decision object with the relevant information.
 */
class AliEmcalTriggerSelection: public TNamed {
public:

  /**
   * @brief Dummy constructor
   *
   * Used by I/O, not to be used by the user
   */
  AliEmcalTriggerSelection();

  /**
   * @brief Main constructor
   *
   * Initialises the trigger selection
   *
   * @param name name of the trigger selection
   * @param cuts(optional) selection cuts to be applied
   * @param alias(optional) trigger alias
   */
  AliEmcalTriggerSelection(const char *name, const AliEmcalTriggerSelectionCuts * const cuts, const AliEmcalTriggerAlias *alias = nullptr);

  /**
   * @brief Destructor
   */
  virtual ~AliEmcalTriggerSelection() {}

  /**
   * @brief Get the trigger alias 
   * @return Trigger alias (nullptr if not set)
   */
  const AliEmcalTriggerAlias *GetTriggerAlias() const { return fTriggerAlias; }

  /**
   * @brief Get the selection cuts used in the trigger selection
   * @return Selection cuts used in this trigger selection
   */
  const AliEmcalTriggerSelectionCuts *GetSelectionCuts() const { return fSelectionCuts; }

  /**
   * @brief Set the selection cuts used in this trigger selection
   * @param[in] cuts Selection cuts used to fire the trigger
   */
  void SetSelectionCuts(const AliEmcalTriggerSelectionCuts * const cuts) { fSelectionCuts = cuts; }

  void SetTriggerAlias(const AliEmcalTriggerAlias *const alias) { fTriggerAlias = alias; }

  /**
   * Perform event selection based on user-defined criteria and create an output trigger decision containing
   * the threshold, the main patch which fired the decision, and all other patches which would have fired the
   * decision as well.
   *
   * @param[in] reconstructedPatches A list of input patches, created by the trigger patch maker and read out from the
   * input event
   * @param[in] rhocontainer Container with rho values for the given event
   * @return the trigger decision (an event is selected when it has a main patch that fired the decision)
   */
  AliEmcalTriggerDecision * MakeDecison(const TClonesArray * const reconstructedPatches, const AliEmcalTriggerSelectionCuts::RhoForTrigger &rhocontainer) const;

  /**
   * @brief Output stream operator
   * 
   * Logging trigger selection configuration (name, cuts) to the output stream
   * 
   * @param stream Stream used for logging
   * @param sel Trigger selection object to be logged
   * @return Output stream after logging
   */
  friend std::ostream& ::operator<<(std::ostream &stream, const AliEmcalTriggerSelection &sel);
protected:
  const AliEmcalTriggerSelectionCuts  *fSelectionCuts;    ///< Cuts used for the trigger patch selection
  const AliEmcalTriggerAlias *fTriggerAlias;              ///< Trigger alias for the trigger selection

  /// \cond CLASSIMP
  ClassDef(AliEmcalTriggerSelection, 1);    // EMCAL trigger selection component
  /// \endcond

private:
  AliEmcalTriggerSelection(const AliEmcalTriggerSelection &ref);
  AliEmcalTriggerSelection &operator=(const AliEmcalTriggerSelection &ref);

  /**
   * @brief Helper function for the stream operator
   * 
   * Performs logging
   * 
   * @param stream Stream used for logging
   */
  void PrintStream(std::ostream &stream) const;
};

}
}



#endif /* ALIEMCALTRIGGERSELECTION_H */
