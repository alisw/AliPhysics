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
   * Called by the user
   */
  AliEmcalTriggerDecisionContainer(const char *name);

  /**
   * @brief Destructor
   */
  virtual ~AliEmcalTriggerDecisionContainer() {}

  /**
   * @brief Clear container with trigger decisions
   */
  void Reset();

  /**
   * @brief Clear function of the trigger decision container
   * 
   * Overwriting the Clear function preventing to reset the name.
   * Onlu clearing the trigger decision objects content.
   * @param option Not used
   */
  virtual void Clear(Option_t *option="") { Reset(); }

  /**
   * @brief Add trigger decision to the container
   *
   * @param[in] decision Trigger decision, created by the trigger selection task
   */
  void AddTriggerDecision(AliEmcalTriggerDecision * const decision);

  /**
   * @brief Find a trigger decision with a given name in the trigger decision container
   *
   * @param[in] decname the name of the trigger decision object
   * @return the trigger decision (NULL if not found)
   */
  const AliEmcalTriggerDecision *FindTriggerDecision(const char *name) const;

  /**
   * @brief Get container with trigger decision results
   * @return List of trigger decision objects
   */
  const TList *GetListOfTriggerDecisions() const { return &fContainer; }

  /**
   * @brief Checks whether the events is selected for a given trigger type
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
