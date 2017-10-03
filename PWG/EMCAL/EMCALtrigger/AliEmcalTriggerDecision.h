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

class AliEmcalTriggerSelectionCuts;

/**
 * @class AliEmcalTriggerDecision
 * @brief Container for trigger decision
 * @ingroup EMCALTRGFW
 * @author Markus Fasel <markus.fasel@cern.ch>
 * @since Dec 17, 2014
 *
 * Object storing the result of the EMCAL trigger decision. The result is appended to the
 * input event and can be read out by consumer tasks.
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
   * task according to the name.
   *
   * @param[in] name Name of the decision object
   * @param[in] title Title of the decision object
   */
  AliEmcalTriggerDecision(const char *name, const char *title = "");

  /**
   * @brief Destructor
   */
  virtual ~AliEmcalTriggerDecision() {}

  const AliEMCALTriggerPatchInfo *GetMainPatch() const { return fMainPatch; }
  const AliEmcalTriggerSelectionCuts *GetSelectionCuts() const { return fSelectionCuts; }
  const TList *GetAcceptedPatches() const { return &fAcceptedPatches; }
  Bool_t IsSelected() const { return fMainPatch != NULL; }

  void SetSelectionCuts(const AliEmcalTriggerSelectionCuts * const cuts) { fSelectionCuts = cuts; }
  void SetMainPatch(const AliEMCALTriggerPatchInfo * const mainpatch) { fMainPatch = mainpatch; }

  /**
   * Add accepted patch to the trigger decision
   *
   * @param[in] patch the accepted patch
   */
  void AddAcceptedPatch(AliEMCALTriggerPatchInfo * const acceptedPatch);

protected:
  const AliEMCALTriggerPatchInfo          *fMainPatch;         ///< Main trigger patch which fires the decision
  const AliEmcalTriggerSelectionCuts      *fSelectionCuts;     ///< Pointer to the cuts used for the trigger selection
  TList                                    fAcceptedPatches;   ///< All trigger patches which are accepted as well

  ClassDef(AliEmcalTriggerDecision, 1);               // Container of the trigger decision information
private:
  AliEmcalTriggerDecision(const AliEmcalTriggerDecision &ref);
  AliEmcalTriggerDecision &operator=(const AliEmcalTriggerDecision &ref);
};

}
}

#endif /* ALIEMCALTRIGGERDECISION_H */
