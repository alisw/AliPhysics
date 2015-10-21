/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
/*
 * Container storing all trigger decisions by the trigger selection task
 *
 * Author: Markus Fasel
 */
#include "AliEmcalTriggerDecision.h"
#include "AliEmcalTriggerDecisionContainer.h"

ClassImp(AliEmcalTriggerDecisionContainer)

//______________________________________________________________________________
AliEmcalTriggerDecisionContainer::AliEmcalTriggerDecisionContainer():
  TNamed(),
  fContainer()
{
  /*
   * Dummy constructor, for I/O, not to be called by the user
   */
  fContainer.SetOwner();
}

//______________________________________________________________________________
AliEmcalTriggerDecisionContainer::AliEmcalTriggerDecisionContainer(const char* name):
  TNamed(name, ""),
  fContainer()
{
  /*
   * Main constructor, called by the user
   */
  fContainer.SetOwner();
}

//______________________________________________________________________________
void AliEmcalTriggerDecisionContainer::Reset() {
  /*
   * Clear container with trigger decisions
   */
  fContainer.Clear();
}

//______________________________________________________________________________
void AliEmcalTriggerDecisionContainer::AddTriggerDecision(AliEmcalTriggerDecision* const decision) {
  /*
   * Add trigger decision to the container
   *
   * @param decision: Trigger decision, created by the trigger selection task
   */
  fContainer.Add(decision);
}

//______________________________________________________________________________
const AliEmcalTriggerDecision* AliEmcalTriggerDecisionContainer::FindTriggerDecision(const char* decname) const {
  /*
   * Find a trigger decision with a given name in the trigger decision container
   *
   * @param decname: the name of the trigger decision object
   * @return: the trigger decision (NULL if not found)
   */
  return dynamic_cast<const AliEmcalTriggerDecision *>(fContainer.FindObject(decname));
}
