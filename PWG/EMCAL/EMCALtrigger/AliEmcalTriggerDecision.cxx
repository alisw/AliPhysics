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
 * Object storing the result of the EMCAL trigger decision. The result is appended to the
 * input event and can be read out by consumer tasks.
 *
 * Author: Markus Fasel
 */

#include "AliEmcalTriggerDecision.h"
#include "AliEMCALTriggerPatchInfo.h"

ClassImp(AliEmcalTriggerDecision)

//______________________________________________________________________________
AliEmcalTriggerDecision::AliEmcalTriggerDecision():
  TNamed(),
  fMainPatch(NULL),
  fSelectionCuts(NULL),
  fAcceptedPatches()
{
  /*
   * Dummy constructor, needed for I/O, not to be used by the user
   */
  fAcceptedPatches.SetOwner(kFALSE);
}

//______________________________________________________________________________
AliEmcalTriggerDecision::AliEmcalTriggerDecision(const char *name, const char *title):
  TNamed(name, title),
  fMainPatch(NULL),
  fSelectionCuts(NULL),
  fAcceptedPatches()
{
  /*
   * The main (named) constructor. The decision object can be read out later by the consumer
   * task according to the name.
   *
   * @param name: Name of the decision object
   * @param title: Title of the decision object
   */
  fAcceptedPatches.SetOwner(kFALSE);
}

//______________________________________________________________________________
void AliEmcalTriggerDecision::AddAcceptedPatch(AliEMCALTriggerPatchInfo * const acceptedPatch){
  /*
   * Add accepted patch to the trigger decision
   *
   * @param patch: the accepted patch
   */
  fAcceptedPatches.Add(acceptedPatch);
}
