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
 * Object performing an offline EMCAL trigger decision based on user defined criterions
 * (trigger patch type, energy threshold,...). The main method MakeTriggerDecision performs
 * an event selection and creates a trigger decision object with the relevant information.
 *
 * Author: Markus Fasel
 */
#include <vector>
#include <TClonesArray.h>

#include "AliEMCALTriggerPatchInfo.h"
#include "AliEmcalTriggerDecision.h"
#include "AliEmcalTriggerSelection.h"
#include "AliEmcalTriggerSelectionCuts.h"

ClassImp(AliEmcalTriggerSelection)

//______________________________________________________________________________
AliEmcalTriggerSelection::AliEmcalTriggerSelection() :
  TNamed(),
  fSelectionCuts(NULL),
  fOutputName("")
{
  /*
   * Dummy constructor, used by I/O, not to be used by the user
   */
}

//______________________________________________________________________________
AliEmcalTriggerSelection::AliEmcalTriggerSelection(const char *name, const AliEmcalTriggerSelectionCuts * const cuts):
  TNamed(name, ""),
  fSelectionCuts(cuts),
  fOutputName("")
{
  /*
   * Main constructor, initialises the trigger selection
   *
   * @param name: name of the trigger selection
   * @param cuts(optional): selection cuts to be applied
   */
}

//______________________________________________________________________________
AliEmcalTriggerDecision* AliEmcalTriggerSelection::MakeDecison(const TClonesArray * const inputPatches) const {
  /*
   * Perform event selection based on user-defined criteria and create an output trigger decision containing
   * the threshold, the main patch which fired the decision, and all other patches which would have fired the
   * decision as well.
   *
   * @param input patches: A list of input patches, created by the trigger patch maker and read out from the
   * input event
   * @return: the trigger decision (an event is selected when it has a main patch that fired the decision)
   */
  AliEmcalTriggerDecision *result = new AliEmcalTriggerDecision(fOutputName.Data());
  TIter patchIter(inputPatches);
  AliEMCALTriggerPatchInfo *patch(NULL);
  std::vector<AliEMCALTriggerPatchInfo *> selectedPatches;
  while((patch = dynamic_cast<AliEMCALTriggerPatchInfo *>(patchIter()))){
    if(fSelectionCuts->IsSelected(patch)){
      selectedPatches.push_back(patch);
    }
  }
  // Find the main patch
  AliEMCALTriggerPatchInfo *mainPatch(NULL), *testpatch(NULL);
  for(std::vector<AliEMCALTriggerPatchInfo *>::iterator it = selectedPatches.begin(); it != selectedPatches.end(); ++it){
    testpatch = *it;
    if(!mainPatch) mainPatch = testpatch;
    else if(fSelectionCuts->CompareTriggerPatches(testpatch, mainPatch) > 0) mainPatch = testpatch;
    result->AddAcceptedPatch(testpatch);
  }
  if(mainPatch) result->SetMainPatch(mainPatch);
  result->SetSelectionCuts(fSelectionCuts);
  return result;
}
