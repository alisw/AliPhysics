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
 * Task providing an event selection for EMCAL-triggered events based on the
 * reconstructed EMCAL trigger patches
 *
 * Author: Markus Fasel
 */

#include "AliEmcalTriggerDecision.h"
#include "AliEmcalTriggerDecisionContainer.h"
#include "AliEmcalTriggerSelection.h"
#include "AliAnalysisTaskEmcalTriggerSelection.h"

ClassImp(AliAnalysisTaskEmcalTriggerSelection)

//______________________________________________________________________________
AliAnalysisTaskEmcalTriggerSelection::AliAnalysisTaskEmcalTriggerSelection():
  AliAnalysisTaskEmcal(),
  fGlobalDecisionContainerName(),
  fTriggerSelections()
{
  /*
   * Dummy constructor, only for I/O, not to be used by the user
   */
  fTriggerSelections.SetOwner(kTRUE);
}

//______________________________________________________________________________
AliAnalysisTaskEmcalTriggerSelection::AliAnalysisTaskEmcalTriggerSelection(const char* name):
  AliAnalysisTaskEmcal(name, kFALSE),
  fGlobalDecisionContainerName(),
  fTriggerSelections()
{
  /*
   * Main constructor, to be called by the users
   */
  fTriggerSelections.SetOwner(kTRUE);
}

//______________________________________________________________________________
void AliAnalysisTaskEmcalTriggerSelection::AddTriggerSelection(AliEmcalTriggerSelection * const selection){
  /*
   * Add trigger selection to the trigger selection task
   *
   * @param selection: the trigger selection to be added
   */
  fTriggerSelections.Add(selection);
}

//______________________________________________________________________________
Bool_t AliAnalysisTaskEmcalTriggerSelection::Run(){
  /*
   * Run over all trigger selections, and append the selection to the global trigger selection container
   */
  AliEmcalTriggerDecisionContainer *cont = GetGlobalTriggerDecisionContainer();
  cont->Reset();
  TClonesArray *triggerPatches(fTriggerPatchInfo);
  TIter selectionIter(&fTriggerSelections);
  AliEmcalTriggerSelection *selection(NULL);
  while((selection = dynamic_cast<AliEmcalTriggerSelection *>(selectionIter()))){
    cont->AddTriggerDecision(selection->MakeDecison(triggerPatches));
  }
  return kTRUE;
}

//______________________________________________________________________________
AliEmcalTriggerDecisionContainer *AliAnalysisTaskEmcalTriggerSelection::GetGlobalTriggerDecisionContainer(){
  /*
   * Find the main trigger container in the input event. If not available, create it and add it to the input event
   */
  AliEmcalTriggerDecisionContainer *cont(dynamic_cast<AliEmcalTriggerDecisionContainer *>(fInputEvent->FindListObject(fGlobalDecisionContainerName.Data())));
  if(!cont){
    cont = new AliEmcalTriggerDecisionContainer(fGlobalDecisionContainerName.Data());
    fInputEvent->AddObject(cont);
  }
  return cont;
}
