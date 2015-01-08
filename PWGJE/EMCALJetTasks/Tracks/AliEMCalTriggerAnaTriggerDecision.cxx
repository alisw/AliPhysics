/**************************************************************************
 * Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
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
 * Class performing the selection of triggered events
 *
 * Author:
 *    Markus Fasel
 */
#include <TClonesArray.h>
#include <TString.h>
#include "AliVEvent.h"
#include "AliEmcalTriggerPatchInfo.h"
#include "AliEMCalTriggerAnaTriggerDecision.h"
#include "AliEMCalTriggerEventData.h"

ClassImp(EMCalTriggerPtAnalysis::AliEMCalTriggerAnaTriggerDecision)

namespace EMCalTriggerPtAnalysis {

//______________________________________________________________________________
AliEMCalTriggerAnaTriggerDecision::AliEMCalTriggerAnaTriggerDecision() :
    fSwapThresholds(kFALSE),
    fIsMinBias(kFALSE)
{
  /*
   * Main constructor
   */
  Reset();
}

//______________________________________________________________________________
void AliEMCalTriggerAnaTriggerDecision::Create(const AliEMCalTriggerEventData* const data) {
  /*
   * Steer creation of the trigger decision
   *
   * @param data: all event information
   */
  Reset();
  MakeDecisionFromPatches(*(data->GetTriggerPatchContainer()));
  MakeDecisionFromString(data->GetRecEvent()->GetFiredTriggerClasses());
}

//______________________________________________________________________________
void AliEMCalTriggerAnaTriggerDecision::Reset() {
  for(int itrg = 0; itrg < 4; itrg++){
    fDecisionFromPatches[itrg] = kFALSE;
    fDecisionFromString[itrg] = kFALSE;
  }
}

//______________________________________________________________________________
void AliEMCalTriggerAnaTriggerDecision::MakeDecisionFromString(const TString& triggerstring) {
  /*
   * Create trigger decision from trigger string
   *
   * @param triggerstring: the trigger string
   */
  if(triggerstring.Contains("EJ1")) fDecisionFromString[kTAEMCJHigh]  = kTRUE;
  if(triggerstring.Contains("EJ2")) fDecisionFromString[kTAEMCJLow]   = kTRUE;
  if(triggerstring.Contains("EG1")) fDecisionFromString[kTAEMCGHigh]  = kTRUE;
  if(triggerstring.Contains("EG2")) fDecisionFromString[kTAEMCGLow]   = kTRUE;
}

//______________________________________________________________________________
void AliEMCalTriggerAnaTriggerDecision::MakeDecisionFromPatches(const TClonesArray& listOfPatches) {
  /*
   * Create trigger decision from trigger patches. In case swap thresholds is requested, the low threshold
   * triggers are replaced by the high threshold triggers and vice versa
   *
   * @param triggerstring: the TClonesArray of the trigger patches, created by the trigger patch maker
   */
  TIter patchIter(&listOfPatches);
  AliEmcalTriggerPatchInfo *mypatch(NULL);
  while((mypatch = dynamic_cast<AliEmcalTriggerPatchInfo *>(patchIter()))){
    if(mypatch->IsJetHigh()) fDecisionFromPatches[fSwapThresholds ? kTAEMCJLow : kTAEMCJHigh] = kTRUE;
    if(mypatch->IsJetLow()) fDecisionFromPatches[fSwapThresholds ? kTAEMCJHigh : kTAEMCJLow] = kTRUE;
    if(mypatch->IsGammaHigh()) fDecisionFromPatches[fSwapThresholds ? kTAEMCGLow : kTAEMCGHigh] = kTRUE;
    if(mypatch->IsGammaLow()) fDecisionFromPatches[fSwapThresholds ? kTAEMCGHigh : kTAEMCJLow] = kTRUE;
  }
}

} /* namespace EMCalTriggerPtAnalysis */
