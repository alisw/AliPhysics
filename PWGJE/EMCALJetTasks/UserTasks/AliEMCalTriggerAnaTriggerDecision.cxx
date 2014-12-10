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
