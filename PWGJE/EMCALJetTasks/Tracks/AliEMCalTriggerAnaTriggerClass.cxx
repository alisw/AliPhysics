/**************************************************************************
 * Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
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
#include "AliLog.h"
#include "AliVEvent.h"
#include "AliEMCalTriggerAnaTriggerClass.h"
#include "AliEMCalTriggerAnaTriggerDecision.h"
#include "AliEMCalTriggerEventData.h"

ClassImp(PWGJE::EMCALJetTasks::AliEMCalTriggerAnaTriggerClass)
ClassImp(PWGJE::EMCALJetTasks::AliEMCalTriggerAnaPatternObject)
ClassImp(PWGJE::EMCALJetTasks::AliEMCalTriggerAnaPatternContainer)
ClassImp(PWGJE::EMCALJetTasks::AliEMCalTriggerAnaTriggerPatchTypeObject)

using namespace PWGJE::EMCALJetTasks;

/**
 * Dummy (I/O) constructor - not to be used
 */
AliEMCalTriggerAnaTriggerClass::AliEMCalTriggerAnaTriggerClass() :
    TNamed(),
    fDecisionFromTriggerBits(kFALSE),
    fDecisionFromTriggerString(kFALSE),
    fDecisionFromTriggerPatches(kFALSE),
    fIsMinBiasTrigger(kFALSE),
    fTriggerBits(0),
    fTriggerStringPattern(),
    fEmcalTriggerHandler(NULL)
{
}

/**
 * Named constructor - Defines the trigger class with a name and a title
 * \param name Name of the trigger class
 * \param title A short description
 */
AliEMCalTriggerAnaTriggerClass::AliEMCalTriggerAnaTriggerClass(const char *name, const char *title) :
    TNamed(name, title),
    fDecisionFromTriggerBits(kFALSE),
    fDecisionFromTriggerString(kFALSE),
    fDecisionFromTriggerPatches(kFALSE),
    fIsMinBiasTrigger(kFALSE),
    fTriggerBits(0),
    fTriggerStringPattern(),
    fTriggerPatchTypes(),
    fEmcalTriggerHandler(NULL)
{
  fTriggerPatchTypes.SetOwner(kTRUE);
}

/**
 * Destructor - nothing to do
 */
AliEMCalTriggerAnaTriggerClass::~AliEMCalTriggerAnaTriggerClass()  {}

/**
 * Selection of events according to the trigger class. In case any condition fails, the others from that time on are not checked anymore.
 * \param triggerevnet The event data to check
 * \return True if the event is selected for this trigger class, false otherwise
 * \throw TriggerMethodUndefinedException in case no method to select events is defined
 * \throw EventCorruptionException if the reconstructed event is missing or the trigger patch container is missing
 * \throw PatchHandlerMissingException if the trigger patch handler is not set
 */
bool AliEMCalTriggerAnaTriggerClass::IsEventTriggered(const AliEMCalTriggerEventData *const triggerevent) const{
  if(!(fDecisionFromTriggerBits || fDecisionFromTriggerString || fDecisionFromTriggerPatches))
    throw TriggerMethodUndefinedException(this->GetName());
  bool result = kTRUE;

  if(fDecisionFromTriggerBits){
    result = result && (triggerevent->GetTriggerBitSelection() & fTriggerBits);
  }
  if(!result) return kFALSE;

  if(fDecisionFromTriggerString){
    if(!triggerevent->GetRecEvent())
      throw EventCorruptionException(this->GetName(), "Reconstructed event missing");
    result = result && fTriggerStringPattern.CheckTriggerString(triggerevent->GetRecEvent()->GetFiredTriggerClasses().Data());
  }
  if(!result) return kFALSE;

  if(fDecisionFromTriggerPatches){
    if(fEmcalTriggerHandler){
      if(triggerevent->GetTriggerPatchContainer())
        throw EventCorruptionException(this->GetName(), "Trigger patch container missing");
      for(TIter typeiter = TIter(&fTriggerPatchTypes).Begin(); typeiter != TIter::End(); ++typeiter){
        AliEMCalTriggerAnaTriggerPatchTypeObject *patchtype = static_cast<AliEMCalTriggerAnaTriggerPatchTypeObject *>(*typeiter);
        result = result && fEmcalTriggerHandler->IsTriggered(patchtype->GetTriggerType(), kTriggerPatches);
      }
    } else
      throw PatchHandlerMissingException(this->GetName());
  }

  return result;
}

/**
 * Match patter in the trigger string
 * \param triggerstring The trigger string to check
 * \return True if the pattern is requested and found or not requested and not found, false otherwise
 */
Bool_t AliEMCalTriggerAnaPatternObject::MatchTriggerString(const char *triggerstring) const {
  Bool_t patternmatch = TString(triggerstring).Contains(fPattern);
  return fInString ? patternmatch : !patternmatch;
}

/**
 * Check triggerstring for all patterns defined in the container.
 * \param triggerstring The triggerstring to check
 * \return True if all patterns are correctly matched, false otherwise
 */
Bool_t AliEMCalTriggerAnaPatternContainer::CheckTriggerString(const char *triggerstring) const{
  Bool_t result = kTRUE;

  for(TIter piter = TIter(&fPatterns).Begin(); piter != TIter::End(); ++piter){
    const AliEMCalTriggerAnaPatternObject *toCheck = static_cast<const AliEMCalTriggerAnaPatternObject *>(*piter);
    result = result && toCheck->MatchTriggerString(triggerstring);
  }
  return result;
}
