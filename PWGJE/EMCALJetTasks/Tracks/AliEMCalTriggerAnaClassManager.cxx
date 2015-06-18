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
#include <TObjArray.h>

#include "AliEMCalTriggerAnaClassManager.h"
#include "AliEMCalTriggerAnaTriggerClass.h"
#include "AliEMCalTriggerAnaTriggerDecision.h"

/// \cond CLASSIMP
ClassImp(EMCalTriggerPtAnalysis::AliEMCalTriggerAnaClassManager)
/// \endcond

namespace EMCalTriggerPtAnalysis {

/**
 * Dummy constructor
 */
AliEMCalTriggerAnaClassManager::AliEMCalTriggerAnaClassManager():
  TNamed(),
  fTriggerClasses(NULL),
  fSelected(NULL)
{
}

/**
 * Named constructor, creating also arrays
 * \param name Name of the handler
 */
AliEMCalTriggerAnaClassManager::AliEMCalTriggerAnaClassManager(const char* name):
  TNamed(name, ""),
  fTriggerClasses(NULL),
  fSelected(NULL)
{
  fTriggerClasses = new TObjArray;
  fTriggerClasses->SetOwner(kTRUE);
  fSelected = new TObjArray;
  fSelected->SetOwner(kFALSE);
}

/**
 * Destructor
 */
AliEMCalTriggerAnaClassManager::~AliEMCalTriggerAnaClassManager() {
  if(fTriggerClasses) delete fTriggerClasses;
  if(fSelected) delete fSelected;
}

/**
 * For each trigger class test whether event is selected for the class and mark as selected
 * \param trgevent The event data to check.
 */
void AliEMCalTriggerAnaClassManager::PerformEventSelection(AliEMCalTriggerEventData* trgevent) {
  fSelected->Clear();
  for(TIter clsiter = TIter(fTriggerClasses).Begin(); clsiter != TIter::End(); ++clsiter){
    AliEMCalTriggerAnaTriggerClass * myclass = static_cast<AliEMCalTriggerAnaTriggerClass *>(*clsiter);
    if(myclass->IsEventTriggered(trgevent)) fSelected->Add(myclass);
  }
}

/**
 * Add new trigger class to the manager
 * \param triggerclass
 */
void AliEMCalTriggerAnaClassManager::AddTriggerClass(AliEMCalTriggerAnaTriggerClass* triggerclass) {
  fTriggerClasses->Add(triggerclass);
}
/**
 * Forward trigger decision handler to all trigger classes
 * \param triggerdecision The trigger decision for the given event
 */
void AliEMCalTriggerAnaClassManager::SetTriggerDecision(AliEMCalTriggerAnaTriggerDecision* triggerdecision) {
  for(TIter clsiter = TIter(fTriggerClasses).Begin(); clsiter != TIter::End(); ++clsiter){
    (static_cast<AliEMCalTriggerAnaTriggerClass *>(*clsiter))->SetTriggerDecisionHandler(triggerdecision);
  }
}

/**
 * Check whether among the selected trigger classes we find at least one minimum bias trigger
 * \return True if among the selected classes we find a minimum bias trigger
 */
bool AliEMCalTriggerAnaClassManager::HasMinBiasTrigger() const {
  bool result = kFALSE;
  for(TIter clsiter = TIter(fSelected).Begin(); clsiter != TIter::End(); ++clsiter){
    if((static_cast<AliEMCalTriggerAnaTriggerClass *>(*clsiter))->IsMinBiasTrigger()){
      result = kTRUE;
      break;
    }
  }
  return kTRUE;
}

} /* namespace EMCalTriggerPtAnalysis */
