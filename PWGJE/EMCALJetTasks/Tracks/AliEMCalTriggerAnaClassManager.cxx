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

ClassImp(PWGJE::EMCALJetTasks::AliEMCalTriggerAnaClassManager)

using namespace PWGJE::EMCALJetTasks;

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
 * Copy constructor
 * @param ref Reference for the copy
 */
AliEMCalTriggerAnaClassManager::AliEMCalTriggerAnaClassManager(const AliEMCalTriggerAnaClassManager &ref):
  TNamed(ref),
  fTriggerClasses(NULL),
  fSelected(NULL)
{
  fTriggerClasses = new TObjArray;
  fTriggerClasses->SetOwner(kFALSE);
  // Copy trigger classes - only pointers
  for(TIter classiter = TIter(ref.fTriggerClasses).Begin(); classiter != TIter::End(); ++classiter){
    fTriggerClasses->Add(*classiter);
  }
  fSelected = new TObjArray;
  fSelected->SetOwner(kFALSE);
}

/**
 * Assignment operator
 * @param ref Reference for assignment
 * @return Trigger class manager after assignment
 */
AliEMCalTriggerAnaClassManager &AliEMCalTriggerAnaClassManager::operator=(const AliEMCalTriggerAnaClassManager &ref){
  TNamed::operator=(ref);
  if(this != &ref){
    fSelected->Clear();
    fTriggerClasses->Clear();
    // Copy trigger classes - only pointers
    for(TIter classiter = TIter(ref.fTriggerClasses).Begin(); classiter != TIter::End(); ++classiter){
      fTriggerClasses->Add(*classiter);
    }
  }
  return *this;
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
 * \throw TriggerManagerEmptyException
 */
void AliEMCalTriggerAnaClassManager::PerformEventSelection(AliEMCalTriggerEventData* trgevent) {
  fSelected->Clear();
  if(!fTriggerClasses->GetEntries()) throw TriggerManagerEmptyException();
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

/**
 * Get the list of selected trigger classes. Event selection has to be performed before.
 * \return The list of selected trigger classes.
 * \throw TriggerManagerEmptyException
 */
TObjArray * AliEMCalTriggerAnaClassManager::GetSelectedTriggerClasses() const {
  if(!fTriggerClasses->GetEntries()) throw TriggerManagerEmptyException();
  return fSelected;
}
/**
 * Get list of all trigger classes
 * \return List of all trigger classes
 * \throw TriggerManagerEmptyException
 */
TObjArray * AliEMCalTriggerAnaClassManager::GetAllTriggerClasses() const {
  if(!fTriggerClasses->GetEntries()) throw TriggerManagerEmptyException();
  return fTriggerClasses;
}
