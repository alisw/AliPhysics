/*
***********************************************************
  Implementation of AliReducedEventCut class.
  Contact: iarsene@cern.ch
  2016/09/07
  *********************************************************
*/

#ifndef ALIREDUCEDEVENTCUT_H
#include "AliReducedEventCut.h"
#endif

#include "AliReducedBaseEvent.h"
#include "AliReducedEventInfo.h"
#include "AliReducedVarManager.h"

ClassImp(AliReducedEventCut)

//____________________________________________________________________________
AliReducedEventCut::AliReducedEventCut() :
  AliReducedVarCut(),
  fEventTagFilterEnabled(kFALSE),
  fEventFilter(0),
  fEventTriggerMaskEnabled(kFALSE), 
  fEventTriggerMask(0)
{
  //
  // default constructor
  //
}

//____________________________________________________________________________
AliReducedEventCut::AliReducedEventCut(const Char_t* name, const Char_t* title) :
  AliReducedVarCut(name, title),
  fEventTagFilterEnabled(kFALSE),
  fEventFilter(0),
  fEventTriggerMaskEnabled(kFALSE), 
  fEventTriggerMask(0)
{
  //
  // named constructor
  //
}

//____________________________________________________________________________
AliReducedEventCut::~AliReducedEventCut() {
  //
  // destructor
  //
}

//____________________________________________________________________________
Bool_t AliReducedEventCut::IsSelected(TObject* obj) {
  //
  // apply cuts
  //
  if(!obj->InheritsFrom(AliReducedBaseEvent::Class())) return kFALSE;
  
  //Fill values
  Float_t values[AliReducedVarManager::kNVars];
  AliReducedVarManager::FillEventInfo((AliReducedBaseEvent*)obj, values);
  
  return IsSelected(obj, values);
}


//____________________________________________________________________________
Bool_t AliReducedEventCut::IsSelected(TObject* obj, Float_t* values) {
   //
   // apply cuts
   //      
   if(!obj->InheritsFrom(AliReducedBaseEvent::Class())) return kFALSE;
   
   AliReducedBaseEvent* event = (AliReducedBaseEvent*)obj;
   if(fEventTagFilterEnabled && !(event->EventTag() & fEventFilter)) return kFALSE;

   if(fEventTriggerMaskEnabled) {
     if(!obj->InheritsFrom(AliReducedEventInfo::Class())) return kFALSE;
     AliReducedEventInfo* eventInfo = (AliReducedEventInfo*)obj;
     if(fEventTriggerMaskEnabled && !(eventInfo->TriggerMask() & fEventTriggerMask)) return kFALSE;
   }
   
   return AliReducedVarCut::IsSelected(values);   
}
