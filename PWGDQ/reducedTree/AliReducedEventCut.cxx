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

#include <vector>
#include "AliReducedBaseEvent.h"
#include "AliReducedEventInfo.h"
#include "AliReducedVarManager.h"

ClassImp(AliReducedEventCut)

//____________________________________________________________________________
AliReducedEventCut::AliReducedEventCut() :
  AliReducedVarCut(),
  fEventTagFilterEnabled(kFALSE),
  fEventTagFilterExclude(0),
  fEventFilter(0),
  fEventTriggerMaskEnabled(kFALSE), 
  fEventTriggerMask(0),
  fEventTriggerClassEnabled(kFALSE),
  fEventTriggerClassLogicalOr(kFALSE),
  fEventTriggerClass(),
  fEventL1MaskEnabled(kFALSE),
  fEventL1Mask(0),
  fEventL0MaskEnabled(kFALSE),
  fEventL0Mask(0)
{
  //
  // default constructor
  //
}

//____________________________________________________________________________
AliReducedEventCut::AliReducedEventCut(const Char_t* name, const Char_t* title) :
  AliReducedVarCut(name, title),
  fEventTagFilterEnabled(kFALSE),
  fEventTagFilterExclude(0),
  fEventFilter(0),
  fEventTriggerMaskEnabled(kFALSE), 
  fEventTriggerMask(0),
  fEventTriggerClassEnabled(kFALSE),
  fEventTriggerClassLogicalOr(kFALSE),
  fEventTriggerClass(),
  fEventL1MaskEnabled(kFALSE),
  fEventL1Mask(0),
  fEventL0MaskEnabled(kFALSE),
  fEventL0Mask(0)
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
   if (fEventTagFilterEnabled && (fEventTagFilterExclude & fEventFilter) && (fEventTagFilterExclude & (event->EventTag() & fEventFilter))) return kFALSE;  // exclusion selection
   if (fEventTagFilterEnabled && !(fEventTagFilterExclude & fEventFilter) && !(event->EventTag() & fEventFilter)) return kFALSE;                           // inclusion selection

   if(fEventTriggerMaskEnabled) {
     if(!obj->InheritsFrom(AliReducedEventInfo::Class())) return kFALSE;
     AliReducedEventInfo* eventInfo = (AliReducedEventInfo*)obj;
     if(!(eventInfo->TriggerMask() & fEventTriggerMask)) return kFALSE;
   }

  if (fEventTriggerClassEnabled) {
    if(!obj->InheritsFrom(AliReducedEventInfo::Class())) return kFALSE;
    AliReducedEventInfo* eventInfo = (AliReducedEventInfo*)obj;
    TString trgClasses = eventInfo->TriggerClass();
    Int_t counter = 0;
    for (Int_t i=0; i<fEventTriggerClass.size(); ++i) {
      if (trgClasses.Contains(fEventTriggerClass[i].Data())) counter++;
    }
    if (fEventTriggerClassLogicalOr && !counter) return kFALSE;
    if (!fEventTriggerClassLogicalOr && counter!=fEventTriggerClass.size()) return kFALSE;
  }

   if(fEventL1MaskEnabled) {
      if(!obj->InheritsFrom(AliReducedEventInfo::Class())) return kFALSE;
      AliReducedEventInfo* eventInfo = (AliReducedEventInfo*)obj;
      if(!(eventInfo->L1TriggerInputs() & fEventL1Mask)) return kFALSE;
   }
   
   if(fEventL0MaskEnabled) {
      if(!obj->InheritsFrom(AliReducedEventInfo::Class())) return kFALSE;
      AliReducedEventInfo* eventInfo = (AliReducedEventInfo*)obj;
      if(!(eventInfo->L0TriggerInputs() & fEventL0Mask)) return kFALSE;
   }
   
   return AliReducedVarCut::IsSelected(values);   
}
