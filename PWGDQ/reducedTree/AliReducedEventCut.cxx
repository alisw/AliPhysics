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
#include "AliReducedVarManager.h"

ClassImp(AliReducedEventCut)

//____________________________________________________________________________
AliReducedEventCut::AliReducedEventCut() :
  AliReducedVarCut()
{
  //
  // default constructor
  //
}

//____________________________________________________________________________
AliReducedEventCut::AliReducedEventCut(const Char_t* name, const Char_t* title) :
  AliReducedVarCut(name, title)
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
   
   return AliReducedVarCut::IsSelected(obj, values);   
}
