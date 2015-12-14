/*
***********************************************************
  Implementation of AliReducedEventCut class.
  Contact: iarsene@cern.ch
  2015/09/10
  *********************************************************
*/

#ifndef ALIREDUCEDEVENTCUT_H
#include "AliReducedEventCut.h"
#endif

#include "AliReducedEventInfo.h"

ClassImp(AliReducedEventCut)

//____________________________________________________________________________
AliReducedEventCut::AliReducedEventCut() :
  AliReducedInfoCut(),
  fEventTagMask(0),
  fCutOnEventTag(kFALSE),  
  fTriggerMask(0),
  fCutOnTriggerMask(kFALSE),
  fVertexZRange(),
  fCutOnVertexZ(kFALSE),
  fCentVZERORange(),
  fCutOnCentralityVZERO(kFALSE),
  fCutOnPhysicsSelection(kFALSE)
{
  //
  // default constructor
  //
}

//____________________________________________________________________________
AliReducedEventCut::AliReducedEventCut(const Char_t* name, const Char_t* title) :
  AliReducedInfoCut(name, title),
  fEventTagMask(0),
  fCutOnEventTag(kFALSE),  
  fTriggerMask(0),
  fCutOnTriggerMask(kFALSE),
  fVertexZRange(),
  fCutOnVertexZ(kFALSE),
  fCentVZERORange(),
  fCutOnCentralityVZERO(kFALSE),
  fCutOnPhysicsSelection(kFALSE)
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
   AliReducedBaseEvent* baseEv = (AliReducedBaseEvent*)obj;
   
  AliReducedEventInfo* ev = NULL;
  if(obj->IsA()==AliReducedEventInfo::Class()) ev = (AliReducedEventInfo*)obj;
    
  if(fCutOnEventTag && !(fEventTagMask&baseEv->EventTag())) return kFALSE;
  if(fCutOnVertexZ && (baseEv->Vertex(2)<fVertexZRange[0] || baseEv->Vertex(2)>fVertexZRange[1])) return kFALSE;
  if(fCutOnCentralityVZERO && (baseEv->CentralityVZERO()<fCentVZERORange[0] || baseEv->CentralityVZERO()>fCentVZERORange[1])) return kFALSE;
  if(ev && fCutOnTriggerMask && !(fTriggerMask&ev->TriggerMask())) return kFALSE;
  if(ev && fCutOnPhysicsSelection && !ev->IsPhysicsSelection()) return kFALSE;
  
  return kTRUE;
}

