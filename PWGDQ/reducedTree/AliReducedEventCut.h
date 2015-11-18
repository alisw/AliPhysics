// Class for cutting on event information
// Author: Ionut-Cristian Arsene (iarsene@cern.ch)
//   10/09/2015

#ifndef ALIREDUCEDEVENTCUT_H
#define ALIREDUCEDEVENTCUT_H

#include "AliReducedInfoCut.h"

//_________________________________________________________________________
class AliReducedEventCut : public AliReducedInfoCut {

 public:
  AliReducedEventCut();
  AliReducedEventCut(const Char_t* name, const Char_t* title);
  virtual ~AliReducedEventCut();
  
  virtual Bool_t IsSelected(TObject* obj);
  
  void SetEventTagMask(ULong64_t mask) {fEventTagMask = mask; fCutOnEventTag=kTRUE;}
  void SetTriggerMask(ULong64_t mask) {fTriggerMask = mask; fCutOnTriggerMask=kTRUE;}
  void SetVertexZRange(Float_t min, Float_t max) {fVertexZRange[0]=min; fVertexZRange[1]=max; fCutOnVertexZ=kTRUE;}
  void SetCentVZERORange(Float_t min, Float_t max) {fCentVZERORange[0]=min; fCentVZERORange[1]=max; fCutOnCentralityVZERO = kTRUE;}
  void SetCutOnPhysicsSelection(Bool_t cut) {fCutOnPhysicsSelection=cut;}
  
 protected: 
   
  ULong64_t fEventTagMask;
  Bool_t    fCutOnEventTag;
  
  ULong64_t fTriggerMask;
  Bool_t    fCutOnTriggerMask;
  
  Float_t   fVertexZRange[2];
  Bool_t    fCutOnVertexZ;
  
  Float_t   fCentVZERORange[2];
  Bool_t    fCutOnCentralityVZERO;
  
  Bool_t    fCutOnPhysicsSelection;
  
  
  AliReducedEventCut(const AliReducedEventCut &c);
  AliReducedEventCut& operator= (const AliReducedEventCut &c);
  
  ClassDef(AliReducedEventCut,1);
};

#endif