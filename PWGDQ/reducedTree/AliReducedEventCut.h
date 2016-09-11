// Class for cutting on ALICE Var manager and other event specific information
// Author: Ionut-Cristian Arsene (iarsene@cern.ch)
//   07/09/2016

#ifndef ALIREDUCEDEVENTCUT_H
#define ALIREDUCEDEVENTCUT_H

#include "AliReducedVarCut.h"
#include "AliReducedVarManager.h"

//_________________________________________________________________________
class AliReducedEventCut : public AliReducedVarCut {

 public:
  AliReducedEventCut();
  AliReducedEventCut(const Char_t* name, const Char_t* title);
  virtual ~AliReducedEventCut();

  virtual Bool_t IsSelected(TObject* obj);
  virtual Bool_t IsSelected(TObject* obj, Float_t* values);
  
 protected: 
      
  // Cuts on event specific quantities
   
  AliReducedEventCut(const AliReducedEventCut &c);
  AliReducedEventCut& operator= (const AliReducedEventCut &c);
  
  ClassDef(AliReducedEventCut,1);
};

#endif
