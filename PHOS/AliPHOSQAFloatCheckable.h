#ifndef ALIPHOSQAFLOATCHECKABLE_H
#define ALIPHOSQAFLOATCHECKABLE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
// Class for a QA checkable that is a Float_t    
//                  
//*-- Author: Yves Schutz (SUBATECH)


// --- ROOT system ---

// --- Standard library ---

// --- AliRoot header files ---

#include "AliPHOSQAVirtualCheckable.h"

class AliPHOSQAFloatCheckable : public AliPHOSQAVirtualCheckable {

public:

  AliPHOSQAFloatCheckable(){}           // default ctor not to be used
  AliPHOSQAFloatCheckable(const char * name) ;          // ctor
  AliPHOSQAFloatCheckable(AliPHOSQAFloatCheckable& obj)
    : AliPHOSQAVirtualCheckable(obj) {assert(0==1);}
  virtual ~AliPHOSQAFloatCheckable() ; // dtor

  virtual Float_t GetValue() const { return fValue ; }
  virtual  void Print() const ; 
  virtual void Reset() { fValue=0.; fChange=kFALSE ; }
  void Set(Float_t value) ; 
  void Update(Float_t value) ;

private:
  
  Float_t fValue ; 

  ClassDef(AliPHOSQAFloatCheckable,1)  // description 

};

#endif // ALIPHOSQAFloatCheckable_H
