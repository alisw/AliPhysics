#ifndef ALIPHOSQAINTCHECKABLE_H
#define ALIPHOSQAINTCHECKABLE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
// Class for a QA checkable that is an Int    
//                  
//*-- Author: Yves Schutz (SUBATECH)


// --- ROOT system ---

// --- Standard library ---

// --- AliRoot header files ---

#include "AliPHOSQAVirtualCheckable.h"

class AliPHOSQAIntCheckable : public AliPHOSQAVirtualCheckable {

public:

  AliPHOSQAIntCheckable(){}           // default ctor not to be used
  AliPHOSQAIntCheckable(const char * name) ;          // ctor
  AliPHOSQAIntCheckable(AliPHOSQAIntCheckable& obj)
    : AliPHOSQAVirtualCheckable(obj) {assert(0==1);}
  virtual ~AliPHOSQAIntCheckable() ; // dtor

  virtual Float_t GetValue() const { return (Float_t)fValue ; }
  virtual  void Print() const ; 
  virtual void Reset() { fValue=0; fChange=kFALSE ; }
  void Set(Int_t value) ; 
  void Update(Int_t value) ;  

private:
  
  Int_t fValue ; 

  ClassDef(AliPHOSQAIntCheckable,1)  // description 

};

#endif // ALIPHOSQAIntCheckable_H
