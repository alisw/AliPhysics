#ifndef ALIPHOSQAOBJECTCHECKABLE_H
#define ALIPHOSQAOBJECTCHECKABLE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
// Abstact Class for a QA checkable that is a TObject    
//                  
//*-- Author: Yves Schutz (SUBATECH)


// --- ROOT system ---

// --- Standard library ---

// --- AliRoot header files ---

#include "AliPHOSQAVirtualCheckable.h"

class AliPHOSQAObjectCheckable : public AliPHOSQAVirtualCheckable {

public:

  AliPHOSQAObjectCheckable(){
    fObject = 0;
  }           // default ctor not to be used
  AliPHOSQAObjectCheckable(const char * name) ;          // ctor
  AliPHOSQAObjectCheckable(AliPHOSQAObjectCheckable& obj) : AliPHOSQAVirtualCheckable(obj)
  {assert(0==1);}
  virtual ~AliPHOSQAObjectCheckable() ; // dtor

  virtual TObject * GetObject() const { return fObject ; }
  virtual Float_t GetValue() const {return 0. ;} 
  virtual void Print() const ; 
  virtual void Reset() { fChange=kFALSE ; }
  virtual void Set(TObject * obj) {fObject = obj ;} 
  virtual void Update(TObject * value) {} ; 

private:
  
  TObject *  fObject ; 

  ClassDef(AliPHOSQAObjectCheckable,1)  // description 

};

#endif // ALIPHOSQAObjectCheckable_H
