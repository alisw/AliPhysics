#ifndef ALIDETECTORRECOPARAM_H
#define ALIDETECTORRECOPARAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Base Class for Detector reconstruction parameters                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "TNamed.h"

class AliDetectorRecoParam : public TNamed
{
  
 public: 
  AliDetectorRecoParam();
  virtual ~AliDetectorRecoParam();
  void  Print(Option_t *option) const {Dump();}
  const Int_t * GetEventType() { return fEventType;}
protected:
  Int_t fEventType[5]; // Reconstruction - event type        
  
  ClassDef(AliDetectorRecoParam, 1)
};


#endif
