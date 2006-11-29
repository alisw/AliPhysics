#ifndef T0V2_H
#define T0V2_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager and hits classes for set:T0     //
////////////////////////////////////////////////
 
#include "AliT0.h"
 
class AliT0v2 : public AliT0 {
  
public:
  AliT0v2() {};
  AliT0v2(const char *name, const char *title);
  virtual       ~AliT0v2() {}
  virtual void   CreateGeometry();
  virtual void   CreateMaterials();
  virtual void   DrawModule() const;
  virtual void   Init();
  virtual Int_t  IsVersion() const {return 0;}
  virtual void   StepManager();
  virtual void SetHitsAddressBranch(TBranch *b1,TBranch *b2)
    {b1->SetAddress(&fHits); b2->SetAddress(&fPhotons);}
   
protected:
   Int_t fIdSens1; // Sensetive volume  in T0
 
  ClassDef(AliT0v2, 1)  //Class for T0 version 0
};

#endif


