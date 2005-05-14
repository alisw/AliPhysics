#ifndef STARTV2_H
#define STARTV2_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager and hits classes for set:START     //
////////////////////////////////////////////////
 
#include "AliSTART.h"
 
class AliSTARTv2 : public AliSTART {
  
public:
  AliSTARTv2() {};
  AliSTARTv2(const char *name, const char *title);
  virtual       ~AliSTARTv2() {}
  virtual void   CreateGeometry();
  virtual void   CreateMaterials();
  virtual void   DrawModule() const;
  virtual void   Init();
  virtual Int_t  IsVersion() const {return 0;}
  virtual void   StepManager();
  virtual void SetHitsAddressBranch(TBranch *b1,TBranch *b2)
    {b1->SetAddress(&fHits); b2->SetAddress(&fPhotons);}
   
protected:
   Int_t fIdSens1; // Sensetive volume  in START
 
  ClassDef(AliSTARTv2, 1)  //Class for START version 0
};

#endif


