#ifndef VZEROv0_H
#define VZEROv0_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


///////////////////////////////////////////////////
//  Manager and hits classes for set : VZERO     //
///////////////////////////////////////////////////
 
#include "AliVZERO.h"
#include "TFile.h"
#include "TH1.h"

class AliVZEROv0 : public AliVZERO {
  
public:
  AliVZEROv0();
  AliVZEROv0(const char *name, const char *title);
  virtual       ~AliVZEROv0() {}
  virtual void   AddHit(Int_t  , Int_t *, Float_t *); 
  virtual void   AddDigits(Int_t* , Int_t* );
  virtual void   CreateGeometry();
  virtual void   BuildGeometry();
  virtual void   CreateMaterials();
  virtual void   DrawModule();
  virtual void   Init();
  virtual void   MakeBranch(Option_t *option);
  virtual Int_t  IsVersion() const {return 0;}
  virtual void   StepManager();
 
public:
   Int_t         fIdSens1;      // Sensitive volume  in VZERO
   Int_t         digits[1];  
   Int_t         tracks[5];   
   
  ClassDef(AliVZEROv0,1)  //Class for VZERO version 0
};

#endif


