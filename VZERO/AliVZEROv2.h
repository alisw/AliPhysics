#ifndef ALIVZEROV2_H
#define ALIVZEROV2_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


///////////////////////////////////////////////////
//                                               //
//  Manager and hits classes for set : VZERO     //
//                                               //
///////////////////////////////////////////////////
 
#include "AliVZERO.h"

class AliVZEROv2 : public AliVZERO {
  
public:
  AliVZEROv2();
  AliVZEROv2(const char *name, const char *title);
  virtual       ~AliVZEROv2() {}
  virtual void   AddHit(Int_t track, Int_t *vol, Float_t *hits); 
  virtual void   AddDigits(Int_t *tracks, Int_t *digits);
  virtual void   CreateGeometry();
  virtual void   CreateMaterials();
  virtual void   DrawModule() const;
  virtual void   Init();
  virtual void   MakeBranch(Option_t *option);
  virtual Int_t  IsVersion() const {return 2;}
  virtual void   StepManager();
   
  ClassDef(AliVZEROv2,1)  //Class for VZERO version 2
};

#endif


