#ifndef AliRICHv0_h
#define AliRICHv0_h

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliRICH.h"

class AliRICHv0 : public AliRICH 
{    
public:
  AliRICHv0():AliRICH()                                            {;}
  AliRICHv0(const char *name, const char *title);
  virtual       ~AliRICHv0()                                       {;}
  virtual void   Init()                                            {;}
  virtual Int_t  IsVersion()                                  const{return 0;}
  virtual void   StepManager();
protected:
  ClassDef(AliRICHv0,1)  //RICH coarse version for material budget study and debuging
};
	
#endif
