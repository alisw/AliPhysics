#ifndef AliRICHv0_h
#define AliRICHv0_h

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliRICH.h"

class AliRICHv0 : public AliRICH 
{    
public:
                 AliRICHv0():AliRICH()                                              {;}       //default ctor
                 AliRICHv0(const char *name, const char *title):AliRICH(name,title) {;}       //named ctor
  virtual       ~AliRICHv0()                                         {;}                      //dtor
  virtual void   Init()                                              {;}                      //interface from AliRICH
  virtual Int_t  IsVersion()                                    const{return 0;}              //interface from AliRICH
  virtual void   StepManager();                                                               //interface from AliRICH
protected:
  ClassDef(AliRICHv0,1)  //RICH coarse version for material budget study and debuging
};
	
#endif
