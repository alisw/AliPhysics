#ifndef AliRICHv0_h
#define AliRICHv0_h

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


#include "AliRICH.h"

class AliRICHv0 : public AliRICH 
{    
public:
  inline                AliRICHv0():AliRICH()                             {;}
  inline                AliRICHv0(const char *name, const char *title);
  inline virtual       ~AliRICHv0()                                       {;}
  inline virtual void   Init()                                            {;}
  inline virtual Int_t  IsVersion()                                  const{return 0;}
         virtual void   StepManager();
protected:
  ClassDef(AliRICHv0,1)  //RICH coarse version for naterial budget study and debug
};
	
#endif
