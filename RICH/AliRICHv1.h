#ifndef AliRICHv1_h
#define AliRICHv1_h

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliRICH.h"

class AliRICHv1 : public AliRICH 
{
public:
  inline                AliRICHv1():AliRICH()                           {;}
                        AliRICHv1(const char *name, const char *title);
  inline virtual       ~AliRICHv1()                                     {;}
         virtual void   Init();
  inline virtual Int_t  IsVersion()                                const{return 1;}
         virtual void StepManager();
private:
  ClassDef(AliRICHv1,1)//RICH full version for simulation
};
		
#endif
