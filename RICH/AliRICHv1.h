#ifndef AliRICHv1_h
#define AliRICHv1_h

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliRICH.h"

class AliRICHv1 : public AliRICH 
{
public:
                 AliRICHv1():AliRICH()                                               {;}
                 AliRICHv1(const char *name, const char *title):AliRICH(name,title)  {;}
  virtual       ~AliRICHv1()                                                         {;}
  virtual void   Init()                                                              {;}
  virtual Int_t  IsVersion()                                                    const{return 1;}
  
  virtual void   StepManager();                                         //full slow step manager
          Bool_t IsLostByFresnel();                                     //checks if the photon lost on Fresnel reflection  
          void   GenerateFeedbacks(Int_t iChamber,Float_t eloss=0);     //generates feedback photons; eloss=0 for photon
private:
  ClassDef(AliRICHv1,1)                                                 //RICH full version for simulation
};
		
#endif
