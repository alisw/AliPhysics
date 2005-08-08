#ifndef AliRICHv1_h
#define AliRICHv1_h

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliRICH.h"             //base class 
#include "AliRICHDigitizer.h"  //CreateDigitizer()

class AliRICHv1 : public AliRICH 
{
public:
                 AliRICHv1():AliRICH()                                               {;}          //default ctor
                 AliRICHv1(const char *name, const char *title):AliRICH(name,title)  {;}          //named ctor
  virtual       ~AliRICHv1()                                                         {;}          //dtor
//framework part  
  void   Init()                                                              {;}
  Int_t          IsVersion        (                          )const{return 1;}  
  void           StepManager      (                          );                                           //called from AliSimulation or AliRun when transport particles
  void           Hits2SDigits     (                          );                                           //called from AliSimulation for Hits->SDigits
  AliDigitizer*  CreateDigitizer  (AliRunDigitizer *pMan     )const{return new AliRICHDigitizer(pMan);}   //called from AliSimulation for SDigits->Digits
  void           Digits2Raw       (                          );                                           //called from AliSimulation for Digits->Raw
//private part  
          Bool_t IsLostByFresnel  ();                                                                     //checks if the photon lost on Fresnel reflection  
          void   GenerateFeedbacks(Int_t iChamber,Float_t eloss=0);                                       //generates feedback photons; eloss=0 for photon
protected:
  ClassDef(AliRICHv1,1)                                                 //RICH full version for simulation
};
		
#endif
