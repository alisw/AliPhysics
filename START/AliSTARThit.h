#ifndef AliSTARThit_H
#define AliSTARThit_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
////////////////////////////////////////////////
//  Manager and hits classes for set:START     //
////////////////////////////////////////////////
 
#include "AliDetector.h"
#include "AliHit.h"
 
class AliSTARThit : public AliHit {
public:
  Int_t      fVolume;   //T0 arm mark
  Int_t      fPmt;      //PMT number in the arm  
  Int_t      fParticle; //Primary particle ID
  Float_t    fEdep;     //Energy deposition
  Float_t    fEtot;     //Energy of primary particle at the entrance to radiator 
  Float_t    fTime;     //Primary particle TOF 
 
public:
   AliSTARThit(){}//Empty ctor
   AliSTARThit(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits);
   virtual ~AliSTARThit(){}//Empty virtual dtor
  
   ClassDef(AliSTARThit,1)  //Hits for detector START
};


#endif//AliSTARThit_H
