#ifndef ALISTARThit_H
#define ALISTARThit_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
////////////////////////////////////////////////
//  Manager and hits classes for set:START     //
////////////////////////////////////////////////
 
#include "AliHit.h"
 
class AliSTARThit : public AliHit {
public:
  AliSTARThit(){}//Empty ctor
  AliSTARThit(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits);
  virtual ~AliSTARThit(){}//Empty virtual dtor
  Int_t Volume() const {return fVolume;}
  Int_t Pmt() const {return fPmt;}
  Float_t Particle() const {return fParticle;} 
  Double_t Etot() const {return fEtot;}
  Float_t Time() const {return fTime;}

private:
  Int_t      fVolume;   //T0 arm mark
  Int_t      fPmt;      //PMT number in the arm  
  Int_t      fParticle; //Primary particle ID
  Double_t    fEtot;     //Energy of primary particle at the entrance to radiator 
  Float_t    fTime;     //Primary particle TOF 
 
  
   ClassDef(AliSTARThit,2)  //Hits for detector START
};


#endif//ALISTARThit_H
