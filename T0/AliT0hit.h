#ifndef ALIT0hit_H
#define ALIT0hit_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
////////////////////////////////////////////////
//  Manager and hits classes for set:T0     //
////////////////////////////////////////////////
 
#include "AliHit.h"
 
class AliT0hit : public AliHit {
public:
  AliT0hit(){}//Empty ctor
  AliT0hit(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits);
  virtual ~AliT0hit(){}//Empty virtual dtor
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
 
  
   ClassDef(AliT0hit,2)  //Hits for detector T0
};

typedef AliT0hit AliSTARThit; // for backward compatibility

#endif//ALIT0hit_H
