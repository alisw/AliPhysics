#ifndef STARTHIT_H
#define STARTHIT_H
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
  Int_t      fVolume;
  Int_t      fPmt;
  Int_t      fParticle;     //Particle identificator
  Float_t    fEdep;    //Energy deposition
  Float_t    fEtot;    //Energy of particle 
  Float_t    fTime;    //Time of flight 
 
public:
  AliSTARThit() {}
  AliSTARThit(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits);
  virtual ~AliSTARThit() {}
  
  ClassDef(AliSTARThit,1)  //Hits for detector START
};
#endif
