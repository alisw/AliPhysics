#ifndef ALIGENBOX_H
#define ALIGENBOX_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


// Generator for particles in a preset
// kinematic range (flat distribution)
// Comments and suggestions: andreas.morsch@cern.ch


#include "AliGenerator.h"
class AliGenBox : public AliGenerator
{
 public:

  AliGenBox();
  AliGenBox(Int_t npart);
  virtual ~AliGenBox() {}
  virtual void Generate();
  virtual void Init();
  virtual void SetEtaRange(Float_t etamin, Float_t etamax)
      {SetBit(kEtaRange);fEtaMin = etamin; fEtaMax = etamax;}
  virtual void SetPart(Int_t part) {fIpart=part;}
  virtual void SetParticleType(Int_t part) {SetPart(part);}
protected:

  Int_t fIpart; // Particle type
  Float_t fEtaMin;  // Minimum eta 
  Float_t fEtaMax;  // Maximum eta
  ClassDef(AliGenBox,2) // Square box random generator
};

#endif
