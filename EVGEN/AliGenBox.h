#ifndef ALIGENBOX_H
#define ALIGENBOX_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliGenerator.h"
class AliGenBox : public AliGenerator
{
 public:

  AliGenBox();
  AliGenBox(Int_t npart);
  virtual ~AliGenBox() {}
  virtual void Generate();
  virtual void Init();
  virtual void SetPart(Int_t part) {fIpart=part;}
protected:

  Int_t fIpart; // Particle type

  ClassDef(AliGenBox,1) // Square box random generator
};

#endif
