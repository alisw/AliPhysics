#ifndef ALIGENTHETASLICE_H
#define ALIGENTHETASLICE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


// Generates n particles with in the same phi angle, varies theta
// In equidistant intervals
// Comments and suggestions: Jiri.Chudoba@cern.ch


#include "AliGenerator.h"
class AliGenThetaSlice : public AliGenerator
{
 public:

  AliGenThetaSlice();
  AliGenThetaSlice(Int_t npart);
  virtual ~AliGenThetaSlice() {}
  virtual void Generate();
  virtual void Init();
  virtual void SetPart(Int_t part) {fIpart=part;}
protected:

  Int_t fIpart; // Particle type

  ClassDef(AliGenThetaSlice,1) // theta slices phi constant random generator
};

#endif
