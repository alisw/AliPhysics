#ifndef ALIGENHIJINGPARA_H
#define ALIGENHIJINGPARA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliGenerator.h"
#include "TF1.h"

class AliGenHIJINGpara : public AliGenerator
{
 public:

  AliGenHIJINGpara();
  AliGenHIJINGpara(Int_t npart);
  AliGenHIJINGpara(const AliGenHIJINGpara &HIJINGpara);
     
  virtual ~AliGenHIJINGpara();
  virtual void Generate();
  virtual void Init();
  AliGenHIJINGpara & operator=(const AliGenHIJINGpara & rhs);
 protected:

  TF1* fPtpi; // Parametrised pt distribution for pi
  TF1* fPtka; // Parametrised pt distribution for ka
  TF1* fETApic; // Parametrised eta distribution for pi
  TF1* fETAkac; // Parametrised eta distribution fro ka

  ClassDef(AliGenHIJINGpara,1) // Hijing parametrisation generator
};
#endif










