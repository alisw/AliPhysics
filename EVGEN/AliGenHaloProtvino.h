#ifndef ALIGENHALOPROTVINO_H
#define ALIGENHALOPROTVINO_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id $ */


#include "AliGenerator.h"
#include <TString.h>

// Read background particles from a boundary source
// Very specialized generator to simulate background from beam halo.
// Author: andreas.morsch@cern.ch

class AliGenHaloProtvino : public AliGenerator
{
public:
    AliGenHaloProtvino();
    AliGenHaloProtvino(Int_t npart);
    AliGenHaloProtvino(const AliGenHaloProtvino &HaloProtvino);
    virtual ~AliGenHaloProtvino();
    virtual void Init();
    virtual void SetFileName(TString filename) {fFileName=TString(filename);}
    virtual void Generate();
    virtual Float_t GassPressureWeight(Float_t zPrimary);
    AliGenHaloProtvino & operator=(const AliGenHaloProtvino & rhs);

protected:
  FILE *fp;                             // ! Pointer to file
  TString  fFileName;                   //   Choose the file
  ClassDef(AliGenHaloProtvino,1)        // LHC background boundary source (Protvino Group results)
};
#endif






