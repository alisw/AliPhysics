#ifndef ALIGENHALO_H
#define ALIGENHALO_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


#include "AliGenerator.h"

// Read background particles from a FLUKA boundary source file

class AliGenHalo : public AliGenerator
{
public:
    AliGenHalo();
    AliGenHalo(Int_t npart);
    AliGenHalo(const AliGenHalo &Halo);
    virtual ~AliGenHalo();
    virtual void Init();
    virtual void SetFileName(const Text_t *filname) {fFileName=filname;}
    virtual void Generate();
    AliGenHalo & operator=(const AliGenHalo & rhs);
protected:
  FILE *fp;                             //   Pointer to file
  const Text_t     *fFileName;          // ! Choose the file
  ClassDef(AliGenHalo,1) // LHC background boundary source (MARS input)
};
#endif






