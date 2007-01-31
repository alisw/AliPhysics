#ifndef ALIGENHALO_H
#define ALIGENHALO_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


#include "AliGenerator.h"
#include <TString.h>

// Read background particles from a boundary source
// Very specialized generator to simulate background from beam halo.
// Author: andreas.morsch@cern.ch

class AliGenHalo : public AliGenerator
{
 public:
    AliGenHalo();
    AliGenHalo(Int_t npart);
    virtual ~AliGenHalo();
    virtual void Init();
    virtual void SetFileName(TString filename) {fFileName=TString(filename);}
    virtual void Generate();
 protected:
    FILE *fp;                             // ! Pointer to file
    TString  fFileName;                   //   Choose the file
 private:
    AliGenHalo(const AliGenHalo &Halo);
    AliGenHalo & operator=(const AliGenHalo & rhs);

    ClassDef(AliGenHalo,1) // LHC background boundary source (MARS input)
};
#endif






