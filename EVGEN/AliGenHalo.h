#ifndef AliGenHalo_H
#define AliGenHalo_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/////////////////////////////////////////////////////////
//  Manager and hits classes for set:MUON version 0    //
/////////////////////////////////////////////////////////
#include "AliGenerator.h"
#include "TNamed.h"
#include "TF1.h"
#include "TArrayF.h"
#include "TTree.h"

// Read background particles from a FLUKA boundary source file

class AliGenHalo : public AliGenerator
{
 
protected:
  FILE *fp;                             //! Pointer to file
  const Text_t     *fFileName;          //!Choose the file
  
public:
    AliGenHalo();
    AliGenHalo(Int_t npart);
    virtual ~AliGenHalo();
    virtual void Init();
    virtual void SetFileName(const Text_t *filname) {fFileName=filname;}
    virtual void Generate();
    ClassDef(AliGenHalo,1)
};
#endif






