#ifndef ALIGENKINE_H
#define ALIGENKINE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */



#include "AliGenerator.h"
#include "AliHeader.h"
#include "TTree.h"
#include "TFile.h"
#include "TArrayI.h"
#include "TParticle.h"

class TParticle;
class AliStack;
class TClonesArray;

// Read background particles from a FLUKA boundary source file

class AliGenKine : public AliGenerator
{
 public:
    AliGenKine();
    AliGenKine(Int_t npart);
    virtual ~AliGenKine();
    // Initialise 
    virtual void Init() {}
    // set file name for data file
    virtual void SetFileName (const Text_t *filname)  {fFileName  = filname;}
    // generate event
    virtual void Generate();
protected:
    const Text_t     *fFileName;          //! Choose the file
    Int_t             fNcurrent;          // points to the next entry
    Int_t             fNp;                // number of particles
    TFile            *fFile;              // ! pointer to file
    TFile            *fBaseFile;          // ! pointer to base file
    AliStack         *fStack;             // ! Particle stack
    
//
    ClassDef(AliGenKine,1)                // Generate particles from external file
};
#endif














