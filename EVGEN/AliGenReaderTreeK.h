#ifndef ALIGENREADERTreeK_H
#define ALIGENREADERTREEK_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliGenReader.h"
class TFile;
class AliStack;
class AliHeader;


class AliGenReaderTreeK : public AliGenReader
{
 public:
    AliGenReaderTreeK();
    
    AliGenReaderTreeK(const AliGenReaderTreeK &reader){;}
    virtual ~AliGenReaderTreeK(){;}
    // Initialise 
    virtual void Init();
    // Read
    virtual Int_t NextEvent();
    virtual TParticle*  NextParticle();
    AliGenReaderTreeK & operator=(const AliGenReader & rhs){return *this;}
 protected:
    Int_t             fNcurrent;          // points to the next entry
    Int_t             fNparticle;         // Next particle in list
    Int_t             fNp;                // number of particles
    TFile            *fFile;              //! pointer to file
    TFile            *fBaseFile;          //! pointer to base file
    AliStack         *fStack;             //! Particle stack
    AliHeader        *fHeader;            //! Pointer to event header
    ClassDef(AliGenReaderTreeK,1) // Read particles from TreeK
};
#endif






