#ifndef ALIGENREADERECALHIJING_H
#define ALIGENREADERECALHIJING_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//
// Realisation of AliGenReader to be used with AliGenExtFile
// It reads Hijing events from a ntuple like event structure.
// Author: andreas.morsch@cern.ch
//
#include "AliGenReader.h"


class AliGenReaderEcalHijing : public AliGenReader
{
 public:
    AliGenReaderEcalHijing();
    
    AliGenReaderEcalHijing(const AliGenReaderEcalHijing &reader);
    virtual ~AliGenReaderEcalHijing(){;}
    // Initialise 
    virtual void Init();
    // Read
    virtual Int_t NextEvent();
    virtual TParticle*  NextParticle();
    virtual void RewindEvent(){;}
    AliGenReaderEcalHijing & operator=(const AliGenReaderEcalHijing & rhs);

 protected:
    Int_t             fNcurrent;      // points to the next entry
    Int_t             fNparticle;     // number of particles
    
    TTree            *fTreeNtuple;    // pointer to the TTree
    //Declaration of leaves types
    Int_t           fNjatt;           // Number of particles
    Int_t           fNahij;           // Number of particles in alice accept. 
    Int_t           fNphij;           // ?
    Int_t           fKhij[10000];     // particle code
    Float_t         fPxhij[10000];    // px
    Float_t         fPyhij[10000];    // py
    Float_t         fPzhij[10000];    // pz
    Float_t         fEhij[10000];     // energy
 private:
    void Copy(TObject&) const;
    
    ClassDef(AliGenReaderEcalHijing,1) // Read particles from cwn-ntuple
};
#endif






