#ifndef ALIGENREADERECALHIJING_H
#define ALIGENREADERECALHIJING_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliGenReader.h"


class AliGenReaderEcalHijing : public AliGenReader
{
 public:
    AliGenReaderEcalHijing();
    
    AliGenReaderEcalHijing(const AliGenReaderEcalHijing &reader):AliGenReader(reader)
	{reader.Copy(*this);}
    virtual ~AliGenReaderEcalHijing(){;}
    // Initialise 
    virtual void Init();
    // Read
    virtual Int_t NextEvent();
    virtual TParticle*  NextParticle();
    AliGenReaderEcalHijing & operator=(const AliGenReaderEcalHijing & rhs);
 private:
    void Copy(AliGenReaderEcalHijing&) const;
 protected:
    Int_t             fNcurrent;      // points to the next entry
    Int_t             fNparticle;     // 
    
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
    ClassDef(AliGenReaderEcalHijing,1) // Read particles from cwn-ntuple
};
#endif






