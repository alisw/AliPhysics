#ifndef ALIGENREADERCWN_H
#define ALIGENREADERCWN_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//
// Realisation of AliGenReader to be used with AliGenExtFile
// It reads events from a ntuple like event structure.
// Author: andreas.morsch@cern.ch
//
#include "AliGenReader.h"
#include <Rtypes.h>


class AliGenReaderCwn : public AliGenReader
{
 public:
    AliGenReaderCwn();
    AliGenReaderCwn(const AliGenReaderCwn &reader);
    virtual ~AliGenReaderCwn();
        // Initialise 
    virtual void Init();
    // Read
    virtual Int_t NextEvent();
    virtual TParticle*  NextParticle();
    virtual void RewindEvent(){;}
    AliGenReaderCwn & operator=(const AliGenReaderCwn & rhs);
    
 protected:
    Int_t             fNcurrent;      // points to the next entry
    Int_t             fNparticle;     // particle number in event
    Int_t             fNparticleMax;  // number of particles in event    
    TTree            *fTreeNtuple;    // pointer to the TTree
    //Declaration of leaves types
    Int_t           fNihead;          // Number of entries in integer header  
    Int_t           fIhead[12];       // Integer header
    Int_t           fNrhead;          // Number of entries in float header
    Float_t         fRhead[6];        // Float header
    UInt_t          fIdpart;          // Particle type
    Float_t         fTheta;           // Theta 
    Float_t         fPhi;             // Phi
    Float_t         fP;               // Total momentum
    Float_t         fE;               // Total energy
 private:
    void Copy(TObject&) const;
    ClassDef(AliGenReaderCwn,1) // Read particles from cwn-ntuple
};
#endif






