#ifndef ALIGENREADER_H
#define ALIGENREADER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "TObject.h"

class TParticle;

class AliGenReader : public TObject
{
 public:
    AliGenReader():fFileName(NULL),fCode(kPDG){;}
    AliGenReader(const AliGenReader &reader):fFileName(NULL),fCode(kPDG){;}
    
    virtual ~AliGenReader(){;}
    // Initialise 
    virtual void Init() {}
    // set file name of data file
    virtual void SetFileName(const Text_t *filname) {fFileName=filname;}
    // Read
    virtual Int_t NextEvent(){return 0;}
    enum Code_t {kPDG, kGEANT3};
    void SetParticleCode(Code_t code) {fCode = code;}
    virtual TParticle* NextParticle(){return NULL;}
    virtual void RewindEvent();
        
    AliGenReader & operator=(const AliGenReader & rhs);
 protected:
    const Text_t *fFileName;      // Name of file to read from
    Code_t        fCode;          // Particle code type
    
    ClassDef(AliGenReader,1) //Generate particles from external file
};
#endif






