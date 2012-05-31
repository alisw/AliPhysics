#ifndef ALIGENREADER_H
#define ALIGENREADER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

// Interface for reading events from files.
// Realisations of this interface have to be used with AliGenExFile.
// Author: andreas.morsch@cern.ch

#include "TObject.h"

class TParticle;
class AliRunLoader;

class AliGenReader : public TObject
{
 public:
    AliGenReader():fFileName(NULL),fCode(kPDG){;}
    AliGenReader(const AliGenReader &reader)
	:TObject(reader), fFileName(NULL), fCode(kPDG){reader.Copy(*this);}
    virtual ~AliGenReader(){;}
    virtual void SetFileName(const Text_t *filname) {fFileName=filname;}
    virtual AliRunLoader * GetRunLoader() const {return 0x0;}
    virtual void Init()                                                    = 0;
    virtual Int_t NextEvent()                                              = 0;
    virtual TParticle* NextParticle()                                      = 0;
    virtual void RewindEvent()                                             = 0;
    typedef enum {kPDG, kGEANT3} Code_t;
    void SetParticleCode(Code_t code) {fCode = code;}
    AliGenReader & operator=(const AliGenReader & rhs);

 protected:
    const Text_t *fFileName;      // Name of file to read from
    Code_t        fCode;          // Particle code type
 private:
    void Copy(TObject&) const;
    
    ClassDef(AliGenReader,1) //Generate particles from external file
};
#endif






