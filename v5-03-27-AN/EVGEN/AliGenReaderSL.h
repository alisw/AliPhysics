#ifndef ALIGENREADERSL_H
#define ALIGENREADERSL_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

// Realisations of the AliGenReader interface to be used with AliGenExFile.
// NextEvent() loops over events 
// and NextParticle() loops over particles. 
// This implementation reads various StarLight output formats 
// Author: andreas.morsch@cern.ch

#include "AliGenReader.h"

class TParticle;

class AliGenReaderSL : public AliGenReader
{
 public:
    AliGenReaderSL():fFile(0), fNParticles(0), fFormat(0) {;}
    AliGenReaderSL(const AliGenReaderSL &reader)
	:AliGenReader(reader), fFile(0), fNParticles(0), fFormat(0)  {reader.Copy(*this);}
    virtual ~AliGenReaderSL(){;}
    virtual void Init();
    virtual Int_t NextEvent();
    virtual TParticle* NextParticle();
    virtual void RewindEvent();
    virtual void SetFormat(Int_t format) {fFormat = format;}
    AliGenReaderSL & operator=(const AliGenReaderSL & rhs);

 protected:
    FILE *fFile;          // pointer to the file
    Int_t fNParticles;    // Number of particles
    Int_t fFormat;        // File format
 private:
    void Copy(TObject&) const;
    
    ClassDef(AliGenReaderSL, 1) //Generate particles from external file
};
#endif






