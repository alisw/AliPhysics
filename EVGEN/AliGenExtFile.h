#ifndef ALIGENEXTFILE_H
#define ALIGENEXTFILE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


// Event generator that can read events from a files.
// The reading is performed by a realisation of AliGenReader specific to the file format.
// Author: andreas.morsch@cern.ch

#include "AliGenMC.h"
#include "AliGenReader.h"

class TTree;

class AliGenExtFile : public AliGenMC
{
 public:
    AliGenExtFile();
    AliGenExtFile(Int_t npart);
    AliGenExtFile(const AliGenExtFile &ext);
    virtual ~AliGenExtFile();
    // Initialise 
    virtual void Init();
    // generate event
    virtual void Generate();
    AliGenExtFile & operator=(const AliGenExtFile & rhs);
    void SetReader(AliGenReader* reader) {fReader = reader;}
 protected:
    void CdEventFile();
    void Copy(TObject&) const;
    const Text_t     *fFileName;      //! File to read from
    AliGenReader     *fReader;        //! Reader to read the file
    
  ClassDef(AliGenExtFile,1) //Generate particles from external file
};
#endif






