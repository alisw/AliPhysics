#ifndef ALIGENEXTFILE_H
#define ALIGENEXTFILE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


// Event generator that can read the old ALICE event format based on CW-ntuples
// http://consult.cern.ch/alice/Internal_Notes/1995/32/abstract
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
    void Copy(AliGenExtFile&) const;
    const Text_t     *fFileName;      //! File to read from
    AliGenReader     *fReader;        //! Reader to read the file
    
  ClassDef(AliGenExtFile,1) //Generate particles from external file
};
#endif






