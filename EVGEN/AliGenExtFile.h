#ifndef ALIGENEXTFILE_H
#define ALIGENEXTFILE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


// Event generator that can read events from a files.
// The reading is performed by a realisation of AliGenReader specific to the file format.
// Author: andreas.morsch@cern.ch

#include "AliGenMC.h"
class AliGenReader;


class TTree;

class AliGenExtFile : public AliGenMC
{
 public:
    AliGenExtFile();
    AliGenExtFile(Int_t npart);
     virtual ~AliGenExtFile();
    // Initialise 
    virtual void Init();
    // generate event
    virtual void Generate();
    void SetReader(AliGenReader* reader) {fReader = reader;}
    void SetStartEvent(Int_t startEvent) {fStartEvent = startEvent;}
    AliGenReader* Reader() const {return fReader;}
 protected:
    void CdEventFile();
    const Text_t     *fFileName;      //! File to read from
    AliGenReader     *fReader;        //! Reader to read the file
    Int_t  fStartEvent; //! Start event number

 private:
    AliGenExtFile(const AliGenExtFile &ext);
    AliGenExtFile & operator=(const AliGenExtFile & rhs);
    
  ClassDef(AliGenExtFile,2) //Generate particles from external file
};
#endif






