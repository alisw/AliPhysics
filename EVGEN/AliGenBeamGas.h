#ifndef ALIGENBEAMGAS_H
#define ALIGENBEAMGAS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


#include "AliGenExtFile.h"
#include "AliGenReader.h"

class TTree;

class AliGenBeamGas : public AliGenExtFile
{
 public:
    AliGenBeamGas();
    AliGenBeamGas(const AliGenBeamGas &beamgas);
    virtual ~AliGenBeamGas();
    //
    virtual void SetNumberOfInteractions(Int_t n) 
	{fInteractions = n;}
    // Initialise 
    virtual void Init();
    // generate event
    virtual void Generate();
 private:
    void Copy(AliGenBeamGas&) const;
 protected:
    Int_t fInteractions;
    
    ClassDef(AliGenBeamGas,1) //Generate for beam gas interactions
	
};
#endif






