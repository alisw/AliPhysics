#ifndef ALIGENBEAMGAS_H
#define ALIGENBEAMGAS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//
// Generator to simulate beam gas interactions.
// At present single interactions are read from an external file. 
// Author: andreas.morsch@cern.ch

#include "AliGenExtFile.h"

class AliGenBeamGas : public AliGenExtFile
{
 public:
    AliGenBeamGas();
    virtual ~AliGenBeamGas();
    //
    virtual void SetNumberOfInteractions(Int_t n) 
	{fInteractions = n;}
    // Initialise 
    virtual void Init();
    // generate event
    virtual void Generate();
 protected:
    Int_t fInteractions;    // Number of interactions
 private:
    AliGenBeamGas(const AliGenBeamGas &beamgas);
    AliGenBeamGas & operator=(const AliGenBeamGas &beamgas);
    
    ClassDef(AliGenBeamGas,1) //Generator for beam gas interactions
	
};
#endif






