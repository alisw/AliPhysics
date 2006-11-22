#ifndef ALIJETMCREADER_H
#define ALIJETMCREADER_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
// Jet MC Reader 
// MC reader for jet analysis
// Author: Mercedes Lopez Noriega (mercedes.lopez.noriega@cern.ch)

#include "AliJetESDReader.h"

class AliJetMCReader : public AliJetESDReader
{
 public: 
    AliJetMCReader();
    virtual ~AliJetMCReader();
    Bool_t FillMomentumArray(Int_t event);
    
 protected:
    Float_t fPdgC;   // Pdg code
    ClassDef(AliJetMCReader,1)
};
 
#endif
