#ifndef ALIJETAODREADERHEADER_H
#define ALIJETAODREADERHEADER_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
// Jet AOD Reader Header
// Header for the AOD reader in the jet analysis
// Author: Davide Perrino (davide.perrino@cern.ch)

#include "AliJetReaderHeader.h"
 
class AliJetAODReaderHeader : public AliJetReaderHeader
{

 public:
  AliJetAODReaderHeader();
  virtual ~AliJetAODReaderHeader();
  
  // Getters
  Int_t   GetNaod()       const {return fNaod;}
	    
  // Setters
  virtual void SetNumberOfAOD(Int_t i=1) {fNaod = i;}
    
 protected:
  Int_t   fNaod;           // number of aods
  
  ClassDef(AliJetAODReaderHeader,1);
};
 
#endif
