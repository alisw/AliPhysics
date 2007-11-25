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
  UInt_t  GetTestFilterMask(){return fTestFilterMask;}	    

  // Setters
  virtual void SetNumberOfAOD(Int_t i=1) {fNaod = i;}    
  virtual void SetTestFilterMask(UInt_t i){fTestFilterMask = i;}

 protected:
  Int_t   fNaod;           // number of aods
  UInt_t  fTestFilterMask; // Filter Mask for jets, not tested if 0
  
  ClassDef(AliJetAODReaderHeader,2);
};
 
#endif
