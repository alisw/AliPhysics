#ifndef ALIJETESDMCREADER_H
#define ALIJETESDMCREADER_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
// Jet ESD Reader 
// ESD reader for jet analysis (it reads the esd and the MC trees)
// Author: Mercedes Lopez Noriega (mercedes.lopez.noriega@cern.ch)

#include "AliJetESDReader.h"
class AliJetESDReaderHeader;
class AliHeader;



class AliJetESDmcReader : public AliJetESDReader
{
 public: 
  AliJetESDmcReader();
  virtual ~AliJetESDmcReader();
  // Setters
  Bool_t FillMomentumArray(Int_t event); 
  void   OpenInputFiles();
   
 protected:
  TChain*    fChainMC;   // chain for mc information
  AliHeader* fAliHeader; //! Event header
  ClassDef(AliJetESDmcReader,1)
};
 
#endif
