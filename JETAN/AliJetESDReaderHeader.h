#ifndef ALIJETESDREADERHEADER_H
#define ALIJETESDREADERHEADER_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
// Jet ESD Reader Header
// Header for the ESD reader in the jet analysis
// Author: Mercedes Lopez Noriega (mercedes.lopez.noriega@cern.ch)

#include "AliJetReaderHeader.h"
 
class AliJetESDReaderHeader : public AliJetReaderHeader
{

 public:
  AliJetESDReaderHeader();
  virtual ~AliJetESDReaderHeader();
  
  // Getters
  Bool_t  ReadSignalOnly() const  {return fReadSignalOnly;}
  Bool_t  ReadBkgdOnly() const  {return fReadBkgdOnly;}
	  
  // Setters
  virtual void SetReadSignalOnly(Bool_t flag = kTRUE) {fReadSignalOnly = flag;}
  virtual void SetReadBkgdOnly(Bool_t flag = kTRUE) {fReadBkgdOnly = flag;}
  
 protected:
  //parameters set by user
  Bool_t  fReadSignalOnly; // read particles from signal event only
  Bool_t  fReadBkgdOnly;   // read particles from bkgd event only
  ClassDef(AliJetESDReaderHeader,2);
};
 
#endif
