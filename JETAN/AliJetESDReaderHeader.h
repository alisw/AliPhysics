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
  Float_t GetDCA()        const {return fDCA;}       
  Float_t GetTLength()    const {return fTLength;}    
  Bool_t  ReadSignalOnly() const  {return fReadSignalOnly;}
  Bool_t  ReadBkgdOnly() const  {return fReadBkgdOnly;}
  Int_t   GetNesd()       const {return fNesd;}
	    
  // Setters
  virtual void SetDCA(Float_t dca = 0.0) {fDCA = dca;}
  virtual void SetTLength(Float_t length = 0.0) {fTLength = length;}
  virtual void SetReadSignalOnly(Bool_t flag = kTRUE) {fReadSignalOnly = flag;}
  virtual void SetReadBkgdOnly(Bool_t flag = kTRUE) {fReadBkgdOnly = flag;}
  virtual void SetNumberOfESD(Int_t i=1) {fNesd = i;}
    
 protected:
  //parameters set by user
  Float_t fDCA;            // dca cut
  Float_t fTLength;        // track length cut
  Bool_t  fReadSignalOnly; // read particles from signal event only
  Bool_t  fReadBkgdOnly;   // read particles from bkgd event only
  Int_t   fNesd;           // number of esds
  
  ClassDef(AliJetESDReaderHeader,2);
};
 
#endif
