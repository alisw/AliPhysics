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
  Float_t GetPtCut()       const  {return fPtCut;}
  Float_t GetDCA()         const  {return fDCA;} // not working so far..(always 0)
  Float_t GetTLength()     const  {return fTLength;} // not working so far.. (always 0)
  Bool_t  ReadSignalOnly() const  {return fReadSignalOnly;}
  
	  
  // Setters
  virtual void SetPtCut(const Float_t par = 2.0) {fPtCut = par;}
  virtual void SetDCA(const Float_t dca = 3.0) {fDCA = dca;}
  virtual void SetTLength(const Float_t length = 0.0) {fTLength = length;}
  virtual void SetReadSignalOnly(Bool_t flag = kTRUE) {fReadSignalOnly = flag;}
  
 protected:
  //parameters set by user
  Float_t fPtCut;          // pt cut 
  Float_t fDCA;            // dca cut
  Float_t fTLength;        // track length cut
  Bool_t  fReadSignalOnly; // read particles from signal event only
  ClassDef(AliJetESDReaderHeader,1);
};
 
#endif
