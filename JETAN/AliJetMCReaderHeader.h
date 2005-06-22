#ifndef ALIJETMCREADERHEADER_H
#define ALIJETMCREADERHEADER_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
// Jet MC Reader Header
// Header for the MC reader in the jet analysis
// Author: Mercedes Lopez Noriega (mercedes.lopez.noriega@cern.ch)

#include "AliJetReaderHeader.h"
 
class AliJetMCReaderHeader : public AliJetReaderHeader
{

 public:
  AliJetMCReaderHeader();
  virtual ~AliJetMCReaderHeader();
  
  // Getters
  Double_t GetPtCut() const {return fPtCut;}
  Float_t GetDCA() const {return fDCA;} // not working so far..(always 0)
  
  // Setters
  virtual void SetPtCut(const Double_t par=2.0) {fPtCut = par;}
  virtual void SetDCA(const Float_t dca = 3.0) {fDCA = dca;} // this is for the track!!! (recommended by P.H.)
  
 protected:
  //parameters set by user
  Double_t fPtCut;
  Float_t fDCA;


  ClassDef(AliJetMCReaderHeader,1);
};
 
#endif
