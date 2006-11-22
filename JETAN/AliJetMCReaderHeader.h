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
 protected:
  ClassDef(AliJetMCReaderHeader,1);
};
 
#endif
