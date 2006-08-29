#ifndef ALIPXCONEJETFINDER_H
#define ALIPXCONEJETFINDER_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
 
//---------------------------------------------------------------------
// Pxcone Jet finder 
// manages the search for jets 
// Author: jgcn@mda.cinvestav.mx
//---------------------------------------------------------------------

#include "AliJetFinder.h"

class AliPxconeJetHeader;

class AliPxconeJetFinder : public AliJetFinder 
{
 public:

  AliPxconeJetFinder();
  ~AliPxconeJetFinder();

  // getters

  // setters
  void SetJetHeader(AliPxconeJetHeader* h) {fHeader= h;}
  // others
  void Reset();
  void FindJets();
  void WriteJHeaderToFile();

 protected:
  AliPxconeJetFinder(const AliPxconeJetFinder& rPxJet);
  AliPxconeJetFinder& operator = (const AliPxconeJetFinder& rhsh);

  AliPxconeJetHeader* fHeader;         // pointer to jet header

  ClassDef(AliPxconeJetFinder,1)
};

#endif
