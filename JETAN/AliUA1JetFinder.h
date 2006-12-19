#ifndef ALIUA1JETFINDER_H
#define ALIUA1JETFINDER_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
 
//---------------------------------------------------------------------
// UA1 Jet finder 
// manages the search for jets 
// Author: jgcn@mda.cinvestav.mx
// (code adapted from EMCAL directory)
//---------------------------------------------------------------------

#include "AliJetFinder.h"
class AliUA1JetHeader;
class TH2F;

class AliUA1JetFinder : public AliJetFinder 
{
 public:

  AliUA1JetFinder();
  ~AliUA1JetFinder();

  // getters

  // setters
  void SetJetHeader(AliUA1JetHeader* h) {fHeader= h;}
  // others
  void FindJetsTPC();
  void FindJets();
  void Reset();
  void Init();
  void WriteJHeaderToFile();

 protected:

  AliUA1JetFinder(const AliUA1JetFinder& rUA1Finder);
  AliUA1JetFinder& operator = (const AliUA1JetFinder& ruaf);

  AliUA1JetHeader* fHeader;         // pointer to jet header
  TH2F           * fLego;           //! Lego Histo

  ClassDef(AliUA1JetFinder,1)
};

#endif
