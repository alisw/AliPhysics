// @(#) $Id$

#ifndef ALIHLTRECONSTRUCTOR_H
#define ALIHLTRECONSTRUCTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#ifdef use_reconstruction
#include "AliReconstructor.h"
#include "AliL3Logging.h"

class AliHLTReconstructor: public AliReconstructor {
public:
  AliHLTReconstructor();
  AliHLTReconstructor(Bool_t doTracker, Bool_t doHough);
  virtual ~AliHLTReconstructor();

  virtual void         Reconstruct(AliRunLoader* runLoader) const;
  virtual void         FillESD(AliRunLoader* runLoader, AliESD* esd) const;
  void SetDoBench(Bool_t b){fDoBench=b;}
  void SetDoCleanup(Bool_t b){fDoCleanUp=b;}

private:
  void ReconstructWithConformalMapping(AliRunLoader* runLoader,Int_t iEvent) const;
  void ReconstructWithHoughTransform(AliRunLoader* runLoader,Int_t iEvent) const;
  void FillESDforConformalMapping(AliESD* esd,Int_t iEvent) const;
  void FillESDforHoughTransform(AliESD* esd,Int_t iEvent) const;


  Bool_t fDoHough;   //do the hough transform
  Bool_t fDoTracker; //do the standard conformal tracker
  Bool_t fDoBench;   //store the benchmark results
  Bool_t fDoCleanUp; //delete tmp tracking files

  ClassDef(AliHLTReconstructor, 0)   // class for the TPC reconstruction
};
#endif

#endif
