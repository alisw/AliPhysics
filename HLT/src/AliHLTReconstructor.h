#ifndef ALIHLTRECONSTRUCTOR_H
#define ALIHLTRECONSTRUCTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#ifdef use_reconstruction
#include "AliReconstructor.h"
#include "AliL3Logging.h"

class AliHLTReconstructor: public AliReconstructor {
public:
  AliHLTReconstructor(): AliReconstructor() {
    AliL3Log::fgLevel=AliL3Log::kWarning;
  };
  virtual ~AliHLTReconstructor() {};
  virtual void         Reconstruct(AliRunLoader* runLoader) const;
  virtual void         FillESD(AliRunLoader* runLoader, AliESD* esd) const;

private:
  void ReconstructWithConformalMapping(AliRunLoader* runLoader,Int_t iEvent) const;
  void ReconstructWithHoughTransform(AliRunLoader* runLoader,Int_t iEvent) const;
  void FillESDforConformalMapping(AliESD* esd,Int_t iEvent) const;
  void FillESDforHoughTransform(AliESD* esd,Int_t iEvent) const;

  ClassDef(AliHLTReconstructor, 0)   // class for the TPC reconstruction
};
#endif

#endif
