#ifndef ALIPHOSPREPROCESSOR_H
#define ALIPHOSPREPROCESSOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
// Class AliPHOSPreprocessor
///////////////////////////////////////////////////////////////////////////////


#include "AliPreprocessor.h"
#include "TFile.h"

class AliPHOSPreprocessor : public AliPreprocessor {
public:

  AliPHOSPreprocessor();
  AliPHOSPreprocessor(AliShuttleInterface* shuttle);

protected:

  virtual UInt_t Process(TMap* valueSet);
  Bool_t ProcessLEDRun();
  Bool_t FindBadChannelsEmc();
  Bool_t CalibrateEmc();
  Float_t HG2LG(Int_t module, Int_t X, Int_t Z, TFile* f);

  ClassDef(AliPHOSPreprocessor,1);

};

#endif
