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

class TList;
class AliPHOSEmcBadChannelsMap;
class AliPHOSEmcCalibData;

class AliPHOSPreprocessor : public AliPreprocessor {
 public:

  AliPHOSPreprocessor();
  AliPHOSPreprocessor(AliShuttleInterface* shuttle);

  Bool_t ProcessDCS() { return kFALSE; }

 protected:

  virtual UInt_t Process(TMap* valueSet);
  Bool_t ProcessLEDRun();
  Bool_t FindBadChannelsEmc();
  Bool_t CalibrateEmc();
  Float_t HG2LG(Int_t module, Int_t X, Int_t Z, TFile* f);

 private:

  Bool_t DoCalibrateEmc(Int_t system, TList* sources, const AliPHOSEmcBadChannelsMap* badMap, AliPHOSEmcCalibData& calibData);
  Bool_t DoFindBadChannelsEmc(Int_t system, TList* sources, AliPHOSEmcBadChannelsMap& badMap);

  ClassDef(AliPHOSPreprocessor,2);

};

#endif
