#ifndef AliMFTDigitizer_H
#define AliMFTDigitizer_H

/* Copyright (c) 1998-2001, ALICE Experiment at CERN, All rights reserved *
 * See cxx source for full Copyright notice                               */

//====================================================================================================================================================
//
//      Digitizer class for the ALICE Muon Forward Tracker
//
//      Contact author: antonio.uras@cern.ch
//
//====================================================================================================================================================

#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliRunDigitizer.h"
#include "AliLoader.h"
#include "AliLog.h"
#include "AliMFTDigit.h"
#include "AliMFTSegmentation.h"
#include "TObjArray.h"
#include "TClonesArray.h"
#include "AliDigitizer.h"

//====================================================================================================================================================

class AliMFTDigitizer : public AliDigitizer {

public:

  AliMFTDigitizer();
  AliMFTDigitizer(AliRunDigitizer *pRunDig);
  virtual ~AliMFTDigitizer() { }

  void Exec(Option_t *option);
  void SDigits2Digits(TClonesArray *pSDigitList, TObjArray *pDigitLst);

  void MergeDigits(AliMFTDigit *mainDig, AliMFTDigit *digToSum);
  
  
protected:
 
  static const Int_t fNMaxPlanes = 20;     // max number of MFT planes
  static const Int_t fNMaxMCTracks = 10;   // max MC tracks sharing a digit  

  Int_t fNPlanes;    
  
  AliMFTSegmentation *fSegmentation;

private:

  AliMFTDigitizer (const AliMFTDigitizer& mftDigitizer);             // dummy copy constructor
  AliMFTDigitizer &operator=(const AliMFTDigitizer& mftDigitizer);   // dummy assignment operator

  ClassDef(AliMFTDigitizer,1)
    
};

//====================================================================================================================================================

#endif


