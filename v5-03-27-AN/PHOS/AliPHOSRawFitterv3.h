#ifndef ALIPHOSRAWFITTERV3_H
#define ALIPHOSRAWFITTERV3_H
/* Copyright(c) 2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                          */

/* $Id: $ */

// This class extracts the PHOS "digits" of current event
// (amplitude,time, position,gain) from the raw stream 
// provided by AliRawReader. 
// Uses FastFitting procedure to evaluate time and estimate sample quality

#include "AliPHOSRawFitterv0.h"

class AliPHOSRawFitterv3 : public AliPHOSRawFitterv0 {

public:

  AliPHOSRawFitterv3();
  AliPHOSRawFitterv3(const AliPHOSRawFitterv3& rawFitter);
  AliPHOSRawFitterv3& operator = (const AliPHOSRawFitterv3& rawFitter);
  virtual ~AliPHOSRawFitterv3();

  virtual Bool_t Eval(const UShort_t *signal, Int_t sigStart, Int_t sigLength);

private:
  
  ClassDef(AliPHOSRawFitterv3,1)
};

#endif
