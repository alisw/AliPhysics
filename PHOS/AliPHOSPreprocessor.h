#ifndef ALIPHOSPREPROCESSOR_H
#define ALIPHOSPREPROCESSOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
// Class AliPHOSPreprocessor
///////////////////////////////////////////////////////////////////////////////


#include "AliPreprocessor.h"

class AliPHOSPreprocessor : public AliPreprocessor {
public:

  AliPHOSPreprocessor();
  AliPHOSPreprocessor(const char* detector, AliShuttleInterface* shuttle);

protected:

  virtual void Initialize(Int_t run, UInt_t startTime, UInt_t endTime);
  virtual UInt_t Process(TMap* valueSet);

  ClassDef(AliPHOSPreprocessor,0);

};

#endif
