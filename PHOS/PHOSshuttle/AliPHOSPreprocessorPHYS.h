#ifndef ALIPHOSPREPROCESSORPHYS_H
#define ALIPHOSPREPROCESSORPHYS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////////////////////
// Class AliPHOSPreprocessorPHYS
///////////////////////////////////////////////////////////////////////////////


#include "AliPreprocessor.h"

class AliPHOSPreprocessorPHYS : public AliPreprocessor {
 public:

  AliPHOSPreprocessorPHYS();
  AliPHOSPreprocessorPHYS(AliShuttleInterface* shuttle);

 protected:

  virtual UInt_t Process(TMap* valueSet);
  Bool_t CalibratePhys();

 private:

  ClassDef(AliPHOSPreprocessorPHYS,1)  // Preprocessor for PHYS event
};
#endif
