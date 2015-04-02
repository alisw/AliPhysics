#ifndef ALIPHOSCPVPREPROCESSOR_H
#define ALIPHOSCPVPREPROCESSOR_H
/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////////////////////
// Class AliPHOSCpvPreprocessor - CPV preprocessor.
///////////////////////////////////////////////////////////////////////////////


#include "AliPreprocessor.h"
#include "TFile.h"

class AliPHOSCpvPreprocessor : public AliPreprocessor {
public:

  AliPHOSCpvPreprocessor();
  AliPHOSCpvPreprocessor(AliShuttleInterface* shuttle);

protected:

  virtual UInt_t Process(TMap* valueSet);

  ClassDef(AliPHOSCpvPreprocessor,2);

};

#endif
