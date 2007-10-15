#ifndef ALIHLTALTROBUNCH_H
#define ALIHLTALTROBUNCH_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "Rtypes.h"
#include <iostream>

using namespace std;


class AliHLTAltroBunch
{
 public:
  AliHLTAltroBunch();
  ~ AliHLTAltroBunch();
  UInt_t *fData;
  Int_t fBunchDataSize;
  Int_t fBunchSize;
  unsigned int fEndTimeBin;
 private:
  unsigned int fStartTimeBin;
};


#endif

