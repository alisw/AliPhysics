//-*- Mode: C++ -*-
// $Id: AliHLTPHOSRcuFFTDataStruct.h 34951 2009-09-23 14:35:38Z phille $

#ifndef ALIHLTCALORCUFFTDATASTRUCT_H
#define ALIHLTCALORCUFFTDATASTRUCT_H
#include "Rtypes.h"
#include "AliHLTCaloConstants.h"


struct AliHLTCaloRcuFFTDataStruct
{
  int fDataLength;
  //   Double_t fTest [TEST];
  Double_t fGlobalAccumulatedPSD[NGAINS][ALTROMAXSAMPLES];
  Double_t fGlobalLastPSD[NGAINS][ALTROMAXSAMPLES];
};

#endif
