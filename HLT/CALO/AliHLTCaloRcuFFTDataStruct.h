//-*- Mode: C++ -*-
// $Id: AliHLTPHOSRcuFFTDataStruct.h 34951 2009-09-23 14:35:38Z phille $

#ifndef ALIHLTCALORCUFFTDATASTRUCT_H
#define ALIHLTCALORCUFFTDATASTRUCT_H
#include "Rtypes.h"

//#include "AliHLTCaloConstants.h"

#include "AliHLTCaloConstants.h"

// using AliHLTCaloConstants::TEST1;
// using AliHLTCaloConstants::TEST2;

using CALO::NGAINS;
using CALO::ALTROMAXSAMPLES;

struct AliHLTCaloRcuFFTDataStruct
{
  int fDataLength;
  
  //  Double_t fTest1[AliHLTCaloConstants::TEST1];
  //  Double_t fTest2[AliHLTCaloConstants::TEST2];

  Double_t fGlobalAccumulatedPSD[NGAINS][ALTROMAXSAMPLES];
  Double_t fGlobalLastPSD[NGAINS][ALTROMAXSAMPLES];
};

#endif
