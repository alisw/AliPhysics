//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTPHOSRCUFFTDATASTRUCT_H
#define ALIHLTPHOSRCUFFTDATASTRUCT_H


#include "Rtypes.h"

//#include "AliHLTPHOSConstant.h"
#include "AliHLTPHOSConstants.h"

//using namespace PhosHLTConst;

struct AliHLTPHOSRcuFFTDataStruct
{
  int fDataLength;
  Double_t fGlobalAccumulatedPSD[NGAINS][ALTROMAXSAMPLES];
  Double_t fGlobalLastPSD[NGAINS][ALTROMAXSAMPLES];
  //  Double_t fDummy[64][64][NGAINS];

};

#endif
