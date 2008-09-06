#ifndef ALIHLTPHOSRCUFFTDATASTRUCT_H
#define ALIHLTPHOSRCUFFTDATASTRUCT_H


#include "Rtypes.h"
#include "AliHLTPHOSConstants.h"

using namespace PhosHLTConst;

struct AliHLTPHOSRcuFFTDataStruct
{
  int fDataLength;
  Double_t fGlobalAccumulatedPSD[N_GAINS][ALTRO_MAX_SAMPLES];
  Double_t fGlobalLastPSD[N_GAINS][ALTRO_MAX_SAMPLES];
  //  Double_t fDummy[64][64][N_GAINS];

};

#endif
