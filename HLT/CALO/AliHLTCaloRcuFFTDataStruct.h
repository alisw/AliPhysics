//-*- Mode: C++ -*-
// $Id: AliHLTPHOSRcuFFTDataStruct.h 34951 2009-09-23 14:35:38Z phille $

#ifndef ALIHLTCALORCUFFTDATASTRUCT_H
#define ALIHLTCALORCUFFTDATASTRUCT_H
#include "Rtypes.h"

//#include "AliHLTCaloConstant.h"
#include "AliHLTCaloConstants.h"

// using namespace CaloHLTConst;

struct AliHLTCaloRcuFFTDataStruct
{
  int fDataLength;

  //  Double_t fGlobalAccumulatedPSD[NGAINS][ALTROMAXSAMPLES];
  //  Double_t fGlobalLastPSD[NGAINS][ALTROMAXSAMPLES];
    

  Double_t fGlobalAccumulatedPSD[fgkNGAINS][fgkALTROMAXSAMPLES];
  Double_t fGlobalLastPSD[fgkNGAINS][fgkALTROMAXSAMPLES];

  //  Double_t fGlobalAccumulatedPSD[AliHLTCaloConstants::fgkNGAINS][AliHLTCaloConstants::fgkALTROMAXSAMPLES];
  // Double_t fGlobalLastPSD[AliHLTCaloConstants::fgkNGAINS][AliHLTCaloConstants::fgkALTROMAXSAMPLES]; 


};

#endif
