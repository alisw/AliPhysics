
/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        *
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors: Oystein Djuvsland                                     *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#ifndef ALIHLTPHOSDIGITDATASTRUCT_H
#define ALIHLTPHOSDIGITDATASTRUCT_H

#include "AliHLTPHOSBase.h"

struct AliHLTPHOSDigitDataStruct
{
  Int_t fX;
  Int_t fZ;
  Int_t fModule;
  Float_t fAmplitude;
  Float_t fTime;
  Float_t fEnergy;
  Int_t fGain;
  
  Int_t fData[512];

  Int_t fCrazyness; 
  Float_t fBaseline;

  void SetRawData(Int_t *data)
  {
    for(Int_t i = 0; i < 512; i++)
      {
	fData[i] = data[i];
      }
  }
};

#endif

