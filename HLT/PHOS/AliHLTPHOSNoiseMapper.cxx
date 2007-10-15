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

#include "AliHLTPHOSNoiseMapper.h"
#include "AliHLTPHOSBase.h"
#include "AliHLTPHOSDigitDataStruct.h"
#include "AliHLTPHOSDigitContainerDataStruct.h"
#include "AliHLTPHOSConstants.h"

using namespace std;

AliHLTPHOSNoiseMapper::AliHLTPHOSNoiseMapper()
 : AliHLTPHOSBase(),
 //fChannelArray(0),
 fNoiseThreshold(10)
{
  for(Int_t x = 0; x < N_XCOLUMNS_MOD; x++)
  {
    for(Int_t z = 0; z < N_ZROWS_MOD; z++)
    {
      for(Int_t gain = 0; gain < N_GAINS; gain++)
      {
	fChannelArray[x][z][gain] = 0;
      }
    }
  }
}


AliHLTPHOSNoiseMapper::~AliHLTPHOSNoiseMapper()
{
}


void
AliHLTPHOSNoiseMapper::MapNoisyChannels(AliHLTPHOSDigitContainerDataStruct *digitContainerPtr)
{
  AliHLTPHOSDigitDataStruct *digitPtr = 0;
  for(Int_t i = 0; i < digitContainerPtr->fNDigits; i++)
  {
    digitPtr = &(digitContainerPtr->fDigitDataStruct[i]);
    if(digitPtr->fAmplitude > fNoiseThreshold)
    {
      fChannelArray[digitPtr->fX][digitPtr->fZ][digitPtr->fGain]++;
    }
  }
}

void
AliHLTPHOSNoiseMapper::GetChannelArray(Int_t channelArray[N_XCOLUMNS_MOD][N_ZROWS_MOD][N_GAINS])
{
  for(Int_t x = 0; x < N_XCOLUMNS_MOD; x++)
  {
    for(Int_t z = 0; z < N_ZROWS_MOD; z++)
    {
      for(Int_t gain = 0; gain < N_GAINS; gain++)
      {
	channelArray[x][z][gain] = fChannelArray[x][z][gain];
      }
    }
  }
}
