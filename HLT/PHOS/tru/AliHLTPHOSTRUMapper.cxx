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

#include "AliHLTPHOSTRUMapper.h"

#define N_CHANNELS_ADC 8
#define N_ADCS_TRU 14
#define N_ADCS_MASK 2
#define N_XCRYSTALS_TRU 16
#define N_ZCRYSTALS_TRU 28
#define N_ZCRYSTALS 56
#define N_XCRYSTALS 64
#define N_XCRYSTALS_ADC 16
#define N_ZCRYSTALS_ADC 2
#define CHANNEL_MASK 1

AliHLTPHOSTRUMapper::AliHLTPHOSTRUMapper()
{
}


AliHLTPHOSTRUMapper::~AliHLTPHOSTRUMapper()
{
}

void
AliHLTPHOSTRUMapper::GetChannelMask(x, z, Short_t *maskArray)
{
  Int_t chNb = -1;
  Int_t ADCNb = -1;
  
  for(Int_t i = 0; i < N_ADCS_TRU/N_ADCS_MASK; i++)
  {
    maskArray[i] = 0;
  }
  
  chNb = fChannelLookupArrayPtr[x];
  
  ADCNb = N_ADCS_TRU - z/N_ZCRYSTALS_ADC;
  
  maskArray[(ADCNb-1)/N_ADCS_MASK] = (CHANNEL_MASK << chNb)  << ((1 - ADCNb % N_ADCS_MASK) * N_CHANNELS_ADC );
  
}

Int_t
AliHLTPHOSTRUMapper::GetTRUMask(x, z, Short_t *maskArray)
{
  Int_t TRUNb = x/N_XCRYSTALS_TRU + z/N_ZCRYSTALS_TRU;
  
  GetChannelMask(x, z, maskArray);
  
  return TRUNb;
}

void 
AliHLTPHOSTRUMapper::GetLookUpTable()
{
  ifstream infile;
  
  inFile.open("channelLookup.txt");
  
  for(Int_t i = 0; i < 8; i++)
  {
    inFile >> fChannelLookupArrayPtr[i];  
  }
    
  inFile.close();
}
  
  
  
  
  
  
  
  