
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


#ifndef ALIHLTPHOSNOISEMAPPER_H
#define ALIHLTPHOSNOISEMAPPER_H

#include <AliHLTPHOSBase.h>

class AliHLTPHOSDigitContainerDataStruct;

/**
		

	@author Oystein Djuvsland <oystein.djuvsland@gmail.com>
*/
class AliHLTPHOSNoiseMapper : public AliHLTPHOSBase
{
public:
    AliHLTPHOSNoiseMapper();

    ~AliHLTPHOSNoiseMapper();
    
    void MapNoisyChannels(AliHLTPHOSDigitContainerDataStruct *fDigitContainer);
    
    //void SetChannelArray ( Int_t* val ) { fChannelArrayPtr = val; }
    void GetChannelArray(Int_t channelArray[N_XCOLUMNS_MOD][N_ZROWS_MOD][N_GAINS]);
    void SetNoiseThreshold ( Float_t val ) { fNoiseThreshold = val; }
    Float_t GetNoiseThreshold() const { return fNoiseThreshold; }
    
	
private:
  
  Int_t fChannelArray[N_XCOLUMNS_MOD][N_ZROWS_MOD][N_GAINS];
  Float_t fNoiseThreshold;
  
};

#endif
