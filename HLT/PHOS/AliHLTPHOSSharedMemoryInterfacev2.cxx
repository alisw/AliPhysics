/**************************************************************************
 * This file is property of and copyright by the Experimental Nuclear     *
 * Physics Group, Dep. of Physics                                         *
 * University of Oslo, Norway, 2007                                       *
 *                                                                        *
 * Author: Per Thomas Hille <perthi@fys.uio.no> for the ALICE HLT Project.*
 * Contributors are mentioned in the code where appropriate.              *
 * Please report bugs to perthi@fys.uio.no                                *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/


#include "AliHLTPHOSSharedMemoryInterfacev2.h"
#include "AliHLTPHOSChannelDataHeaderStruct.h"
#include "AliHLTPHOSChannelDataStruct.h"
#include "AliHLTLogging.h"


AliHLTPHOSSharedMemoryInterfacev2::AliHLTPHOSSharedMemoryInterfacev2(): 
  fCurrentChannel(0),
  fChannelDataPtr(0),
  fIsSetMemory(false),
  fHasRawData(false),
  fMaxCnt(0),
  fCurrentCnt(0)
{
  
}


AliHLTPHOSSharedMemoryInterfacev2::~AliHLTPHOSSharedMemoryInterfacev2()
{

}


AliHLTPHOSChannelDataStruct*   
AliHLTPHOSSharedMemoryInterfacev2::NextChannel()
{
  // Comment
  AliHLTPHOSChannelDataStruct* tmpChannelPtr = 0;
  if(fCurrentCnt < fMaxCnt)
    {
      tmpChannelPtr = reinterpret_cast<AliHLTPHOSChannelDataStruct*>(fChannelDataPtr);
      fCurrentCnt++;
      fChannelDataPtr += sizeof(AliHLTPHOSChannelDataStruct);
      if(fHasRawData == true)
	{
	  //	  HLTDebug("Raw data interface not yet implemented, ignoring raw data");
	}
      return tmpChannelPtr;
    }
  else
    {
      Reset();
      return 0;
    }
  return 0;
}

void
AliHLTPHOSSharedMemoryInterfacev2::SetMemory(AliHLTPHOSChannelDataHeaderStruct* channelDataHeaderPtr)
{
  
  //Shutting up rule checker
  fHasRawData = channelDataHeaderPtr->fHasRawData; 
  fMaxCnt = channelDataHeaderPtr->fNChannels;
  fChannelDataPtr = reinterpret_cast<AliHLTUInt8_t*>(channelDataHeaderPtr) + sizeof(AliHLTPHOSChannelDataHeaderStruct);
  if(fHasRawData == true)
    {
      //  HLTWarning("Raw data interface not yet implemented, ignoring raw data");
    }
  fIsSetMemory = true;
}


void
AliHLTPHOSSharedMemoryInterfacev2::Reset()
{
  //Shutting up rule checker
  fCurrentCnt = 0;
  fIsSetMemory = false;
  fHasRawData = false;
}
