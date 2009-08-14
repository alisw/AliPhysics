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

#include "AliHLTPHOSCaloClusterReader.h"
#include "AliHLTPHOSCaloClusterDataStruct.h"
#include "AliHLTPHOSCaloClusterHeaderStruct.h"

AliHLTPHOSCaloClusterReader::AliHLTPHOSCaloClusterReader(): 
  fCurrentClusterPtr(0),
  fIsSetMemory(false),
  fMaxCnt(0),
  fCurrentCnt(0)
{
  //See header file for documentation
}


AliHLTPHOSCaloClusterReader::~AliHLTPHOSCaloClusterReader()
{
  //See header file for documentation
}


AliHLTPHOSCaloClusterDataStruct*   
AliHLTPHOSCaloClusterReader::NextCluster()
{
  // See header file for documentation
  AliHLTPHOSCaloClusterDataStruct* tmpChannelPtr = 0;
  if(fIsSetMemory == false) 
    {
      return 0;
    }
  if(fCurrentCnt < fMaxCnt)
    {
      //cout << "Reader: E = " << fCurrentClusterPtr->fEnergy << endl;
      tmpChannelPtr = reinterpret_cast<AliHLTPHOSCaloClusterDataStruct*>(fCurrentClusterPtr);
      fCurrentCnt++;
      fCurrentClusterPtr = (AliHLTPHOSCaloClusterDataStruct*)((UChar_t*)fCurrentClusterPtr + sizeof(AliHLTPHOSCaloClusterDataStruct) + (fCurrentClusterPtr->fNCells-1)*(sizeof(Float_t) + sizeof(Short_t)));
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
AliHLTPHOSCaloClusterReader::SetMemory(AliHLTPHOSCaloClusterHeaderStruct* clusterHeaderPtr)
{
  //See header file for documentation
  fMaxCnt = clusterHeaderPtr->fNClusters;
  fCurrentClusterPtr = reinterpret_cast<AliHLTPHOSCaloClusterDataStruct*>((UChar_t*)(clusterHeaderPtr) + sizeof(AliHLTPHOSCaloClusterHeaderStruct));
  fIsSetMemory = true;
}


void
AliHLTPHOSCaloClusterReader::Reset()
{
  //See header file for documentation
  fCurrentCnt = 0;
  fIsSetMemory = false;
}
