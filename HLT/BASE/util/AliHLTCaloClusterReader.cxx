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

#include "AliHLTCaloClusterReader.h"
#include "AliHLTCaloClusterDataStruct.h"

AliHLTCaloClusterReader::AliHLTCaloClusterReader(): 
  fCurrentClusterPtr(0),
  fIsSetMemory(false),
  fMaxCnt(0),
  fCurrentCnt(0)
{
  //See header file for documentation
}


AliHLTCaloClusterReader::~AliHLTCaloClusterReader()
{
  //See header file for documentation
}


AliHLTCaloClusterDataStruct*   
AliHLTCaloClusterReader::NextCluster()
{
  // See header file for documentation
  AliHLTCaloClusterDataStruct* tmpChannelPtr = 0;

  // Check if we have set our memory buffer
  if(fIsSetMemory == false) 
    {
      return 0;
    }

  // Check if we don't read more clusters than we have
  if(fCurrentCnt < fMaxCnt)
    {
      // So, we get our cluster
      tmpChannelPtr = reinterpret_cast<AliHLTCaloClusterDataStruct*>(fCurrentClusterPtr);
      
      // increment the number of clusters read
      fCurrentCnt++;
      
      // Move our cluster pointer to give us a cluster next time
      // The increment is defined as the size of the cluster data struct 
      // + the number of cells minus the one already included in the struct
      // times the amount of data for each cell (absolute ID (short) and amplitude fraction (float))
      fCurrentClusterPtr = (AliHLTCaloClusterDataStruct*)((UChar_t*)fCurrentClusterPtr 
							      + sizeof(AliHLTCaloClusterDataStruct) 
							      + (fCurrentClusterPtr->fNCells-1)*(sizeof(Float_t) + sizeof(Short_t)));
      // return the cluster
      return tmpChannelPtr;
    }
  // if we have read our clusters we reset our memory pointer and returns 0;
  else
    {
      Reset();
      return 0;
    }

  // will never happen, but anyway...
  return 0;
}

bool
AliHLTCaloClusterReader::GetCell(AliHLTCaloClusterDataStruct *clusterPtr, UShort_t &cellId, Double32_t &cellAmp, UInt_t index)
{
  // See header file for documentation

  // check if the index is within bounds
  if(index < clusterPtr->fNCells)
    {
      // the absolute ID is located in the memory address of the first ID plus the size of the pair of cell properties times the index
      cellId = *(UShort_t*)((UChar_t*)(&(clusterPtr->fCellsAbsId)) + index * (sizeof(Short_t) + sizeof(Float_t)));
      // the amplitude fraction is located in the memory address of the first ID plus the size of the pair of cell properties times the index
      cellAmp = *(Float_t*)((UChar_t*)(&(clusterPtr->fCellsAmpFraction)) + index * (sizeof(Short_t) + sizeof(Float_t)));
      return true;
    }
  else return false;
}

void
AliHLTCaloClusterReader::SetMemory(const AliHLTCaloClusterHeaderStruct* clusterHeaderPtr)
{
  //See header file for documentation
  fMaxCnt = clusterHeaderPtr->fNClusters;
  fCurrentClusterPtr = reinterpret_cast<AliHLTCaloClusterDataStruct*>((UChar_t*)(clusterHeaderPtr) + sizeof(AliHLTCaloClusterHeaderStruct));
  fIsSetMemory = true;
}


void
AliHLTCaloClusterReader::Reset()
{
  //See header file for documentation
  fCurrentCnt = 0;
  fIsSetMemory = false;
}
