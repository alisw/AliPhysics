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
#include "AliHLTCaloDigitDataStruct.h"

AliHLTCaloClusterReader::AliHLTCaloClusterReader(): 
  fCurrentClusterPtr(0),
  fIsSetMemory(false),
  fMaxCnt(0),
  fCurrentCnt(0),
  fDigitsPointer(0),
  fNDigits(0)
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
      
//      printf("CR: Energy: %f, cluster pointer: %x\n", tmpChannelPtr->fEnergy, tmpChannelPtr);
      
      // increment the number of clusters read
      fCurrentCnt++;
      
      // Move our cluster pointer to give us a cluster next time
      // The increment is defined as the size of the cluster data struct 
      // + the number of cells minus the one already included in the struct
      // times the amount of data for each cell (absolute ID (short) and amplitude fraction (float))
//      fCurrentClusterPtr = (AliHLTCaloClusterDataStruct*)((UChar_t*)fCurrentClusterPtr 
//							      + sizeof(AliHLTCaloClusterDataStruct) 
//							      + (fCurrentClusterPtr->fNCells-1)*(sizeof(Float_t) + sizeof(Short_t)));
      fCurrentClusterPtr = reinterpret_cast<AliHLTCaloClusterDataStruct*>(reinterpret_cast<UChar_t*>(fCurrentClusterPtr)
							      + sizeof(AliHLTCaloClusterDataStruct) 
							      + sizeof(AliHLTCaloCellDataStruct)*(fCurrentClusterPtr->fNCells - 1));;
							      //+ (fCurrentClusterPtr->fNCells-1)*(sizeof(Float_t) + sizeof(Short_t))
							      //- sizeof(Short_t)); //TODO: Why?;
						
							      
      // return the cluster
//      printf("CR: Energy: %f, number of cells: %d, cluster pointer: %x, next cluster pointer: %x\n", tmpChannelPtr->fEnergy, tmpChannelPtr->fNCells, tmpChannelPtr, fCurrentClusterPtr);
//      if(fCurrentCnt < fMaxCnt-1) printf("CR: Next cluster Energy: %f, number of cells: %d, cluster pointer: %x\n", fCurrentClusterPtr->fEnergy, fCurrentClusterPtr->fNCells,  fCurrentClusterPtr);

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
       cellId = (&(clusterPtr->fCaloCells))[index].fCellsAbsId;
       cellAmp = (&(clusterPtr->fCaloCells))[index].fCellsAmpFraction;
       return true;
    }
  else return false;
}

void
AliHLTCaloClusterReader::SetMemory(const AliHLTCaloClusterHeaderStruct* clusterHeaderPtr)
{
  //See header file for documentation
  //printf("CR: Number of clusters in event: %d, Number of digits in event: %d\n", clusterHeaderPtr->fNClusters, clusterHeaderPtr->fNDigits);
  fMaxCnt = clusterHeaderPtr->fNClusters;
  fCurrentClusterPtr = reinterpret_cast<AliHLTCaloClusterDataStruct*>((UChar_t*)(clusterHeaderPtr) + sizeof(AliHLTCaloClusterHeaderStruct) + sizeof(AliHLTCaloDigitDataStruct)*clusterHeaderPtr->fNDigits);
  fNDigits = clusterHeaderPtr->fNDigits;
  fDigitsPointer = reinterpret_cast<AliHLTCaloDigitDataStruct*>((UChar_t*)(clusterHeaderPtr) + sizeof(AliHLTCaloClusterHeaderStruct));
  fIsSetMemory = true;
}


void
AliHLTCaloClusterReader::Reset()
{
  //See header file for documentation
  fCurrentCnt = 0;
  fIsSetMemory = false;
}
