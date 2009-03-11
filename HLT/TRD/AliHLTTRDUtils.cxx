/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// HLT TRD Utillities Class                                                  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliHLTTRDUtils.h"
#include <TClonesArray.h>
#include "AliHLTTRDTrack.h"
#include "AliHLTTRDCluster.h"

AliHLTUInt32_t AliHLTTRDUtils::AddClustersToOutput(TClonesArray* inClusterArray, AliHLTUInt8_t* outBlockPtr)
{
  AliTRDcluster* cluster = 0;
  AliHLTUInt32_t addedSize = 0;
  //  == OUTdatatype pointer
  AliHLTTRDCluster * outPtr = (AliHLTTRDCluster*)outBlockPtr;
  
  if (inClusterArray){
    Int_t nbEntries  = inClusterArray->GetEntries();
    for (Int_t iCluster = 0; iCluster<nbEntries; iCluster++){
      //cout << "Geting cluster #" << iCluster << endl;
      UInt_t blockSize=0;
      
      cluster = dynamic_cast<AliTRDcluster*>(inClusterArray->At(iCluster));

      AliHLTTRDCluster *hltCluster = new (outPtr) AliHLTTRDCluster(cluster);
      //cout << Form("cluster added at 0x%x (%i)\n",outPtr,outPtr);

      blockSize = sizeof(*hltCluster);

      addedSize += blockSize;
      outBlockPtr += blockSize;
      outPtr = (AliHLTTRDCluster*)outBlockPtr;
    }
  }
  return addedSize;
  
}

AliHLTUInt32_t AliHLTTRDUtils::AddTracksToOutput(TClonesArray* inTrackArray, AliHLTUInt8_t* output)
{
  cout << "\nWriting tracks to the Memory\n ============= \n";
  AliTRDtrackV1* track = 0;
  AliHLTUInt32_t addedSize = 0;
  AliHLTUInt8_t *iterPtr = output;
  AliHLTTRDTrack * outPtr = (AliHLTTRDTrack*)iterPtr;
  
  if (inTrackArray){
    Int_t nbTracks  = inTrackArray->GetEntries();
    for (Int_t iTrack = 0; iTrack<nbTracks; iTrack++){
      AliHLTUInt32_t trackSize=0;
      
      track = dynamic_cast<AliTRDtrackV1*>(inTrackArray->At(iTrack));
      //track->Print();
      
      AliHLTTRDTrack *hltTrack = new (outPtr) AliHLTTRDTrack(track);
      trackSize = hltTrack->GetSize();
      addedSize += trackSize;
      HLTDebug("addedSize %i, trackSize %i", addedSize, trackSize);
      
      iterPtr += trackSize;
      outPtr = (AliHLTTRDTrack*)iterPtr;
    }
  }
  return addedSize;
  
}

/**
 * Read cluster to the TClonesArray from the memory 
 */
//============================================================================
AliHLTUInt32_t AliHLTTRDUtils::ReadClusters(TClonesArray *outArray, void* inputPtr, AliHLTUInt32_t size)
{
  //HLTDebug("\nReading clusters from the Memory\n ============= \n");
  AliHLTTRDCluster * curCluster;
  UInt_t clusterSize = sizeof(AliHLTTRDCluster), curSize = 0;
  Int_t i=0;
  
  curCluster = (AliHLTTRDCluster*) inputPtr;
  while (curSize + clusterSize <= size)
    {
      //  HLTDebug(" fX = %f; fY = %f; fZ = %f", curCluster->fX, curCluster->fY, curCluster->fZ);

      AliTRDcluster* curTRDCluster = new((*outArray)[i]) AliTRDcluster();
      curCluster->ExportTRDCluster(curTRDCluster);
      //      HLTDebug(" fX = %f; fY = %f; fZ = %f", curTRDCluster->GetX(), curTRDCluster->GetY(), curTRDCluster->GetZ());
      curSize += clusterSize; 
      i++;
      curCluster++;
      //cout << " current readed size is " << curSize << "/" << size << endl;
    }
  
  return i;
}

AliHLTUInt32_t AliHLTTRDUtils::ReadTracks(TClonesArray *outArray, void* inputPtr, AliHLTUInt32_t size)
{
  AliHLTUInt8_t* iterPtr = (AliHLTUInt8_t* )inputPtr;
  
  //cout << "\nReading tracks from the Memory\n ============= \n";
  HLTDebug ("\nReading tracks from the Memory\n ============= \n");
  AliHLTTRDTrack * hltTrack;
  AliHLTUInt32_t trackSize = 0, curSize = 0;
  Int_t counter=0;
  
  while (curSize < size)
    {
      hltTrack = (AliHLTTRDTrack*) iterPtr;
      HLTDebug("curSize %i, size %i",curSize, size);
      
      trackSize = hltTrack->GetSize();
      HLTDebug("GetSize() %i", trackSize);

      hltTrack->ReadTrackletsFromMemory(iterPtr + sizeof(AliHLTTRDTrack));

      AliTRDtrackV1* curTRDTrack = new((*outArray)[counter]) AliTRDtrackV1();
      hltTrack->ExportTRDTrack(curTRDTrack);
      
      curSize += trackSize; 
      iterPtr += trackSize;
      counter++;
    }

  //CheckTrackArray(outArray);
  
  return counter;
}

