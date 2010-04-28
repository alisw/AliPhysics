// $Id$
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
#include "AliHLTTRDTracklet.h"
#include "AliHLTTRDCluster.h"
#include "AliHLTExternalTrackParam.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliPID.h"

ClassImp(AliHLTTRDUtils)

AliHLTUInt32_t AliHLTTRDUtils::AddClustersToOutput(TClonesArray* inClusterArray, AliHLTUInt8_t* outBlockPtr, Int_t nTimeBins)
{
  AliTRDcluster* cluster = 0;
  AliHLTUInt32_t addedSize = 0;
  //  == OUTdatatype pointer
  AliHLTTRDCluster* outPtr = (AliHLTTRDCluster*)outBlockPtr;
  
  if (inClusterArray){
    Int_t nbEntries  = inClusterArray->GetEntries();
    for (Int_t iCluster = 0; iCluster<nbEntries; iCluster++){
      //cout << "Geting cluster #" << iCluster << endl;
      
      cluster = (AliTRDcluster*)(inClusterArray->At(iCluster));

      AliHLTTRDCluster *hltCluster = new (outPtr) AliHLTTRDCluster(cluster);
      //cout << Form("cluster added at 0x%x (%i)\n",outPtr,outPtr);

      addedSize += sizeof(*hltCluster);
      outBlockPtr += sizeof(*hltCluster);
      outPtr = (AliHLTTRDCluster*)outBlockPtr;
    }
  }

  Int_t *TBptr = (Int_t*)outPtr;
  *TBptr = nTimeBins;
  
  addedSize += sizeof(*TBptr);
  outBlockPtr += sizeof(*TBptr);

  return addedSize;
  
}

AliHLTUInt32_t AliHLTTRDUtils::AddTracksToOutput(TClonesArray* inTrackArray, AliHLTUInt8_t* output, Int_t nTimeBins)
{
  //cout << "\nWriting tracks to the Memory\n ============= \n";

  Int_t *TBptr = (Int_t*)output;
  *TBptr = nTimeBins;

  AliTRDtrackV1* track = 0;
  AliHLTUInt32_t addedSize = sizeof(*TBptr);
  AliHLTUInt8_t* iterPtr = output+addedSize;
  AliHLTTRDTrack* outPtr = (AliHLTTRDTrack*)iterPtr;

  if (inTrackArray){
    Int_t nbTracks  = inTrackArray->GetEntries();
    for (Int_t iTrack = 0; iTrack<nbTracks; iTrack++){
      AliHLTUInt32_t trackSize=0;
      
      track = (AliTRDtrackV1*)(inTrackArray->At(iTrack));
      //track->Print();
      
      AliHLTTRDTrack *hltTrack = new (outPtr) AliHLTTRDTrack(track);
      trackSize = hltTrack->GetSize();
      addedSize += trackSize;
      //HLTDebug("addedSize %i, trackSize %i", addedSize, trackSize);
      
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
AliHLTUInt32_t AliHLTTRDUtils::ReadClusters(TClonesArray *outArray, void* inputPtr, AliHLTUInt32_t size, Int_t* nTimeBins)
{
  //HLTDebug("\nReading clusters from the Memory\n ============= \n");
  AliHLTTRDCluster * curCluster;
  UInt_t clusterSize = sizeof(AliHLTTRDCluster), curSize = 0;
  Int_t counter=outArray->GetEntriesFast();

  if(nTimeBins){
    *nTimeBins=*(Int_t*)(((AliHLTUInt8_t*)inputPtr)+size-sizeof(Int_t));
    //HLTDebug("Reading number of time bins from input block: %d", *nTimeBins);
  }
  size-=sizeof(*nTimeBins);

  curCluster = (AliHLTTRDCluster*) inputPtr;
  while (curSize < size)
    {
      //HLTDebug(" fX = %f; fY = %f; fZ = %f", curCluster->fX, curCluster->fY, curCluster->fZ);

      AliTRDcluster* curTRDCluster = new((*outArray)[counter]) AliTRDcluster();
      curCluster->ExportTRDCluster(curTRDCluster);
      //HLTDebug(" fX = %f; fY = %f; fZ = %f", curTRDCluster->GetX(), curTRDCluster->GetY(), curTRDCluster->GetZ());
      curSize += clusterSize; 
      counter++;
      curCluster++;
      //cout << " current readed size is " << curSize << "/" << size << endl;
    }
  
  return counter;
}

AliHLTUInt32_t AliHLTTRDUtils::ReadTracks(TClonesArray *outArray, void* inputPtr, AliHLTUInt32_t size, Int_t* nTimeBins)
{
  if(nTimeBins){
    *nTimeBins=*(Int_t*)inputPtr;
    //HLTDebug("Reading number of time bins from input block: %d", *nTimeBins);
  }

  AliHLTUInt8_t* iterPtr = ((AliHLTUInt8_t*)inputPtr)+sizeof(*nTimeBins);
  
  //cout << "\nReading tracks from the Memory\n ============= \n";
  //HLTDebug ("\nReading tracks from the Memory\n ============= \n");
  AliHLTTRDTrack * hltTrack;
  AliHLTUInt32_t trackSize = 0, curSize = sizeof(*nTimeBins);
  Int_t counter=outArray->GetEntriesFast();
  
  while (curSize < size)
    {
      hltTrack = (AliHLTTRDTrack*) iterPtr;
      //HLTDebug("curSize %i, size %i",curSize, size);
      
      trackSize = hltTrack->GetSize();
      //HLTDebug("GetSize() %i", trackSize);

      // hltTrack->ReadTrackletsFromMemory(iterPtr + sizeof(AliHLTTRDTrack));

      AliTRDtrackV1* curTRDTrack = new((*outArray)[counter]) AliTRDtrackV1();
      hltTrack->ExportTRDTrack(curTRDTrack);
      
      curSize += trackSize; 
      iterPtr += trackSize;
      counter++;
    }

  //CheckTrackArray(outArray);
  
  return counter;
}

AliHLTUInt32_t AliHLTTRDUtils::AddESDToOutput(const AliESDEvent* const esd, AliHLTUInt8_t* const outBlockPtr)
{
  AliESDtrack* esdTrack = 0;
  AliHLTUInt8_t* iterPtr = outBlockPtr;

  AliHLTTracksData* trksData = new(iterPtr) AliHLTTracksData;
  iterPtr += sizeof(AliHLTTracksData);
  trksData->fCount=0;
  
  if(esd){
    Double_t pid[5];
    for(Int_t i=0; i<esd->GetNumberOfTracks(); i++){
      esdTrack=esd->GetTrack(i);
      if(!esdTrack)continue;
      AliHLTExternalTrackParam* trk = new(iterPtr) AliHLTExternalTrackParam;
      iterPtr += sizeof(AliHLTExternalTrackParam);
      trk->fAlpha = esdTrack->GetAlpha();
      trk->fX = esdTrack->GetX();
      trk->fY = esdTrack->GetY();
      trk->fZ = esdTrack->GetZ();
      trk->fSinPsi = esdTrack->GetSnp();
      trk->fTgl = esdTrack->GetTgl();
      trk->fq1Pt = esdTrack->GetSigned1Pt();
      trk->fC[0] = esdTrack->GetSigmaY2();
      trk->fC[1] = esdTrack->GetSigmaZY();
      trk->fC[2] = esdTrack->GetSigmaZ2();
      trk->fC[3] = esdTrack->GetSigmaSnpY();
      trk->fC[4] = esdTrack->GetSigmaSnpZ();
      trk->fC[5] = esdTrack->GetSigmaSnp2();
      trk->fC[6] = esdTrack->GetSigmaTglY();
      trk->fC[7] = esdTrack->GetSigmaTglZ();
      trk->fC[8] = esdTrack->GetSigmaTglSnp();
      trk->fC[9] = esdTrack->GetSigmaTgl2();
      trk->fC[10] = esdTrack->GetSigma1PtY();
      trk->fC[11] = esdTrack->GetSigma1PtZ();
      trk->fC[12] = esdTrack->GetSigma1PtSnp();
      trk->fC[13] = esdTrack->GetSigma1PtTgl();
      trk->fC[14] = esdTrack->GetSigma1Pt2();
      esdTrack->GetTRDpid(pid);
      //trk->fTRDpid = pid[AliPID::kElectron]; ...
      trk->fNPoints = 0;
      trksData->fCount++;
    }
  }
  return iterPtr - outBlockPtr;
}

void AliHLTTRDUtils::EmulateHLTClusters(TClonesArray* clusterArray)
{
  AliHLTUInt32_t estimatedSize = (clusterArray->GetEntriesFast()+1)*sizeof(AliHLTTRDCluster);
  AliHLTUInt8_t* pBlock = (AliHLTUInt8_t*)malloc(estimatedSize);
  AliHLTUInt32_t size = AddClustersToOutput(clusterArray, pBlock);
  clusterArray->Delete();
  ReadClusters(clusterArray, pBlock, size);
  free(pBlock);
}

void AliHLTTRDUtils::EmulateHLTTracks(TClonesArray* trackArray)
{
  AliHLTUInt32_t estimatedSize = (trackArray->GetEntriesFast()+1)*(sizeof(AliHLTTRDTrack)+6*(sizeof(AliHLTTRDTracklet)+30*sizeof(AliHLTTRDCluster)));
  AliHLTUInt8_t* pBlock = (AliHLTUInt8_t*)malloc(estimatedSize);
  AliHLTUInt32_t size = AddTracksToOutput(trackArray, pBlock);
  trackArray->Delete();
  ReadTracks(trackArray, pBlock, size);
  free(pBlock);
}

AliHLTUInt32_t AliHLTTRDUtils::GetSM(AliHLTUInt32_t spec)
{
  spec = (spec&-spec); // use only least significant bit
  spec -= 1;           // as spec is now power of 2, this creates ones..
  int count = 0;       // .. which are counted
  while (spec) {
    count++;
    spec &= spec - 1;
  }
  return count;
}
