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

//  @file   AliHLTTRDUtils.cxx
//  @author Theodor Rascanu
//  @date   
//  @brief  Utilities needed the HLT::TRD code.
// 

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// HLT TRD Utilities Class                                                   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliHLTTRDUtils.h"
#include <TClonesArray.h>
#include "AliHLTTRDTrack.h"
#include "AliHLTTRDTracklet.h"
#include "AliHLTTRDCluster.h"
#include "AliHLTExternalTrackParam.h"
#include "AliTRDtransform.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliPID.h"

ClassImp(AliHLTTRDUtils)

AliHLTUInt32_t AliHLTTRDUtils::AddClustersToOutput(const TClonesArray *const inClusterArray, AliHLTUInt8_t *const outBlockPtr, Int_t nTimeBins)
{
  AliTRDcluster* cluster = 0;
  AliHLTUInt32_t addedSize = 0;
  Int_t lastDet = -1;

  if (inClusterArray){
    AliHLTTRDClustersArray* clsArr = NULL;
    Int_t nbEntries = inClusterArray->GetEntries();
    for (Int_t iCluster = 0; iCluster<nbEntries; iCluster++){
      cluster = (AliTRDcluster*)(inClusterArray->At(iCluster));
      if(lastDet!=cluster->GetDetector()){
	lastDet=cluster->GetDetector();
	clsArr = new(outBlockPtr+addedSize) AliHLTTRDClustersArray(lastDet);
	addedSize += sizeof(AliHLTTRDClustersArray);
      }
      
      new (&clsArr->fCluster[clsArr->fCount]) AliHLTTRDClustersArray::cluster_type(cluster);
      clsArr->fCount++;
      addedSize += sizeof(AliHLTTRDClustersArray::cluster_type);
    }
  }

  Int_t *TBptr = (Int_t*)(outBlockPtr+addedSize);
  *TBptr = nTimeBins;
  
  addedSize += sizeof(*TBptr);

  return addedSize;
  
}

AliHLTUInt32_t AliHLTTRDUtils::AddTracksToOutput(const TClonesArray *const inTrackArray, AliHLTUInt8_t *const output, Int_t nTimeBins)
{

  Int_t *TBptr = (Int_t*)output;
  *TBptr = nTimeBins;

  AliTRDtrackV1* track = 0;
  AliHLTUInt32_t addedSize = sizeof(*TBptr);

  if (inTrackArray){
    Int_t nbTracks  = inTrackArray->GetEntries();
    for (Int_t iTrack = 0; iTrack<nbTracks; iTrack++){
      AliHLTUInt32_t trackSize=0;
      
      track = (AliTRDtrackV1*)(inTrackArray->At(iTrack));
      
      AliHLTTRDTrack *hltTrack = new (output+addedSize) AliHLTTRDTrack(track);
      trackSize = hltTrack->GetSize();
      addedSize += trackSize;
    }
  }
  return addedSize;
  
}

AliHLTUInt32_t AliHLTTRDUtils::AddTracksToOutputAlt(const TClonesArray *const inTrackArray, AliHLTUInt8_t *const block, Int_t nTimeBins)
{
  AliHLTUInt32_t addedSize = 0;

  AliHLTUInt64_t *TBptr = (AliHLTUInt64_t*)block;
  *TBptr = nTimeBins;

  addedSize += sizeof(AliHLTUInt64_t);

  if(!inTrackArray) return addedSize;

  Int_t nbTracks  = inTrackArray->GetEntriesFast();
  for (Int_t i = 0; i<nbTracks; i++){
    AliTRDtrackV1* inTrack = (AliTRDtrackV1*)(inTrackArray->At(i));
    if(inTrack)addedSize+=AliHLTTRDTrack::SaveAt(block+addedSize, inTrack);
  }

  return addedSize;
  
}

/**
 * Read cluster to the TClonesArray from the memory 
 */
//============================================================================
AliHLTUInt32_t AliHLTTRDUtils::ReadClusters(TClonesArray *const outArray, const void *const inputPtr, AliHLTUInt32_t size, Int_t* nTimeBins)
{
  const AliHLTUInt8_t* inPtr = (AliHLTUInt8_t*)inputPtr;
  UInt_t curSize = 0;
  Int_t counter = outArray->GetEntriesFast();

  if(nTimeBins){
    *nTimeBins=*(Int_t*)(inPtr+size-sizeof(*nTimeBins));
  }
  size-=sizeof(*nTimeBins);
#ifndef HAVE_NOT_ALITRD_CLUSTERIZER_r42837
  AliTRDtransform trans;
#endif
  while (curSize < size)
    {
      AliHLTTRDClustersArray* clsArr = (AliHLTTRDClustersArray*)(inPtr+curSize);
      curSize+=sizeof(AliHLTTRDClustersArray);
#ifndef HAVE_NOT_ALITRD_CLUSTERIZER_r42837
      trans.SetDetector(clsArr->fDetector);
#endif
      for(Int_t iCluster = 0; iCluster<clsArr->fCount; iCluster++){
	AliTRDcluster* curTRDCluster = new((*outArray)[counter]) AliTRDcluster();
	clsArr->fCluster[iCluster].ExportTRDCluster(curTRDCluster);
	curTRDCluster->SetDetector(clsArr->fDetector);
#ifndef HAVE_NOT_ALITRD_CLUSTERIZER_r42837
	trans.Transform(curTRDCluster);
#endif
	curSize += sizeof(AliHLTTRDClustersArray::cluster_type); 
	counter++;
      }
    }
  
  return counter;
}

AliHLTUInt32_t AliHLTTRDUtils::ReadTracks(TClonesArray *const outArray, const void *const inputPtr, AliHLTUInt32_t size, Int_t* nTimeBins)
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

 AliHLTUInt32_t AliHLTTRDUtils::ReadTracksAlt(TClonesArray *const outArray, const void *const inputPtr, AliHLTUInt32_t size, Int_t* nTimeBins)
 {
   const AliHLTUInt8_t *const block = ((AliHLTUInt8_t*)inputPtr);
   AliHLTUInt32_t readSize = 0;

   if(nTimeBins){
     *nTimeBins=*(AliHLTUInt64_t*)block;
     //HLTDebug("Reading number of time bins from input block: %d", *nTimeBins);
   }

   readSize += sizeof(AliHLTUInt64_t);

   if(!outArray) return readSize;

   Int_t counter=outArray->GetEntriesFast();
   while(readSize<size){
     AliTRDtrackV1 *const outTrack = new((*outArray)[counter]) AliTRDtrackV1;
     readSize+=AliHLTTRDTrack::LoadFrom(outTrack, block+readSize);
     counter++;
   }

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
  AliHLTUInt32_t estimatedSize = (clusterArray->GetEntriesFast()+1)*sizeof(AliHLTTRDClustersArray::cluster_type);
  AliHLTUInt8_t* pBlock = (AliHLTUInt8_t*)malloc(estimatedSize);
  AliHLTUInt32_t size = AddClustersToOutput(clusterArray, pBlock);
  clusterArray->Delete();
  ReadClusters(clusterArray, pBlock, size);
  free(pBlock);
}

void AliHLTTRDUtils::EmulateHLTTracks(TClonesArray* trackArray)
{
  AliHLTUInt32_t estimatedSize = (trackArray->GetEntriesFast()+1)*(sizeof(AliHLTTRDTrack)+6*(sizeof(AliHLTTRDTracklet)+30*sizeof(AliHLTTRDClustersArray::cluster_type)));
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
