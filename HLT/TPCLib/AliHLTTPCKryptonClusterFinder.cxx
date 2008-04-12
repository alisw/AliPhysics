// $Id$
/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Kenneth Aamodt <kenneth.aamodt@student.uib.no>        *
 *                  for The ALICE HLT Project.                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** @file    AliHLTTPCKryptonClusterFinder.cxx
    @author  Kenneth Aamodt, Kalliopi Kanaki
    @date   
    @brief  Krypton Cluster Finder for the TPC
*/

#include "AliHLTTPCDigitReader.h"
#include "AliHLTTPCKryptonClusterFinder.h"
#include "AliHLTTPCTransform.h"

#if __GNUC__ >= 3
using namespace std;
#endif

ClassImp(AliHLTTPCKryptonClusterFinder)

AliHLTTPCKryptonClusterFinder::AliHLTTPCKryptonClusterFinder()
  : AliHLTTPCClusterFinder(),
    fTimebinsInBunch(),
    fIndexOfBunchStart()
{
  //constructor  
}

void AliHLTTPCKryptonClusterFinder::ReBunch(const UInt_t *bunchData , Int_t bunchSize){
  fTimebinsInBunch.clear();
  fIndexOfBunchStart.clear();
  Bool_t newBunch=kTRUE;
  Int_t currentBunchNumber=-1;
  for(Int_t i=0;i<bunchSize;i++){
    if(newBunch){
      if(bunchData[i]>5){
	fIndexOfBunchStart.push_back(i);
	fTimebinsInBunch.push_back(1);
	currentBunchNumber++;
	newBunch=kFALSE;
      }
    }
    else{
      if(bunchData[i]>0){
	fTimebinsInBunch[currentBunchNumber]++;
      }
      else{
	newBunch=kTRUE;
      }
    }
  }
}

void AliHLTTPCKryptonClusterFinder::ReadDataUnsorted(void* ptr,unsigned long size)
{
  //set input pointer
  fPtr = (UChar_t*)ptr;
  fSize = size;

  if(!fVectorInitialized){
    InitializePadArray();
   }

  fDigitReader->InitBlock(fPtr,fSize,fFirstRow,fLastRow,fCurrentPatch,fCurrentSlice);
  
  while(fDigitReader->NextChannel()){
    UInt_t row=fDigitReader->GetRow();
    UInt_t pad=fDigitReader->GetPad();

    fRowPadVector[row][pad]->ClearCandidates();

    while(fDigitReader->NextBunch()){
      if(fDigitReader->GetBunchSize()>1){
	const UInt_t *bunchData= fDigitReader->GetSignals();

	ReBunch(bunchData,fDigitReader->GetBunchSize());
	Int_t rebunchCount=fIndexOfBunchStart.size();
	
	for(Int_t r=0;r<rebunchCount;r++){
	  UInt_t time = fDigitReader->GetTime()+fIndexOfBunchStart[r];
	  AliHLTTPCClusters candidate;
	  candidate.fTotalCharge=0;
	  if(fTimebinsInBunch[r]>2){
	    for(Int_t i=0;i<fTimebinsInBunch[r];i++){
	      candidate.fTotalCharge+=bunchData[i + fIndexOfBunchStart[r]];	
	      candidate.fTime += time*bunchData[i + fIndexOfBunchStart[r]];
	      candidate.fTime2 += time*time*bunchData[i + fIndexOfBunchStart[r]];
	      time++;
	    }
	    if(candidate.fTotalCharge>0){
	      candidate.fMean=candidate.fTime/candidate.fTotalCharge;
	      candidate.fPad=candidate.fTotalCharge*pad;
	      candidate.fPad2=candidate.fPad*pad;
	      candidate.fLastMergedPad=pad;
	      candidate.fRowNumber=row+fDigitReader->GetRowOffset();
	    }
	    if(candidate.fTotalCharge>10 && candidate.fMean<924){
	      fRowPadVector[row][pad]->AddClusterCandidate(candidate);
	    }
	  }
	}
      }
    }
  }
}

void AliHLTTPCKryptonClusterFinder::FindRowClusters()
{
  // see header file for function documentation

  AliHLTTPCClusters* tmpCandidate=NULL;
  for(UInt_t row=5;row<fNumberOfRows-5;row++){
    fRowOfFirstCandidate=row;
    for(UInt_t pad=5;pad<fNumberOfPadsInRow[row]-1-5;pad++){
      AliHLTTPCPad *tmpPad=fRowPadVector[row][pad];
      for(size_t candidate=0;candidate<tmpPad->fClusterCandidates.size();candidate++){
	if(tmpPad->fUsedClusterCandidates[candidate]){
	  continue;
	}
	if((Int_t)tmpPad->fClusterCandidates[candidate].fMean<100 || (Int_t)tmpPad->fClusterCandidates[candidate].fMean>AliHLTTPCTransform::GetNTimeBins()-100){
	  continue;
	}
	tmpCandidate=&tmpPad->fClusterCandidates[candidate];
	UInt_t tmpTotalCharge=tmpCandidate->fTotalCharge;
	ComparePads(fRowPadVector[row][pad+1],tmpCandidate,pad+1);
	if(tmpCandidate->fTotalCharge>tmpTotalCharge){
	  //we have a cluster
	  tmpCandidate->fPad=tmpCandidate->fPad/tmpCandidate->fTotalCharge;
	  fClusters.push_back(*tmpCandidate);
	}
      }
    }
  }

  HLTDebug("Found %d normal clusters.",fClusters.size());
}

void AliHLTTPCKryptonClusterFinder::FindKryptonClusters()
{
  if(fClusters.size()<2){
    return;
  }

  for(UInt_t i=0; i<fClusters.size();i++){
    AliHLTTPCClusters * tmpCluster=&fClusters[i];
    if(tmpCluster->fFlags==99){//quickfix to check if a cluster is used already
      continue;
    }
    
    //adds "normal" clusters belonging to the krypton cluster
    for(UInt_t j=i+1;j<fClusters.size();j++){
      AliHLTTPCClusters * nextCluster=&fClusters[j];

      if(nextCluster->fFlags==99){//quickfix to check if a cluster is used already
	continue;
      }

      if(tmpCluster->fRowNumber==nextCluster->fRowNumber-1){//Checks if the row numbers are ok (next to eachother)
	if(abs((Int_t)(tmpCluster->fPad) - (Int_t)(nextCluster->fPad))<3){ // checks if the pad numbers are ok
	  if(abs((Int_t)tmpCluster->fMean-(Int_t)nextCluster->fMean)<2){
	    tmpCluster->fMean=nextCluster->fMean;
	    tmpCluster->fTotalCharge+=nextCluster->fTotalCharge;
	    tmpCluster->fPad=nextCluster->fPad;
	    if(tmpCluster->fFlags!=99){//means that this is the first time normal clusters match
	      CheckForCandidateOnPreviousPad(tmpCluster);
	    }
	    tmpCluster->fRowNumber=nextCluster->fRowNumber;
	    nextCluster->fFlags=99;
	    tmpCluster->fFlags=99;
	    if(j!=fClusters.size()-1){
	      continue;
	    }
	    
	  }
	}
      }
      
      if(tmpCluster->fFlags==99){
	//adds a clustercandidate on next row if present TODO
	/*	if(tmpCluster->fFlags==99){
	  for(Int_t p=-1;p<2;p++){
	  AliHLTTPCPad *padAfter=fRowPadVector[tmpCluster->fRowNumber+1][tmpCluster->fPad+p];
	  for(UInt_t c=0;c<padAfter->fClusterCandidates.size();c++)
	  if(abs((Int_t)tmpCluster->fMean - (Int_t)padAfter->fClusterCandidates[c].fMean)<2){
	  tmpCluster->fTotalCharge+=padAfter->fClusterCandidates[c].fTotalCharge;
	  }
	  }//end add clustercandidate if present
	  }*/
	HLTInfo("Krypton cluster found, charge: %d   in patch number: %d",tmpCluster->fTotalCharge,fCurrentPatch);
	break;
      }


    }
  }//end add "normal" clusters belonging to the krypton cluster

}

void AliHLTTPCKryptonClusterFinder::CheckForCandidateOnPreviousPad(AliHLTTPCClusters *tmpCluster){
  if(tmpCluster->fRowNumber>1){
    for(Int_t p=-1;p<2;p++){
      if(tmpCluster->fPad+p>0 && tmpCluster->fPad+p<fNumberOfPadsInRow[tmpCluster->fRowNumber-1]){
	AliHLTTPCPad *prevPad=fRowPadVector[tmpCluster->fRowNumber-1-AliHLTTPCTransform::GetFirstRow(fCurrentPatch)][tmpCluster->fPad+p];
	for(UInt_t i=0;i<prevPad->fClusterCandidates.size();i++){
	  if(abs((Int_t)prevPad->fClusterCandidates[i].fMean - (Int_t)tmpCluster->fMean)<2 && prevPad->fUsedClusterCandidates[i]==0){
	    tmpCluster->fTotalCharge += prevPad->fClusterCandidates[i].fTotalCharge;
	  }
	}
      }
    }
  }
}
