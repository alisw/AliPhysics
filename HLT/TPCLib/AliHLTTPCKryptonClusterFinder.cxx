// $Id$
//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Kenneth Aamodt <kenneth.aamodt@student.uib.no>        *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/** @file    AliHLTTPCKryptonClusterFinder.cxx
    @author  Kenneth Aamodt, Kalliopi Kanaki
    @date   
    @brief  Krypton Cluster Finder for the TPC
*/

#include "AliHLTTPCDigitReader.h"
#include "AliHLTTPCKryptonClusterFinder.h"
#include "AliHLTTPCTransform.h"
#include "AliHLTTPCPad.h"
#include "AliHLTTPCClusters.h"
#include "TFile.h"
#include "AliHLTTPCSpacePointData.h"


#if __GNUC__ >= 3
using namespace std;
#endif

ClassImp(AliHLTTPCKryptonClusterFinder)

AliHLTTPCKryptonClusterFinder::AliHLTTPCKryptonClusterFinder()
  : AliHLTTPCClusterFinder(),
    fHWAddressVector(),
    fTimebinsInBunch(),
    fIndexOfBunchStart(),
    fMaxQOfCluster(0),
    fSelectionMinRowNumber(0),
    fSelectionMaxRowNumber(0),
    fNKryptonClusters(0),
    fMaxOutputSize(0)
{
  //constructor  
}

AliHLTTPCKryptonClusterFinder::~AliHLTTPCKryptonClusterFinder(){
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
  
  fHWAddressVector.clear();
  fNKryptonClusters=0;
  
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
	  candidate.fQMax=0;
	  if(fTimebinsInBunch[r]>2){
	    for(Int_t i=0;i<fTimebinsInBunch[r];i++){
	      candidate.fTotalCharge+=bunchData[i + fIndexOfBunchStart[r]];	
	      candidate.fTime += time*bunchData[i + fIndexOfBunchStart[r]];
	      candidate.fTime2 += time*time*bunchData[i + fIndexOfBunchStart[r]];
	      if(bunchData[i + fIndexOfBunchStart[r]]>candidate.fQMax){
		candidate.fQMax = bunchData[i + fIndexOfBunchStart[r]];
	      }
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

Bool_t AliHLTTPCKryptonClusterFinder::ComparePads(AliHLTTPCPad *nextPad,AliHLTTPCClusters* cluster,Int_t nextPadToRead){
  //Checking if we have a match on the next pad
  for(UInt_t candidateNumber=0;candidateNumber<nextPad->fClusterCandidates.size();candidateNumber++){
    AliHLTTPCClusters *candidate =&nextPad->fClusterCandidates[candidateNumber]; 
    if(cluster->fMean-candidate->fMean==1 || candidate->fMean-cluster->fMean==1 || cluster->fMean-candidate->fMean==0){
      cluster->fMean=candidate->fMean;
      cluster->fTotalCharge+=candidate->fTotalCharge;
      cluster->fTime += candidate->fTime;
      cluster->fTime2 += candidate->fTime2;
      cluster->fPad+=candidate->fPad;
      cluster->fPad2=candidate->fPad2;
      cluster->fLastMergedPad=nextPad->GetPadNumber();
      if(candidate->fQMax>cluster->fQMax){
	cluster->fQMax = candidate->fQMax;
      }

      //setting the matched pad to used
      nextPad->fUsedClusterCandidates[candidateNumber]=1;
      nextPadToRead++;
      if(nextPadToRead<(Int_t)fNumberOfPadsInRow[fRowOfFirstCandidate]){
	nextPad=fRowPadVector[fRowOfFirstCandidate][nextPadToRead];
	ComparePads(nextPad,cluster,nextPadToRead);
      }
      else{
	return kFALSE;
      }
    }
    else{
      return kFALSE;
    }
  }
  return kFALSE;
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
	tmpCandidate->fFirstPad=pad;
	UInt_t tmpTotalCharge=tmpCandidate->fTotalCharge;
	
	ComparePads(fRowPadVector[row][pad+1],tmpCandidate,pad+1);
	if(tmpCandidate->fTotalCharge>tmpTotalCharge){
	  //we have a cluster
	  tmpCandidate->fRowNumber= tmpCandidate->fRowNumber*tmpCandidate->fTotalCharge;
	  //	  tmpCandidate->fPad = tmpCandidate->fPad*tmpCandidate->fTotalCharge;
	  fClusters.push_back(*tmpCandidate);
	}
      }
    }
  }

  HLTDebug("Found %d normal clusters.",fClusters.size());
}

void AliHLTTPCKryptonClusterFinder::FindKryptonClusters()
{
  fMaxQOfCluster=0;
  if(fClusters.size()<2){
    return;
  }
  for(UInt_t i=0; i<fClusters.size();i++){
    AliHLTTPCClusters * tmpCluster=&fClusters[i];
    if(tmpCluster->fFlags==99){//quickfix to check if a cluster is used already
      continue;
    }

    UInt_t prevRow=tmpCluster->fRowNumber/tmpCluster->fTotalCharge;

    if((Int_t)tmpCluster->fQMax>fMaxQOfCluster){
      fMaxQOfCluster = tmpCluster->fQMax;
    }
    
    //adds "normal" clusters belonging to the krypton cluster
    for(UInt_t j=i+1;j<fClusters.size();j++){
      AliHLTTPCClusters * nextCluster=&fClusters[j];

      if(nextCluster->fFlags==99){//quickfix to check if a cluster is used already
	continue;
      }
      if(prevRow == (UInt_t)(nextCluster->fRowNumber/nextCluster->fTotalCharge)-1){//Checks if the row numbers are ok (next to eachother)
	if(abs((Int_t)(tmpCluster->fPad/tmpCluster->fTotalCharge) - (Int_t)(nextCluster->fPad/nextCluster->fTotalCharge))<3){ // checks if the pad numbers are ok
	  if(abs((Int_t)tmpCluster->fMean-(Int_t)nextCluster->fMean)<2){
	    tmpCluster->fMean=nextCluster->fMean;
	    tmpCluster->fTotalCharge+=nextCluster->fTotalCharge;
	    tmpCluster->fRowNumber+=nextCluster->fRowNumber;
	    tmpCluster->fPad+=nextCluster->fPad;
	    tmpCluster->fTime+=nextCluster->fTime;
	    if((Int_t)nextCluster->fQMax>fMaxQOfCluster){
	      fMaxQOfCluster = nextCluster->fQMax;
	    }
	    	    


	    if(tmpCluster->fFlags!=99){//means that this is the first time normal clusters match
	      CheckForCandidateOnPreviousRow(tmpCluster);
	      Int_t minFirst=0;
	      Int_t maxFirst=0;
	      if(tmpCluster->fFirstPad>1){
		minFirst=2;
	      }
	      if(tmpCluster->fLastMergedPad+2<(UInt_t)AliHLTTPCTransform::GetNPads(prevRow)){
		maxFirst=2;
	      }
		
	      for(UInt_t ap = tmpCluster->fFirstPad -minFirst; ap<=tmpCluster->fLastMergedPad+maxFirst; ap++){
		fHWAddressVector.push_back((AliHLTUInt16_t)fDigitReader->GetAltroBlockHWaddr(tmpCluster->fRowNumber/tmpCluster->fTotalCharge-AliHLTTPCTransform::GetFirstRow(fCurrentPatch),ap));
	      }
	    }
	      
	    UInt_t minNext=0;
	    UInt_t maxNext=0;
	    if(nextCluster->fFirstPad>1){
	      minNext=2;
	    }
	    if(nextCluster->fLastMergedPad+2<(UInt_t)AliHLTTPCTransform::GetNPads(prevRow+1)){
	      maxNext=2;
	    }
	    for(UInt_t ap = nextCluster->fFirstPad-minNext; ap<=nextCluster->fLastMergedPad+maxNext; ap++){
	      fHWAddressVector.push_back((AliHLTUInt16_t)fDigitReader->GetAltroBlockHWaddr(nextCluster->fRowNumber/nextCluster->fTotalCharge-AliHLTTPCTransform::GetFirstRow(fCurrentPatch),ap));
	      HLTDebug("Pushing back hw address %d from row: %d and Pad: %d",fDigitReader->GetAltroBlockHWaddr(nextCluster->fRowNumber/nextCluster->fTotalCharge-AliHLTTPCTransform::GetFirstRow(fCurrentPatch),ap),nextCluster->fRowNumber/nextCluster->fTotalCharge,ap);
	    }
	      
	    prevRow=nextCluster->fRowNumber/nextCluster->fTotalCharge;
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


	//convert the (AliHLTClusters) cluster to spacepointdata and add it to the output array.
	if (fNKryptonClusters*sizeof(AliHLTTPCSpacePointData)>=fMaxOutputSize*sizeof(AliHLTUInt32_t)) {
	  HLTWarning("Buffer too small too add more spacepoints: %d of %d byte(s) already used",fNKryptonClusters*sizeof(AliHLTTPCSpacePointData) , fMaxOutputSize*sizeof(AliHLTUInt32_t));
	  return;
	}
	if(tmpCluster->fTotalCharge>10){
	  Float_t xyz[3]={0,0,0};
	  Int_t thisrow=-1;
	  Int_t thissector=-1;
	  AliHLTTPCTransform::Slice2Sector(fCurrentSlice, (Int_t)(tmpCluster->fRowNumber/tmpCluster->fTotalCharge), thissector, thisrow);
	  AliHLTTPCTransform::Raw2Local(xyz, thissector, thisrow,(Float_t)(tmpCluster->fPad/tmpCluster->fTotalCharge),(Float_t)(tmpCluster->fTime/tmpCluster->fTotalCharge));
	  fSpacePointData[fNKryptonClusters].fID = fCurrentSlice*10 +fCurrentPatch;
	  fSpacePointData[fNKryptonClusters].fX = xyz[0];
	  fSpacePointData[fNKryptonClusters].fY = xyz[1];
	  fSpacePointData[fNKryptonClusters].fZ = xyz[2];
	  fSpacePointData[fNKryptonClusters].fCharge = tmpCluster->fTotalCharge;
	  fSpacePointData[fNKryptonClusters].fQMax = tmpCluster->fQMax;
	  fSpacePointData[fNKryptonClusters].fPadRow = tmpCluster->fRowNumber/tmpCluster->fTotalCharge;
	  HLTDebug("Krypton cluster found");
	  HLTDebug("xyz=[%f,%f,%f]",fSpacePointData[fNKryptonClusters].fX,fSpacePointData[fNKryptonClusters].fY,fSpacePointData[fNKryptonClusters].fZ);
	  HLTDebug("TotalCharge = %d    and   QMax = %d",fSpacePointData[fNKryptonClusters].fCharge,fSpacePointData[fNKryptonClusters].fQMax);
	  fNKryptonClusters++;
	  break;
	}
      }
    }
  }//end add "normal" clusters belonging to the krypton cluster

  //resets the candidates for every pad and the fClusters(row clusters)
  for(UInt_t row=0;row<fNumberOfRows;row++){
    for(UInt_t pad=0;pad<fNumberOfPadsInRow[row];pad++){
      fRowPadVector[row][pad]->fClusterCandidates.clear();
    }
  }
  fClusters.clear();
}

void AliHLTTPCKryptonClusterFinder::CheckForCandidateOnPreviousRow(AliHLTTPCClusters *tmpCluster){
  if(tmpCluster->fRowNumber>1){
    for(Int_t p=-1;p<2;p++){
      if(tmpCluster->fPad+p>0 && tmpCluster->fPad+p<fNumberOfPadsInRow[tmpCluster->fRowNumber/tmpCluster->fTotalCharge-1]){
	if(tmpCluster->fTotalCharge==0){
	  HLTDebug("Charge of tmpCluster in AliHLTTPCKryptonClusterFinder::CheckForCandidateOnPreviousRow is 0");
	  return;
	}
	if((Int_t)(tmpCluster->fRowNumber/tmpCluster->fTotalCharge-1-AliHLTTPCTransform::GetFirstRow(fCurrentPatch))<0){
	  HLTDebug("AliHLTTPCKryptonClusterFinder::CheckForCandidateOnPreviousRow:    Rownumber is below 0");
	  return;
	}
	if(tmpCluster->fRowNumber/tmpCluster->fTotalCharge-1-AliHLTTPCTransform::GetFirstRow(fCurrentPatch)>fNumberOfRows){
	  HLTDebug("AliHLTTPCKryptonClusterFinder::CheckForCandidateOnPreviousRow:    Rownumber is too high");
	  return;
	}
	AliHLTTPCPad *prevPad=fRowPadVector[tmpCluster->fRowNumber/tmpCluster->fTotalCharge-1-AliHLTTPCTransform::GetFirstRow(fCurrentPatch)][tmpCluster->fPad/tmpCluster->fTotalCharge+p];
	for(UInt_t i=0;i<prevPad->fClusterCandidates.size();i++){
	  if(abs((Int_t)prevPad->fClusterCandidates[i].fMean - (Int_t)tmpCluster->fMean)<2 && prevPad->fUsedClusterCandidates[i]==0){
	    tmpCluster->fTotalCharge += prevPad->fClusterCandidates[i].fTotalCharge;
	    fHWAddressVector.push_back((AliHLTUInt16_t)fDigitReader->GetAltroBlockHWaddr(tmpCluster->fRowNumber/tmpCluster->fTotalCharge-1-AliHLTTPCTransform::GetFirstRow(fCurrentPatch),tmpCluster->fPad/tmpCluster->fTotalCharge+p));
	      //	    fHWAddressVector.push_back((AliHLTUInt16_t)fDigitReader->GetAltroBlockHWaddr(prevPad->GetRowNumber(),prevPad->GetPadNumber()));
	    HLTDebug("Pushing back hw address %d",fDigitReader->GetAltroBlockHWaddr(tmpCluster->fRowNumber/tmpCluster->fTotalCharge-1-AliHLTTPCTransform::GetFirstRow(fCurrentPatch),tmpCluster->fPad/tmpCluster->fTotalCharge+p));
	  }
	}
      }
    }
  }
}

void AliHLTTPCKryptonClusterFinder::SetSelection(Int_t minRow, Int_t maxRow){
  fSelectionMinRowNumber=minRow;
  fSelectionMaxRowNumber=maxRow;
}

