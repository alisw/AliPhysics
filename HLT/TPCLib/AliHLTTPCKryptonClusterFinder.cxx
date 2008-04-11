// $Id$
/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Anders Vestbo, Constantin Loizides, Jochen Thaeder    *
 *                  Kenneth Aamodt <kenneth.aamodt@student.uib.no>        *
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
    @author  Kenneth Aamodt kenneth.aamodt@student.uib.no
    @date   
    @brief  Krypton Cluster Finder for the TPC
*/

#include "AliHLTTPCDigitReader.h"
#include "AliHLTTPCRootTypes.h"
#include "AliHLTTPCLogging.h"
#include "AliHLTTPCKryptonClusterFinder.h"
#include "AliHLTTPCDigitData.h"
#include "AliHLTTPCTransform.h"
#include "AliHLTTPCSpacePointData.h"
#include "AliHLTTPCMemHandler.h"
#include "AliHLTTPCPad.h"
#include <sys/time.h>

#if __GNUC__ >= 3
using namespace std;
#endif

ClassImp(AliHLTTPCKryptonClusterFinder)

AliHLTTPCKryptonClusterFinder::AliHLTTPCKryptonClusterFinder()
  :
  fSpacePointData(NULL),
  fDigitReader(NULL),
  fPtr(NULL),
  fSize(0),
  fFirstRow(0),
  fLastRow(0),
  fCurrentRow(0),
  fCurrentSlice(0),
  fCurrentPatch(0),
  fMatch(1),
  fThreshold(10),
  fNClusters(0),
  fMaxNClusters(0),
  fXYErr(0.2),
  fZErr(0.3),
  fVectorInitialized(kFALSE),
  fRowPadVector(),
  fClusters(),
  fTimebinsInBunch(),
  fIndexOfBunchStart(),
  fNumberOfPadsInRow(NULL),
  fNumberOfRows(0),
  fRowOfFirstCandidate(0)
{
  //constructor  
}

AliHLTTPCKryptonClusterFinder::~AliHLTTPCKryptonClusterFinder()
{
  //destructor
  if(fVectorInitialized){
    DeInitializePadArray();
  }
  if(fNumberOfPadsInRow){
    delete [] fNumberOfPadsInRow;
    fNumberOfPadsInRow=NULL;
  }
}
 
void AliHLTTPCKryptonClusterFinder::InitializePadArray(){
  // see header file for class documentation
  
  if(fCurrentPatch>5||fCurrentPatch<0){
    HLTFatal("Patch is not set");
    return;
  }

  HLTDebug("Patch number=%d",fCurrentPatch);

  fFirstRow = AliHLTTPCTransform::GetFirstRow(fCurrentPatch);
  fLastRow = AliHLTTPCTransform::GetLastRow(fCurrentPatch);

  fNumberOfRows=fLastRow-fFirstRow+1;
  fNumberOfPadsInRow= new UInt_t[fNumberOfRows];

  memset( fNumberOfPadsInRow, 0, sizeof(Int_t)*(fNumberOfRows));

  for(UInt_t i=0;i<fNumberOfRows;i++){
    fNumberOfPadsInRow[i]=AliHLTTPCTransform::GetNPads(i+fFirstRow);
    AliHLTTPCPadVector tmpRow;
    for(UInt_t j=0;j<fNumberOfPadsInRow[i];j++){
      AliHLTTPCPad *tmpPad = new AliHLTTPCPad(2);
      tmpPad->SetID(i,j);
      tmpRow.push_back(tmpPad);
    }
    fRowPadVector.push_back(tmpRow);
  }
  fVectorInitialized=kTRUE;
}

Int_t AliHLTTPCKryptonClusterFinder::DeInitializePadArray()
{
  // see header file for class documentation
  for(UInt_t i=0;i<fNumberOfRows;i++){
    for(UInt_t j=0;j<fNumberOfPadsInRow[i];j++){
      delete fRowPadVector[i][j];
      fRowPadVector[i][j]=NULL;
    }
    fRowPadVector[i].clear();
  }
  fRowPadVector.clear();
  return 1;
} 

void AliHLTTPCKryptonClusterFinder::InitSlice(Int_t slice,Int_t patch,Int_t firstrow, Int_t lastrow,Int_t nmaxpoints)
{
  //init slice
  fNClusters = 0;
  fMaxNClusters = nmaxpoints;
  fCurrentSlice = slice;
  fCurrentPatch = patch;
  fFirstRow = firstrow;
  fLastRow = lastrow;
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
      cluster->fLastMergedPad=candidate->fPad;

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

void AliHLTTPCKryptonClusterFinder::FindNormalClusters()
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
	//	if(abs((Int_t)(tmpCluster->fPad/tmpCluster->fTotalCharge) - (Int_t)(nextCluster->fPad/nextCluster->fTotalCharge))<3){ // checks if the pad numbers are ok
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

void AliHLTTPCKryptonClusterFinder::PrintClusters()
{
  // see header file for class documentation
  for(size_t i=0;i<fClusters.size();i++){
    HLTInfo("Cluster number: %d",i);
    HLTInfo("Row: %d \t Pad: %d",fClusters[i].fRowNumber,fClusters[i].fFirstPad);
    HLTInfo("Total Charge:   %d",fClusters[i].fTotalCharge);
    HLTInfo("fPad:           %d",fClusters[i].fPad);
    HLTInfo("PadError:       %d",fClusters[i].fPad2);
    HLTInfo("TimeMean:       %d",fClusters[i].fTime);
    HLTInfo("TimeError:      %d",fClusters[i].fTime2);
    HLTInfo("EndOfCluster:");
  }
}
